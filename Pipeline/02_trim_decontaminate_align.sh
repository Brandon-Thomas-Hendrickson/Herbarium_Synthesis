#!/bin/bash
# =============================================================================
# 02_trim_decontaminate_align.sh
# Complete read-processing pipeline:
#   fastp (QC + trim) → bowtie2 (decontamination) → bwa mem (alignment)
#   → samtools sort/fixmate/markdup → coverage-based downsampling
#
# Prerequisites:
#   - config.sh edited with correct paths
#   - Reference genome indexed (bwa index $GENOME)
#   - Contamination index built (01_build_database.sh)
#   - Raw FASTQ files in $RAW named <sample>_1.fq.gz / <sample>_2.fq.gz
#
# Outputs (per sample in $ALIGNDIR):
#   <sample>_downsampled.bam + .bai   – final analysis-ready BAM
#   <sample>_bamstats.txt             – bamtools statistics
#   <sample>_qualimap/                – qualimap report
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

mkdir -p "$FASTQCDIR" "$TRIMDIR" "$CONTDIR" "$DECONDIR" "$ALIGNDIR"

# ---------------------------------------------------------------------------
# 1. Index reference genome (skip if index already exists)
# ---------------------------------------------------------------------------
if [[ ! -f "${GENOME}.bwt" ]]; then
    echo "[bwa] Indexing reference genome..."
    bwa index "$GENOME"
fi

# ---------------------------------------------------------------------------
# 2. Build sample list from raw FASTQ file names
# ---------------------------------------------------------------------------
echo "[list] Building sample list from $RAW..."
ls "$RAW" | grep "${S1}$" | sed "s/${S1}//" | sort -u > "$SAMPLE_LIST"
echo "  Found $(wc -l < "$SAMPLE_LIST") samples."

# ---------------------------------------------------------------------------
# 3. FastQC on raw reads
# ---------------------------------------------------------------------------
echo "[fastqc] Running FastQC on raw reads..."
for SAMPLE in $(cat "$SAMPLE_LIST"); do
    fastqc "$RAW/${SAMPLE}${S1}" "$RAW/${SAMPLE}${S2}" -o "$FASTQCDIR" -t 2
done

# ---------------------------------------------------------------------------
# 4. Trim with fastp
# ---------------------------------------------------------------------------
echo "[fastp] Trimming reads..."
for SAMPLE in $(cat "$SAMPLE_LIST"); do
    fastp \
        -i "$RAW/${SAMPLE}${S1}" \
        -I "$RAW/${SAMPLE}${S2}" \
        -o "$TRIMDIR/${SAMPLE}${S1}" \
        -O "$TRIMDIR/${SAMPLE}${S2}" \
        --cut_right \
        --dedup \
        -h "$TRIMDIR/${SAMPLE}.html" \
        -g \
        -w "$THREADS_ALIGN"
done

# ---------------------------------------------------------------------------
# 5. Decontamination with bowtie2
#    Reads aligning to the bacterial/fungal database are saved as contaminated.
#    Unaligned read pairs (clean) are written to $DECONDIR.
# ---------------------------------------------------------------------------
echo "[bowtie2] Decontaminating reads..."
for SAMPLE in $(cat "$SAMPLE_LIST"); do
    bowtie2 -p "$THREADS_ALIGN" \
        -x "$CONTGENOME" \
        -1 "$TRIMDIR/${SAMPLE}${S1}" \
        -2 "$TRIMDIR/${SAMPLE}${S2}" \
        --un-conc-gz "$DECONDIR/${SAMPLE}_CLEAN" \
        > "${SAMPLE}_FUN_BAC_REMOVED.sam"

    # Save contaminated reads for reference
    samtools view -hbS -F 4 "${SAMPLE}_FUN_BAC_REMOVED.sam" \
        > "$CONTDIR/${SAMPLE}_CONTAMINATED.bam"
    rm "${SAMPLE}_FUN_BAC_REMOVED.sam"

    # Rename bowtie2 --un-conc-gz output to standard naming convention
    mv "$DECONDIR/${SAMPLE}_CLEAN.1" "$DECONDIR/${SAMPLE}${S1}"
    mv "$DECONDIR/${SAMPLE}_CLEAN.2" "$DECONDIR/${SAMPLE}${S2}"
done

# ---------------------------------------------------------------------------
# 6. Align clean reads with bwa mem
# ---------------------------------------------------------------------------
echo "[bwa] Aligning clean reads..."
BAM_LIST_TMP=$(mktemp)
for SAMPLE in $(cat "$SAMPLE_LIST"); do
    bwa mem -t "$THREADS_ALIGN" "$GENOME" \
        "$DECONDIR/${SAMPLE}${S1}" \
        "$DECONDIR/${SAMPLE}${S2}" \
        | samtools view -bS \
        | samtools sort -@ "$THREADS_SORT" \
        -o "$ALIGNDIR/${SAMPLE}.bam"
    echo "$ALIGNDIR/${SAMPLE}_downsampled.bam" >> "$BAM_LIST_TMP"
done
mv "$BAM_LIST_TMP" "$BAM_LIST"

# ---------------------------------------------------------------------------
# 7. samtools fixmate pipeline (requires name-sorted input)
# ---------------------------------------------------------------------------
echo "[samtools] Running fixmate → coordinate sort → markdup..."
for SAMPLE in $(cat "$SAMPLE_LIST"); do
    # 7a. Sort by query name (required by fixmate)
    samtools sort -n -@ "$THREADS_SORT" \
        "$ALIGNDIR/${SAMPLE}.bam" \
        -o "$ALIGNDIR/${SAMPLE}_grouped.bam"

    # 7b. Fix mate pair information
    samtools fixmate -@ "$THREADS_SORT" -m \
        "$ALIGNDIR/${SAMPLE}_grouped.bam" \
        "$ALIGNDIR/${SAMPLE}_fixmate.bam"

    # 7c. Coordinate sort (required by markdup)
    samtools sort -@ "$THREADS_SORT" \
        "$ALIGNDIR/${SAMPLE}_fixmate.bam" \
        -o "$ALIGNDIR/${SAMPLE}_fixmate_sorted.bam"

    # 7d. Mark duplicates
    samtools markdup -@ "$THREADS_SORT" \
        "$ALIGNDIR/${SAMPLE}_fixmate_sorted.bam" \
        "$ALIGNDIR/${SAMPLE}_marked.bam"

    # 7e. Index
    samtools index -@ "$THREADS_SORT" "$ALIGNDIR/${SAMPLE}_marked.bam"

    # Clean up intermediates
    rm -f "$ALIGNDIR/${SAMPLE}_grouped.bam" \
          "$ALIGNDIR/${SAMPLE}_fixmate.bam" \
          "$ALIGNDIR/${SAMPLE}_fixmate_sorted.bam"
done

# ---------------------------------------------------------------------------
# 8. QC stats: bamtools + qualimap
# ---------------------------------------------------------------------------
echo "[stats] Generating BAM statistics..."
for SAMPLE in $(cat "$SAMPLE_LIST"); do
    bamtools stats -in "$ALIGNDIR/${SAMPLE}_marked.bam" \
        > "$ALIGNDIR/${SAMPLE}_bamstats.txt"
    qualimap bamqc \
        -bam "$ALIGNDIR/${SAMPLE}_marked.bam" \
        -outdir "$ALIGNDIR/${SAMPLE}_qualimap" \
        -outfile "${SAMPLE}_qualimap_report.txt" \
        -outformat TXT
done

# ---------------------------------------------------------------------------
# 9. Coverage-based downsampling
#    Target: 2–4X mean coverage.
#    Samples already ≤ 4X are renamed without modification.
#    Samples with < 2X coverage are kept as-is (downsampling not possible).
# ---------------------------------------------------------------------------
echo "[downsample] Normalising coverage to 2–4X..."
for SAMPLE in $(cat "$SAMPLE_LIST"); do
    MARKED="$ALIGNDIR/${SAMPLE}_marked.bam"
    FINAL="$ALIGNDIR/${SAMPLE}_downsampled.bam"

    COVERAGE=$(samtools depth "$MARKED" \
        | awk '{sum+=$3; n++} END {if(n>0) print sum/n; else print 0}')

    if (( $(echo "$COVERAGE > 4" | bc -l) )); then
        # Iteratively subsample until coverage drops into 2–4X window
        cp "$MARKED" "$FINAL"
        while (( $(echo "$COVERAGE > 4" | bc -l) )); do
            FRAC=$(echo "3 / $COVERAGE" | bc -l)
            # bc -l returns a fraction; clamp to (0,1)
            FRAC=$(echo "if ($FRAC > 1) 1 else $FRAC" | bc -l)
            samtools view -s "$FRAC" -b "$FINAL" \
                -o "$ALIGNDIR/${SAMPLE}_tmp.bam"
            mv "$ALIGNDIR/${SAMPLE}_tmp.bam" "$FINAL"
            COVERAGE=$(samtools depth "$FINAL" \
                | awk '{sum+=$3; n++} END {if(n>0) print sum/n; else print 0}')
        done
        echo "  $SAMPLE: downsampled to ${COVERAGE}X"
    else
        cp "$MARKED" "$FINAL"
        echo "  $SAMPLE: coverage ${COVERAGE}X – no downsampling needed"
    fi
    samtools index "$FINAL"
done

echo "Done. Final BAM files are in $ALIGNDIR (*_downsampled.bam)."
