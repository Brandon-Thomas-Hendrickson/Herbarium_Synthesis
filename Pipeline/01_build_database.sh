#!/bin/bash
# =============================================================================
# 01_build_database.sh
# Build a bowtie2-indexed contamination database from NCBI bacterial and
# fungal complete genomes.
#
# Run ONCE before processing any samples.
# Outputs:
#   $GENOMEDIR/GB_BAC_FUN.fasta        – concatenated reference FASTA
#   $GENOMEDIR/all_genera_index.*      – bowtie2 index files
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

mkdir -p "$GENOMEDIR" BAC_FUN_GENOMES

# ---------------------------------------------------------------------------
# 1. Download NCBI assembly summaries
# ---------------------------------------------------------------------------
echo "[1/5] Downloading NCBI assembly summaries..."
wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt \
    -O assembly_summary_bac.txt
wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt \
    -O assembly_summary_fun.txt

# ---------------------------------------------------------------------------
# 2. Filter for complete genomes only (field 12 = assembly_level,
#    field 20 = ftp_path in the NCBI assembly_summary format)
# ---------------------------------------------------------------------------
echo "[2/5] Filtering for complete genomes..."
awk -F '\t' '$12=="Complete Genome" {print $8, $20}' assembly_summary_bac.txt \
    > assembly_complete_bac.txt
awk -F '\t' '$12=="Complete Genome" {print $8, $20}' assembly_summary_fun.txt \
    > assembly_complete_fun.txt

# Combine and get one representative per genus (first word of organism name)
cat assembly_complete_bac.txt assembly_complete_fun.txt \
    | sort -u -k1,1 \
    | awk '{print $2}' \
    > ftp_list.txt

# ---------------------------------------------------------------------------
# 3. Download genomes
# ---------------------------------------------------------------------------
echo "[3/5] Downloading complete genomes (this may take several hours)..."
while IFS= read -r ftp_path; do
    # Skip lines that are not valid ftp/https paths
    [[ "$ftp_path" =~ ^ftp ]] || continue
    fname="${ftp_path##*/}"
    wget -q "${ftp_path}/${fname}_genomic.fna.gz" -P BAC_FUN_GENOMES/ || true
done < ftp_list.txt

# ---------------------------------------------------------------------------
# 4. Decompress and concatenate
# ---------------------------------------------------------------------------
echo "[4/5] Decompressing and concatenating genomes..."
gunzip -f BAC_FUN_GENOMES/*.gz
cat BAC_FUN_GENOMES/*.fna > "$GENOMEDIR/GB_BAC_FUN.fasta"
echo "Combined FASTA written to $GENOMEDIR/GB_BAC_FUN.fasta"

# ---------------------------------------------------------------------------
# 5. Build bowtie2 index
#    NOTE: NCBI recommends bowtie2 >= 2.2.4 for large databases.
#    The index is built with --packed to reduce memory usage.
# ---------------------------------------------------------------------------
echo "[5/5] Building bowtie2 index (this may take several hours)..."
bowtie2-build --packed "$GENOMEDIR/GB_BAC_FUN.fasta" "$GENOMEDIR/all_genera_index"
echo "Bowtie2 index written to $GENOMEDIR/all_genera_index"

echo "Done. Contamination database is ready."
