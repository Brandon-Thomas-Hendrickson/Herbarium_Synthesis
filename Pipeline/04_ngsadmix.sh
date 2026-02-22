#!/bin/bash
# =============================================================================
# 04_ngsadmix.sh
# Population admixture analysis pipeline:
#   1. ngsRelate  – identify and remove related individuals (r > 0.25)
#   2. degenotate – identify 4-fold degenerate (neutral) sites
#   3. NGSadmix   – K=1..K_MAX, 5 random seeds each
#   4. evalAdmix  – cross-validation for each run
#   5. Summarise log-likelihoods for best-K selection
#
# Input:  $ANGSDDIR/Herbarium.beagle.gz  (from step 03)
#         $ANGSDDIR/Herbarium.mafs.gz    (from step 03)
# Output: $ADMIXDIR/clumppak_ready.txt   (K + log-likelihood table)
#         $ADMIXDIR/ADMIX{K}.{seed}.Q    (ancestry proportions per run)
#
# SLURM directives below; adjust account/partition in config.sh.
# =============================================================================

#SBATCH --job-name ngsadmix
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH -t 3-00:00:00
#SBATCH --account ${SLURM_ACCOUNT}
#SBATCH --partition ${SLURM_PARTITION}

set -euo pipefail
source "$(dirname "$0")/config.sh"

mkdir -p "$ADMIXDIR" "$DEGENDIR"

BEAGLE_RAW="$ANGSDDIR/Herbarium.beagle.gz"
BEAGLE_FILTERED="$ANGSDDIR/Herbarium_filtered.beagle.gz"
BEAGLE_4FOLD="$ANGSDDIR/Herbarium_4fold.beagle.gz"

# ---------------------------------------------------------------------------
# 1. Build ngsRelate (compile from source if not present)
# ---------------------------------------------------------------------------
if [[ ! -x "$NGSRELATE" ]]; then
    echo "[ngsRelate] Building from source..."
    SOFTDIR=$(dirname "$NGSRELATE")
    mkdir -p "$(dirname "$SOFTDIR")"
    cd "$(dirname "$SOFTDIR")"
    git clone https://github.com/samtools/htslib
    git clone https://github.com/ANGSD/ngsRelate
    cd htslib && make -j4 && cd ../ngsRelate && make HTSSRC=../htslib/
    cd "$OLDPWD"
fi

# ---------------------------------------------------------------------------
# 2. Run ngsRelate to detect related individuals
# ---------------------------------------------------------------------------
echo "[ngsRelate] Computing pairwise relatedness..."
zcat "$ANGSDDIR/Herbarium.mafs.gz" | cut -f5 | tail -n +2 \
    > "$ANGSDDIR/freq.txt"

N_SAMPLES=$(zcat "$BEAGLE_RAW" | head -1 | awk '{print (NF-3)/3}')

"$NGSRELATE" \
    -g "$BEAGLE_RAW" \
    -n "$N_SAMPLES" \
    -f "$ANGSDDIR/freq.txt" \
    -O "$ANGSDDIR/ngsrelate_out.txt" \
    -m 1 \
    -p "$THREADS_ADMIX" \
    -i 5000

# ---------------------------------------------------------------------------
# 3. Identify individuals to remove (relatedness > 0.25)
#    Outputs: $ANGSDDIR/individuals_to_remove.txt (0-based sample indices)
# ---------------------------------------------------------------------------
echo "[relatedness] Filtering related individuals..."
python3 - <<'PYEOF'
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd

input_file  = os.environ["ANGSDDIR"] + "/ngsrelate_out.txt"
output_file = os.environ["ANGSDDIR"] + "/individuals_to_remove.txt"

df = pd.read_csv(input_file, sep="\t")
to_remove = set()
for _, row in df.iterrows():
    if row["rab"] > 0.25:
        a, b = int(row["a"]), int(row["b"])
        if a not in to_remove:
            to_remove.add(b)

print(f"  Flagging {len(to_remove)} related individuals for removal.")
with open(output_file, "w") as f:
    for idx in sorted(to_remove):
        f.write(f"{idx}\n")
PYEOF

# ---------------------------------------------------------------------------
# 4. Filter related individuals from Beagle (column-based removal)
#    Beagle format: marker | allele1 | allele2 | [3 cols per individual ...]
#    Individual i (0-based) occupies columns: 3+3i, 3+3i+1, 3+3i+2
# ---------------------------------------------------------------------------
echo "[beagle] Removing related individuals from Beagle file..."
python3 - <<'PYEOF'
import gzip, sys, os

beagle_in  = os.environ["BEAGLE_RAW"]
beagle_out = os.environ["BEAGLE_FILTERED"]
remove_file = os.environ["ANGSDDIR"] + "/individuals_to_remove.txt"

with open(remove_file) as f:
    remove_set = set(int(line.strip()) for line in f if line.strip())

with gzip.open(beagle_in, "rt") as fin, gzip.open(beagle_out, "wt") as fout:
    for line_num, line in enumerate(fin):
        fields = line.rstrip("\n").split("\t")
        # Columns 0,1,2 are marker/allele1/allele2; then triplets per sample
        keep = [0, 1, 2]
        n_samples = (len(fields) - 3) // 3
        for i in range(n_samples):
            if i not in remove_set:
                base = 3 + 3 * i
                keep += [base, base + 1, base + 2]
        fout.write("\t".join(fields[c] for c in keep) + "\n")
PYEOF

# Export for embedded Python scripts
export BEAGLE_RAW BEAGLE_FILTERED BEAGLE_4FOLD ANGSDDIR ADMIXDIR

# ---------------------------------------------------------------------------
# 5. Get 4-fold degenerate sites with degenotate
# ---------------------------------------------------------------------------
echo "[degenotate] Identifying 4-fold degenerate sites..."
"$DEGENOTATE" -a "$ANNOTATION" -g "$GENOME" -o "$DEGENDIR" -d " "

# Build list of 4-fold site IDs in "Chr_Pos" format (matching Beagle col 1)
awk '$5==4 {print $1"_"$2}' "$DEGENDIR/degeneracy-all-sites.bed" \
    > "$DEGENDIR/degen_4fold.txt"

# ---------------------------------------------------------------------------
# 6. Extract 4-fold sites from filtered Beagle
# ---------------------------------------------------------------------------
echo "[beagle] Extracting 4-fold degenerate sites..."
python3 - <<'PYEOF'
import gzip, os

beagle_in  = os.environ["BEAGLE_FILTERED"]
beagle_out = os.environ["BEAGLE_4FOLD"]
sites_file = os.environ["DEGENDIR"] + "/degen_4fold.txt"

with open(sites_file) as f:
    keep_sites = set(line.strip() for line in f if line.strip())

with gzip.open(beagle_in, "rt") as fin, gzip.open(beagle_out, "wt") as fout:
    for line in fin:
        marker = line.split("\t", 1)[0]
        if marker == "marker" or marker in keep_sites:
            fout.write(line)
PYEOF

# ---------------------------------------------------------------------------
# 7. Run NGSadmix K=1..K_MAX, multiple seeds
# ---------------------------------------------------------------------------
echo "[NGSadmix] Running admixture analysis (K=1..$K_MAX, ${#SEEDS[@]} seeds each)..."
for K in $(seq 1 "$K_MAX"); do
    for SEED in "${SEEDS[@]}"; do
        echo "  K=${K} seed=${SEED}"
        "$NGSADMIX" \
            -likes "$BEAGLE_4FOLD" \
            -K "$K" \
            -P "$THREADS_ADMIX" \
            -o "$ADMIXDIR/ADMIX${K}.${SEED}" \
            -minMaf "$MIN_MAF" \
            -seed "$SEED"
    done
done

# ---------------------------------------------------------------------------
# 8. Build evalAdmix (compile from source if not present)
# ---------------------------------------------------------------------------
if [[ ! -x "$EVALADMIX" ]]; then
    echo "[evalAdmix] Building from source..."
    cd "$(dirname "$EVALADMIX")"
    git clone https://github.com/GenisGE/evalAdmix.git "$(basename "$(dirname "$EVALADMIX")")"
    cd "$(basename "$(dirname "$EVALADMIX")")"
    make
    cd "$OLDPWD"
fi

# ---------------------------------------------------------------------------
# 9. Run evalAdmix cross-validation
# ---------------------------------------------------------------------------
echo "[evalAdmix] Running cross-validation..."
cd "$ADMIXDIR"
for K in $(seq 1 "$K_MAX"); do
    for SEED in "${SEEDS[@]}"; do
        "$EVALADMIX" \
            -beagle "$BEAGLE_4FOLD" \
            -q "ADMIX${K}.${SEED}.Q" \
            -P "$THREADS_ADMIX" \
            -o "evalADMIX${K}.${SEED}"
    done
done
cd "$OLDPWD"

# ---------------------------------------------------------------------------
# 10. Extract log-likelihoods for best-K selection
# ---------------------------------------------------------------------------
echo "[summary] Extracting log-likelihoods..."
> "$ADMIXDIR/K_likelihood.txt"
for K in $(seq 1 "$K_MAX"); do
    for SEED in "${SEEDS[@]}"; do
        LOG="$ADMIXDIR/ADMIX${K}.${SEED}.log"
        grep -A1 "best like=" "$LOG" \
            | awk -v K="$K" -v S="$SEED" '{print K, S, $2}' \
            >> "$ADMIXDIR/K_likelihood.txt"
    done
done
sed -i 's/like=//g' "$ADMIXDIR/K_likelihood.txt"
awk '{print $1, $3}' "$ADMIXDIR/K_likelihood.txt" \
    > "$ADMIXDIR/clumppak_ready.txt"

echo "Done."
echo "  K log-likelihoods : $ADMIXDIR/K_likelihood.txt"
echo "  CLUMPPAK input    : $ADMIXDIR/clumppak_ready.txt"
echo "  Q-matrices        : $ADMIXDIR/ADMIX*.Q"
