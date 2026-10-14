
#!/usr/bin/env bash
#SBATCH --job-name=bakta_contigs
#SBATCH --error=JOB_OUT/bakta_%j.err
#SBATCH --output=JOB_OUT/bakta_%j.out
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

set -Eeuo pipefail

# Usage:
# sbatch workflow/04_functional_analysis/09_bakta_annotate_contigs.slurm.sh SAMPLE_ID ASSEMBLY_FASTA
#
# Required:
# BAKTA_ENV=/path/to/conda_BAKTA
# BAKTA_DB=/path/to/BAKTA_DB/db
#
# Optional:
# OUTDIR=/path/to/analysis_output

SAMPLE_ID="${1:?Usage: sbatch $0 SAMPLE_ID ASSEMBLY_FASTA}"
ASSEMBLY="${2:?Usage: sbatch $0 SAMPLE_ID ASSEMBLY_FASTA}"

BAKTA_ENV="${BAKTA_ENV:?Set BAKTA_ENV}"
BAKTA_DB="${BAKTA_DB:?Set BAKTA_DB}"
OUTDIR="${OUTDIR:-.}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

[[ -s "$ASSEMBLY" ]] || {
    echo "ERROR: Assembly FASTA not found or empty: $ASSEMBLY" >&2
    exit 1
}
ASSEMBLY="$(realpath "$ASSEMBLY")"

[[ -d "$BAKTA_DB" ]] || {
    echo "ERROR: BAKTA database directory not found: $BAKTA_DB" >&2
    exit 1
}

# Output directory is unique to the sample.
RESULT_DIR="${OUTDIR}/${SAMPLE_ID}/bakta_contigs"

# Use temporary storage without overwriting the system TMPDIR variable.
SCRATCH_ROOT="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
WORKDIR="$(mktemp -d "${SCRATCH_ROOT%/}/${SAMPLE_ID}_bakta_XXXXXX")"

cleanup() {
    rm -rf "$WORKDIR"
}
trap cleanup EXIT

mkdir -p "$WORKDIR/tmp" "$WORKDIR/output" "$RESULT_DIR"

module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$BAKTA_ENV"

command -v bakta >/dev/null 2>&1 || {
    echo "ERROR: bakta not found in the active environment." >&2
    exit 1
}

echo "Annotating assembly contigs for sample: $SAMPLE_ID"
echo "Assembly: $ASSEMBLY"

bakta \
    --db "$BAKTA_DB" \
    --meta \
    --skip-plot \
    --skip-trna \
    --skip-tmrna \
    --skip-rrna \
    --skip-ncrna \
    --skip-ncrna-region \
    --skip-crispr \
    --skip-pseudo \
    --skip-gap \
    --skip-ori \
    --threads "$THREADS" \
    --keep-contig-headers \
    --min-contig-length 200 \
    --output "$WORKDIR/output" \
    --tmp-dir "$WORKDIR/tmp" \
    "$ASSEMBLY"

# Check that annotation output exists.
[[ -s "$WORKDIR/output/${SAMPLE_ID}.ffn" ]] || {
    echo "ERROR: Expected nucleotide feature FASTA was not produced." >&2
    echo "Inspect BAKTA output and logs in the job log." >&2
    exit 1
}

# Preserve all generated outputs for downstream gene-centric analysis.
cp -a "$WORKDIR/output/." "$RESULT_DIR/"

echo "BAKTA completed for $SAMPLE_ID"
echo "Results saved to: $RESULT_DIR"

