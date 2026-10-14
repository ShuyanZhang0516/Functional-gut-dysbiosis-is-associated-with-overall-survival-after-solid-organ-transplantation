#!/usr/bin/env bash
#SBATCH --job-name=bin_quantification
#SBATCH --error=JOB_OUT/bin_quantification_%j.err
#SBATCH --output=JOB_OUT/bin_quantification_%j.out
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

set -Eeuo pipefail

# Usage:
# sbatch workflow/03_quality_control/06_metawrap_quantify_bins.slurm.sh SAMPLE_ID
#
# Required environment variable:
# METAWRAP_ENV=/path/to/conda_metawrap
#
# Optional:
# OUTDIR=/path/to/analysis_output

SAMPLE_ID="${1:?Usage: sbatch $0 SAMPLE_ID}"
OUTDIR="${OUTDIR:-.}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
METAWRAP_ENV="${METAWRAP_ENV:?Set METAWRAP_ENV to your MetaWRAP conda environment}"

SAMPLE_DIR="${OUTDIR}/${SAMPLE_ID}"
BIN_DIR="${SAMPLE_DIR}/bins_metawrap_refined"
READ1="${SAMPLE_DIR}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq"
READ2="${SAMPLE_DIR}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq"
FINAL_DIR="${SAMPLE_DIR}/bins_quantification"

# Check required inputs.
[[ -d "$BIN_DIR" ]] || {
    echo "ERROR: Refined bin directory not found: $BIN_DIR" >&2
    exit 1
}
[[ -f "$READ1" ]] || {
    echo "ERROR: Read 1 not found: $READ1" >&2
    exit 1
}
[[ -f "$READ2" ]] || {
    echo "ERROR: Read 2 not found: $READ2" >&2
    exit 1
}

shopt -s nullglob
BINS=("$BIN_DIR"/*.fa)
shopt -u nullglob

(( ${#BINS[@]} > 0 )) || {
    echo "ERROR: No .fa bin files found in $BIN_DIR" >&2
    exit 1
}

# Load the software environment.
module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$METAWRAP_ENV"

command -v metawrap >/dev/null 2>&1 || {
    echo "ERROR: metawrap not found in the active environment." >&2
    exit 1
}

# Create a temporary working directory.
SCRATCH_ROOT="${TMPDIR:-/tmp}"
WORKDIR="$(mktemp -d "${SCRATCH_ROOT%/}/${SAMPLE_ID}_quant_bins_XXXXXX")"

cleanup() {
    rm -rf "$WORKDIR"
}
trap cleanup EXIT

mkdir -p "$WORKDIR/bins" "$WORKDIR/quantification" "$FINAL_DIR"

# Stage refined bins on scratch.
cp "${BINS[@]}" "$WORKDIR/bins/"

echo "=========================================="
echo "Quantifying bins for sample: $SAMPLE_ID"
echo "=========================================="

metawrap quant_bins \
    -t "$THREADS" \
    -b "$WORKDIR/bins" \
    -o "$WORKDIR/quantification" \
    "$READ1" \
    "$READ2"

ABUNDANCE_TABLE="$WORKDIR/quantification/bin_abundance_table.tab"

[[ -s "$ABUNDANCE_TABLE" ]] || {
    echo "ERROR: MetaWRAP abundance table missing or empty: $ABUNDANCE_TABLE" >&2
    echo "Check the MetaWRAP logs and input read/bin compatibility." >&2
    exit 1
}

cp "$ABUNDANCE_TABLE" \
   "$FINAL_DIR/${SAMPLE_ID}_bin_abundances.tab"

echo "Quantification complete."
echo "Abundance table: $FINAL_DIR/${SAMPLE_ID}_bin_abundances.tab"

