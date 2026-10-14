#!/usr/bin/env bash
#SBATCH --job-name=checkm
#SBATCH --output=JOB_OUT/checkm_%j.out
#SBATCH --error=JOB_OUT/checkm_%j.err
#SBATCH --mem=50G
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

set -Eeuo pipefail
shopt -s nullglob

# Usage:
# sbatch workflow/03_quality_control/05_checkm_quality_assessment.slurm.sh SAMPLE_ID
#
# Required environment variable:
# METAWRAP_ENV

SAMPLE_ID="${1:?Usage: sbatch script SAMPLE_ID}"

OUTDIR="${OUTDIR:-.}"
SCRATCH_BASE="${TMPDIR:-/tmp}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

: "${METAWRAP_ENV:?Set METAWRAP_ENV to the Conda environment containing CheckM}"

SAMPLE_DIR="${OUTDIR}/${SAMPLE_ID}"
BIN_DIR="${SAMPLE_DIR}/bins_metawrap_refined"
RESULT_DIR="${SAMPLE_DIR}/bins_checkM"

[[ -d "$BIN_DIR" ]] || {
    echo "ERROR: Refined-bin directory not found: $BIN_DIR" >&2
    exit 1
}

BIN_FILES=("$BIN_DIR"/*.fa)

[[ "${#BIN_FILES[@]}" -gt 0 ]] || {
    echo "ERROR: No refined .fa bins found in $BIN_DIR" >&2
    exit 1
}

mkdir -p "$RESULT_DIR"

WORKDIR="$(mktemp -d "${SCRATCH_BASE%/}/${SAMPLE_ID}_checkm.XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT

mkdir -p "$WORKDIR/input_bins" "$WORKDIR/checkm_output"

# Stage only the refined genome FASTA files.
cp "${BIN_FILES[@]}" "$WORKDIR/input_bins/"

# Activate the environment.
module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$METAWRAP_ENV"

command -v checkm >/dev/null 2>&1 || {
    echo "ERROR: CheckM not found in the activated environment" >&2
    exit 1
}

echo "Starting CheckM lineage workflow for ${SAMPLE_ID}"
echo "Input bins: ${#BIN_FILES[@]}"
echo "Threads: ${THREADS}"

# Estimate genome completeness and contamination.
checkm lineage_wf \
    -t "$THREADS" \
    --pplacer_threads "$THREADS" \
    -x fa \
    --tab_table \
    -f "$WORKDIR/${SAMPLE_ID}_checkm_results.tsv" \
    "$WORKDIR/input_bins" \
    "$WORKDIR/checkm_output"

[[ -s "$WORKDIR/${SAMPLE_ID}_checkm_results.tsv" ]] || {
    echo "ERROR: CheckM summary table was not generated" >&2
    exit 1
}

# Save the summary table and full CheckM output for traceability.
cp "$WORKDIR/${SAMPLE_ID}_checkm_results.tsv" \
   "$RESULT_DIR/${SAMPLE_ID}_checkm_results.tsv"

cp -a "$WORKDIR/checkm_output" \
   "$RESULT_DIR/${SAMPLE_ID}_checkm_output"

echo "CheckM completed for ${SAMPLE_ID}"
echo "Summary: ${RESULT_DIR}/${SAMPLE_ID}_checkm_results.tsv"
echo "Detailed output: ${RESULT_DIR}/${SAMPLE_ID}_checkm_output"

