#!/usr/bin/env bash
#SBATCH --job-name=gtdbtk_classify
#SBATCH --error=JOB_OUT/gtdbtk_%j.err
#SBATCH --output=JOB_OUT/gtdbtk_%j.out
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

set -Eeuo pipefail

# Usage:
# sbatch workflow/03_quality_control/07_gtdbtk_classify_bins.slurm.sh SAMPLE_ID
#
# Required environment variables:
# GTDBTK_ENV=/path/to/conda_GTDB_TK
# GTDBTK_DATA_PATH=/path/to/GTDB-Tk/reference_database
#
# Optional:
# OUTDIR=/path/to/analysis_output

SAMPLE_ID="${1:?Usage: sbatch $0 SAMPLE_ID}"
OUTDIR="${OUTDIR:-.}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

GTDBTK_ENV="${GTDBTK_ENV:?Set GTDBTK_ENV to your GTDB-Tk conda environment}"
export GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:?Set GTDBTK_DATA_PATH to your GTDB-Tk database}"

SAMPLE_DIR="${OUTDIR}/${SAMPLE_ID}"
BIN_DIR="${SAMPLE_DIR}/bins_metawrap_refined"
FINAL_DIR="${SAMPLE_DIR}/bins_GTDBK"

# Check input bins.
[[ -d "$BIN_DIR" ]] || {
    echo "ERROR: Refined bin directory not found: $BIN_DIR" >&2
    exit 1
}

shopt -s nullglob
BINS=("$BIN_DIR"/*.fa)
shopt -u nullglob

(( ${#BINS[@]} > 0 )) || {
    echo "ERROR: No .fa MAG files found in $BIN_DIR" >&2
    exit 1
}

# Load software environment.
module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$GTDBTK_ENV"

command -v gtdbtk >/dev/null 2>&1 || {
    echo "ERROR: gtdbtk not found in the active environment." >&2
    exit 1
}

# Confirm the configured database path exists.
[[ -d "$GTDBTK_DATA_PATH" ]] || {
    echo "ERROR: GTDB-Tk database directory not found: $GTDBTK_DATA_PATH" >&2
    exit 1
}

# Create a separate scratch workspace.
SCRATCH_ROOT="${TMPDIR:-/tmp}"
WORKDIR="$(mktemp -d "${SCRATCH_ROOT%/}/${SAMPLE_ID}_gtdbtk_XXXXXX")"

cleanup() {
    rm -rf "$WORKDIR"
}
trap cleanup EXIT

mkdir -p "$WORKDIR/genomes" "$WORKDIR/gtdbtk_out" "$FINAL_DIR"

# Stage MAGs separately from GTDB-Tk output.
cp "${BINS[@]}" "$WORKDIR/genomes/"

echo "=========================================="
echo "GTDB-Tk classification"
echo "Sample: $SAMPLE_ID"
echo "Number of MAGs: ${#BINS[@]}"
echo "=========================================="

gtdbtk classify_wf \
    --skip_ani_screen \
    --genome_dir "$WORKDIR/genomes" \
    --out_dir "$WORKDIR/gtdbtk_out" \
    --extension fa \
    --cpus "$THREADS"

# Confirm that GTDB-Tk produced at least one domain summary.
SUMMARY_FILES=(
    "$WORKDIR/gtdbtk_out"/gtdbtk.bac120.summary.tsv
    "$WORKDIR/gtdbtk_out"/gtdbtk.ar53.summary.tsv
)

SUMMARY_FOUND=0
for FILE in "${SUMMARY_FILES[@]}"; do
    if [[ -s "$FILE" ]]; then
        SUMMARY_FOUND=1
    fi
done

if [[ "$SUMMARY_FOUND" -ne 1 ]]; then
    echo "ERROR: No non-empty GTDB-Tk domain summary table was produced." >&2
    echo "Inspect GTDB-Tk logs and warnings in the scratch output." >&2
    exit 1
fi

# Preserve all GTDB-Tk results, logs, summaries, and intermediate outputs.
cp -a "$WORKDIR/gtdbtk_out/." "$FINAL_DIR/"

echo "GTDB-Tk classification complete."
echo "Results saved in: $FINAL_DIR"

