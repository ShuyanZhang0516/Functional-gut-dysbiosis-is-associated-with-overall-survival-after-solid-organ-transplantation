#!/usr/bin/env bash
#SBATCH --job-name=metawrap_refinement
#SBATCH --output=JOB_OUT/metawrap_refinement_%j.out
#SBATCH --error=JOB_OUT/metawrap_refinement_%j.err
#SBATCH --mem=48G
#SBATCH --time=06:30:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

set -Eeuo pipefail
shopt -s nullglob

# Usage:
# sbatch workflow/02_binning/04_metawrap_bin_refinement.sh SAMPLE_ID
#
# Required environment variable:
# METAWRAP_ENV

SAMPLE_ID="${1:?Usage: sbatch script SAMPLE_ID}"

OUTDIR="${OUTDIR:-.}"
SCRATCH_BASE="${TMPDIR:-/tmp}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
MEM_GB="${METAWRAP_MEM_GB:-48}"

: "${METAWRAP_ENV:?Set METAWRAP_ENV to the MetaWRAP Conda environment path}"

SAMPLE_DIR="${OUTDIR}/${SAMPLE_ID}"
BIN_ROOT="${SAMPLE_DIR}/bins_metawrap"
FINAL_DIR="${SAMPLE_DIR}/bins_metawrap_refined"

INPUT_METABAT="${BIN_ROOT}/metabat2_bins"
INPUT_MAXBIN="${BIN_ROOT}/maxbin2_bins"
INPUT_CONCOCT="${BIN_ROOT}/concoct_bins"

# Avoid overwriting a previous result.
if [[ -e "$FINAL_DIR" ]]; then
    echo "ERROR: Output directory already exists: $FINAL_DIR" >&2
    echo "Move it aside or remove it deliberately before rerunning." >&2
    exit 1
fi

for DIR in "$INPUT_METABAT" "$INPUT_MAXBIN" "$INPUT_CONCOCT"; do
    [[ -d "$DIR" ]] || {
        echo "ERROR: Missing binner output directory: $DIR" >&2
        exit 1
    }

    FILES=("$DIR"/*.fa)
    if [[ "${#FILES[@]}" -eq 0 ]]; then
        echo "ERROR: No .fa bin files found in $DIR" >&2
        exit 1
    fi
done

WORKDIR="$(mktemp -d "${SCRATCH_BASE%/}/${SAMPLE_ID}_refinement.XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT

mkdir -p \
    "$WORKDIR/bins_maxbin2" \
    "$WORKDIR/bins_metabat2" \
    "$WORKDIR/bins_concoct"

# Stage bin FASTA files on scratch.
cp "$INPUT_MAXBIN"/*.fa "$WORKDIR/bins_maxbin2/"
cp "$INPUT_METABAT"/*.fa "$WORKDIR/bins_metabat2/"
cp "$INPUT_CONCOCT"/*.fa "$WORKDIR/bins_concoct/"

# Load MetaWRAP.
module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$METAWRAP_ENV"

command -v metawrap >/dev/null 2>&1 || {
    echo "ERROR: metawrap not found in the activated environment" >&2
    exit 1
}

echo "Starting MetaWRAP bin refinement for ${SAMPLE_ID}"
echo "Minimum completeness: 50%"
echo "Maximum contamination: 10%"
echo "Quick mode enabled"

metawrap bin_refinement \
    -t "$THREADS" \
    -m "$MEM_GB" \
    --quick \
    -c 50 \
    -x 10 \
    -o "$WORKDIR/bins_metawrap_refined" \
    -A "$WORKDIR/bins_maxbin2" \
    -B "$WORKDIR/bins_metabat2" \
    -C "$WORKDIR/bins_concoct"

REFINED_DIR="$WORKDIR/bins_metawrap_refined/metawrap_50_10_bins"

[[ -d "$REFINED_DIR" ]] || {
    echo "ERROR: Expected refined-bin directory not found: $REFINED_DIR" >&2
    exit 1
}

REFINED_BINS=("$REFINED_DIR"/*.fa)

if [[ "${#REFINED_BINS[@]}" -eq 0 ]]; then
    echo "ERROR: No refined FASTA bins were produced" >&2
    exit 1
fi

# Copy and rename bins only after successful refinement.
mkdir -p "$FINAL_DIR"

for BIN in "${REFINED_BINS[@]}"; do
    BIN_NAME="$(basename "$BIN")"
    RENAMED="${BIN_NAME/bin./${SAMPLE_ID}.bin.}"
    cp "$BIN" "$FINAL_DIR/$RENAMED"
done

# Preserve other refinement reports when available.
for REPORT in "$WORKDIR/bins_metawrap_refined"/*.stats \
              "$WORKDIR/bins_metawrap_refined"/*.txt \
              "$WORKDIR/bins_metawrap_refined"/*.tsv; do
    [[ -f "$REPORT" ]] && cp "$REPORT" "$FINAL_DIR/"
done

echo "Refinement completed for ${SAMPLE_ID}"
echo "Number of refined bins: ${#REFINED_BINS[@]}"
echo "Refined bins saved in: ${FINAL_DIR}"

# Downstream CheckM, GTDB-Tk and quantification jobs should be
# submitted only after verifying their scripts and input paths.

