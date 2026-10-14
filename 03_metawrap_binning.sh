#!/usr/bin/env bash
#SBATCH --job-name=metawrap_binning
#SBATCH --output=JOB_OUT/metawrap_binning_%j.out
#SBATCH --error=JOB_OUT/metawrap_binning_%j.err
#SBATCH --mem=16G
#SBATCH --time=06:30:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

set -Eeuo pipefail

# Usage:
# sbatch workflow/02_binning/03_metawrap_binning.sh SAMPLE_ID
#
# Required environment variable:
# MET AWRAP_ENV (see configuration below; variable name is MET AWRAP_ENV
# without the space: MET AWRAP_ENV is documented here for readability).

SAMPLE_ID="${1:?Missing sample ID}"

OUTDIR="${OUTDIR:-.}"
SCRATCH_BASE="${TMPDIR:-/tmp}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
MEM_GB="${METAWRAP_MEM_GB:-16}"

: "${METAWRAP_ENV:?Set the MetaWRAP Conda environment path}"

SAMPLE_DIR="${OUTDIR}/${SAMPLE_ID}"
READ1="${SAMPLE_DIR}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq"
READ2="${SAMPLE_DIR}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq"
ASSEMBLY="${SAMPLE_DIR}/assembly_megahit/${SAMPLE_ID}_megahit_contigs.fa"
FINAL_DIR="${SAMPLE_DIR}/bins_metawrap"

[[ -s "$READ1" ]] || { echo "ERROR: R1 not found: $READ1" >&2; exit 1; }
[[ -s "$READ2" ]] || { echo "ERROR: R2 not found: $READ2" >&2; exit 1; }
[[ -s "$ASSEMBLY" ]] || { echo "ERROR: Assembly not found: $ASSEMBLY" >&2; exit 1; }

# Avoid accidentally overwriting a previous run.
if [[ -e "$FINAL_DIR" ]]; then
    echo "ERROR: Output directory already exists: $FINAL_DIR" >&2
    echo "Move it aside or remove it deliberately before rerunning." >&2
    exit 1
fi

mkdir -p "$(dirname "$FINAL_DIR")"

WORKDIR="$(mktemp -d "${SCRATCH_BASE%/}/${SAMPLE_ID}_metawrap.XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT

module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$MET AWRAP_ENV"

command -v repair.sh >/dev/null 2>&1 || {
    echo "ERROR: repair.sh is not available in the activated environment" >&2
    exit 1
}
command -v metawrap >/dev/null 2>&1 || {
    echo "ERROR: metawrap is not available in the activated environment" >&2
    exit 1
}

# Repair paired reads for the binning tools.
mkdir -p "$WORKDIR/reads_fixed"

repair.sh \
    in="$READ1" \
    in2="$READ2" \
    out="$WORKDIR/reads_fixed/${SAMPLE_ID}_rdy_1.fastq" \
    out2="$WORKDIR/reads_fixed/${SAMPLE_ID}_rdy_2.fastq" \
    outs="$WORKDIR/reads_fixed/${SAMPLE_ID}_singletons.fastq"

FIXED_R1="$WORKDIR/reads_fixed/${SAMPLE_ID}_rdy_1.fastq"
FIXED_R2="$WORKDIR/reads_fixed/${SAMPLE_ID}_rdy_2.fastq"

[[ -s "$FIXED_R1" && -s "$FIXED_R2" ]] || {
    echo "ERROR: BBTools repair did not produce both paired-read files" >&2
    exit 1
}

# Run MetaWRAP binning with MetaBAT2, MaxBin2 and CONCOCT.
metawrap binning \
    -o "$WORKDIR/bins_metawrap" \
    -t "$THREADS" \
    -m "$MEM_GB" \
    -a "$ASSEMBLY" \
    --metabat2 \
    --maxbin2 \
    --concoct \
    -l 1000 \
    "$FIXED_R1" \
    "$FIXED_R2"

# Verify expected binner output directories.
for BIN_DIR in metabat2_bins maxbin2_bins concoct_bins; do
    if [[ ! -d "$WORKDIR/bins_metawrap/$BIN_DIR" ]]; then
        echo "ERROR: Expected MetaWRAP output missing: $BIN_DIR" >&2
        exit 1
    fi
done

# Preserve the complete binner outputs, including unbinned contigs
# and tool-generated files, for traceability.
mkdir -p "$FINAL_DIR"
cp -a "$WORKDIR/bins_metawrap/." "$FINAL_DIR/"

echo "MetaWRAP binning completed for ${SAMPLE_ID}"
echo "Binner outputs: ${FINAL_DIR}"

