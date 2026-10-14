#!/usr/bin/env bash
#SBATCH --job-name=megahit
#SBATCH --output=JOB_OUT/megahit_%j.out
#SBATCH --error=JOB_OUT/megahit_%j.err
#SBATCH --mem=40G
#SBATCH --time=09:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

set -Eeuo pipefail

# Usage:
# sbatch workflow/01_preprocessing/02_megahit_assembly.sh SAMPLE_ID
#
# Required environment variable:
# MEGAHIT_ENV

SAMPLE_ID="${1:?Missing sample ID}"

OUTDIR="${OUTDIR:-.}"
SCRATCH_BASE="${TMPDIR:-/tmp}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
MEGAHIT_MEM_FRACTION="${MEGAHIT_MEM_FRACTION:-0.85}"

: "${MEGAHIT_ENV:?Set MEGAHIT_ENV to the MEGAHIT Conda environment}"

CLEAN_DIR="${OUTDIR}/${SAMPLE_ID}/clean_reads"
ASSEMBLY_DIR="${OUTDIR}/${SAMPLE_ID}/assembly_megahit"

READ1="${CLEAN_DIR}/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq"
READ2="${CLEAN_DIR}/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq"

[[ -s "$READ1" ]] || { echo "ERROR: R1 not found: $READ1" >&2; exit 1; }
[[ -s "$READ2" ]] || { echo "ERROR: R2 not found: $READ2" >&2; exit 1; }
[[ -d "$MEGAHIT_ENV" ]] || { echo "ERROR: MEGAHIT environment not found" >&2; exit 1; }

mkdir -p "$ASSEMBLY_DIR"

WORKDIR="$(mktemp -d "${SCRATCH_BASE%/}/${SAMPLE_ID}_megahit.XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT

module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$MEGAHIT_ENV"

command -v megahit >/dev/null 2>&1 || {
    echo "ERROR: megahit not found in MEGAHIT_ENV" >&2
    exit 1
}

echo "Starting MEGAHIT for ${SAMPLE_ID}"
megahit \
    -1 "$READ1" \
    -2 "$READ2" \
    -o "$WORKDIR/assembly" \
    -m "$MEGAHIT_MEM_FRACTION" \
    -t "$THREADS" \
    --tmp-dir "$WORKDIR/tmp"

CONTIGS="$WORKDIR/assembly/final.contigs.fa"
ASSEMBLY_LOG="$WORKDIR/assembly/log"

[[ -s "$CONTIGS" ]] || {
    echo "ERROR: MEGAHIT final contigs were not produced" >&2
    exit 1
}

cp "$CONTIGS" \
    "$ASSEMBLY_DIR/${SAMPLE_ID}_megahit_contigs.fa"

if [[ -f "$ASSEMBLY_LOG" ]]; then
    cp "$ASSEMBLY_LOG" \
        "$ASSEMBLY_DIR/${SAMPLE_ID}_megahit_log.log"
fi

echo "MEGAHIT assembly completed for ${SAMPLE_ID}"
echo "Contigs: ${ASSEMBLY_DIR}/${SAMPLE_ID}_megahit_contigs.fa"

