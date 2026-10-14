#!/usr/bin/env bash
#SBATCH --job-name=preprocess
#SBATCH --output=JOB_OUT/preprocess_%j.out
#SBATCH --error=JOB_OUT/preprocess_%j.err
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

set -Eeuo pipefail

# Usage:
# sbatch workflow/01_preprocessing/01_read_preprocessing.sh \
#   SAMPLE_ID READ1.fastq.gz READ2.fastq.gz
#
# Required environment variables:
# ADAPTERS, HUMAN_DB, BBDUK_ENV, KNEADDATA_ENV

SAMPLE_ID="${1:?Missing sample ID}"
READ1="${2:?Missing R1 FASTQ.gz}"
READ2="${3:?Missing R2 FASTQ.gz}"

OUTDIR="${OUTDIR:-.}"
SCRATCH_BASE="${TMPDIR:-/tmp}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

: "${ADAPTERS:?Set ADAPTERS to the adapter FASTA file}"
: "${HUMAN_DB:?Set HUMAN_DB to the KneadData human reference database}"
: "${BBDUK_ENV:?Set BBDUK_ENV to the BBDuk Conda environment}"
: "${KNEADDATA_ENV:?Set KNEADDATA_ENV to the KneadData Conda environment}"

SAMPLE_OUT="${OUTDIR}/${SAMPLE_ID}"
CLEAN_DIR="${SAMPLE_OUT}/clean_reads"
QC_PRE="${SAMPLE_OUT}/qc_preclean"
QC_POST="${SAMPLE_OUT}/qc_postclean"

[[ -f "$READ1" ]] || { echo "ERROR: R1 not found: $READ1" >&2; exit 1; }
[[ -f "$READ2" ]] || { echo "ERROR: R2 not found: $READ2" >&2; exit 1; }
[[ -f "$ADAPTERS" ]] || { echo "ERROR: Adapter file not found: $ADAPTERS" >&2; exit 1; }
[[ -d "$HUMAN_DB" ]] || { echo "ERROR: Human database not found: $HUMAN_DB" >&2; exit 1; }

mkdir -p "$CLEAN_DIR" "$QC_PRE" "$QC_POST"

WORKDIR="$(mktemp -d "${SCRATCH_BASE%/}/${SAMPLE_ID}_preprocess.XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT

module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"

# Stage input reads on scratch and decompress them.
cp "$READ1" "$WORKDIR/R1.fastq.gz"
cp "$READ2" "$WORKDIR/R2.fastq.gz"
gunzip "$WORKDIR/R1.fastq.gz" "$WORKDIR/R2.fastq.gz"

R1="$WORKDIR/R1.fastq"
R2="$WORKDIR/R2.fastq"

# Quality control before cleaning.
conda activate "$KNEADDATA_ENV"
fastqc -t "$THREADS" -o "$QC_PRE" "$R1" "$R2"

# Adapter trimming with BBDuk.
conda activate "$BBDUK_ENV"
command -v bbduk.sh >/dev/null 2>&1 || {
    echo "ERROR: bbduk.sh not found in BBDUK_ENV" >&2
    exit 1
}

bbduk.sh \
    in1="$R1" \
    in2="$R2" \
    out1="$WORKDIR/R1_trimmed.fastq" \
    out2="$WORKDIR/R2_trimmed.fastq" \
    ref="$ADAPTERS" \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    threads="$THREADS" \
    2>&1 | tee "$SAMPLE_OUT/bbduk.log"

# Host-read removal and quality trimming with KneadData.
conda activate "$KNEADDATA_ENV"
command -v kneaddata >/dev/null 2>&1 || {
    echo "ERROR: kneaddata not found in KNEADDATA_ENV" >&2
    exit 1
}

mkdir -p "$WORKDIR/kneaddata_output"

kneaddata \
    --input "$WORKDIR/R1_trimmed.fastq" \
    --input "$WORKDIR/R2_trimmed.fastq" \
    --output "$WORKDIR/kneaddata_output" \
    --output-prefix "${SAMPLE_ID}" \
    --reference-db "$HUMAN_DB" \
    --threads "$THREADS" \
    --trimmomatic-options "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" \
    --bypass-trf \
    --log "$SAMPLE_OUT/kneaddata.log"

# Locate the paired, host-filtered output files.
mapfile -t R1_FILES < <(
    find "$WORKDIR/kneaddata_output" -type f -name '*_paired_1.fastq'
)

if [[ "${#R1_FILES[@]}" -ne 1 ]]; then
    echo "ERROR: Expected one KneadData paired R1 output; found ${#R1_FILES[@]}" >&2
    find "$WORKDIR/kneaddata_output" -maxdepth 1 -type f -print >&2
    exit 1
fi

KNEAD_R1="${R1_FILES[0]}"
KNEAD_R2="${KNEAD_R1%_1.fastq}_2.fastq"

[[ -s "$KNEAD_R1" && -s "$KNEAD_R2" ]] || {
    echo "ERROR: KneadData paired output files are missing or empty" >&2
    exit 1
}

# Standardized filenames consumed by the assembly script.
cp "$KNEAD_R1" "$CLEAN_DIR/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq"
cp "$KNEAD_R2" "$CLEAN_DIR/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq"

# Quality control after cleaning.
fastqc -t "$THREADS" -o "$QC_POST" \
    "$CLEAN_DIR/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq" \
    "$CLEAN_DIR/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq"

echo "Read preprocessing completed for ${SAMPLE_ID}"
echo "Clean reads: ${CLEAN_DIR}"

