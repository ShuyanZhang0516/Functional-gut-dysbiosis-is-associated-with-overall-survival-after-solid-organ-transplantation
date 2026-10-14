#!/usr/bin/env bash
#SBATCH --job-name=gene_coverage
#SBATCH --output=gene_coverage_%j.out
#SBATCH --error=gene_coverage_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G

set -Eeuo pipefail

# Usage:
# sbatch 01_gene_coverage.slurm.sh \
#   SAMPLE_ID READ1 READ2 INDEX_PREFIX GENE_CATALOGUE OUTPUT_DIR

SAMPLE_ID="${1:?Provide sample ID}"
READ1="${2:?Provide forward cleaned reads}"
READ2="${3:?Provide reverse cleaned reads}"
INDEX_PREFIX="${4:?Provide Bowtie2 index prefix}"
GENE_CATALOGUE="${5:?Provide representative gene catalogue FASTA}"
OUTPUT_DIR="${6:?Provide output directory}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Load modules appropriate for your HPC system.
module load Bowtie2
module load SAMtools
module load CoverM

mkdir -p "$OUTPUT_DIR"

for FILE in "$READ1" "$READ2" "$GENE_CATALOGUE"; do
    [[ -s "$FILE" ]] || {
        echo "ERROR: File missing or empty: $FILE" >&2
        exit 1
    }
done

[[ -s "${INDEX_PREFIX}.1.bt2" ||
   -s "${INDEX_PREFIX}.1.bt2l" ]] || {
    echo "ERROR: Bowtie2 index not found: $INDEX_PREFIX" >&2
    exit 1
}

BAM="${OUTPUT_DIR}/${SAMPLE_ID}.genes.sorted.bam"
COVERAGE="${OUTPUT_DIR}/${SAMPLE_ID}.gene_coverage.tsv"
LOG="${OUTPUT_DIR}/${SAMPLE_ID}.bowtie2.log"

echo "Sample: $SAMPLE_ID"
echo "Mapping reads to representative gene catalogue..."

# Bowtie2 alignment, piped directly into coordinate-sorted BAM.
bowtie2 \
    --very-sensitive \
    --threads "$THREADS" \
    -x "$INDEX_PREFIX" \
    -1 "$READ1" \
    -2 "$READ2" \
    2> "$LOG" \
| samtools sort \
    -@ "$THREADS" \
    -o "$BAM" -

samtools quickcheck "$BAM"

echo "Calculating per-gene coverage with CoverM..."

coverm contig \
    --bam-files "$BAM" \
    --methods mean covered_bases count length \
    --contig-end-exclusion 0 \
    --threads "$THREADS" \
    --output-file "$COVERAGE"

echo "Finished sample: $SAMPLE_ID"
echo "BAM: $BAM"
echo "Coverage table: $COVERAGE"

