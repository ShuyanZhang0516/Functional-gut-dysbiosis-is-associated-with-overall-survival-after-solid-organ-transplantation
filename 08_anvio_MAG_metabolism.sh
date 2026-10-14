#!/usr/bin/env bash
#SBATCH --job-name=anvio_metabolism
#SBATCH --error=JOB_OUT/anvio_%j.err
#SBATCH --output=JOB_OUT/anvio_%j.out
#SBATCH --mem=100G
#SBATCH --time=25:00:00
#SBATCH --cpus-per-task=16
#SBATCH --open-mode=truncate

set -Eeuo pipefail

# Usage:
# sbatch workflow/04_functional_analysis/08_anvio_MAG_metabolism.sh SAMPLE_ID /path/to/MAG.fa
#
# Required environment variables:
# ANVIO_ENV=/path/to/anvio8
# KEGG_DATA_DIR=/path/to/kegg_data_anvio
#
# Optional:
# OUTDIR=/path/to/analysis_output

SAMPLE_ID="${1:?Usage: sbatch $0 SAMPLE_ID MAG_FASTA}"
MAG_FASTA="${2:?Usage: sbatch $0 SAMPLE_ID MAG_FASTA}"

ANVIO_ENV="${ANVIO_ENV:?Set ANVIO_ENV to your anvi'o conda environment}"
KEGG_DATA_DIR="${KEGG_DATA_DIR:?Set KEGG_DATA_DIR to your anvi'o KEGG data directory}"
OUTDIR="${OUTDIR:-.}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Resolve the input FASTA to an absolute path before changing directories.
[[ -f "$MAG_FASTA" ]] || {
    echo "ERROR: MAG FASTA not found: $MAG_FASTA" >&2
    exit 1
}
MAG_FASTA="$(realpath "$MAG_FASTA")"

[[ -d "$KEGG_DATA_DIR" ]] || {
    echo "ERROR: KEGG data directory not found: $KEGG_DATA_DIR" >&2
    exit 1
}

MAG_FILE="$(basename "$MAG_FASTA")"
MAG_ID="${MAG_FILE%.fa}"
MAG_ID="${MAG_ID%.fasta}"
MAG_ID="${MAG_ID%.fna}"

RESULT_DIR="${OUTDIR}/${SAMPLE_ID}/${MAG_ID}"

# Load anvi'o.
module purge
module load Anaconda3/2022.05
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ANVIO_ENV"

for CMD in anvi-script-reformat-fasta \
           anvi-gen-contigs-database \
           anvi-run-kegg-kofams \
           anvi-estimate-metabolism; do
    command -v "$CMD" >/dev/null 2>&1 || {
        echo "ERROR: Required command not found: $CMD" >&2
        exit 1
    }
done

mkdir -p "$RESULT_DIR"
cd "$RESULT_DIR"

echo "=========================================="
echo "anvi'o MAG metabolism analysis"
echo "Sample: $SAMPLE_ID"
echo "MAG: $MAG_ID"
echo "Input: $MAG_FASTA"
echo "=========================================="

# 1. Reformat FASTA headers for anvi'o.
REFORMATTED="${MAG_ID}.reformat.fa"

anvi-script-reformat-fasta \
    "$MAG_FASTA" \
    --simplify-names \
    --output-file "$REFORMATTED"

# 2. Create the anvi'o contigs database.
CONTIGS_DB="${MAG_ID}-CONTIGS.db"

anvi-gen-contigs-database \
    --contigs-fasta "$REFORMATTED" \
    --project-name "$MAG_ID" \
    --output-db-path "$CONTIGS_DB" \
    --num-threads "$THREADS"

# 3. Annotate genes using KEGG KOfam profiles.
anvi-run-kegg-kofams \
    -c "$CONTIGS_DB" \
    -T "$THREADS" \
    --kegg-data-dir "$KEGG_DATA_DIR"

# 4. Estimate metabolic modules and export annotation hits.
anvi-estimate-metabolism \
    -c "$CONTIGS_DB" \
    -O "$MAG_ID" \
    --kegg-data-dir "$KEGG_DATA_DIR" \
    --output-modes hits,modules

echo "=========================================="
echo "Completed: $SAMPLE_ID / $MAG_ID"
echo "Results: $RESULT_DIR"
echo "=========================================="

