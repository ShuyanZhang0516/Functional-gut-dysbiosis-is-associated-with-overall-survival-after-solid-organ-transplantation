```bash
#!/usr/bin/env bash
#
# Gene catalogue functional annotation with KofamScan
#
# Steps:
#   1. Split the representative gene FASTA into chunks.
#   2. Translate nucleotide CDS sequences to proteins.
#   3. Submit KofamScan annotation jobs to SLURM.
#   4. Filter KofamScan results with R.
#   5. Merge filtered annotations into one table.
#
# Usage:
#   bash gene_catalogue_annotation.sh \
#     /path/to/cohort_genes_rep_seq.fasta \
#     /path/to/output_directory \
#     [sequences_per_chunk]
#
# Required software:
#   Bash, seqtk (not required by this version),
#   EMBOSS transeq, KofamScan, R, SLURM
#
# Configure these paths for your cluster before running.

set -Eeuo pipefail

# ============================================================
# 0. Configuration
# ============================================================

INPUT_FASTA="${1:?Usage: $0 input.fasta output_directory [chunk_size]}"
OUTDIR="${2:?Please provide an output directory}"
CHUNK_SIZE="${3:-50000}"

THREADS="${SLURM_CPUS_PER_TASK:-56}"
KOFAM_ENV="${KOFAM_ENV:-/home3/p306886/.conda/envs/kofamscan}"
KO_LIST="${KO_LIST:-/scratch/p306886/kofamscan/ko_list}"
KO_PROFILES="${KO_PROFILES:-/scratch/p306886/kofamscan/profiles}"

# This script is expected to be stored in the repository.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/simplifying_KofamScan_output.R"

# Output directories
SPLIT_DIR="${OUTDIR}/split"
PROTEIN_DIR="${OUTDIR}/proteins"
ANNOTATION_DIR="${OUTDIR}/kofamscan"
FILTERED_DIR="${OUTDIR}/filtered"
LOG_DIR="${OUTDIR}/logs"

mkdir -p "$SPLIT_DIR" "$PROTEIN_DIR" \
         "$ANNOTATION_DIR" "$FILTERED_DIR" "$LOG_DIR" JOB_OUT

if [[ ! -s "$INPUT_FASTA" ]]; then
    echo "ERROR: Input FASTA is missing or empty: $INPUT_FASTA" >&2
    exit 1
fi

if ! [[ "$CHUNK_SIZE" =~ ^[0-9]+$ ]] || (( CHUNK_SIZE < 1 )); then
    echo "ERROR: Chunk size must be a positive integer." >&2
    exit 1
fi

# ============================================================
# Step 1. Split the representative gene FASTA
# ============================================================
#
# Each output file contains at most CHUNK_SIZE sequences.
# Headers are preserved so gene IDs can be traced throughout
# the annotation workflow.
#
# A fresh output directory is recommended for every new run.

echo "STEP 1: Splitting gene catalogue"

if compgen -G "${SPLIT_DIR}/split*.fa" >/dev/null; then
    echo "ERROR: Existing split files found in ${SPLIT_DIR}."
    echo "Use a fresh output directory to avoid mixing old and new files."
    exit 1
fi

awk -v chunk_size="$CHUNK_SIZE" -v outdir="$SPLIT_DIR" '
    /^>/ {
        if (n % chunk_size == 0) {
            file = sprintf("%s/split%d.fa", outdir,
                           int(n / chunk_size) + 1)
        }
        print > file
        n++
        next
    }
    {
        if (n == 0) {
            print "ERROR: Invalid FASTA: sequence before first header." \
                > "/dev/stderr"
            exit 1
        }
        print > file
    }
    END {
        if (n == 0) {
            print "ERROR: No FASTA sequences found." > "/dev/stderr"
            exit 1
        }
        printf "Split %d sequences into %d chunks.\n",
               n, int((n + chunk_size - 1) / chunk_size) \
               > "/dev/stderr"
    }
' "$INPUT_FASTA"

# ============================================================
# Step 2. Convert nucleotide CDSs to protein sequences
# ============================================================
#
# IMPORTANT:
# Translate only if the nucleotide sequences are valid CDSs
# whose reading frame and orientation are known.
#
# If you have matching BAKTA-predicted proteins for the
# representative gene IDs, use those instead.
#
# The example below uses frame 1 and bacterial/archaeal
# translation table 11. Validate this for your data first.

echo "STEP 2: Translating gene sequences"

module load EMBOSS

shopt -s nullglob
SPLIT_FILES=("${SPLIT_DIR}"/split*.fa)

for fasta in "${SPLIT_FILES[@]}"; do
    base="$(basename "$fasta" .fa)"
    protein="${PROTEIN_DIR}/${base}.faa"

    if [[ -s "$protein" ]]; then
        echo "Keeping existing protein file: $protein"
        continue
    fi

    transeq \
        -sequence "$fasta" \
        -outseq "$protein" \
        -frame 1 \
        -table 11

    echo "Translated $base"
done

PROTEIN_FILES=("${PROTEIN_DIR}"/split*.faa)

if (( ${#PROTEIN_FILES[@]} == 0 )); then
    echo "ERROR: No protein FASTA files were generated." >&2
    exit 1
fi

# ============================================================
# Step 3. Submit KofamScan jobs
# ============================================================
#
# Each protein chunk is annotated independently.
# The SLURM worker is generated below so this workflow can
# be maintained as one main file.

echo "STEP 3: Submitting KofamScan jobs"

if [[ ! -s "$KO_LIST" || ! -d "$KO_PROFILES" ]]; then
    echo "ERROR: Check KO_LIST and KO_PROFILES database paths." >&2
    exit 1
fi

KOFAM_WORKER="${OUTDIR}/kofamscan_worker.slurm.sh"

cat > "$KOFAM_WORKER" <<'WORKER'
#!/usr/bin/env bash
#SBATCH --job-name=kofamscan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=23:00:00
#SBATCH --mem=18G

set -Eeuo pipefail

FAA="$1"
OUT="$2"
KOFAM_ENV="$3"
KO_LIST="$4"
KO_PROFILES="$5"

module load Anaconda3/2024.02-1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$KOFAM_ENV"

mkdir -p "$OUT"

name="$(basename "$FAA" .faa)"
tmpbase="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
tmpdir="$(mktemp -d "${tmpbase%/}/kofam_${name}_XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT

exec_annotation \
    -f detail-tsv \
    --cpu "${SLURM_CPUS_PER_TASK:-56}" \
    --tmp-dir "$tmpdir" \
    -o "${OUT}/${name}_annotated.tsv" \
    "$FAA" \
    --ko-list "$KO_LIST" \
    --profile "$KO_PROFILES"
WORKER

chmod +x "$KOFAM_WORKER"

KOFAM_JOB_IDS=()

for faa in "${PROTEIN_FILES[@]}"; do
    name="$(basename "$faa" .faa)"
    output="${ANNOTATION_DIR}/${name}_annotated.tsv"

    # Skip completed outputs only if they are nonempty.
    if [[ -s "$output" ]]; then
        echo "Annotation output already exists: $output"
        continue
    fi

    job_id="$(sbatch --parsable \
        --output="${LOG_DIR}/kofam_%j.out" \
        --error="${LOG_DIR}/kofam_%j.err" \
        "$KOFAM_WORKER" "$faa" "$ANNOTATION_DIR" \
        "$KOFAM_ENV" "$KO_LIST" "$KO_PROFILES")"

    KOFAM_JOB_IDS+=("$job_id")
    echo "Submitted $name: SLURM job $job_id"
done

# Wait for submitted jobs to finish.
# Requires the SLURM accounting command sacct to be available.
for job_id in "${KOFAM_JOB_IDS[@]}"; do
    while true; do
        state="$(sacct -n -X -j "$job_id" --format=State 2>/dev/null \
                 | awk 'NF {print $1; exit}')"

        case "$state" in
            COMPLETED)
                break
                ;;
            FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)
                echo "ERROR: KofamScan job $job_id ended with $state" >&2
                exit 1
                ;;
            *)
                sleep 30
                ;;
        esac
    done
done

# Verify that each protein chunk has an annotation result.
for faa in "${PROTEIN_FILES[@]}"; do
    name="$(basename "$faa" .faa)"
    if [[ ! -s "${ANNOTATION_DIR}/${name}_annotated.tsv" ]]; then
        echo "ERROR: Missing KofamScan output for $name" >&2
        exit 1
    fi
done

# ============================================================
# Step 4. Filter KofamScan results with R
# ============================================================
#
# The R script selects a best KO hit per gene after applying
# the configured E-value and score filters.
#
# Verify its expected input columns against your installed
# KofamScan version before processing all chunks.

echo "STEP 4: Filtering KofamScan results"

FILTER_WORKER="${OUTDIR}/filter_worker.slurm.sh"

cat > "$FILTER_WORKER" <<'WORKER'
#!/usr/bin/env bash
#SBATCH --job-name=kofam_filter
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=4G

set -Eeuo pipefail

INPUT="$1"
OUTPUT="$2"
R_SCRIPT="$3"

module load R

Rscript "$R_SCRIPT" "$INPUT" "$OUTPUT"
WORKER

FILTER_JOB_IDS=()

for input in "${ANNOTATION_DIR}"/split*_annotated.tsv; do
    [[ -e "$input" ]] || continue

    name="$(basename "$input" .tsv)"
    output="${FILTERED_DIR}/${name}_filtered.tsv"

    if [[ -s "$output" ]]; then
        echo "Filtered output already exists: $output"
        continue
    fi

    job_id="$(sbatch --parsable \
        --output="${LOG_DIR}/filter_%j.out" \
        --error="${LOG_DIR}/filter_%j.err" \
        "$FILTER_WORKER" "$input" "$output" "$R_SCRIPT")"

    FILTER_JOB_IDS+=("$job_id")
    echo "Submitted filtering job $job_id for $name"
done

for job_id in "${FILTER_JOB_IDS[@]}"; do
    while true; do
        state="$(sacct -n -X -j "$job_id" --format=State 2>/dev/null \
                 | awk 'NF {print $1; exit}')"

        case "$state" in
            COMPLETED)
                break
                ;;
            FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)
                echo "ERROR: Filtering job $job_id ended with $state" >&2
                exit 1
                ;;
            *)
                sleep 30
                ;;
        esac
    done
done

# ============================================================
# Step 5. Merge filtered annotation tables
# ============================================================
#
# Keep one header and append data rows from every filtered TSV.

echo "STEP 5: Merging filtered annotations"

FILTERED_FILES=("${FILTERED_DIR}"/split*_annotated_filtered.tsv)

if (( ${#FILTERED_FILES[@]} == 0 )); then
    echo "ERROR: No filtered annotation files found." >&2
    exit 1
fi

MERGED="${OUTDIR}/merged_annotation_filtered.tsv"

head -n 1 "${FILTERED_FILES[0]}" > "$MERGED"

for file in "${FILTERED_FILES[@]}"; do
    tail -n +2 "$file" >> "$MERGED"
done

echo "=========================================="
echo "Gene catalogue KO annotation completed."
echo "Merged annotation: $MERGED"
echo "=========================================="
```
