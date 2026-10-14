#!/bin/bash

#SBATCH --job-name=MMSeqs2_clustering
#SBATCH --output=MMSeqs2_clustering_%j.log
#SBATCH --error=MMSeqs2_clustering_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=23:00:00
#SBATCH --mem=180G

# Load the cd-hit module (assuming it is available as a module on your system)
module load MMseqs2

# Get the input file name
input_file=$1

mkdir -p ${input_file%.ffn}_tmp/

# Output file name
output_file="${input_file%.ffn}.fa.MMSeqs2_cluster"

# Run cd-hit with specified parameters
mmseqs easy-linclust --threads 8 --split-memory-limit 100G "$input_file" "$output_file" ${input_file%.ffn}_tmp/ --min-seq-id 0.8 -c 0.9

