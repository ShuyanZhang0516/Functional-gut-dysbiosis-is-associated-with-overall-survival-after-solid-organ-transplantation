#!/bin/bash

#SBATCH --job-name=seqtk
#SBATCH --output=seqtk_%j.log
#SBATCH --error=seqtk_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=120G

module load seqtk

input_file=$1
output_file=$2

seqtk seq -A "$input_file" | \
awk '
  /^>/ { 
    if (length(sequence) >= 500) 
      print header "\n" sequence; 
    header = $0; 
    sequence = ""; 
    next 
  } 
  { 
    sequence = sequence "" $0 
  } 
  END { 
    if (length(sequence) >= 500) 
      print header "\n" sequence 
  }' >  "$output_file"
