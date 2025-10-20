#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --cpus-per-task=3
#SBATCH --qos=cpu_normal
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=sequences_concatination
#SBATCH --output=/home/pgsb/sepideh.jafarian/bmw_project/logs/%x.%j.%a.out
#SBATCH --error=/home/pgsb/sepideh.jafarian/bmw_project/logs/%x.%j.%a.err     
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sepideh.jafarian@helmholtz-munich.de



INPUT_DIR=/lustre/groups/pgsb/workspaces/sepideh.jafarian/bmw_project/wheat_raw_data
OUTPUT_DIR=/lustre/groups/pgsb/workspaces/sepideh.jafarian/bmw_project/wheat_raw_data_merged

# Change to the input directory
cd "$INPUT_DIR" || exit


for file in BMW_*_1.fq.gz; do
  sample=$(basename "$file" _1.fq.gz)
  sample_number=$(echo "$sample" | cut -d'_' -f2)
  echo "$sample"
  echo "${sample_number}"
  cat "$file" >> "$OUTPUT_DIR/BMW_$sample_number"_forward_concatenated.fq.gz
done


for file in BMW_*_2.fq.gz; do
  sample=$(basename "$file" _2.fq.gz)
  sample_number=$(echo "$sample" | cut -d'_' -f2)
  echo "$sample"
  echo "${sample_number}"
  cat "$file" >> "$OUTPUT_DIR/BMW_$sample_number"_reverse_concatenated.fq.gz
done

