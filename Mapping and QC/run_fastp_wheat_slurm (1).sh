#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --cpus-per-task=3
#SBATCH --qos=cpu_normal
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=fastp_bmw_wheat
#SBATCH --output=/home/pgsb/sepideh.jafarian/bmw_project/logs/%x.%j.%a.out
#SBATCH --error=/home/pgsb/sepideh.jafarian/bmw_project/logs/%x.%j.%a.err
#SBATCH --array=1-85     
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sepideh.jafarian@helmholtz-munich.de

#set -x

source ~/.bash_profile
source ~/.bashrc

PROJECT_TABLE=/lustre/groups/pgsb/workspaces/sepideh.jafarian/bmw_project/project_table.txt

name=$(cat $PROJECT_TABLE | head -n${SLURM_ARRAY_TASK_ID} | tail -n1)

TOOL=/home/pgsb/sepideh.jafarian/miniconda3/envs/QC/bin/fastp # version 0.23.4
OUTDIR=/lustre/groups/pgsb/workspaces/sepideh.jafarian/bmw_project/wheat_trimmed
INDIR=/lustre/groups/pgsb/workspaces/sepideh.jafarian/bmw_project/wheat_raw_data_merged


fq1_in="$INDIR/${name}_forward_concatenated.fq.gz"
fq2_in="$INDIR/${name}_reverse_concatenated.fq.gz"
fq1_out="$OUTDIR/${name}_1_trimmed.fq.gz"
fq2_out="$OUTDIR/${name}_2_trimmed.fq.gz"

date
echo "Forward read: $fq1_in"
echo "Reverse read: $fq2_in"
echo "Processing sample $name"

$TOOL -i $fq1_in -o $fq1_out -I $fq2_in -O $fq2_out \
	--adapter_sequence AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --adapter_sequence_r2 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG --length_required 50 \
       	--html "$OUTDIR/${name}_fastp.html" --json "$OUTDIR/${name}_fastp.json" --thread 3 

echo "Fastp has run successfully."

# sbatch -N 1 /lustre/groups/pgsb/workspaces/sepideh.jafarian/bmw_project/run_fastp_wheat_slurm.sh
