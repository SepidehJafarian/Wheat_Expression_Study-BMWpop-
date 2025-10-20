#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --cpus-per-task=3
#SBATCH --qos=cpu_normal
#SBATCH --reservation=pgsb_users
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=sequences_concatination
#SBATCH --output=/home/pgsb/sepideh.jafarian/bmw_project/log_bmw_mapping/%x.%j.%a.out
#SBATCH --error=/home/pgsb/sepideh.jafarian/bmw_project/log_bmw_mapping/%x.%j.%a.err     
#SBATCH --ntasks=1
#SBATCH --array=1-85
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sepideh.jafarian@helmholtz-munich.de

#set -x #verbose debugging from bash

source ~/.bash_profile
source ~/.bashrc


PROJECT_TABLE=/lustre/groups/pgsb/workspaces/sepideh.jafarian/bmw_project/project_table.txt



name=$(cat $PROJECT_TABLE | head -n${SLURM_ARRAY_TASK_ID} | tail -n1)

TOOL=/home/pgsb/sepideh.jafarian/miniconda3/envs/kallisto/bin/kallisto # version 3.27.0
OUTDIR=/lustre/groups/pgsb/workspaces/sepideh.jafarian/bmw_project/kallisto_bmw_results
INDIR=/lustre/groups/pgsb/workspaces/sepideh.jafarian/bmw_project/wheat_trimmed
INDEX=/lustre/groups/pgsb/workspaces/sepideh.jafarian/bmw_project/bmw_index.idx



fq1="$INDIR/${name}_1_trimmed.fq.gz"
fq2="$INDIR/${name}_2_trimmed.fq.gz"

mkdir $OUTDIR/$name


date
echo "Forward read: $fq1"
echo "Reverse read: $fq2"
echo "Mapping sample $name"


$TOOL quant --index $INDEX --output $OUTDIR/$name $fq1 $fq2 




echo "Kallito has run successfully."
