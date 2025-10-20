#!/bin/bash

#SBATCH -p cpu_p
#SBATCH -n 16 -N1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=kallisto_index_wheat_BMW
#SBATCH -t 24:00:00
#SBATCH --qos=cpu_normal

source ~/.bash_profile
source ~/.bashrc

kallisto=/home/pgsb/sepideh.jafarian/miniconda3/envs/kallisto/bin/kallisto

export TMPDIR=/home/pgsb/sepideh.jafarian/slurm_temp

DIR=/lustre/groups/pgsb/workspaces/sepideh.jafarian/bmw_project

cd $DIR

$kallisto index -i bmw_index.idx /lustre/groups/pgsb/datasets/annotations/wheatCS_IWGSCv1.1/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706_cds.fasta 
