#!/bin/bash

#SBATCH --time=03-00:00:00
#SBATCH --array [1-56]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=2G
#SBATCH --job-name=GLIMPSEphasing
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr

set -x
set -e

#Load the modules
module load gcc bcftools parallel htslib curl boost

#source the functions
source functions_GLIMPSE_1.0_3Kowls.sh

#Extract SS from task id
ss=$(sed -n "$SLURM_ARRAY_TASK_ID"p ./Chromosome.list)

#Launch the function
CHR_phasing_impute /path/to/lcWGS_${ss}.vcf.gz ./RP_${ss}.vcf.gz ./GLIMPSE_chunks/Chunks_${ss}.txt ./GLIMPSE_phasing 10000 30
