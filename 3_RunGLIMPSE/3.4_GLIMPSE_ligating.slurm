#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --array [1-56]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=GLIMPSEligating
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr

set -x
set -e

#Load the modules
module load gcc bcftools parallel htslib curl boost

#source the functions
source ./functions_GLIMPSE_1.0_3Kowls.sh

#Extract SS from task id
ss=$(sed -n "$SLURM_ARRAY_TASK_ID"p ./Chromosome.list)

#Launch the function
CHR_ligating_chunks ${ss} ./GLIMPSE_phasing ./GLIMPSE_ligating 1
