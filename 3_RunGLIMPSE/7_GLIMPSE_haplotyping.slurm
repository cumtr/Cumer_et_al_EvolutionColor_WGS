#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --array [1-56]
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=GLIMPSEhaplotyping
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr
#SBATCH --account jgoudet_barn_owl

set -x
set -e

cd /users/elavanc1/3KOWLS/2_GLIMPSE/data

#Load the modules
module load gcc bcftools parallel htslib curl boost

#source the functions
source /users/elavanc1/3KOWLS/2_GLIMPSE/scripts/functions_GLIMPSE_1.0_3Kowls.sh

#Extract SS from task id
ss=$(sed -n "$SLURM_ARRAY_TASK_ID"p ./inputFILES/GLIMPSESuperScaffoldsManySNPs.list)

#Launch the function
CHR_sample_haplotypes ${ss} ./5_GLIMPSE_ligating ./6_GLIMPSE_haplotypes
