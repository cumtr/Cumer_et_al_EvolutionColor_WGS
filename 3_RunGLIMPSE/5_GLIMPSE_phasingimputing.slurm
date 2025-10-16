#!/bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time=03-00:00:00
#SBATCH --array [1-56]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=2G
#SBATCH --job-name=GLIMPSEphasing
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr
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
CHR_phasing_impute ./2_merged3KowlsVCFs/${ss}.vcf.gz /work/FAC/FBM/DEE/jgoudet/barn_owl/Common/3000owls/4_REF_PANEL/data/Out_ShapeIt/RPall502/PerSuperScaffolds/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_${ss}.vcf.gz \
./3_GLIMPSE_chunks/Chunks_${ss}.txt ./4_GLIMPSE_phasing 10000 30
