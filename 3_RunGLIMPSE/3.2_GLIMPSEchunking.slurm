#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=GLIMPSEchunking
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr

set -x
set -e

#Load the modules
module load gcc bcftools parallel htslib curl boost

#source the functions
source functions_GLIMPSE_1.0_3Kowls.sh

#First step is to re-create the list of ss for GLIMPSE (becasue we splitted SS42 into PAR and unPAR regions)
find ./2_merged3KowlsVCFs -name "*.vcf.gz" | xargs -I {} bash -c 'basename -s ".vcf.gz" {}' > Chromosome.list

#Then simply launch the function
CHR_chunking Chromosome.list RP_Dir GLIMPSE_chunks 2
