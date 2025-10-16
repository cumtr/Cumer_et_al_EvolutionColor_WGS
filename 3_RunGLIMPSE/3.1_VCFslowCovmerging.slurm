#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=500G
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr
#SBATCH --account jgoudet_barn_owl

set -x
set -e

#Load the modules
module load gcc bcftools parallel htslib curl boost

#source the functions
source functions_GLIMPSE_1.0_3Kowls.sh

#Launch the function
INDVs_VCF_merging ./Chromosome.list ./OutDir 30


### For GLIMPSE we need to PAR and UNPAR region of SS42 separately !!! So we'll just split it there
#PAR region
bcftools view -t Super-Scaffold_42:1-4415000,Super-Scaffold_42:7364001-7389000 -O z -o ./Super-Scaffold_42PAR.vcf.gz ./Super-Scaffold_42.vcf.gz
bcftools index Super-Scaffold_42PAR.vcf.gz
#Non PAR region
bcftools view -t Super-Scaffold_42:4415001-7364000,Super-Scaffold_42:7389001-41426968 -O z -o ./Super-Scaffold_42nonPAR.vcf.gz ./Super-Scaffold_42.vcf.gz
bcftools index Super-Scaffold_42nonPAR.vcf.gz

#rm the UNSPLIT SS42
rm Super-Scaffold_42.vcf.gz*

### Process the RP if needed
# make sure the ref pannel files are also split !
#Par region
bcftools view -t Super-Scaffold_42:1-4415000,Super-Scaffold_42:7364001-7389000 -O z -o RP_Super-Scaffold_42PAR.vcf.gz RP_Super-Scaffold_42.vcf.gz
bcftools index RP_Super-Scaffold_42PAR.vcf.gz
#NonPAR region
bcftools view -t Super-Scaffold_42:4415001-7364000,Super-Scaffold_42:7389001-41426968 -O z -o RP_Super-Scaffold_42nonPAR.vcf.gz RP_Super-Scaffold_42.vcf.gz
bcftools index RP_Super-Scaffold_42nonPAR.vcf.gz
