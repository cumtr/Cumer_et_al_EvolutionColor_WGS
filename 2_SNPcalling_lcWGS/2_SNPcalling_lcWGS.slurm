#!/bin/bash

#SBATCH --time=03-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=4G
#SBATCH --job-name=SNPsCalling3Kowls
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr

set -x
set -e


#Load the modules
module load gcc bcftools parallel htslib curl boost

# define the fonction 
#This function requires BCFTools & parallel --> on curnagl: module load gcc/9.3.0 bcftools/1.12 parallel/20200822
REF_SNPs_Calling_perINDV(){

	#First argument is the path to bam files
	BAMpath=$1

	#Second argument passed to the function is the CHRlist from REF panel
	CHRlist=$2

	#Second argument is the fasta REF genome
	fastaREF=$3

	#Third argument is output directory
	outputDIR=$4

	#Forth argument is threads number
	threads=$5

        #Creating the TMP sample specific samples list file (for ploidy)
	parallel --rpl '{/..} s:^.*/::;s:_markdup\.bam::;' -j${threads} 'grep "{1/..}" /Path/to/INDVSex.tsv > ./inputFILES/TMPlist_{1/..}.tsv' ::: $(find "${BAMpath}" -maxdepth 1 -name "*.bam")

	#Calling the SNPs
	parallel --rpl '{/..} s:^.*/::;s:_markdup\.bam::;' -j${threads} "bcftools mpileup -f ${fastaREF} -I -E -a 'FORMAT/DP' -T \
	/Path/to/RefPannel_{2}.vcf.gz -r {2} {1} -O u | bcftools call --ploidy-file /Path/to/ploidyFile.txt \
	-S ./inputFILES/TMPlist_{1/..}.tsv -Aim -C alleles -T /Path/to/SNPs_{2}.tsv.gz -O z -o ${outputDIR}/{1/..}_{2}.vcf.gz" \
	::: $(find "${BAMpath}" -maxdepth 1 -name "*.bam") ::: $(grep -v '#' ${CHRlist} | awk '{print $0}')

	#Indexing
	parallel --rpl '{/..} s:^.*/::;s:_markdup\.bam::;' -j${threads} bcftools index ${outputDIR}/{1/..}_{2}.vcf.gz ::: $(find "${BAMpath}" -maxdepth 1 -name "*.bam") ::: $(grep -v '#' ${CHRlist} | awk '{print $0}')

	#rm tmp sample list files
	rm ./inputFILES/TMPlist_*.tsv
}


# Run the fuction

REF_SNPs_Calling_perINDV /Path/to/BamFileDir/ /Path/to/Chromosome.list /Path/to/Reference.fasta /Path/to/OtuputDir/ 30

