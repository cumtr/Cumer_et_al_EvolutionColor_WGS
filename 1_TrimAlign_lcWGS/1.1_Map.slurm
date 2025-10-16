#!/bin/bash

#SBATCH --partition cpu
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=12G
#SBATCH --job-name=Mapping
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr

# Load the modules
module load gcc
module load bwa
module load samtools

# short path
myPATH=/My/Path/

# define the function
Mapping_and_Sorting_READS(){

	#First argument passed to the script is the sample with path and extension
	sample=$1

	#Extract sample name as RingID_LANE_PLATE (without extension and path)
	sampname=$(basename -s "_R1_paired.fq.gz" ${sample})

	#Save short path
	myPATH=/My/Path/1_TrimmAlign/data

	##MAPPING
	in_forward=${myPATH}/0_trimm/${sampname}_R1_paired.fq.gz
	in_reverse=${myPATH}/0_trimm/${sampname}_R2_paired.fq.gz
	ref='/My/Path/to/my/Reference.fasta'
	mapped=${myPATH}/1_align/${sampname}_mapped.sam

	bwa mem -M -t 5 $ref $in_forward $in_reverse > $mapped

	## add RG and make bam

	#Extract info for BAM file
	ID=${sampname} #UNIQUE identifier --> here: RingID_LaneFlowcell
	BC=$(zcat ${sample}| head | awk -F: 'NR==1 {print $10}') #Adapter, extarcted from fq file
	SM=$(echo ${sampname} | cut -d'_' -f1) #RingID

	#Add ReadGroup and directly sort the BAM
	samtools addreplacerg -r "ID:$ID" \
		-r "BC:$BC" \
		-r "CN:GTF" \
		-r "LB:$SM" \
		-r "PL:ILLUMINA" \
		-r "PM:NovaSeq6000" \
		-r "SM:$SM" \
		$mapped | samtools sort -o ${myPATH}/1_align/${sampname}_sorted.bam -O bam

	#rm the mapped SAM
	rm $mapped
}

export -f Mapping_and_Sorting_READS

#Extracting the read to map and sort
runFILE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${myPATH}/samples_NAMES_for_ARRAYS_mapping_1.list)

#launching the function
Mapping_and_Sorting_READS ${runFILE}
