#!/bin/bash

#SBATCH --partition cpu
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=Trimming
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr

# load the modules
module load gcc
module load trimmomatic/0.39

# sort path
myPATH=/My/Path/


# define the function
Trimming(){

	#save the path
	myPATHOUT=/My/Path/Out/

	#First argument passed to the function is sample name WITH PATH (in case one sample has several files)
	sample=$1

	#get sample flowcell (1, 2 ou 3) (if an indv has been sequenced on two different flowcells but in the same lane)
	flowcell=$(dirname ${sample} | cut -d'/' -f11 | cut -d'-' -f2 | cut -c1)
	#Create a new sample name WITHOUT PATH for output
	NEWsamplename=$(echo $(basename ${sample}))P${flowcell}

	#input, no path needed (in sample name)
	forward_in=${sample}_R1_001.fastq.gz
	reverse_in=${sample}_R2_001.fastq.gz
	#outputs paths neeed
	forward_out_paired=${myPATHOUT}/1_TrimmAlign/data/0_trimm/${NEWsamplename}_R1_paired.fq.gz
	forward_out_unpaired=${myPATHOUT}/1_TrimmAlign/data/0_trimm/${NEWsamplename}_R1_unpaired.fq.gz
	reverse_out_paired=${myPATHOUT}/1_TrimmAlign/data/0_trimm/${NEWsamplename}_R2_paired.fq.gz
	reverse_out_unpaired=${myPATHOUT}/1_TrimmAlign/data/0_trimm/${NEWsamplename}_R2_unpaired.fq.gz

	#launch trommomatic with primers
	trimmomatic PE -threads 4 -phred33 ${forward_in} ${reverse_in} \
	${forward_out_paired} ${forward_out_unpaired} \
	${reverse_out_paired} ${reverse_out_unpaired} \
	MINLEN:70 \
	ILLUMINACLIP:/My/Path/to/trimmomatic/0.36/adapters/NexteraPE-PE.fa:2:30:10:3:"true"

	#rm unpared (won't be used)
	rm $forward_out_unpaired $reverse_out_unpaired
}

export -f Trimming

#get input file (from array ID)
run_FILE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${myPATH}/1_TrimmAlign/data/Sample_Names.list)

#launch the function
Trimming ${run_FILE}

