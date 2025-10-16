#!/bin/bash

#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --job-name=StatBams
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr

#Load the modules
module load gcc
module load samtools

#short path
myPATH=/My/Path/1_TrimmAlign/data

# define function
StatBamFiles(){

	#First argument passed to the script is the sample with path and extension
	sample=$1

        #Extract sample name as RingID_LANE_PLATE (without extension and path)
        sampname=$(basename -s "" ${sample})

	#Perpare output files
	input="$sampname"_markdup.bam
	output="$sampname"_markdup.stat

	#Save short path
	myPATH=/scratch/elavanc1/3Kowls/1_TrimmAlign/data

	### Process the files ####

	# Extract mean coverage
	samtools depth -Q 30 $sample |  awk -v r=$sampname '{sum+=$3} END { print r,sum/1249867532}' > ${myPATH}/3_mergedBAMs/$output

}

export -f StatBamFiles

#Extracting the read to map and sort
runFILE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${myPATH}/samples_NAMES_for_STATS_Bams.list)

#launching the function
StatBamFiles ${runFILE}


