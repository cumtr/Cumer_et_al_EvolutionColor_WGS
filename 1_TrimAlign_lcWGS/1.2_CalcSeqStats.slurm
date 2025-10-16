#!/bin/bash

#SBATCH --partition cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=16G
#SBATCH --job-name=SeqStats
#SBATCH -o STD/%x_%A_%a.stdout
#SBTACH -e STD/%x_%A_%a.stderr

# load the modules
module load gcc samtools/1.12 parallel/20200822

# short path
myPATH=/My/Path/

# define function
CalcSeqStats(){

	#save path shortcut
	myPATH=myPATH=/My/Path/1_TrimmAlign/data

	#first argument passed to the script is file (with path, extensiom everything)
	file=$1
	#Then sample name (with lane and plate and everything)
	Sample=`basename -s "_sorted.bam" ${file}`
	NbReads=`samtools view ${file} | wc -l`
	NbReadsUnmapped=`samtools view ${file} | awk '$7=="*"' | wc -l`
	NbReadsWellMapped=`samtools view ${file} | awk '$5==60' | wc -l`
	echo -e "${NbReads}\t${NbReadsUnmapped}\t${NbReadsWellMapped}" > ${myPATH}/2_SeqStats/sampFiles/ReadsNB/${Sample}_ReadNBs.txt
	samtools view ${file} | awk '$7=="=" && $9>0 {print $9}' | sort -n | uniq -c > ${myPATH}/2_SeqStats/sampFiles/HistDistanceReads/${Sample}_hist.txt
}

export -f CalcSeqStats

#Apply the function
parallel -j 40 ::: CalcSeqStats ::: $(find ${myPATH}/1_align/ -type f)
