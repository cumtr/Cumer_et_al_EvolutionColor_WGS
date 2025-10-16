#This scripts contains all the fonctions used in the GLIMPSE 1.0 3K Owls pipeline

### Splitting REF panel per CHRs ###

#This function requires BCFTools & parallel --> on curnagl: module load gcc/9.3.0 bcftools/1.12 parallel/20200822

CHR_Splitting(){

	#First argument passed to the function is VCF
	REF_panel=$1

	#Second argument passed to the function is CHR list with one line per CHR name (header starting with # authorized)
	CHRlist=$2

	#Third argument passed to the function is output directory
	outputDIR=$3

	#Forth argument is threads number
	threads=$4

	#Launching the subsetting in parallel
	parallel -j${threads} bcftools view -r {1} -O z -o ${outputDIR}/{1}.vcf.gz {2} ::: $(grep -v '#' ${CHRlist} | awk '{print $0}') ::: ${REF_panel}

	#Indexing the subset VCF
	parallel -j${threads} bcftools index ${outputDIR}/{}.vcf.gz ::: $(grep -v '#' ${CHRlist} | awk '{print $0}')

}


### Extracting Variable positions from REF PANEL (or any VCF) ###

#This function requires BCFTools & parallel --> on curnagl: module load gcc/9.3.0 bcftools/1.12 parallel/20200822 and htslib

REF_SNPs_Extraction(){

	#First argument passed to the function is CHR list with one line per CHR name (header starting with # authorized)
	VCFlist=$1

	#Second argument passed to the function is output directory
	outputDIR=$2

	#Third argument is threads number
	threads=$3

	#Pass this to TSV
	parallel --rpl '{/..} s:^.*/::;s:\.[^.]+$::;s:\.[^.]+$::;' -j${threads} "bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {} | bgzip -c > ${outputDIR}/SNPs_{/..}.tsv.gz" ::: $(grep -v '#' ${VCFlist} | awk '{print $0}')

	#Index again
	parallel --rpl '{/..} s:^.*/::;s:\.[^.]+$::;s:\.[^.]+$::;' -j${threads} tabix -s1 -b2 -e2 ${outputDIR}/SNPs_{/..}.tsv.gz ::: $(grep -v '#' ${VCFlist} | awk '{print $0}')

}


### Merging the INDIVIDUALS BAM (per CHR) ###

#This function requires BCFTools & parallel --> on curnagl: module load gcc/9.3.0 bcftools/1.12 parallel/20200822

INDVs_VCF_merging(){

	#First argument passed to the function is the CHRlist from REF panel
	CHRlist=$1

	#Second argument is path the output directory
	outputDIR=$2

	#Third argument is threads numbers (in TOTAL, I divide it between number of jobs and cores per job for merging)
	threads=$3

	#Generating the files with one individual VCF per line for each CHR
	parallel -j${threads} "ls ./SNPcalling/*_{}.vcf.gz > ./inputFILES/TMP_{}_INDVs.list" ::: $(grep -v '#' ${CHRlist} | awk '{print $0}')

	#increase max number of open files allowed for one process (default 1024)
	ulimit -n 3000

	#Merging the VCF
	parallel -j $(echo $((${threads}/3))) "bcftools merge --threads $(echo $((${threads}/10))) -m none -r {} -l ./inputFILES/TMP_{}_INDVs.list -O z -o ${outputDIR}/{}.vcf.gz" ::: $(grep -v '#' ${CHRlist} | awk '{print $0}')

	#Indexing
	parallel -j${threads} bcftools index -f ${outputDIR}/{}.vcf.gz ::: $(grep -v '#' ${CHRlist} | awk '{print $0}')

	#Rm the per INDV vcfs
	#Rm the files with one individual VCF per line for each CHR
	rm ./inputFILES/TMP_*_INDVs.list

}

### Separate CHR into chunks with the GLIMPSE_chunk tool ###

#This function requires GLIMPSE & parallel --> on curnagl: module load gcc/9.3.0 xz bzip2 curl htslib boost parallel/20200822

CHR_chunking(){

	#First argument passed to the function is the CHRlist from REF panel
	CHRlist=$1

	#Second argument is path to REF VCFs
	REFpath=$2

	#Third argument is output directory
	outputDIR=$3

	#Forth argument is threads number (per INDV for the different CHRs)
	threads=$4

	#Save shortcut for GLIMPSE_chunk executable
	GLIMPSE_chunk=/path/to/GLIMPSE_chunk

	#Chunking the CHR
	parallel --rpl '{/..} s:[nP].*::;' -j${threads} "${GLIMPSE_chunk} --input ${REFpath}/RPAll502_Libnames_Filtered_WhatshapPhased_ShapeitPhased_{}.vcf.gz \
	--region {/..} --window-size 2000000 --buffer-size 200000 --output ${outputDIR}/Chunks_{}.txt" ::: $(grep -v '#' ${CHRlist} | awk '{print $0}')

}

### Phase & Impute CHR (per chunks) with the GLIMPSE_phase tool ###

#This function requires GLIMPSE & BCFTools --> on curnagl: module load gcc/9.3.0 xz bzip2 curl htslib boost bcftools/1.12

#This function is designed to be used in an array with on SS per ARRAY task !!!

CHR_phasing_impute(){

	#First argument passed to the function is the target VCF (low cov)
	VCF=$1

	#Extract CHR Name from VCF
	CHR=$(basename -s ".vcf.gz" ${VCF})

	#Second argument is the REF VCFs
	REF=$2

	#Third argument is the the chunk file (generated in CHR_chunking)
	CHUNKfile=$3

	#Forth argument is output directory
	outputDIR=$4

	#Fifth argument is ne estimated
	Ne=$5

	#Sixth argument is threads number
	threads=$6

	#Save shortcut for GLIMPSE_phase executable
	GLIMPSE_phase=/path/to/GLIMPSE_phase

	#the phasing is slightyl different for autosaumes and sexual CHR (actually only one parameter different bu anyway)

	if [[ ${CHR} == "Super-Scaffold_13" || ${CHR} == "Super-Scaffold_42nonPAR" ]]; then
	    param=3
	else
	    param=2
	fi

	#Phasing + Impute, we don't need parallel here because we can use several threads which will parallelize per INDV
	while IFS="" read -r LINE || [ -n "$LINE" ];
	do
		printf -v ID "%0${param}d" $(echo $LINE | cut -d" " -f1)
		IRG=$(echo $LINE | cut -d" " -f3)
		ORG=$(echo $LINE | cut -d" " -f4)
		OUT=${outputDIR}/${CHR}_${ID}.bcf
		${GLIMPSE_phase} --input ${VCF} --reference ${REF} --thread ${threads} --ne ${Ne} --burnin 100 --main 15 --input-region ${IRG} --map ./RecombinationMap/${CHR}.map.txt \
		--samples-file ${CHR}.ploidy --output-region ${ORG} --output ${OUT}
		bcftools index -f ${OUT}
	done < ${CHUNKfile}

}

### Ligate the different chunks together with the GLIMPSE_ligate tool ###

#This function requires GLIMPSE & parallel & BCFTools --> on curnagl: module load gcc/9.3.0 xz bzip2 curl htslib boost parallel/20200822 bcftools/1.12

CHR_ligating_chunks(){

	#First argument passed to the function is the Super Scaffold
	CHR=$1

	#Second argument is the path to phased & Imputed BCFs to ligate
	VCFpath=$2

	#Third argument is output directory
	outputDIR=$3

	#Forth argument is threads number
	threads=$4

	#Save shortcut for GLIMPSE_ligate executable
	GLIMPSE_ligate=/path/to/GLIMPSE_ligate

	#Create the lists per VCF
	find ${VCFpath}/ -name "${CHR}_*.bcf" > ${outputDIR}/TMP_${CHR}_chunks.list

	#If it's the nonPAR region of the ZCHR, we need to ligate males and females separately then merge the final VCF
	if [[ ${CHR} == "Super-Scaffold_13" || ${CHR} == "Super-Scaffold_42nonPAR" ]]; then

	    #Loop through chunks file and extract only females AND change FPLOIDY line in VCF header, repass to BCF
	    for file in $(cat ${outputDIR}/TMP_${CHR}_chunks.list); do chunk=$(basename -s ".bcf" ${file} | cut -d"_" -f3); bcftools view -S ./inputFILES/3Kowls_FEMALES.list -O v -o ${VCFpath}/${CHR}_${chunk}_FEMALES.vcf ./${VCFpath}/${CHR}_${chunk}.bcf;\
	    sed -i 's/##FPLOIDY=-2/##FPLOIDY=1/' ${VCFpath}/${CHR}_${chunk}_FEMALES.vcf; bcftools view -O u -o ${VCFpath}/${CHR}_${chunk}_FEMALES.bcf ${VCFpath}/${CHR}_${chunk}_FEMALES.vcf; rm ${VCFpath}/${CHR}_${chunk}_FEMALES.vcf;\
	    bcftools index ${VCFpath}/${CHR}_${chunk}_FEMALES.bcf; done

	    #same for males
	    for file in $(cat ${outputDIR}/TMP_${CHR}_chunks.list); do chunk=$(basename -s ".bcf" ${file} | cut -d"_" -f3); bcftools view -S ./inputFILES/3Kowls_MALES.list -O v -o ${VCFpath}/${CHR}_${chunk}_MALES.vcf ./${VCFpath}/${CHR}_${chunk}.bcf;\
	    sed -i 's/##FPLOIDY=-2/##FPLOIDY=2/' ${VCFpath}/${CHR}_${chunk}_MALES.vcf; bcftools view -O u -o ${VCFpath}/${CHR}_${chunk}_MALES.bcf ${VCFpath}/${CHR}_${chunk}_MALES.vcf; rm ${VCFpath}/${CHR}_${chunk}_MALES.vcf;\
	    bcftools index ${VCFpath}/${CHR}_${chunk}_MALES.bcf; done

            #Recreate TMP file with FEMALES only
            find ${VCFpath}/ -name "${CHR}_*_FEMALES.bcf" > ${outputDIR}/TMP_${CHR}_chunks_FEMALES.list

            #Recreate TMP file with MALES only
            find ${VCFpath}/ -name "${CHR}_*_MALES.bcf" > ${outputDIR}/TMP_${CHR}_chunks_MALES.list

	    #LIGATION FEMALES
	    ${GLIMPSE_ligate} --input ${outputDIR}/TMP_${CHR}_chunks_FEMALES.list --output ${outputDIR}/${CHR}_phased_imputed_FEMALES.bcf

            #LIGATION FEMALES
            ${GLIMPSE_ligate} --input ${outputDIR}/TMP_${CHR}_chunks_MALES.list --output ${outputDIR}/${CHR}_phased_imputed_MALES.bcf

	    #rm the chunk list
	    rm ${outputDIR}/TMP_${CHR}_chunks*.list

	#else we just do regular ligation
	else

	    #Ligating the CHR
	    ${GLIMPSE_ligate} --input ${outputDIR}/TMP_${CHR}_chunks.list --output ${outputDIR}/${CHR}_phased_imputed.bcf

	    #rm the chunk list
            rm ${outputDIR}/TMP_${CHR}_chunks.list

	fi

}

### Get the most probabl haplotype  with the GLIMPSE_sample tool ###

CHR_sample_haplotypes(){

	#First argument passed to the function contains Super-Scaffold
	CHR=$1

        #Second argument passed to the function is the directory containing the phased & Imputed & ligated VCF (one per ss)
        VCFpath=$2

        #Third argument is output directory
        outputDIR=$3

        #Save shortcut for GLIMPSE_ligate executable
        GLIMPSE_sample=/path/to/GLIMPSE_sample

	#If it's the nonPAR region of the ZCHR, we need to ligate males and females separately then merge the final VCF
        if [[ ${CHR} == "Super-Scaffold_13" || ${CHR} == "Super-Scaffold_42nonPAR" ]]; then

	    #FEMALES first
	    ${GLIMPSE_sample} --input ${VCFpath}/${CHR}_phased_imputed_FEMALES.bcf --solve --output ${outputDIR}/${CHR}_phased_imputed_haplotypes_FEMALES.bcf

	    #then MALES
	    ${GLIMPSE_sample} --input ${VCFpath}/${CHR}_phased_imputed_MALES.bcf --solve --output ${outputDIR}/${CHR}_phased_imputed_haplotypes_MALES.bcf

	    #The we merge both into a VCF (which we can edit with sed for FPLOIDY)
	    bcftools merge -O v -o ${outputDIR}/${CHR}_phased_imputed_haplotypes_nogoodorder.vcf -l <(find ${outputDIR}/${CHR}_phased_imputed_haplotypes_*ALES.bcf)

	    #we change FPLOIDY
	    sed -i 's/##FPLOIDY=1/##FPLOIDY=-2/' ${outputDIR}/${CHR}_phased_imputed_haplotypes_nogoodorder.vcf

	    #Change order of the samples to match all other VCFs AND pass to bcf (for consistency with other SS)
	    bcftools view -S ./inputFILES/3KowlsINDVs.list -O u -o ${outputDIR}/${CHR}_phased_imputed_haplotypes.bcf ${outputDIR}/${CHR}_phased_imputed_haplotypes_nogoodorder.vcf
	    bcftools index ${outputDIR}/${CHR}_phased_imputed_haplotypes.bcf

	    #rm no good order bcf
            rm ${outputDIR}/${CHR}_phased_imputed_haplotypes_nogoodorder.vcf*

	    #rm FEMALES and MALES only BCFs
	    rm ${outputDIR}/${CHR}_phased_imputed_haplotypes_*ALES.bcf*

	#else regular haplotyping
	else

            ${GLIMPSE_sample} --input ${VCFpath}/${CHR}_phased_imputed.bcf --solve --output ${outputDIR}/${CHR}_phased_imputed_haplotypes.bcf

	fi

}
