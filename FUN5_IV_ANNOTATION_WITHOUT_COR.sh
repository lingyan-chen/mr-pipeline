#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=annovar
#SBATCH --ntasks=1
#SBATCH --time=2:59:59
#SBATCH -p skylake-himem
# #SBATCH --array=2-369
#SBATCH --output=slurm/annovar_%A.out
#SBATCH --error=slurm/annovar_%A.err
#SBATCH --mail-type=FAIL


#! Optionally modify the environment seen by the application
## load modules for analyses ##
module purge
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)odule
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:
module load slurm
module load use.own
module load gcc-5.4.0-gcc-4.8.5-fis24gg
module load htslib-1.9-gcc-5.4.0-p2taavl
module load mpfr-4.0.1-gcc-5.4.0-bernlyg
module load r-3.5.1-gcc-5.4.0-lk6xb2i

## load modules for genetics analysis from ceuadmin ##
module load ceuadmin/qctool/v1.4-linux-x86_64
module load ceuadmin/LDstore/1.1
module load ceuadmin/plink-bgi/archive/1.90
module load ceuadmin/tabix/0.2.6

## Make R aware of the new library location
echo "R_LIBS_USER=/home/lc753/privatemodules/RLibs/RLibs-r-3.5.1"    > ~/.Renviron

## set working directory ##
DIR="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/MR/olink_vs_stroke"
DIR_CODE="${DIR}/code"
DIR_slurm="${DIR_CODE}/slurm"
DIR_DATA="${DIR}/data"

## GWAS summary stats ##
## exposure: olink proteom GWAS ##
DIR_IN_proteinGWAS="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/pQTL/data/olink"
DIR_IN_proteinGWAS_tabix=${DIR_IN_proteinGWAS}
DIR_OUT_proteinGWAS="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/MR/olink_vs_stroke/data/olink_summary_stats"

# ## exposure: metabolon GWAS ##
# DIR_IN_metabGWAS="/rds/project/jmmh2/rds-jmmh2-projects/metabolon_metabolomics/interval/gwas/interval_epic_meta_analysis"
# DIR_IN_metabGWAS_tabix="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/gwas/metabolon/METABOLON"
# DIR_OUT_metabGWAS="${DIR_DATA}/metabolon_summary_stats"

## Primary outcome: stroke GWAS ##
DIR_IN_strokeGWAS="/rds/project/jmmh2/rds-jmmh2-projects/proteomics_twas/prediction/data_sets/gwas_summary_statistics/phenoscanner_format_tabixed"
DIR_OUT_strokeGWAS="${DIR_DATA}/megastroke_summary_stats"
DIR_IN_strokeGWAS_SNPlist="/rds/user/lc753/hpc-work/gwas/megastroke"

# ## Secondary outcome: stroke risk factors GWAS ##
# DIR_IN_riskFactorGWAS="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/gwas/StrokeRiskFactors"
# DIR_OUT_riskFactorGWAS="${DIR_DATA}/StrokeRiskFactors_summary_stats"
# DIR_IN_riskFactorGWAS_SNPlist="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/gwas/StrokeRiskFactors"

# ## ukbbGWAS ##
# DIR_IN_ukbbGWAS="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/gwas/UKBB_SAIGE_HRC/ALL"
# DIR_OUT_ukbbGWAS="${DIR_DATA}/ukbb_summary_stats"
# DIR_IN_ukbbGWAS_SNPlist="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/gwas/UKBB_SAIGE_HRC"

## set directories for ANNOVAR
DIR_SNPLIST="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/MR/olink_vs_stroke/data/olink_pqtl"
DIR_ANNOVAR="/home/lc753/rds/hpc-work/apps/annovar"
DIR_ANNOVAR_OUT="${DIR_SNPLIST}/annovar"

## creat directory for output if not exist ##
if [ ! -d ${DIR_ANNOVAR_OUT} ]
then
	mkdir -p ${DIR_ANNOVAR_OUT}
fi


# ## read gwas data for my metabolite -- lineID is the ${SLURM_ARRAY_TASK_ID}
# # SLURM_ARRAY_TASK_ID=2
# # myMetabID="M37181"
# # myMetabID="M46356"
# myMetabID=$(awk -F'\t' -v metabolite_id_line=${SLURM_ARRAY_TASK_ID} 'NR==metabolite_id_line{print $1}' ${DIR}/List_Metabolites_ALL.txt)
# echo ${myMetabID}


##### ---------------------------------- #####
##### step 1. prepare data for "annovar" #####
##### ---------------------------------- #####
# for mySNPset in p1_5e8_r2_0.01 p1_5e8_r2_0.1 p1_5e8_r2_0.2; do
	# # mySNPset="p1_5e8_r2_0.01"
	# echo ${mySNPset}
	
	# ## combine exposure summary stats for all Metabolites:
	# cat ${DIR_SNPLIST}/${mySNPset}/exposure_sumstats_* >>${DIR_SNPLIST}/exposure_sumstats_${mySNPset}.txt
	# egrep -w "effect_allele.exposure" ${DIR_SNPLIST}/exposure_sumstats_${mySNPset}.txt | head -1 >${DIR_SNPLIST}/exposure_sumstats_${mySNPset}_header
	# egrep -vw "effect_allele.exposure" ${DIR_SNPLIST}/exposure_sumstats_${mySNPset}.txt >${DIR_SNPLIST}/exposure_sumstats_${mySNPset}
	# cat ${DIR_SNPLIST}/exposure_sumstats_${mySNPset}_header  ${DIR_SNPLIST}/exposure_sumstats_${mySNPset} >${DIR_SNPLIST}/exposure_sumstats_${mySNPset}.txt

	# ## subset columns into annovar input format:
	# cat ${DIR_SNPLIST}/IVs_LDclumped_at_${mySNPset}.txt | awk '{print $1}' >${DIR_SNPLIST}/SNPID_for_${mySNPset}.txt
	# sed -i -e 's/chr//g'  ${DIR_SNPLIST}/SNPID_for_${mySNPset}.txt 
	# awk -F":"  '{print $1, $2, $2}' ${DIR_SNPLIST}/SNPID_for_${mySNPset}.txt >${DIR_SNPLIST}/SNPID_for_${mySNPset}_update.txt
	# cat ${DIR_SNPLIST}/IVs_LDclumped_at_${mySNPset}.txt | awk '{print $2, $3}' >${DIR_SNPLIST}/SNPallele_for_${mySNPset}.txt

	# ## combine snp position and allele columns:
	# paste -d'\t' ${DIR_SNPLIST}/SNPID_for_${mySNPset}_update.txt  ${DIR_SNPLIST}/SNPallele_for_${mySNPset}.txt | sort | uniq >${DIR_SNPLIST}/annovar_snplist_from_${mySNPset}.txt

	# ## tidy up space:
	# # rm ${DIR_SNPLIST}/exposure_sumstats_${mySNPset}_header
	# # rm ${DIR_SNPLIST}/exposure_sumstats_${mySNPset}
	# rm ${DIR_SNPLIST}/SNPID_for_${mySNPset}.txt
	# rm ${DIR_SNPLIST}/SNPID_for_${mySNPset}_update.txt
	# rm ${DIR_SNPLIST}/SNPallele_for_${mySNPset}.txt
	
	# echo "format ${mySNPset} snplist for annovar done."
	
# done


# R code:
# ##################################################################################
# #####  combine all IVs for each LD clumped setting - for all olink proteins  #####
# ##################################################################################
# ## set working directory ##
# DIR="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/MR/olink_vs_stroke"
# List_SNPsets <- c("p1_5e8_r2_0.01", "p1_5e8_r2_0.1", "p1_5e8_r2_0.2")

# for(mySNPset in List_SNPsets){
  # print(mySNPset)
  
  # dir="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/MR/olink_vs_stroke/data/olink_pqtl"
  # file_list <- list.files(path = paste(dir, "/", mySNPset, sep = ""), full.names = FALSE, ignore.case = FALSE)
  
  # ## an example: i = 1 ##
  # myProtein <- file_list[1]
  # myResult_ALL <- read.delim(file = paste(dir, "/", mySNPset, "/", file_list[1], sep = ""), header = T, sep = "")
  # myResult_ALL <- cbind(myResult_ALL[, c(1:26)], myProtein)
  # colnames(myResult_ALL)[1] <- "alternate_ids"
  # myResult_ALL$FreqB <- (myResult_ALL$all_AB + myResult_ALL$all_BB *2)/((myResult_ALL$all_AA +myResult_ALL$all_AB + myResult_ALL$all_BB) *2 )
  
  # ## loop through all files ##
  # for (i in 2: length(file_list)){
    # myProtein <- file_list[i]
    # myResult <- read.delim(file = paste(dir, "/", mySNPset, "/", file_list[i], sep = ""), header = T, sep = "")

    # if(dim(myResult)[1] == 0){
      # print(paste(myProtein, " is not available for MR", sep = ""))
    # }else{
      # myResult <- cbind(myResult[, c(1:26)], myProtein)
      # colnames(myResult)[1] <- "alternate_ids"
      # myResult$FreqB <- (myResult$all_AB + myResult$all_BB *2)/((myResult$all_AA +myResult$all_AB + myResult$all_BB) *2 )
      
      # ## combine all IVs
      # myResult_ALL <- rbind(myResult_ALL, myResult)
    # }
  # }
  
  # ## check if the allele T is recoded as "TRUE"
  # for(k in 1: dim(myResult_ALL)[1]){
    # if(myResult_ALL$alleleA[k] == TRUE){ myResult_ALL$alleleA[k] <- as.character("T") }else if(is.na(myResult_ALL$alleleA[k]) == TRUE){ myResult_ALL$alleleA[k] <- as.character("NA")}else{myResult_ALL$alleleA[k] <- myResult_ALL$alleleA[k]}
    # # if(myResult_ALL$alleleB[k] == TRUE){ myResult_ALL$alleleB[k] <- as.character("T") }else if(is.na(myResult_ALL$alleleB[k]) == TRUE){ myResult_ALL$alleleB[k] <- as.character("NA")}else{myResult_ALL$alleleB[k] <- myResult_ALL$alleleB[k]}
  # }
  
  # write.table(myResult_ALL, file = paste(dir, "/IVs_LDclumped_at_", mySNPset, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
# }


for mySNPset in p1_5e8_r2_0.01 p1_5e8_r2_0.1 p1_5e8_r2_0.2; do
	# mySNPset="p1_5e8_r2_0.01"
	echo ${mySNPset}
	## subset columns into annovar input format:
	cat ${DIR_SNPLIST}/IVs_LDclumped_at_${mySNPset}.txt | awk '{print $3,$4, $4, $5, $6}' | sort | uniq >${DIR_SNPLIST}/annovar_snplist_from_${mySNPset}.txt
done


##### ---------------------------------------- #####
##### step 2. Run "annovar" for IVs annotation #####
##### ---------------------------------------- #####
cd ${DIR_ANNOVAR}
for mySNPset in p1_5e8_r2_0.01 p1_5e8_r2_0.1 p1_5e8_r2_0.2; do
	# mySNPset="p1_5e8_r2_0.01"
	echo ${mySNPset}	

	## Run annovar to annotate the target SNPlist ##
	${DIR_ANNOVAR}/table_annovar.pl  ${DIR_SNPLIST}/annovar_snplist_from_${mySNPset}.txt ${DIR_ANNOVAR}/humandb/ -buildver hg19 -out ${DIR_ANNOVAR_OUT}/annovar_snplist_from_${mySNPset} -remove -protocol refGene,cytoBand,exac03,avsnp147,kaviar_20150923,hrcr1 -operation gx,r,f,f,f,f  -nastring . -csvout -polish 

done


# ##### ------------------------------------------------ #####
# ##### step 3. annotate IVs for significant metabolites #####
# ##### ------------------------------------------------ #####
# cd ${DIR_SUMMARY}/withoutCorMatrix
# myDate="2019-12-04"
# mySigPval="5e-05"


# for mySNPset in p1_5e8_r2_0.01 p1_5e8_r2_0.1 p1_5e8_r2_0.2; do
	# # mySNPset="p1_5e8_r2_0.01"
	# echo ${mySNPset}
	
	# cat ${DIR_SUMMARY}/withoutCorMatrix/MRresults_sig_${mySigPval}_${mySNPset}_${myDate}.txt | awk '{print "harmonised_sumstats_"$3"_vs_"$2"_"$1".txt" }' >${DIR_SUMMARY}/withoutCorMatrix/${mySNPset}_sigPairs.txt
	# sed -i -e 's/_Raw//g'  ${DIR_SUMMARY}/withoutCorMatrix/${mySNPset}_sigPairs.txt
	# sed -i -e 's/_Outlier-corrected/_outlierCorrected/g'  ${DIR_SUMMARY}/withoutCorMatrix/${mySNPset}_sigPairs.txt

	# ## combine harmonized summary stats for all exposure-outcome pairs ##
	# for myFile in `cat ${DIR_SUMMARY}/withoutCorMatrix/${mySNPset}_sigPairs.txt `; do
		# echo ${myFile}
		# cat ${DIR_DATA}/ForMR/${mySNPset}/${myFile} >>${DIR_SUMMARY}/withoutCorMatrix/IVs_${mySNPset}.txt
	# done
	
	# ## remove duplicated headerline ##
	# head -1 ${DIR_SUMMARY}/withoutCorMatrix/IVs_${mySNPset}.txt >${DIR_SUMMARY}/withoutCorMatrix/harmonised_sumstats_for_MRsigPairs_${mySNPset}.txt
	# egrep -vw 'SNP' ${DIR_SUMMARY}/withoutCorMatrix/IVs_${mySNPset}.txt  >>${DIR_SUMMARY}/withoutCorMatrix/harmonised_sumstats_for_MRsigPairs_${mySNPset}.txt
	
	# ## subset annotation for all the IVs ##
	# cat ${DIR_ANNOVAR_OUT}/annovar_snplist_from_${mySNPset}.hg19_multianno.csv | awk -F',' '{print "chr"$1":"$2","$0}'  >${DIR_ANNOVAR_OUT}/annovar_snplist_from_${mySNPset}.hg19_multianno_tempt.csv
	
	# head -1 ${DIR_ANNOVAR_OUT}/annovar_snplist_from_${mySNPset}.hg19_multianno_tempt.csv | awk -F',' '{print $1, $2, $3, $7, $8, $10}' >${DIR_SUMMARY}/withoutCorMatrix/annovar_ivs_for_MRsigPairs_${mySNPset}.txt
	
	# for mySNP in `cat ${DIR_SUMMARY}/withoutCorMatrix/harmonised_sumstats_for_MRsigPairs_${mySNPset}.txt | awk -F' ' '{print $1}' `; do
		# echo ${mySNP}
		# egrep -w ${mySNP} ${DIR_ANNOVAR_OUT}/annovar_snplist_from_${mySNPset}.hg19_multianno_tempt.csv | awk -F',' '{print $1, $2, $3, $7, $8, $10}' >>${DIR_SUMMARY}/withoutCorMatrix/annovar_ivs_for_MRsigPairs_${mySNPset}.txt
	# done
	
	# ## combine harmonized summary stats and annotation into one file ##
	# paste -d'\t' ${DIR_SUMMARY}/withoutCorMatrix/annovar_ivs_for_MRsigPairs_${mySNPset}.txt ${DIR_SUMMARY}/withoutCorMatrix/harmonised_sumstats_for_MRsigPairs_${mySNPset}.txt  >${DIR_SUMMARY}/withoutCorMatrix/harmonised_sumstats_with_annotation_for_MRsigPairs_${mySNPset}.txt

	# rm ${DIR_SUMMARY}/withoutCorMatrix/IVs_${mySNPset}.txt
	# rm ${DIR_ANNOVAR_OUT}/annovar_snplist_from_${mySNPset}.hg19_multianno_tempt.csv
	# rm ${DIR_SUMMARY}/withoutCorMatrix/harmonised_sumstats_for_MRsigPairs_${mySNPset}.txt
	# rm ${DIR_SUMMARY}/withoutCorMatrix/annovar_ivs_for_MRsigPairs_${mySNPset}.txt
	
# done


##### THE END #####