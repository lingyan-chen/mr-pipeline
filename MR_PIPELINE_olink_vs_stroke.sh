#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=mr
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH -p skylake-himem
#SBATCH --array=2-369
#SBATCH --output=slurm/mr_%A_%a.out
#SBATCH --error=slurm/mr_%A_%a.err
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

## set working directory ##
DIR="~/MR/olink_vs_stroke"
DIR_CODE="${DIR}/code"
DIR_DATA="${DIR}/data"

## GWAS summary stats ##
DIR_IN_proteinGWAS="~/pQTL/data/olink"
DIR_OUT_proteinGWAS="~/MR/olink_vs_stroke/data/olink_summary_stats"
DIR_IN_strokeGWAS="~/gwas/megastroke"
DIR_OUT_strokeGWAS="~/MR/olink_vs_stroke/data/megastroke_summary_stats"

## plink bfile for genome-wide from BP's clean INTERVAL data ##
mybfile="~/interval_subset_olink/genotype_files/unrelated_4994_pihat_0.1875_autosomal_imputed_info_0.4_phwe_1e-4_filtered/output/plink_format/interval.imputed.olink.all_chrs.locuszoom_chr_pos_varids"

## LD clumped SNPs list based on pQTLs p-value ##
DIR_OUT_CLUMPED="${DIR_DATA}/LD_clumped_snps"

## olink pQTL for LD clumped SNPs ##
DIR_OUT_FOR_MR_proteinGWAS="${DIR_DATA}/olink_pqtl"

## stroke GWAS for LD clumped SNPs ##
DIR_OUT_FOR_MR_strokeGWAS="${DIR_DATA}/stroke_gwas"

## harmonized summary stats for LD clumped SNPs ##
DIR_OUT_FOR_MR_harmonised_sumstats="${DIR_DATA}/ForMR"

## correlation matrix for MR analysis ##
DIR_OUT_FOR_MR_corMatrix="${DIR_DATA}/cor_matrix"

## creat directory for output if not exist ##
if [ ! -d ${DIR_DATA} ]; then mkdir -p ${DIR_DATA}; fi
if [ ! -d ${DIR_OUT_proteinGWAS} ]; then  mkdir -p ${DIR_OUT_proteinGWAS}; fi
if [ ! -d ${DIR_OUT_strokeGWAS} ]; then  mkdir -p ${DIR_OUT_strokeGWAS}; fi
if [ ! -d ${DIR_OUT_CLUMPED} ]; then mkdir -p ${DIR_OUT_CLUMPED}; fi
if [ ! -d ${DIR_OUT_FOR_MR_proteinGWAS} ]; then mkdir -p ${DIR_OUT_FOR_MR_proteinGWAS}; fi
if [ ! -d ${DIR_OUT_FOR_MR_strokeGWAS} ]; then mkdir -p ${DIR_OUT_FOR_MR_strokeGWAS}; fi
if [ ! -d ${DIR_OUT_FOR_MR_harmonised_sumstats} ]; then mkdir -p ${DIR_OUT_FOR_MR_harmonised_sumstats}; fi
if [ ! -d ${DIR_OUT_FOR_MR_corMatrix} ]; then mkdir -p ${DIR_OUT_FOR_MR_corMatrix}; fi


########################
#####  MR analysis #####
########################
cd ${DIR_CODE}

## target Protein list ##
targetProteinList="${DIR}/List_targetProteins.txt"

##### read my exposure and outcome #####
myPanel=$(awk -v protein_id_line=${SLURM_ARRAY_TASK_ID} 'NR==protein_id_line{print $2}' ${targetProteinList} )
myProteinID=$(awk -v protein_id_line=${SLURM_ARRAY_TASK_ID} 'NR==protein_id_line{print $3}' ${targetProteinList})
myOlinkPID="${myPanel}_${myProteinID}"

echo ${DIR}
echo ${myPanel}
echo ${myProteinID}
echo ${myOlinkPID}

##### --------------------------------------#####
##### STEP 1. format data for Two Sample MR #####
##### --------------------------------------#####
Rscript --slave --vanilla ${DIR_CODE}/MR_PIPELINE/MR_PIPELINE_FUN1_FORMAT_DATA.R \
${myOlinkPID} \
${DIR} \
${DIR_CODE} \
${DIR_DATA}

##### ---------------------------------------------- #####
##### STEP 2. MR analysis without correlation matrix #####
##### ---------------------------------------------- #####
Rscript --slave --vanilla ${DIR_CODE}/MR_PIPELINE_olink_vs_stroke_2020-04-15.R  \
${myOlinkPID} \
${DIR}


##### ----------------------------------------------------- #####
##### STEP 4. combine MR results without correlation matrix #####
##### ----------------------------------------------------- #####
Rscript --slave --vanilla ${DIR_CODE}/FUN4_COMBINE_MR_RESULTS.R  \
${DIR} \
${DIR_CODE} \
${DIR_DATA} \
${myDate}


# ##### ---------------------------------------------------- #####
# ##### STEP 2. calculate correlation matrix for target SNPs #####
# ##### ---------------------------------------------------- #####
# # sh ${DIR_CODE}/MR_PIPELINE/MR_PIPELINE_FUN2_COR_MATRIX.sh

# ##### generate correlation matrix for LD clumping SNPsets at different R2 threshold #####
# # for R2 in 0.01 0.1 0.2 0.3 0.4 0.5; do
# for R2 in 0.3 0.4 0.5; do
# 	echo ${R2}
# 	DIR_OUT_CLUMPED_R2="${DIR_OUT_CLUMPED}/p1_5e8_r2_${R2}"
# 	DIR_OUT_FOR_MR_proteinGWAS_R2="${DIR_OUT_FOR_MR_proteinGWAS}/p1_5e8_r2_${R2}"
# 	DIR_OUT_FOR_MR_strokeGWAS_R2="${DIR_OUT_FOR_MR_strokeGWAS}/p1_5e8_r2_${R2}"
# 	DIR_OUT_FOR_MR_harmonised_sumstats_R2="${DIR_OUT_FOR_MR_harmonised_sumstats}/p1_5e8_r2_${R2}"
# 	DIR_OUT_FOR_MR_corMatrix_R2="${DIR_OUT_FOR_MR_corMatrix}/p1_5e8_r2_${R2}"

# 	if [ ! -d ${DIR_OUT_CLUMPED_R2} ]; then
# 		mkdir -p ${DIR_OUT_CLUMPED_R2}
# 	fi

# 	if [ ! -d ${DIR_OUT_FOR_MR_proteinGWAS_R2} ]; then
# 		mkdir -p ${DIR_OUT_FOR_MR_proteinGWAS_R2}
# 	fi

# 	if [ ! -d ${DIR_OUT_FOR_MR_strokeGWAS_R2} ]; then
# 		mkdir -p ${DIR_OUT_FOR_MR_strokeGWAS_R2}
# 	fi

# 	if [ ! -d ${DIR_OUT_FOR_MR_harmonised_sumstats_R2} ]; then
# 		mkdir -p ${DIR_OUT_FOR_MR_harmonised_sumstats_R2}
# 	fi

# 	if [ ! -d ${DIR_OUT_FOR_MR_corMatrix_R2} ]; then
# 		mkdir -p ${DIR_OUT_FOR_MR_corMatrix_R2}
# 	fi

# 	##### --------------------------------------------------------------------------------------------- #####
# 	##### subset list of target IVs with effect allele for extracting dosage force on the effect allele #####
# 	##### --------------------------------------------------------------------------------------------- #####
# 	cat ${DIR_OUT_FOR_MR_harmonised_sumstats_R2}/exposure_sumstats_${myOlinkPID}_for_*.txt  >>${DIR_OUT_FOR_MR_harmonised_sumstats_R2}/exposure_sumstats_${myOlinkPID}_all.txt
# 	cat ${DIR_OUT_FOR_MR_harmonised_sumstats_R2}/exposure_sumstats_${myOlinkPID}_all.txt | awk '{print $8, $2}' | sort | uniq > ${DIR_OUT_FOR_MR_corMatrix_R2}/snplist_for_${myOlinkPID}.txt
# 	rm "${DIR_OUT_FOR_MR_harmonised_sumstats_R2}/exposure_sumstats_${myOlinkPID}_all.txt"
# 	echo "subset target IVs for ${myOlinkPID} with effect allele for PLINK done."


# 	##### ------------------------------------------------------------------------------------- #####
# 	##### subset genotypes of IVs for all proteins from INTERVAL reference genome file (N=4994) #####
# 	##### ------------------------------------------------------------------------------------- #####
# 	## step 1. subset genotype data from bfile and recode into .ped + .map file ##
# 	echo "Use plink to subset genotype data for target variants"
# 	plink \
# 	--bfile  ${mybfile} \
# 	--extract ${DIR_OUT_FOR_MR_corMatrix_R2}/snplist_for_${myOlinkPID}.txt \
# 	--recode --out ${DIR_OUT_FOR_MR_corMatrix_R2}/snplist_for_${myOlinkPID}

# 	## step 2. recode .ped + .map file to .raw (dosage) file and specify the effect allele (exposure) as reference allele (A1 in PLINK) ##
# 	## dosage are counted on the effect allele in an additive model.
# 	echo "recode .ped + .map file to .raw (dosage) file and specify the effect allele (exposure) as reference allele"
# 	plink \
# 	--file ${DIR_OUT_FOR_MR_corMatrix_R2}/snplist_for_${myOlinkPID} \
# 	--reference-allele ${DIR_OUT_FOR_MR_corMatrix_R2}/snplist_for_${myOlinkPID}.txt \
# 	--recode A --out ${DIR_OUT_FOR_MR_corMatrix_R2}/snplist_for_${myOlinkPID}

# 	## calculate correlation matrix for target snplist based on the specifiy dosage allele ##
# 	echo "calculate correlation matrix for target snplist based on the specifiy dosage allele"
# 	plink \
# 	--file ${DIR_OUT_FOR_MR_corMatrix_R2}/snplist_for_${myOlinkPID} \
# 	--reference-allele ${DIR_OUT_FOR_MR_corMatrix_R2}/snplist_for_${myOlinkPID}.txt \
# 	--r --matrix --out ${DIR_OUT_FOR_MR_corMatrix_R2}/snplist_for_${myOlinkPID}

# done


# ##### ------------------------------------------- #####
# ##### STEP 3. MR analysis with correlation matrix #####
# ##### ------------------------------------------- #####
# Rscript --slave --vanilla ${DIR_CODE}/MR_PIPELINE/MR_PIPELINE_FUN3_MR_ANALYSES_WITH_COR.R  \
# ${myOlinkPID} \
# ${DIR} \
# ${DIR_CODE} \
# ${DIR_DATA}

# echo ${myOlinkPID}
# echo ${DIR}
# echo ${DIR_CODE}
# echo ${DIR_DATA}


# ##### -------------------------------------------------- #####
# ##### STEP 4. combine MR results with correlation matrix #####
# ##### -------------------------------------------------- #####
# Rscript --slave --vanilla ${DIR_CODE}/MR_PIPELINE/MR_PIPELINE_FUN4_COMBINE_MR_RESULTS_WITH_COR.R  \
# ${myOlinkPID} \
# ${DIR} \
# ${DIR_CODE} \
# ${DIR_DATA}


##### THE END #####
