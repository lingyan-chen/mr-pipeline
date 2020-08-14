#-------------------------------------------------------------------------------------------------------------------------------------------------
# 0. Load libraries
#-------------------------------------------------------------------------------------------------------------------------------------------------
.libPaths("/home/lc753/privatemodules/RLibs/RLibs-r-3.5.1")
## load libraries ##
library(TwoSampleMR)
library("plyr")
library(MendelianRandomization)
library(MRPRESSO)
library("psych")   #install.packages("psych")

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Read in command line arguments
#-------------------------------------------------------------------------------------------------------------------------------------------------
# ##### ---------------------------------------------- #####
# ##### STEP 3. MR analysis without correlation matrix #####
# ##### ---------------------------------------------- #####
# Rscript --slave --vanilla ${DIR_CODE}/MR_PIPELINE/MR_PIPELINE_olink_vs_stroke_2020-04-11.R  \
# ${myOlinkPID} \
# ${DIR} 

# echo ${myOlinkPID}
# echo ${DIR}

# Reading in the command line arguments
args<-commandArgs(trailingOnly=TRUE)

# Print input arguments
print(args[1])
print(args[2])
print(getwd())

## load parameters from bash ##
myExposure <- as.character(args[1])
DIR <-  as.character(args[2])


################################################################################
################################################################################
#####           RUN MR analysis - withouth correlation matrix              #####
################################################################################
################################################################################
## set working directory ##
# DIR="C:/Users/lc753/Dropbox/Cambridge/MR/olink_vs_stroke"
# DIR="/rds/project/jmmh2/rds-jmmh2-projects/blood_pressure_genetics/bioinformatics/mr/MR/olink_vs_stroke"
DIR_CODE=paste(DIR, "/code", sep = "")
DIR_DATA=paste(DIR, "/data", sep = "")

## create directories for data/plots/results ##
dir_data = DIR_DATA
dir_plots =  paste(DIR, "/plots", sep = "")
dir_cor_matrix = paste(DIR, "/ForMR/cor_matrix", sep = "")
dir_results_all = paste(DIR, "/results", sep = "")
dir_summary_all = paste(DIR, "/summary", sep = "")

dir.create(dir_data, showWarnings = FALSE)
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_cor_matrix, showWarnings = FALSE)
dir.create(dir_results_all, showWarnings = FALSE)
dir.create(dir_summary_all, showWarnings = FALSE)

## set variables ##
DIR_SNPset = paste(dir_data, "/LD_clumped_snps", sep = "")
DIR_EXPOSURE = paste(dir_data, "/olink_pqtl", sep = "")
DIR_OUTCOME = paste(dir_data, "/stroke_gwas", sep = "")

#####--------------------#####
##### List of parameters #####
#####--------------------#####
## List SNPsets ##
List_SNPsets <- c("p1_5e8_r2_0.01", "p1_5e8_r2_0.1", "p1_5e8_r2_0.2")
# mySNPset <- "p1_5e8_r2_0.1"

# ## List exposures ##
# myExposure <- "cvd3_TFPI___P10646"
# myExposure <- "neuro_TMPRSS5___Q9H3S3"
# myExposure <- "cvd3_EPHB4___P54760"
# myExposure <- "inf1_CD6___Q8WWJ7"
# myExposure <- "cvd2_GDF.2___Q9UK05"
# myExposure <- "cvd3_ICAM.2___P13598"
# myExposure <- "cvd2_MERTK___Q12866"

## List outcomes ##
# List_Outcomes <- c("pad", "vte")
List_Outcomes <- c("Stroke", "Ischemic-stroke", "Large-artery-stroke", "Cardioembolic-stroke", "Small-vessel-stroke")
# myOutcome <- "Cardioembolic-stroke"

## creat directories for correlation matrix and results accounted for correlation ##
dir_results = paste(dir_results_all, "/withoutCorMatrix", sep = "")
# dir_results = paste(dir_results, "/NoExonicSNP", sep = "")
# dir.create(dir_results, showWarnings = FALSE)
# dir_results = paste(dir_results, "/cisSNP", sep = "")
# dir.create(dir_results, showWarnings = FALSE)
dir_summary = paste(dir_summary_all, "/withoutCorMatrix", sep = "")
dir.create(dir_summary, showWarnings = FALSE)


## loop through all (SNPset) - exposure - outcome pairs to do univariable MR analysis ##
for(k in 1: length(List_SNPsets)){
  mySNPset <- List_SNPsets[k]
  print(mySNPset)

  ## creat directory for data/plots/results for my SNPset if not existed yet ##
  dir.create(paste(dir_data, "/ForMR/", mySNPset, sep = ""), showWarnings = FALSE)
  dir.create(paste(dir_plots, "/", mySNPset, sep = ""), showWarnings = FALSE)
  dir.create(paste(dir_results, "/", mySNPset, sep = ""), showWarnings = FALSE)


  for(j in 1: length(List_Outcomes)){
    myOutcome <- as.character( List_Outcomes[j])
    print(myOutcome)

    # myFile = paste(DIR_OUTCOME,"/", mySNPset, "/", myOutcome, "_transethnic_gwas.txt", sep = "")
    myFile = paste(DIR_OUTCOME,"/",mySNPset, "/",myExposure, "_vs_", myOutcome, sep = "")

    if(file.exists(myFile)){
      myOutcome_data <- read.delim(file = myFile, header = T, sep = "", stringsAsFactors = F)

      if(dim(myOutcome_data)[1] == 0){
        print(paste("No available outcome data for ", myOutcome, sep = ""))
      }else{

        # ##### step 1 - format data #####
        # ## load source R code for function "FUN1_FORMAT_DATA.R" ##
        # source(paste(DIR_CODE, "/MR/FUN1_FORMAT_DATA.R", sep = ""))
        # FORMAT_DATA(DIR_SNPset, DIR_EXPOSURE, DIR_OUTCOME, mySNPset, myExposure, myOutcome, dir_data)

  	  	##### step 2 - MR analysis without correlation matrix #####
  	  	## the difference betweeen funtion "FUN3_MR_ANALYSES_WITHOUT_COR.R" & "FUN3_MR_ANALYSES_WITHOUT_COR_2020-04-11.R" is that:
  	  	## the second function force the MR analysis to include the palindromic alleles.  Almost all the GWAS sumamry stats are now coded the allele in the + strand, so alleles can be compared directly.
  	  	
        ## the difference betweeen funtion "FUN3_MR_ANALYSES_WITHOUT_COR_2020-04-11.R" & "FUN3_MR_ANALYSES_WITHOUT_COR_2020-04-15.R" is that:
        ## the second function perform MR-contamination mixture model for IVs = 3.  
        ## If the contamination mixture model give multiple range (95%CI), keep the first one. [ IVs may have been clustered into different groups]

        source(paste(DIR_CODE, "/MR_conMix_pval.R", sep = ""))
  	  	source(paste(DIR_CODE, "/FUN3_MR_ANALYSES_WITHOUT_COR_2020-04-15.R", sep = ""))
  	  	MR_ANALYSES_WITHOUT_COR(mySNPset, myExposure, myOutcome, dir_data, dir_plots, dir_results)

      }

    }else{
      print(paste("No data for ", myOutcome, sep = ""))
    }


  } ## END of List_Outcomes


} ## END of List_SNPsets


##### THE END #####
