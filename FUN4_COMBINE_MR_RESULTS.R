#-------------------------------------------------------------------------------------------------------------------------------------------------
# 0. Load libraries
#-------------------------------------------------------------------------------------------------------------------------------------------------
.libPaths("path-to-R-libs")
## load libraries ##
library(TwoSampleMR)
library("plyr")
library(MendelianRandomization)
library(MRPRESSO)
library("psych")   #install.packages("psych")

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 1. Read in command line arguments
#-------------------------------------------------------------------------------------------------------------------------------------------------

# ##### ----------------------------------------------------- #####
# ##### STEP 4. combine MR results without correlation matrix #####
# ##### ----------------------------------------------------- #####
# Rscript --slave --vanilla ${DIR_CODE}/FUN4_COMBINE_MR_RESULTS.R  \
# ${DIR} \
# ${DIR_CODE} \
# ${DIR_DATA} \
# ${myDate}

# Reading in the command line arguments
args<-commandArgs(trailingOnly=TRUE)

# Print input arguments
print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(getwd())

## load parameters from bash ##
DIR <-  as.character(args[1])
DIR_CODE <- as.character(args[2])
DIR_DATA <- as.character(args[3])
myDate <- as.character(args[4])


########################################
#### Function4 - Combine MRresults #####
########################################
COMBINE_RESULTS <- function(mySNPset, DIR_RESULT, DIR_SUMMARY, DATE){

  file_list <- list.files(path = paste(DIR_RESULT,  "/", mySNPset, sep = ""), pattern = "MRresults_", full.names = FALSE, ignore.case = FALSE)
  myResult_ALL <- read.delim(file = paste(DIR_RESULT, "/", mySNPset, "/", file_list[1], sep = ""), header = T, sep = "\t")

  for (i in 2: length(file_list)){
    myResult <- read.delim(file = paste(DIR_RESULT,  "/", mySNPset, "/", file_list[i], sep = ""), header = T, sep = "\t", stringsAsFactors = F)
    myResult$GlobalTest_pval_PRESSO <- as.character(myResult$GlobalTest_pval_PRESSO)
    myResult_ALL <- rbind(myResult_ALL, myResult)
  }

  ## write out summary to files
  head(myResult_ALL)
  dim(myResult_ALL)
  write.table(myResult_ALL, file = paste(DIR_SUMMARY, "/MRresults_",mySNPset, "_", DATE, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

}


##########################################################################################################
##### DEFINE FUNTION: "CLEAN_MR_RESULTS" - keep one MR results (with outlier-corrected if available) #####
##########################################################################################################
CLEAN_MR_RESULTS <- function(dir_summary, myDate, List_SNPsets){

  for(k in 1:length(List_SNPsets)){
    mySNPset <- List_SNPsets[k]
    
    MRresults_all <- read.delim(file = paste(dir_summary, "/MRresults_", mySNPset, "_", myDate, ".txt", sep = ""), header = T, sep = "\t")
    dim(MRresults_all)
    
    ## keep unique summary stats for each exposure-outcome pair: either raw or outlier-corrected ##
    List_allExposures <- unique(MRresults_all$exposure)
    List_Outcomes <- unique(MRresults_all$outcome)
    length(List_allExposures)
    length(List_Outcomes)
    head(List_allExposures)
    
    MRresults_all_update <- matrix(ncol = ncol(MRresults_all))
    colnames(MRresults_all_update) <- colnames(MRresults_all)
    
    for(i in 1:length(List_allExposures)){
      myExposure <- as.character(List_allExposures[i])
      # print(myExposure)
      
      for(j in 1: length(List_Outcomes)){
        myOutcome <- as.character(List_Outcomes[j])
        nrow = dim(MRresults_all[MRresults_all$exposure == myExposure & MRresults_all$outcome == myOutcome, ])[1]
        nrow_raw = dim(MRresults_all[MRresults_all$exposure == myExposure & MRresults_all$outcome == myOutcome & MRresults_all$MRanalysis_SNPlist == "Raw", ])[1]
        
        ## Use if conditions to choose which line to use as the final results ##
        if(nrow == 1){
          myMRresults_update <- MRresults_all[MRresults_all$exposure == myExposure & MRresults_all$outcome == myOutcome & MRresults_all$MRanalysis_SNPlist == "Raw", ]
          MRresults_all_update <- rbind(MRresults_all_update, myMRresults_update)

        }else if(nrow_raw == 2){
          nsnps_raw = unique(MRresults_all[MRresults_all$exposure == myExposure & MRresults_all$outcome == myOutcome & MRresults_all$MRanalysis_SNPlist == "Raw" , "nsnp"])
          myMRresults_update <- MRresults_all[MRresults_all$exposure == myExposure & MRresults_all$outcome == myOutcome & MRresults_all$MRanalysis_SNPlist == "Raw" & MRresults_all$nsnp == which.max(nsnps_raw),  ]
          MRresults_all_update <- rbind(MRresults_all_update, myMRresults_update)

        }else{
          nsnps_raw = unique(MRresults_all[MRresults_all$exposure == myExposure & MRresults_all$outcome == myOutcome & MRresults_all$MRanalysis_SNPlist == "Raw" , "nsnp"])
          nsnps_oc = MRresults_all[MRresults_all$exposure == myExposure & MRresults_all$outcome == myOutcome & MRresults_all$MRanalysis_SNPlist == "Outlier-corrected" , "nsnp"]
          
          if(nsnps_raw == nsnps_oc){
            myMRresults_update <- MRresults_all[MRresults_all$exposure == myExposure & MRresults_all$outcome == myOutcome & MRresults_all$MRanalysis_SNPlist == "Raw", ]
            MRresults_all_update <- rbind(MRresults_all_update, myMRresults_update)
          }else{
            myMRresults_update <- MRresults_all[MRresults_all$exposure == myExposure & MRresults_all$outcome == myOutcome & MRresults_all$MRanalysis_SNPlist == "Outlier-corrected", ]
            MRresults_all_update <- rbind(MRresults_all_update, myMRresults_update)
          }
        }
        
      } ## END of loop through "List_Outcomes"
      
    } ## END of loop through "List_Exposures"
    
    dim(MRresults_all_update)
    head(MRresults_all_update)
    write.table(MRresults_all_update[-1, ], file = paste(dir_summary, "/MRresults_clean_", mySNPset, "_", myDate, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
 
   } ## END of loop through "List_SNPsets"

} ## END of function "CLEAN_MR_RESULTS"



################################################################################
################################################################################
#####                   RUN FUNCTION 4: COMBINE_RESULTS                    #####
################################################################################
################################################################################
## create directories for data/plots/results ##
dir_data = DIR_DATA
dir_plots =  paste(DIR, "/plots", sep = "")
dir_cor_matrix = paste(dir_data, "/ForMR/cor_matrix", sep = "")
dir_results_all = paste(DIR, "/results", sep = "")
dir_summary_all = paste(DIR, "/summary", sep = "")

dir.create(dir_data, showWarnings = FALSE)
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_cor_matrix, showWarnings = FALSE)
dir.create(dir_results_all, showWarnings = FALSE)
dir.create(dir_summary_all, showWarnings = FALSE)

#####--------------------#####
##### List of parameters #####
#####--------------------#####
## List SNPsets ##
List_SNPsets <- c("p1_5e8_r2_0.01", "p1_5e8_r2_0.1", "p1_5e8_r2_0.2")
List_SNPsets

## List outcomes ##
List_Outcomes <- c("Stroke", "Ischemic-stroke", "Large-artery-stroke", "Cardioembolic-stroke", "Small-vessel-stroke")
List_Outcomes

## creat directories for correlation matrix and results accounted for correlation ##
dir_results = paste(dir_results_all, "/withoutCorMatrix", sep = "")
dir.create(dir_results, showWarnings = FALSE)
dir_summary = paste(dir_summary_all, "/withoutCorMatrix", sep = "")
dir.create(dir_summary, showWarnings = FALSE)


##### step4 - combine all MRresults ######
for(k in 1: length(List_SNPsets)){
  mySNPset <- List_SNPsets[k]
  print(mySNPset)

  # print today's date
  # today <- Sys.Date()
  COMBINE_RESULTS(mySNPset, dir_results, dir_summary, myDate)

}

## clean MR results:
# today <- Sys.Date()
CLEAN_MR_RESULTS(dir_summary, myDate, List_SNPsets)


##### THE END #####
