##########################################################
##### Combine harmonized data as supplementary table #####
##########################################################
rm(list = ls())
## set up directories ##
# DIR = "working-directory"

myDIR = paste(DIR, "/olink_vs_stroke", sep = "")
dir_code = paste(myDIR, "/code/MR", sep = "")
dir_data = paste(myDIR, "/data", sep = "")
dir_results = paste(myDIR, "/results", sep = "")
dir_plots =  paste(myDIR, "/plots", sep = "")
dir_summary = paste(myDIR, "/summary", sep = "")

## set variables ##
DIR_SNPset = paste(dir_data, "/LD_clumped_snps", sep = "")
DIR_EXPOSURE = paste(dir_data, "/olink_pqtl", sep = "")
DIR_OUTCOME = paste(dir_data, "/stroke_gwas", sep = "")

#####--------------------#####
##### List of parameters #####
#####--------------------#####
## List SNPsets ##
List_SNPsets <- c("p1_5e8_r2_0.01", "p1_5e8_r2_0.1", "p1_5e8_r2_0.2")
mySNPset <- List_SNPsets[1]

## List exposures_all ##
List_Exposures <- read.delim(paste(myDIR,"/List_olinkProteins_withGeneInfo.txt", sep = ""), header = T, stringsAsFactors = F, sep = "\t")
head(List_Exposures)
dim(List_Exposures)

## List exposures_sig for stroke (15) ##
List_Exposures_sig <- read.delim(paste(myDIR,"/List_sigProteins_stroke_15.txt", sep = ""), header = T, stringsAsFactors = F, sep = "\t")
List_Exposures_sig <- List_Exposures[ match(List_Exposures_sig$OlinkPID, List_Exposures$OlinkPID), ]
head(List_Exposures_sig)
dim(List_Exposures_sig)

## List outcomes ##
List_Outcomes <- c("Stroke", "Ischemic-stroke", "Large-artery-stroke", "Cardioembolic-stroke", "Small-vessel-stroke")
List_Outcomes


#####--------------#####
##### combine data #####
#####--------------#####
colnanmes_harmonised_sumstats <- c("SNP", "chromosome",  "position", 
                                   "exposure", "effect_allele.exposure",  "other_allele.exposure", "eaf.exposure", 
                                   "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "Fstat",
                                   "outcome", "effect_allele.outcome", "other_allele.outcome",  "eaf.outcome", 
                                   "beta.outcome", "se.outcome", "pval.outcome")

## combine all exposures and all outcomes ##
for(k in 1: length(List_SNPsets)){

  mySNPset <- List_SNPsets[k]
  myExposure <- List_Exposures$OlinkPID[1]
  myOutcome <- List_Outcomes[1]
  mydata_all <- read.delim(paste(dir_data, "/ForMR/", mySNPset, "/harmonised_sumstats_", myExposure,"_vs_", myOutcome, ".txt", sep =  ""), header = T, sep = "\t")
  mydata_all <- mydata_all[, colnanmes_harmonised_sumstats]
  dim(mydata_all)

  ## For all proteins ##
  for(i in 1: dim(List_Exposures)[1]){
    myExposure <- List_Exposures$OlinkPID[i]
    
    for(j in 1: length(List_Outcomes)){
      myOutcome <- List_Outcomes[j]
      
      ## read my data ##
      myfile_raw = paste(dir_data, "/ForMR/", mySNPset, "/harmonised_sumstats_", myExposure,"_vs_", myOutcome, ".txt", sep =  "")
      myfile_oc = paste(dir_data, "/ForMR/", mySNPset, "/harmonised_sumstats_", myExposure,"_vs_", myOutcome, "_outlierCorrected.txt", sep =  "")
      
      if(file.exists(myfile_oc)){
        mydata <- read.delim(myfile_oc, header = T, sep = "\t")
        mydata_update <- mydata[, colnanmes_harmonised_sumstats]

      }else if(file.exists(myfile_raw)){
        mydata <- read.delim(myfile_raw, header = T, sep = "\t")
        mydata_update <- mydata[, colnanmes_harmonised_sumstats]
        
      }else{
        print(paste("There is no available IVs for ", myExposure, sep = ""))
      }
      
      ## combine data ##
      mydata_all <- rbind(mydata_all, mydata_update)
    }
  }
  ## save data into files for each SNPset from LDclumping ##
  head(mydata_all)
  dim(mydata_all)
  mydata_all$R2 <- (mydata_all$beta.exposure * ((2 * mydata_all$eaf.exposure * (1 - mydata_all$eaf.exposure)) ^ 0.5))^2
  
  
  ## if the allele was coded into "TRUE", change into "T" ##
  for(index in 1: dim(mydata_all)[1]){
    
    if(is.na(mydata_all$effect_allele.exposure[index])  | is.na(mydata_all$other_allele.exposure[index])){
      
    }else{
      if(mydata_all$effect_allele.exposure[index] == TRUE){ mydata_all$effect_allele.exposure[index] <- as.character("T") }
      if(mydata_all$other_allele.exposure[index] == TRUE){ mydata_all$other_allele.exposure[index] <- as.character("T") }
      if(mydata_all$effect_allele.outcome[index] == TRUE){ mydata_all$effect_allele.outcome[index] <- as.character("T") }
      if(mydata_all$other_allele.outcome[index] == TRUE){ mydata_all$other_allele.outcome[index] <- as.character("T") }
    }
    
  }

  ## remove duplicated lines when appropriate.
  mydata_all$INDEX = paste(mydata_all$SNP, mydata_all$exposure, mydata_all$outcome, sep = "_")
  mydata_all_update <- mydata_all[which(duplicated(mydata_all$INDEX) == FALSE), ]
  
  # # For stroke significant proteins only:
  # write.table(mydata_all_update, file = paste(dir_summary, "/harmonised_sumstats_for_StrokeProteins_", mySNPset, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

  ## For all proteins:
  write.table(mydata_all_update, file = paste(dir_summary, "/harmonised_sumstats_for_allProteins_", mySNPset, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

}


###################
##### THE END #####
###################
