#######################################################################
#### Function 1 - format exposure and outcome data for MR analyses #####
#######################################################################
FORMAT_DATA <- function(DIR_SNPset, DIR_EXPOSURE, DIR_OUTCOME, mySNPset, myExposure, myOutcome, DIR_DATA){
  
  # DIR_DATA = dir_data
  
  ## load libraries for format data ##
  # .libPaths("/home/lc753/privatemodules/RLibs/RLibs-r-3.5.1")
  ## load libraries ##
  library(TwoSampleMR)
  library("plyr")
  library(MendelianRandomization)
  library(MRPRESSO)
  library("psych")   #install.packages("psych")
  

  ##### 1. FORMAT EXPOSURE DATA #####
  myExposure_data <- read.delim(file = paste(DIR_EXPOSURE,"/", mySNPset, "/", myExposure, sep = ""), header = T, sep = "")
  # myExposure_data$SNPID_NEW <- paste("chr", myExposure_data$chromosome, ":", myExposure_data$position, sep = "")
  
  ## sometimes there are some problems when reading allele "T" in a matrix:
  for(k in 1: dim(myExposure_data)[1]){
    if(myExposure_data$alleleA[k] == TRUE){ myExposure_data$alleleA[k] <- as.character("T") }
    if(myExposure_data$alleleB[k] == TRUE){ myExposure_data$alleleB[k] <- as.character("T") }
  }
  
  ## calculate alleleB (effect allele) frequency
  myExposure_data$FreqB <- (myExposure_data$all_BB * 2 + myExposure_data$all_AB)/(myExposure_data$all_AA + myExposure_data$all_AB + myExposure_data$all_BB)/2
  dim(myExposure_data)
  head(myExposure_data)
  
  ## format exposure data for MRbase TwoSampleMR analysis: ##
  exposure_data<-format_data(myExposure_data,
                             type = "exposure",
                             header = TRUE,
                             snp_col = "rsid",
                             beta_col = "frequentist_add_beta_1",
                             se_col = "frequentist_add_se_1",
                             eaf_col = "FreqB",
                             effect_allele_col = "alleleB",
                             other_allele_col = "alleleA",
                             pval_col = "frequentist_add_pvalue",
                             samplesize_col = "all_total")
  exposure_data$exposure <- myExposure
  exposure_data$units.exposure <- "SD"
  # head(exposure_data)
  # dim(exposure_data)
  
  ## Save exposure IVs info
  myExposure_SNPinfo <-  myExposure_data[match(exposure_data$SNP, myExposure_data$rsid), c(3:4)]
  exposure_data <- cbind(exposure_data, myExposure_SNPinfo)
  # head(exposure_data)
  # write.table(exposure_data, file = paste(DIR_DATA_FOR_MR, "/", mySNPset, "/exposure_sumstats_", myExposure, "_for_", myOutcome, ".txt", sep =  ""), col.names = T, row.names = F, quote = F, sep = "\t")
  
  
  ##### 2. FORMAT OUTCOME DATA #####
  # myOutcome_data <- read.delim(file = paste(DIR_OUTCOME,"/",mySNPset, "/",myOutcome,  "_transethnic_gwas.txt", sep = ""), header = T, sep = "", stringsAsFactors = F)
  
  ## subset outcome data for myExposure ##
  # outcome_data <- format_data(myOutcome_data,
  #                             type = "outcome",
  #                             header = TRUE,
  #                             snps = exposure_data$SNP,
  #                             snp_col = "snp",
  #                             beta_col = "beta",
  #                             se_col = "se",
  #                             effect_allele_col = "effect_allele",
  #                             other_allele_col = "other_allele",
  #                             eaf_col = "eaf",
  #                             pval_col = "pval")
  
  ## read outcome data ##
  myOutcome_data <- read.delim(file = paste(DIR_OUTCOME,"/",mySNPset, "/",myExposure, "_vs_", myOutcome, sep = ""), header = T, sep = "", stringsAsFactors = F)
  
  ## subset outcome data for myExposure ##
  outcome_data <- format_data(myOutcome_data,
                              type = "outcome",
                              header = TRUE,
                              snps = exposure_data$SNP,
                              snp_col = "rsid",
                              beta_col = "beta",
                              se_col = "se",
                              effect_allele_col = "a1",
                              other_allele_col = "a2",
                              eaf_col = "eaf",
                              pval_col = "p")
  outcome_data$outcome <- myOutcome
  outcome_data$units.outcome <- "log odds"
  # head(outcome_data)
  
  ##### 3. HARMONISE EXPOSURE AND OUTCOME DATA #####
  mydat <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data)
  mydat <- mydat[is.na(mydat$beta.outcome) == F, ]           ## remove na column in exposure and outcome 
  mydat <- mydat[is.na(mydat$beta.exposure) == F, ]
  dim(mydat)
  
  ##### 4. Calulate power of Instrumental variants #####
  ## (1) Quantify the strength of the selected instruments using F-statistic ##
  mydat$Fstat <- (mydat$beta.exposure / mydat$se.exposure)^2
  
  ## (2) MR-steiger filtering to detect the correlation of SNPs with exposure and outcom -> remove SNPs whose effects are larger in outcome than in exposure ##
  # mydat <- steiger_filtering(mydat)
  # mydat$steiger_dir
  # dim(mydat)
  # head(mydat)
  
  ## save harmonized summary statst ##
  write.table(mydat, file = paste(DIR_DATA,"/ForMR/", mySNPset,  "/harmonised_sumstats_", myExposure,"_vs_", myOutcome, ".txt", sep =  ""), col.names = T, row.names = F, quote = F, sep = "\t")

  ## save harmonized summary stats without exonic variant ##
  ExonicSNP = "rs7110738"
  mydat_update <- mydat[-which(mydat$SNP == ExonicSNP), ]
  write.table(mydat_update, file = paste(DIR_DATA,"/ForMR/", mySNPset,  "/harmonised_sumstats_", myExposure,"_vs_", myOutcome, "_noExonicSNP.txt", sep =  ""), col.names = T, row.names = F, quote = F, sep = "\t")
  
  
} ## End of function: FORMAT_DATA(DIR_SNPset, DIR_EXPOSURE, DIR_OUTCOME, mySNPset, myExposure, myOutcome, DIR_DATA)


##### THE END #####