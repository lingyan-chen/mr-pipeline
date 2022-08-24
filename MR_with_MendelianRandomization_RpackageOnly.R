##### ---------------------------------------------------------------------------------- #####
##### MR IVW with fixed-effect / random-effect with MendelianRandomization Rpackage ONLY #####
##### ---------------------------------------------------------------------------------- #####
rm(list = ls())

## load R pacakges 
.libPaths("~/RLibs/RLibs-r-3.5.1")
library(MendelianRandomization)

## List Exposures & Outcoems ##
List_Exposures <- c("cvd3_TFPI___P10646", 
                    "neuro_TMPRSS5___Q9H3S3", 
                    "inf1_CD40___P25942",
                    "inf1_CD6___Q8WWJ7",
                    "cvd2_MMP.12___P39900",
                    "cvd3_IL.6RA___P08887")

List_Outcomes <- c("Stroke",
                   "Ischemic-stroke", 
                   "Large-artery-stroke",
                   "Cardioembolic-stroke", 
                   "Small-vessel-stroke")

List_Pth <- c("p1_5e8_r2_0.01", "p1_5e8_r2_0.1", "p1_5e8_r2_0.2")

myExposure <- List_Exposures[6]
myOutcome <- List_Outcomes[4]
myPth <- List_Pth[2]

##### data directory #####
IV_method="LDclump"

DIR="~/MR/olink_vs_stroke"

##### save all data as the same format with the following columns #####
mydat_colnanmes<- c("SNPList", "Pth", "SNP", "chromosome",  "position",
                    "exposure", "effect_allele.exposure",  "other_allele.exposure", "eaf.exposure",
                    "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "Fstat",
                    "outcome", "effect_allele.outcome", "other_allele.outcome",  "eaf.outcome",
                    "beta.outcome", "se.outcome", "pval.outcome")

##### save all results as the same format with the following columns #####
res_colnames <- c("SNPList", "Pth", "exposure", "outcome", "nsnp",
                  "beta_IVW_fixedEffect", "se_IVW_fixedEffect", "pval_IVW_fixedEffect", 
                  "beta_IVW_randomEffect", "se_IVW_randomEffect", "pval_IVW_randomEffect",
                  "Q_IVW", "Q_df_IVW", "Q_pval_IVW",
                  "beta_egger", "se_egger", "pval_egger", 
                  "intercept_egger", "intercept_se_egger", "intercept_pval_egger",
                  "beta_conMix", "CILower_conMix", "CIUpper_conMix", "pval_conMix")

##### loop through all parameters #####
for (myPth in List_Pth){
  print(myPth)
  
  DIR_DATA_FOR_MR = paste(DIR, "/data/ForMR/", myPth, sep = "")
  DIR_PLOT = paste(DIR, "/plots/", myPth, sep = "")
  DIR_RESULT = paste(DIR, "/results/withoutCorMatrix/", myPth, sep = "")
  DIR_SUMMARY = paste(DIR, "/summary/withoutCorMatrix", sep = "")
  
  for(myExposure in List_Exposures){
    print(myExposure)
    
    for(myOutcome in List_Outcomes){
      print(myOutcome)
      
      ## load harmonized data for exposure and outcome ##
      mydat_original <- read.delim(file = paste(DIR_DATA_FOR_MR, "/harmonised_sumstats_", myExposure,"_vs_", myOutcome, ".txt", sep =  ""),  header = T, sep = "\t")
      dim(mydat_original)
      
      mydat <- cbind(IV_method, myPth, mydat_original)
      colnames(mydat)[1] <- "SNPList"
      colnames(mydat)[2] <- "Pth"
      
      ## sometimes there are some problems when reading allele "T" in a matrix:
      for(k in 1: dim(mydat)[1]){
        if(mydat$effect_allele.exposure[k] == TRUE){ mydat$effect_allele.exposure[k] <- as.character("T") }
        if(mydat$other_allele.exposure[k] == TRUE){ mydat$other_allele.exposure[k] <- as.character("T") }
        if(mydat$effect_allele.outcome[k] == TRUE){ mydat$effect_allele.outcome[k] <- as.character("T") }
        if(mydat$other_allele.outcome[k] == TRUE){ mydat$other_allele.outcome[k] <- as.character("T") }
        if(mydat$palindromic[k] == TRUE) { mydat$mr_keep[k] <- as.logical("TRUE") }
      }
      
      mydat_update <- mydat[, mydat_colnanmes]
      
      ## combine all results ##
      if( exists("mydat_ALL") == FALSE ){
        mydat_ALL <- mydat_update
      }else{
        mydat_ALL <- rbind(mydat_ALL, mydat_update)
      }
      
      
      ##### ----------------------------------- #####
      ##### Use MendelianRandomization Rpackage #####
      ##### ----------------------------------- #####
      if(dim(mydat)[1] == 1){
        ## MR wald's ratio ##
        res_ivw <- mr_ivw(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome))
        NSNPs <- res_ivw$SNPs
        beta_ivw <- res_ivw$Estimate
        se_ivw <- res_ivw$StdError
        pval_ivw <- res_ivw$Pvalue
        
        myResults <- cbind(IV_method, myPth, myExposure, myOutcome, NSNPs,
                           beta_ivw, se_ivw, pval_ivw, 
                           NA, NA, NA, NA, NA, NA, 
                           NA, NA, NA, NA, NA, NA, 
                           NA, NA, NA, NA)
        colnames(myResults) <- res_colnames
        
      }else if(dim(mydat)[1] == 2){
        ## MR-IVW with fixed-effect ##
        res_ivw_fixed <- mr_ivw(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome), model = "fixed")
        NSNPs <- res_ivw_fixed$SNPs
        beta_ivw_fixed <- res_ivw_fixed$Estimate
        se_ivw_fixed <- res_ivw_fixed$StdError
        pval_ivw_fixed <- res_ivw_fixed$Pvalue
        
        ## MR-IVW with random-effect ##
        res_ivw_random <- mr_ivw(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome), model = "random")
        beta_ivw_random <- res_ivw_random$Estimate
        se_ivw_random <- res_ivw_random$StdError
        pval_ivw_random <- res_ivw_random$Pvalue
        
        ## Heteogeneity test: IVW  ##
        Q_ivw <- res_ivw_fixed$Heter.Stat[1]
        Q_df_ivw <- NSNPs - 1
        Q_pval_ivw <- res_ivw_fixed$Heter.Stat[2]
        
        ## combine all results ##
        myResults <- cbind(IV_method, myPth, myExposure, myOutcome, NSNPs,
                           beta_ivw_fixed, se_ivw_fixed, pval_ivw_fixed, 
                           beta_ivw_random, se_ivw_random, pval_ivw_random,
                           Q_ivw, Q_df_ivw, Q_pval_ivw,
                           NA, NA, NA, NA, NA, NA, 
                           NA, NA, NA, NA)
        colnames(myResults) <- res_colnames
        
      }else if(dim(mydat)[1] >= 3){
        ## MR-IVW with fixed-effect ##
        res_ivw_fixed <- mr_ivw(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome), model = "fixed")
        NSNPs <- res_ivw_fixed$SNPs
        beta_ivw_fixed <- res_ivw_fixed$Estimate
        se_ivw_fixed <- res_ivw_fixed$StdError
        pval_ivw_fixed <- res_ivw_fixed$Pvalue
        
        ## MR-IVW with random-effect ##
        res_ivw_random <- mr_ivw(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome), model = "random")
        beta_ivw_random <- res_ivw_random$Estimate
        se_ivw_random <- res_ivw_random$StdError
        pval_ivw_random <- res_ivw_random$Pvalue
        
        ## Heteogeneity test: IVW  ##
        Q_ivw <- res_ivw_fixed$Heter.Stat[1]
        Q_df_ivw <- NSNPs - 1
        Q_pval_ivw <- res_ivw_fixed$Heter.Stat[2]
        
        ## MR_Egger ##
        res_egger <- mr_egger(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome))
        beta_egger <- res_egger$Estimate
        se_egger <- res_egger$StdError.Est
        pval_egger <- res_egger$Pvalue.Est
        
        inter_egger <- res_egger$Intercept
        inter_se_egger <- res_egger$StdError.Int 
        inter_pval_egger <- res_egger$Pvalue.Int
        
        ## MR contamination mixture model ##
        # res_conmix <- mr_conmix(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome),  psi = 3, CIMin = NA, CIMax = NA, CIStep = 0.0001)
        res_conmix <- mr_conmix(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome),  psi = 3, CIMin = -1, CIMax = 1, CIStep = 0.0001)
        beta_conmix <- res_conmix$Estimate
        CIlower_conmix <- res_conmix$CILower
        CIupper_conmix <- res_conmix$CIUpper
        pval_conmix <- res_conmix$Pvalue
        
        
        ## combine all results ##
        myResults <- cbind(IV_method, myPth, myExposure, myOutcome, NSNPs,
                           beta_ivw_fixed, se_ivw_fixed, pval_ivw_fixed, 
                           beta_ivw_random, se_ivw_random, pval_ivw_random,
                           Q_ivw, Q_df_ivw, Q_pval_ivw,
                           beta_egger, se_egger, pval_egger, inter_egger, inter_se_egger, inter_pval_egger,
                           beta_conmix, CIlower_conmix, CIupper_conmix, pval_conmix)
        colnames(myResults) <- res_colnames
        
        ## save MR plots ##
        pdf(file=paste(DIR_PLOT, "/",myExposure, "_vs_", myOutcome, "_scatterlot.pdf",sep = ""), width=7, height=7)    
        mr_plot(mr_allmethods(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome), method="all", iterations = 100))
        dev.off()
        
      }
      
      ## combine all results ##
      if( exists("Results_ALL") == FALSE ){
        Results_ALL <- myResults
      }else{
        Results_ALL <- rbind(Results_ALL, myResults)
      }
      
    } ## END of "List_Outcomes"
    
  } ## END of "List_Exposures"
  
} ## END of "List_Pth"

dim(mydat_ALL)
head(mydat_ALL)
write.table(mydat_ALL, file = paste(DIR_SUMMARY, "/harmonised_sumstats_for_StrokeProteins_12July.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)

dim(Results_ALL)
head(Results_ALL)
write.table(Results_ALL, file = paste(DIR_SUMMARY, "/MRresults_with_MendelianRandomization_Rpackage_for_StrokeProteins_12July.txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)


###################
##### THE END #####
###################
