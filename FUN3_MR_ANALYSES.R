#-------------------------------------------------------------------------------------------------------------------------------------------------
# 0. Load libraries
#-------------------------------------------------------------------------------------------------------------------------------------------------
.libPaths("path-to-R-Libs")

# ## load libraries for MR analysis ##
library("varhandle")   
library(MendelianRandomization)
library(MRPRESSO)
library(TwoSampleMR)
library("plyr")
library("psych") 
library("Hmisc")
library("corrplot")


####################################################
#####      DEFINE FUNCTION: "MR_ANALYSIS"      #####
####################################################
MR_ANALYSIS <- function(myExposure, myOutcome, DIR_DATA_FOR_MR, DIR_RESULT, DIR_PLOT){
  
  ## save all results as the same format with the following columns ##
  res_colnames <- c("MRanalysis_SNPlist", "outcome", "exposure", "nsnp",
                    "beta_IVW", "se_IVW", "pval_IVW", "Q_IVW", "Q_df_IVW", "Q_pval_IVW",
                    "beta_egger", "se_egger", "pval_egger", "intercept_egger", "intercept_se_egger", "intercept_pval_egger",
                    "beta_PRESSO", "se_PRESSO", "Tstat_PRESSO", "pval_PRESSO", "GlobalTest_pval_PRESSO", "beta_conMix", "CILower_conMix", "CIUpper_conMix", "pval_conMix")
  
  ##### ----------------------------------------------------------------- #####
  ##### load data for MR analysis: summary stats for exposure and outcome #####
  ##### ----------------------------------------------------------------- #####
  ## load harmonized data ##
  ## check if file exist:
  FileName <- paste(DIR_DATA_FOR_MR, "/harmonised_sumstats_", myExposure,"_vs_", myOutcome, ".txt", sep =  "")
  
  if(!file.exists(FileName)){
    print("exit FUN2_MR_ANALYSES.R")
    
  }else{
    mydat <- read.delim(file = paste(DIR_DATA_FOR_MR, "/harmonised_sumstats_", myExposure,"_vs_", myOutcome, ".txt", sep =  ""),  header = T, sep = "\t")
    dim(mydat)
    
    ##### ----------------------------------------------------------------------------------------------- #####
    ##### Performed MR analysis used IVW, Egger, PRESSO, and Contamination mixture model when IVs allowed #####
    ##### ----------------------------------------------------------------------------------------------- #####
    if(dim(mydat)[1] == 0){
      print(paste("NOT AVAILABLE IVs for ", myExposure, "_vs_", myOutcome, sep = ""))
    }else{
      
      ## sometimes there are some problems when reading allele "T" in a matrix:
      for(k in 1: dim(mydat)[1]){
        if(mydat$effect_allele.exposure[k] == TRUE){ mydat$effect_allele.exposure[k] <- as.character("T") }
        if(mydat$other_allele.exposure[k] == TRUE){ mydat$other_allele.exposure[k] <- as.character("T") }
        if(mydat$effect_allele.outcome[k] == TRUE){ mydat$effect_allele.outcome[k] <- as.character("T") }
        if(mydat$other_allele.outcome[k] == TRUE){ mydat$other_allele.outcome[k] <- as.character("T") }
        if(mydat$palindromic[k] == TRUE) { mydat$mr_keep[k] <- as.logical("TRUE") }
      }
      
      
      if(dim(mydat)[1] == 1){
        ## MR wald's ratio ##
        res_ivw <-  MendelianRandomization::mr_ivw(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome))
        NSNPs <- res_ivw$SNPs
        beta_ivw <- res_ivw$Estimate
        se_ivw <- res_ivw$StdError
        pval_ivw <- res_ivw$Pvalue
        
        res_all <- cbind("Raw", myOutcome, myExposure, NSNPs, beta_ivw, se_ivw, pval_ivw, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        colnames(res_all) <- res_colnames
        write.table(res_all, file = paste(DIR_RESULT, "/MRresults_", myExposure, "_vs_", myOutcome, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
        
        
      }else if(dim(mydat)[1] == 2){
        ## MR-IVW with fixed-effect ##
        res_ivw <- MendelianRandomization::mr_ivw(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome))
        NSNPs <- res_ivw$SNPs
        beta_ivw <- res_ivw$Estimate
        se_ivw <- res_ivw$StdError
        pval_ivw <- res_ivw$Pvalue
        
        ## Heteogeneity test: IVW  ##
        Q_ivw <- res_ivw$Heter.Stat[1]
        Q_df_ivw <- NSNPs - 1
        Q_pval_ivw <- res_ivw$Heter.Stat[2]
        
        ## combine all results ##
        res_all <- cbind("Raw", myOutcome, myExposure, NSNPs,
                         beta_ivw, se_ivw, pval_ivw, Q_ivw, Q_df_ivw, Q_pval_ivw,
                         NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        colnames(res_all) <- res_colnames
        write.table(res_all, file = paste(DIR_RESULT, "/MRresults_",myExposure, "_vs_", myOutcome, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
        
        ## scatter plot using MendelianRandomization package ##
        pdf(file=paste(DIR_PLOT, "/",myExposure, "_vs_", myOutcome, "_scatterPlot_MendelianRandomization.pdf",sep = ""), width=7, height=7)
        myplot <- MendelianRandomization::mr_plot(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome), line="ivw", orientate = TRUE, interactive = FALSE)
        print(myplot)
        dev.off()
        
        ## scatter plot  using TwoSampleMR package #
        pdf(file=paste(DIR_PLOT, "/", myExposure, "_vs_", myOutcome, "_scatterPlot_TwoSampleMR.pdf",sep = ""), width=7, height=7)
        res <- TwoSampleMR::mr(mydat, method="mr_ivw")
        p1 <- TwoSampleMR::mr_scatter_plot(res, mydat)
        print(p1)
        dev.off()
        
        ## forest plot for each IV using TwoSampleMR package. #
        pdf(file=paste(DIR_PLOT, "/", myExposure, "_vs_", myOutcome, "_forestPlot_TwoSampleMR.pdf",sep = ""), width=7, height=7)       ## Forest plot
        res_single <- TwoSampleMR::mr_singlesnp(mydat, single_method="mr_meta_fixed", all_method=c("mr_ivw"))
        p2 <- TwoSampleMR::mr_forest_plot(res_single)
        print(p2)
        dev.off()
        
        
      }else if(dim(mydat)[1] == 3){
        ## MR-IVW with fixed-effect ##
        res_ivw <- MendelianRandomization::mr_ivw(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome))
        NSNPs <- res_ivw$SNPs
        beta_ivw <- res_ivw$Estimate
        se_ivw <- res_ivw$StdError
        pval_ivw <- res_ivw$Pvalue
        
        ## Heteogeneity test: IVW  ##
        Q_ivw <- res_ivw$Heter.Stat[1]
        Q_df_ivw <- NSNPs - 1
        Q_pval_ivw <- res_ivw$Heter.Stat[2]
        
        ## MR_Egger ##
        res_egger <- MendelianRandomization::mr_egger(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome))
        beta_egger <- res_egger$Estimate
        se_egger <- res_egger$StdError.Est
        pval_egger <- res_egger$Pvalue.Est
        
        inter_egger <- res_egger$Intercept
        inter_se_egger <- res_egger$StdError.Int 
        inter_pval_egger <- res_egger$Pvalue.Int
        
        ## MR contamination mixture model ##
        # res_conmix <- mr_conmix(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome),  psi = 3, CIMin = NA, CIMax = NA, CIStep = 0.0001)
        res_conmix <- MendelianRandomization::mr_conmix(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome),  psi = 3, CIMin = -1, CIMax = 1, CIStep = 0.0001)
        beta_conmix <- res_conmix$Estimate
        CIlower_conmix <- res_conmix$CILower
        CIupper_conmix <- res_conmix$CIUpper
        pval_conmix <- res_conmix$Pvalue
        
        
        ## save MR results ##
        res_all <- cbind("Raw", myOutcome, myExposure, NSNPs,
                         beta_ivw, se_ivw, pval_ivw, Q_ivw, Q_df_ivw, Q_pval_ivw,
                         beta_egger, se_egger, pval_egger, inter_egger, inter_se_egger, inter_pval_egger,
                         NA, NA, NA, NA, NA,
                         beta_conmix, CIlower_conmix, CIupper_conmix, pval_conmix)
        colnames(res_all) <- res_colnames
        write.table(res_all, file = paste(DIR_RESULT,  "/MRresults_",myExposure, "_vs_", myOutcome, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
        
        
        ## scatter plot using MendelianRandomization package ##
        pdf(file=paste(DIR_PLOT, "/",myExposure, "_vs_", myOutcome, "_scatterPlot_MendelianRandomization.pdf",sep = ""), width=7, height=7)
        myplot <- MendelianRandomization::mr_plot(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome), line="ivw", orientate = TRUE, interactive = FALSE)
        print(myplot)
        dev.off()
        
        ## scatter plot and forest plot using TwoSampleMR package. ##
        ## scatter plot ##
        pdf(file=paste(DIR_PLOT, "/", myExposure, "_vs_", myOutcome, "_scatterPlot_TwoSampleMR.pdf",sep = ""), width=7, height=7)
        res <- TwoSampleMR::mr(mydat, method = c("mr_ivw", "mr_egger_regression", "mr_simple_median",  "mr_weighted_median"))
        p1 <- TwoSampleMR::mr_scatter_plot(res, mydat)
        print(p1)
        dev.off()
        
        ## forest plot for each IV ##
        pdf(file=paste(DIR_PLOT, "/", myExposure, "_vs_", myOutcome, "_forestPlot_TwoSampleMR.pdf",sep = ""), width=7, height=7)       ## Forest plot
        res_single <- TwoSampleMR::mr_singlesnp(mydat, single_method="mr_meta_fixed", all_method=c("mr_ivw", "mr_egger_regression", "mr_simple_median",  "mr_weighted_median"))
        p2 <- TwoSampleMR::mr_forest_plot(res_single)
        print(p2)
        dev.off()
        
      }else{
        
        tryCatch({
          # write your intended code here
          ##### MR_PRESSO: global test + outlier test + distortiont test #####
          ## update the settings for MR-PRESSO:
          ## NbDistribution  = 5000
          ## SignifThreshold = 0.10  --> remove SNPs that with "Outlier Test P value < 0.10"
          results_presso <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                                      SdOutcome = "se.outcome", SdExposure = "se.exposure",
                                      OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                                      data = mydat, NbDistribution = 5000, SignifThreshold = 0.10)
          capture.output(results_presso, file = paste(DIR_RESULT, "/MRPRESSO_", myExposure, "_VS_", myOutcome,  ".txt", sep=""))
          MR_PRESSO_main_results <- results_presso$`Main MR results`
          
          ## Save Global test P value - pleiotropy test
          if(is.null(results_presso$`MR-PRESSO results`$Pvalue)){
            GlobalTest <- results_presso$`MR-PRESSO results`$`Global Test`$Pvalue
          }else{
            GlobalTest <- results_presso$`MR-PRESSO results`$Pvalue
          }
          
          MR_PRESSO_main_results <- cbind(myOutcome, myExposure, MR_PRESSO_main_results, c(GlobalTest, NA))
          MR_PRESSO_main_results <- unfactor(MR_PRESSO_main_results)
          colnames(MR_PRESSO_main_results)[dim(MR_PRESSO_main_results)[2]] <- "GlobalTest_Pval"
          
          
          ############################
          ##### MR analyses: Raw #####
          ############################
          ## MR-IVW with fixed-effect ##
          res_ivw <- MendelianRandomization::mr_ivw(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome))
          NSNPs <- res_ivw$SNPs
          beta_ivw <- res_ivw$Estimate
          se_ivw <- res_ivw$StdError
          pval_ivw <- res_ivw$Pvalue
          
          ## Heteogeneity test: IVW  ##
          Q_ivw <- res_ivw$Heter.Stat[1]
          Q_df_ivw <- NSNPs - 1
          Q_pval_ivw <- res_ivw$Heter.Stat[2]
          
          ## MR_Egger ##
          res_egger <- MendelianRandomization::mr_egger(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome))
          beta_egger <- res_egger$Estimate
          se_egger <- res_egger$StdError.Est
          pval_egger <- res_egger$Pvalue.Est
          
          inter_egger <- res_egger$Intercept
          inter_se_egger <- res_egger$StdError.Int 
          inter_pval_egger <- res_egger$Pvalue.Int
          
          ## MR contamination mixture model ##
          # res_conmix <- mr_conmix(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome),  psi = 3, CIMin = NA, CIMax = NA, CIStep = 0.0001)
          res_conmix <- MendelianRandomization::mr_conmix(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome),  psi = 3, CIMin = -1, CIMax = 1, CIStep = 0.0001)
          beta_conmix <- res_conmix$Estimate
          CIlower_conmix <- res_conmix$CILower
          CIupper_conmix <- res_conmix$CIUpper
          pval_conmix <- res_conmix$Pvalue
          
          ## save MR results ##
          res_all <- cbind("Raw", myOutcome, myExposure, NSNPs,
                           beta_ivw, se_ivw, pval_ivw, Q_ivw, Q_df_ivw, Q_pval_ivw,
                           beta_egger, se_egger, pval_egger, inter_egger, inter_se_egger, inter_pval_egger,
                           MR_PRESSO_main_results[MR_PRESSO_main_results$`MR Analysis` == "Raw", -c(1:4)],
                           beta_conmix, CIlower_conmix, CIupper_conmix, pval_conmix)
          colnames(res_all) <- res_colnames
          
          ## scatter plot using MendelianRandomization package ##
          pdf(file=paste(DIR_PLOT, "/",myExposure, "_vs_", myOutcome, "_scatterPlot_MendelianRandomization.pdf",sep = ""), width=7, height=7)
          myplot <- MendelianRandomization::mr_plot(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome), line="ivw", orientate = TRUE, interactive = FALSE)
          print(myplot)
          dev.off()
          
          ## scatter plot and forest plot using TwoSampleMR package. ##
          ## scatter plot ##
          pdf(file=paste(DIR_PLOT, "/", myExposure, "_vs_", myOutcome, "_scatterPlot_TwoSampleMR.pdf",sep = ""), width=7, height=7)
          res <- TwoSampleMR::mr(mydat, method = c("mr_ivw", "mr_egger_regression", "mr_simple_median",  "mr_weighted_median"))
          p1 <- TwoSampleMR::mr_scatter_plot(res, mydat)
          print(p1)
          dev.off()
          
          ## forest plot for each IV ##
          pdf(file=paste(DIR_PLOT, "/", myExposure, "_vs_", myOutcome, "_forestPlot_TwoSampleMR.pdf",sep = ""), width=7, height=7)       ## Forest plot
          res_single <- TwoSampleMR::mr_singlesnp(mydat, single_method="mr_meta_fixed", all_method=c("mr_ivw", "mr_egger_regression", "mr_simple_median",  "mr_weighted_median"))
          p2 <- TwoSampleMR::mr_forest_plot(res_single)
          print(p2)
          dev.off()
          
          
          ###########################################
          ##### MR analysese: outlier-corrected #####
          ###########################################
          ## (1) MR-PRESSO outlier = TRUE
          ## data after outlier-corrected ##
          if(is.null(results_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)){
            mydat_outlier_corrected <- mydat
          }else if(results_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`[1] == "No significant outliers" ){
            mydat_outlier_corrected <- mydat
          }else{
            mydat_outlier_corrected <- mydat[- results_presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, ]
          }
          
          ## (2) Fstat < 10
          mydat_outlier_corrected <- mydat_outlier_corrected[mydat_outlier_corrected$Fstat >= 10, ]
          dim(mydat_outlier_corrected)
          write.table(mydat_outlier_corrected, file = paste(DIR_DATA_FOR_MR, "/harmonised_sumstats_", myExposure,"_vs_", myOutcome, "_outlierCorrected.txt", sep=""), col.names = T, row.names = F, sep = "\t", quote = F)
          
          
          ## make sure to include the allele with frequency close to 0.5 ## 
          for(k in 1:dim(mydat_outlier_corrected)[1]){
            if(mydat_outlier_corrected$palindromic[k] == TRUE) { mydat_outlier_corrected$mr_keep[k] <- as.logical("TRUE") }
          }
          
          ##### MR analyses: Outlier-corrected #####
          if(dim(mydat_outlier_corrected)[1] == 1){
            ## MR wald's ratio ##
            res_oc <- MendelianRandomization::mr_ivw(mr_input(bx = mydat_outlier_corrected$beta.exposure, bxse = mydat_outlier_corrected$se.exposure, by = mydat_outlier_corrected$beta.outcome, byse = mydat_outlier_corrected$se.outcome))
            NSNPs_oc <- res_oc$SNPs
            beta_ivw_oc <- res_oc$Estimate
            se_ivw_oc <- res_oc$StdError
            pval_ivw_oc <- res_oc$Pvalue
            
            res_all_oc <- cbind("Outlier-corrected", myOutcome, myExposure, NSNPs_oc, beta_ivw_oc, se_ivw_oc, pval_ivw_oc, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
            colnames(res_all_oc) <- res_colnames
            
          }else if(dim(mydat_outlier_corrected)[1] == 2){
            ## MR_IVW without correlation matrix ##
            res_ivw_oc <-  MendelianRandomization::mr_ivw(mr_input(bx = mydat_outlier_corrected$beta.exposure, bxse = mydat_outlier_corrected$se.exposure, by = mydat_outlier_corrected$beta.outcome, byse = mydat_outlier_corrected$se.outcome))
            NSNPs_oc <- res_ivw_oc$SNPs
            beta_ivw_oc <- res_ivw_oc$Estimate
            se_ivw_oc <- res_ivw_oc$StdError
            pval_ivw_oc <- res_ivw_oc$Pvalue
            
            ## Heteogeneity test: IVW  ##
            Q_ivw_oc <- res_ivw_oc$Heter.Stat[1]
            Q_df_ivw_oc <- NSNPs_oc - 1
            Q_pval_ivw_oc <- res_ivw_oc$Heter.Stat[2]
            
            ## combine all results ##
            res_all_oc <- cbind("Outlier-corrected", myOutcome, myExposure, NSNPs_oc,
                                beta_ivw_oc, se_ivw_oc, pval_ivw_oc, Q_ivw_oc, Q_df_ivw_oc, Q_pval_ivw_oc,
                                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
            colnames(res_all_oc) <- res_colnames
            
          }else{
            ## MR_IVW without correlation matrix ##
            res_ivw_oc <-  MendelianRandomization::mr_ivw(mr_input(bx = mydat_outlier_corrected$beta.exposure, bxse = mydat_outlier_corrected$se.exposure, by = mydat_outlier_corrected$beta.outcome, byse = mydat_outlier_corrected$se.outcome))
            NSNPs_oc <- res_ivw_oc$SNPs
            beta_ivw_oc <- res_ivw_oc$Estimate
            se_ivw_oc <- res_ivw_oc$StdError
            pval_ivw_oc <- res_ivw_oc$Pvalue
            
            ## Heteogeneity test: IVW  ##
            Q_ivw_oc <- res_ivw_oc$Heter.Stat[1]
            Q_df_ivw_oc <- NSNPs_oc - 1
            Q_pval_ivw_oc <- res_ivw_oc$Heter.Stat[2]
            
            ## MR_Egger without correlation matrix ##
            res_egger_oc <-  MendelianRandomization::mr_egger(mr_input(bx = mydat_outlier_corrected$beta.exposure, bxse = mydat_outlier_corrected$se.exposure, by = mydat_outlier_corrected$beta.outcome, byse = mydat_outlier_corrected$se.outcome))
            beta_egger_oc <- res_egger_oc$Estimate
            se_egger_oc <- res_egger_oc$StdError.Est
            pval_egger_oc <- res_egger_oc$Pvalue.Est
            
            inter_egger_oc <- res_egger_oc$Intercept
            inter_se_egger_oc <- res_egger_oc$StdError.Int 
            inter_pval_egger_oc <- res_egger_oc$Pvalue.Int
            
            ## MR contamination mixture model ##
            res_conmix_oc <-  MendelianRandomization::mr_conmix(mr_input(bx = mydat_outlier_corrected$beta.exposure, 
                                                                bxse = mydat_outlier_corrected$se.exposure,
                                                                by = mydat_outlier_corrected$beta.outcome,
                                                                byse = mydat_outlier_corrected$se.outcome),
                                                       psi = 3, CIMin = -1, CIMax = 1, CIStep = 0.0001)
            
            beta_conmix_oc <- res_conmix_oc$Estimate
            CIlower_conmix_oc <- res_conmix_oc$CILower
            CIupper_conmix_oc <- res_conmix_oc$CIUpper
            pval_conmix_oc <- res_conmix_oc$Pvalue
            
            ## save MR results ##
            res_all_oc <- cbind("Outlier-corrected", myOutcome, myExposure, NSNPs_oc,
                                beta_ivw_oc, se_ivw_oc, pval_ivw_oc, Q_ivw_oc, Q_df_ivw_oc, Q_pval_ivw_oc,
                                beta_egger_oc, se_egger_oc, pval_egger_oc, inter_egger_oc, inter_se_egger_oc, inter_pval_egger_oc,
                                MR_PRESSO_main_results[MR_PRESSO_main_results$`MR Analysis` == "Outlier-corrected", -c(1:4)],
                                beta_conmix_oc, CIlower_conmix_oc, CIupper_conmix_oc, pval_conmix_oc)
            colnames(res_all_oc) <- res_colnames
            
            # ## save MR plots ##
            ## scatter plot using MendelianRandomization package ##
            pdf(file=paste(DIR_PLOT, "/",myExposure, "_vs_", myOutcome, "_scatterPlot_outlierCorrected_MendelianRandomization.pdf",sep = ""), width=7, height=7)
            myplot <- MendelianRandomization::mr_plot(mr_input(bx = mydat_outlier_corrected$beta.exposure, bxse = mydat_outlier_corrected$se.exposure, 
                                                               by = mydat_outlier_corrected$beta.outcome, byse = mydat_outlier_corrected$se.outcome), line="ivw", orientate = TRUE, interactive = FALSE)
            print(myplot)
            dev.off()
            
            ## scatter plot and forest plot using TwoSampleMR package. ##
            pdf(file=paste(DIR_PLOT, "/", myExposure, "_vs_", myOutcome, "_scatterPlot_outlierCorrected_TwoSampleMR.pdf",sep = ""), width=7, height=7)      ## Scatter plot
            res_oc <- TwoSampleMR::mr(mydat_outlier_corrected, method=  c("mr_ivw","mr_egger_regression",  "mr_simple_median",  "mr_weighted_median"))
            p1 <- TwoSampleMR::mr_scatter_plot(res_oc, mydat_outlier_corrected)
            print(p1)
            dev.off()

            pdf(file=paste(DIR_PLOT, "/", myExposure, "_vs_", myOutcome, "_forestPlot_outlierCorrected_TwoSampleMR.pdf",sep = ""), width=7, height=7)       ## Forest plot
            res_single_oc <- TwoSampleMR::mr_singlesnp(mydat_outlier_corrected, single_method="mr_meta_fixed", all_method=c("mr_ivw", "mr_egger_regression", "mr_simple_median",  "mr_weighted_median"))
            p2 <- TwoSampleMR::mr_forest_plot(res_single_oc)
            print(p2)
            dev.off()
          
            
          }
          
          ## write out results into files ##
          res_all_combined <- rbind(res_all, res_all_oc)
          write.table(res_all_combined, file = paste(DIR_RESULT, "/MRresults_", myExposure, "_vs_", myOutcome, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
          
          
        }, error = function(e) {
          # if the MR_PRESSO returns error ##
          print("Error")
          
          ## MR-IVW with fixed-effect ##
          res_ivw <- MendelianRandomization::mr_ivw(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome))
          NSNPs <- res_ivw$SNPs
          beta_ivw <- res_ivw$Estimate
          se_ivw <- res_ivw$StdError
          pval_ivw <- res_ivw$Pvalue
          
          ## Heteogeneity test: IVW  ##
          Q_ivw <- res_ivw$Heter.Stat[1]
          Q_df_ivw <- NSNPs - 1
          Q_pval_ivw <- res_ivw$Heter.Stat[2]
          
          ## MR_Egger ##
          res_egger <- MendelianRandomization::mr_egger(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome))
          beta_egger <- res_egger$Estimate
          se_egger <- res_egger$StdError.Est
          pval_egger <- res_egger$Pvalue.Est
          
          inter_egger <- res_egger$Intercept
          inter_se_egger <- res_egger$StdError.Int 
          inter_pval_egger <- res_egger$Pvalue.Int
          
          ## MR contamination mixture model ##
          # res_conmix <- mr_conmix(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome),  psi = 3, CIMin = NA, CIMax = NA, CIStep = 0.0001)
          res_conmix <- MendelianRandomization::mr_conmix(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome),  psi = 3, CIMin = -1, CIMax = 1, CIStep = 0.0001)
          beta_conmix <- res_conmix$Estimate
          CIlower_conmix <- res_conmix$CILower
          CIupper_conmix <- res_conmix$CIUpper
          pval_conmix <- res_conmix$Pvalue
          
          ## save MR results ##
          res_all <- cbind("Raw", myOutcome, myExposure, NSNPs,
                           beta_ivw, se_ivw, pval_ivw, Q_ivw, Q_df_ivw, Q_pval_ivw,
                           beta_egger, se_egger, pval_egger, inter_egger, inter_se_egger, inter_pval_egger,
                           NA, NA, NA, NA, NA,
                           beta_conmix, CIlower_conmix, CIupper_conmix, pval_conmix)
          colnames(res_all) <- res_colnames
          write.table(res_all, file = paste(DIR_RESULT,  "/MRresults_", myExposure, "_vs_", myOutcome, ".txt", sep = ""), col.names = T, row.names = F, sep = "\t", quote = F)
          
          ## scatter plot using MendelianRandomization package ##
          pdf(file=paste(DIR_PLOT, "/",myExposure, "_vs_", myOutcome, "_scatterPlot_MendelianRandomization.pdf",sep = ""), width=7, height=7)
          myplot <- MendelianRandomization::mr_plot(mr_input(bx = mydat$beta.exposure, bxse = mydat$se.exposure, by = mydat$beta.outcome, byse = mydat$se.outcome), line="ivw", orientate = TRUE, interactive = FALSE)
          print(myplot)
          dev.off()

          ## scatter plot and forest plot using TwoSampleMR package. ##
          res <- TwoSampleMR::mr(mydat, method=  c("mr_ivw","mr_egger_regression",  "mr_simple_median",  "mr_weighted_median"))
          ## scatter plot ##
          pdf(file=paste(DIR_PLOT, "/", myExposure, "_vs_", myOutcome, "_scatterPlot_TwoSampleMR.pdf",sep = ""), width=7, height=7)      ## Scatter plot
          p1 <- TwoSampleMR::mr_scatter_plot(res, mydat)
          print(p1)
          dev.off()

          ## forest plot for each IV ##
          pdf(file=paste(DIR_PLOT, "/", myExposure, "_vs_", myOutcome, "_forestPlot_TwoSampleMR.pdf",sep = ""), width=7, height=7)       ## Forest plot
          res_single <- TwoSampleMR::mr_singlesnp(mydat, single_method="mr_meta_fixed", all_method=c("mr_ivw", "mr_egger_regression", "mr_simple_median",  "mr_weighted_median"))
          p2 <- TwoSampleMR::mr_forest_plot(res_single)
          print(p2)
          dev.off()
          
        }, finally = {
          # this will execute no matter what else happened
          print("IVs >= 4")
        }) # End of tryCatch Function
        
      } ## End of if conditions
      
    } ## End of if condition: dim(dat)[1] !=0
    
  } ## END of if condition: file.exists()
  
} ## End of function "MR_ANALYSIS"

##### ! NOT RUN #####
