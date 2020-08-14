#########################################################################
##### MR Contamination mixture method - function return with Pvalue #####
#########################################################################
library(MendelianRandomization)
setClass("MRConMix",
         representation(Exposure = "character",
                        Outcome  = "character",
                        Psi      = "numeric",
                        Estimate = "numeric",
                        CIRange  = "numeric",
                        CILower  = "numeric",
                        CIUpper  = "numeric",
                        CIMin    = "numeric",
                        CIMax    = "numeric",
                        CIStep   = "numeric",
                        Valid    = "numeric",
                        ValidSNPs= "character",
                        Pvalue   = "numeric",
                        Alpha    = "numeric",
                        SNPs = "numeric")
)

#' @docType methods
#' @rdname mr_conmix

setMethod("mr_conmix",
          "MRInput",
          function(object,
                   psi    = 0,
                   CIMin  = NA,
                   CIMax  = NA,
                   CIStep = 0.01,
                   alpha = 0.05){
            
            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse
            
            nsnps = length(Bx)
            
            ratio = By/Bx; ratio.se = abs(Byse/Bx);
            if (is.na(CIMin)) { CIMin = min((By-2*Byse)/Bx) }
            if (is.na(CIMax)) { CIMax = max((By+2*Byse)/Bx) }
            if (psi < 0 | psi == 0) {   psi = 1.5*sd(ratio)  }
            theta = seq(from = CIMin, to = CIMax, by = CIStep)
            iters = length(theta)
            lik=NULL
            for (j1 in 1:iters) {
              lik.inc = -(theta[j1]-ratio)^2/2/ratio.se^2 - log(sqrt(2*pi*ratio.se^2))
              lik.exc = -ratio^2/2/(psi^2+ratio.se^2) - log(sqrt(2*pi*(psi^2+ratio.se^2)))
              valid = (lik.inc>lik.exc)*1
              lik[j1] = sum(c(lik.inc[valid==1], lik.exc[valid==0]))
              if (which.max(lik)==length(lik)) { valid.best = valid }
            }
            phi = ifelse(sum(valid.best)<1.5, 1,
                         max(sqrt(sum(((ratio[valid.best==1]-weighted.mean(ratio[valid.best==1],
                                                                           ratio.se[valid.best==1]^-2))^2*
                                         ratio.se[valid.best==1]^-2))/(sum(valid.best)-1)), 1))
            loglik = lik
            
            lik.inc0 = -ratio^2/2/ratio.se^2 - log(sqrt(2*pi*ratio.se^2))
            lik.exc0 = -ratio^2/2/(psi^2+ratio.se^2) - log(sqrt(2*pi*(psi^2+ratio.se^2)))
            valid = (lik.inc0>lik.exc0)*1
            loglik0 = sum(c(lik.inc0[valid==1], lik.exc0[valid==0]))
            
            
            whichin = which(2*loglik>(2*max(loglik)-qchisq(1-alpha, df=1)*phi^2))
            # provides an index of estimate values in the 95% confidence interval
            betaConMix = CIMin+CIStep*(which.max(loglik)-1)
            # modal estimate
            CIRange    = CIMin+CIStep*(whichin-1);
            CILower <- c(min(CIRange), CIRange[which(diff(CIRange)>1.01*CIStep)+1])
            CIUpper <- c(CIRange[which(diff(CIRange)>1.01*CIStep)], max(CIRange))
            
            Pvalue = pchisq(2*(max(loglik)-loglik0)*phi^2, df=1, lower.tail=FALSE)
            
            return(new("MRConMix",
                       Exposure = object@exposure,
                       Outcome = object@outcome,
                       Psi = as.numeric(psi),
                       
                       Estimate = as.numeric(betaConMix),
                       CIRange  = as.numeric(CIRange),
                       CILower  = as.numeric(CILower),
                       CIUpper  = as.numeric(CIUpper),
                       
                       CIMin    = as.numeric(CIMin),
                       CIMax    = as.numeric(CIMax),
                       CIStep   = as.numeric(CIStep),
                       Valid    = as.numeric(which(valid.best==1)),
                       ValidSNPs= as.character(object$snps[which(valid.best==1)]),
                       Pvalue   = as.numeric(Pvalue),
                       
                       SNPs = nsnps,
                       Alpha = alpha))
          }
)


setMethod("show",
          "MRConMix",
          function(object){
            
            if (object@CIMax%in%object@CIRange & object@CIMin%in%object@CIRange) {
              cat("Confidence interval range too narrow. Please decrease CIMin and increase CIMax and try again.") }
            else if (object@CIMax>max(object@CIRange) & object@CIMin%in%object@CIRange) {
              cat("Lower bound of confidence interval range too high. Please decrease CIMin and try again.") }
            if (object@CIMax%in%object@CIRange & object@CIMin<min(object@CIRange)) {
              cat("Upper bound of confidence interval range too low. Please increase CIMax and try again.") }
            if (object@CIMax>max(object@CIRange) & object@CIMin<min(object@CIRange)) {
              Interval_type <- paste(100*(1-object@Alpha), "% CI", sep = "")
              
              dps = max(ceiling(-log10(object@CIStep)), 1)
              Ranges <- ifelse(sum(diff(object@CIRange)>1.01*object@CIStep)==0, "Single range", "Multiple ranges");
              
              if (Ranges == "Single range") {
                Statistic <- c("Method", "Estimate", Interval_type, "", "p-value")
                Value <- c("ConMix", decimals(object@Estimate, dps),
                           paste(decimals(min(object@CIRange), dps), ",", sep = ""), decimals(max(object@CIRange), dps), signif(object@Pvalue, 3))
                output.table <- data.frame(matrix(Value, nrow = 1))
                colnames(output.table) <- Statistic
                Ranges.text <- "Note: confidence interval is a single range of values.\n"
              }
              
              if (Ranges == "Multiple ranges") {
                Statistic <- c("Method", "Estimate", Interval_type, "")
                Value <- c("ConMix", rep("", length(object@CILower)-1),
                           decimals(object@Estimate, dps), rep("", length(object@CILower)-1),
                           paste(decimals(object@CILower, dps), ",", sep = ""), decimals(object@CIUpper, dps))
                output.table <- data.frame(matrix(Value, nrow = length(object@CILower), byrow=FALSE))
                colnames(output.table) <- Statistic
                Ranges.text <- "Note: confidence interval contains multiple ranges of values.\n"
              }
              
              
              cat("\nContamination mixture method\n")
              cat("(Standard deviation of invalid estimands = ", object@Psi, ")\n\n" , sep = "")
              cat("Number of Variants :", object@SNPs, "\n")
              
              cat("------------------------------------------------------------------\n")
              print(output.table, quote = F, row.names = FALSE, justify = "left")
              cat("------------------------------------------------------------------\n")
              cat(Ranges.text)
            }
          }
)  ## End of mr_conmix function


