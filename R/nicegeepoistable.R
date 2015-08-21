
#' Lisa's GEE Poisson Regression Table Function - fit must be geeglm type model 
#'
#' @param df Dataframe. 
#' @param fit Poisson geeglm (geepack) model object.
#' @param ref If TRUE, the reference category gets its own line (left blank). Default is FALSE.
#' @param labels Covariate labels - default is NA.
#' @param blanks If TRUE, blank lines will be inserted separating covariates - default is FALSE.
#' @param overallp is FALSE always (LR p-value not supported in geeglm)
#' @param est.dec Number of decimal places for OR estimates - default is 4.
#' @param ci.dec Number of decimal places for 95% CI - default is 4.
#' @param pval.dec Number of decimal places for pvalues - default is 4.
#' @keywords pretty table gee poisson regression geeglm geepack 
#' @importFrom xtable xtable
#' @export 
nicegeepoistable <- function(df, 
                        fit, 
                        intercept = FALSE,
                        ref = FALSE,
                        labels = NA, 
                        blanks = FALSE,
                        overallp = FALSE,
                        est.dec = 4,
                        ci.dec = 4,
                        pval.dec = 3){
    library(xtable)
    
    overallp = FALSE ## does not support LR p-value
    
    ciformat <- paste("%.", ci.dec, "f", sep="")
    
    expconf <- function(x){
        paste("[", 
              sprintf(ciformat, round(exp(x), ci.dec)[1]), ", ",
              sprintf(ciformat, round(exp(x), ci.dec)[2]) , "]", sep="")
    }
    
    trim <- function(x) {
        gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
    }
    
    esformat <- paste("%.", est.dec, "f", sep="")
    
    orig <- data.frame(summary(fit)$coef)
    coef_tbl <- data.frame(summary(fit)$coef)
    coef_tbl$Estimate <- exp(summary(fit)$coef[,"Estimate"])
    coef_tbl$Estimate <- sprintf(esformat, round(coef_tbl$Estimate, est.dec))
    
    names(coef_tbl)[grepl("Est", names(coef_tbl))] <- "R_R"
    names(coef_tbl)[grepl("Pr", names(coef_tbl))] <- "p_value"
    
    cimat <- matrix(data = NA, ncol = 2, nrow = nrow(coef_tbl))
    cimat[,1] <- orig$Estimate - (qnorm(0.975)*orig$Std.err)
    cimat[,2] <- orig$Estimate + (qnorm(0.975)*orig$Std.err)
    
    cimat <- data.frame(cimat)
    rownames(cimat) <- rownames(coef_tbl)
    names(cimat) <- c("2.5 %", "97.5 %")
    
    coef_tbl$CI <- apply(cimat,1,expconf)
    
    sformat <- paste("%.", pval.dec, "f", sep="")
    
    p_value2 <- sprintf(sformat, round(coef_tbl$p_value, pval.dec))
    if (pval.dec == 4) p_value2[coef_tbl$p_value < 0.0001] <- "< 0.0001"
    if (pval.dec == 3) p_value2[coef_tbl$p_value < 0.001] <- "< 0.001"
    if (pval.dec == 2) p_value2[coef_tbl$p_value < 0.01] <- "< 0.01"
    coef_tbl$p_value <- p_value2
    
    covs <- strsplit(as.character(fit$formula), "~")[[3]]
    covs <- trim(unlist(strsplit(covs, "+", fixed=T)))
    
    coef_tbl <- coef_tbl[,c("R_R", "CI", "p_value")]
    
    tbl <- NULL
    
    if (intercept == TRUE){
        tbl <- coef_tbl["(Intercept)",]
        if (overallp == TRUE) tbl$Overall_pvalue <- NA
        tbl$Variable <- "(Intercept)"
    }
    
    for (i in 1:length(covs)){
        
        if (attr(fit$terms, "dataClass")[i+1] == "numeric"){
            tmp <- coef_tbl[grepl(covs[i], rownames(coef_tbl)),]
            if (overallp == TRUE) {
                op <- drop1(fit, test = "Chisq")[covs[i],"Pr(>Chi)"]
                op2 <- sprintf(sformat, round(op, pval.dec))
                if (pval.dec == 4) op2[op < 0.0001] <- "< 0.0001"
                if (pval.dec == 3) op2[op < 0.001] <- "< 0.001"
                if (pval.dec == 2) op2[op < 0.01] <- "< 0.01"
                tmp$Overall_pvalue <- op2
            }
            if (is.na(labels[1])) tmp$Variable <- covs[i]
            if (!is.na(labels[1])) tmp$Variable <- labels[i]
            
            if (blanks == TRUE) tbl <- rbind(tbl, blank, tmp)
            if (blanks == FALSE) tbl <- rbind(tbl, tmp)
            
        }
        
        if (attr(fit$terms, "dataClass")[i+1] == "factor" |
                attr(fit$terms, "dataClass")[i+1] == "character"){
            
            df[,covs[i]] <- as.factor(df[,covs[i]] )
            tmp <- coef_tbl[grepl(covs[i], rownames(coef_tbl)),]
            if (overallp == TRUE) tmp$Overall_pvalue <- NA
            title <- data.frame(tmp[1,])
            title[1,] <- NA
            if (overallp == TRUE) {
                op <- drop1(fit, test = "Chisq")[covs[i],"Pr(>Chi)"]
                op2 <- sprintf(sformat, round(op, pval.dec))
                if (pval.dec == 4) op2[op < 0.0001] <- "< 0.0001"
                if (pval.dec == 3) op2[op < 0.001] <- "< 0.001"
                if (pval.dec == 2) op2[op < 0.01] <- "< 0.01"
                title$Overall_pvalue <- op2
            }
            if (is.na(labels[1])) title$Variable <- covs[i]
            if (!is.na(labels[1])) title$Variable <- labels[i]
            blank <- data.frame(tmp[1,])
            blank <- NA
            reference <- data.frame(tmp[1,])
            reference[1,] <- NA
            reference$Variable <- paste("*", levels(df[,covs[i]])[1])
            if (ref == FALSE){
                tmp$Variable <- 
                    paste("*", levels(df[,covs[i]])[2:nlevels(df[,covs[i]])], "vs.", 
                          levels(df[,covs[i]])[1])
            }
            if (ref == TRUE){
                tmp$Variable <- 
                    paste("*", levels(df[,covs[i]])[2:nlevels(df[,covs[i]])])
            }
            
            if (blanks == TRUE){
                if (ref == TRUE) tbl <- rbind(tbl, blank, title, reference, tmp)
                if (ref == FALSE) tbl <- rbind(tbl, blank, title, tmp)
            }
            if (blanks == FALSE){
                if (ref == TRUE) tbl <- rbind(tbl, title, reference, tmp)
                if (ref == FALSE) tbl <- rbind(tbl, title, tmp)
            }
            
        }    
    }
    
    
    tbl <- tbl[,c(ncol(tbl), 2:ncol(tbl)-1)]
    if (overallp == TRUE){
        names(tbl) <- c("Variable", "Rate Ratio", "95% CI", "Wald p-value", "LR p-value")
    }
    if (overallp == FALSE){
        names(tbl) <- c("Variable", "Rate Ratio", "95% CI", "p-value")
    }
    
    if (overallp == TRUE){
        print(xtable(tbl, align="llccrr"), type='html', 
              include.rownames=F)
    }
    if (overallp == FALSE){
        print(xtable(tbl, align = "llccr"), type='html', 
              include.rownames=F)
    }
    return(tbl)
}