
#' Lisa's GEE Regression Table Function - fit must be geeglm type model 
#'
#' @param df Dataframe. 
#' @param fit geeglm (geepack) model object.
#' @param ref If TRUE, the reference category gets its own line (left blank). 
#' Default is FALSE.
#' @param family Model family in quotes ("guassian", "binomial", "poisson"). 
#' @param id id variable.
#' @param corstr corstr type ("exch" is default).
#' @param labels Covariate labels - default is NA.
#' @param blanks If TRUE, blank lines will be inserted separating covariates - default is FALSE.
#' @param overallp is whether to do overall Chisq test per variable
#' @param est.dec Number of decimal places for OR estimates - default is 4.
#' @param ci.dec Number of decimal places for 95% CI - default is 4.
#' @param pval.dec Number of decimal places for pvalues - default is 4.
#' @keywords pretty table gee logistic regression geeglm geepack 
#' @importFrom xtable xtable
#' @export 
nicegeetable <- function(df, 
                        fit, 
                        family = "gaussian",
                        id,
                        corstr="exch",
                        intercept = FALSE,
                        ref = FALSE,
                        labels = NA, 
                        blanks = FALSE,
                        overallp = FALSE,
                        est.dec = 2,
                        ci.dec = 2,
                        pval.dec = 3){
    library(xtable)
    
    exp <- FALSE
    if (family %in% c("binomial", "poisson")){
        exp <- TRUE
    }
    
    ciformat <- paste("%.", ci.dec, "f", sep="")
    
    expconf <- function(x){
        paste("[", 
              sprintf(ciformat, round(exp(x), ci.dec)[1]), ", ",
              sprintf(ciformat, round(exp(x), ci.dec)[2]) , "]", sep="")
    }
    
    conf <- function(x){
        paste("[", 
              sprintf(ciformat, round(x, ci.dec)[1]), ", ",
              sprintf(ciformat, round(x, ci.dec)[2]) , "]", sep="")
    }
    
    trim <- function(x) {
        gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
    }
    
    esformat <- paste("%.", est.dec, "f", sep="")
    
    orig <- data.frame(summary(fit)$coef)
    coef_tbl <- data.frame(summary(fit)$coef)
    
    coef_tbl$Estimate <- summary(fit)$coef[,"Estimate"]
    if (exp == TRUE) coef_tbl$Estimate <- exp(coef_tbl$Estimate)
    coef_tbl$Estimate <- sprintf(esformat, round(coef_tbl$Estimate, est.dec))
    
    names(coef_tbl)[grepl("Est", names(coef_tbl))] <- "R_R"
    names(coef_tbl)[grepl("Pr", names(coef_tbl))] <- "p_value"
    
    cimat <- matrix(data = NA, ncol = 2, nrow = nrow(coef_tbl))
    cimat[,1] <- orig$Estimate - (qnorm(0.975)*orig$Std.err)
    cimat[,2] <- orig$Estimate + (qnorm(0.975)*orig$Std.err)
    
    cimat <- data.frame(cimat)
    rownames(cimat) <- rownames(coef_tbl)
    names(cimat) <- c("2.5 %", "97.5 %")
    
    if (exp == TRUE)  coef_tbl$CI <- apply(cimat,1,expconf)
    if (exp == FALSE) coef_tbl$CI <- apply(cimat,1,conf)
    
    sformat <- paste("%.", pval.dec, "f", sep="")
    
    p_value2 <- sprintf(sformat, round(coef_tbl$p_value, pval.dec))
    if (pval.dec == 4) p_value2[coef_tbl$p_value < 0.0001] <- "< 0.0001"
    if (pval.dec == 3) p_value2[coef_tbl$p_value < 0.001] <- "< 0.001"
    if (pval.dec == 2) p_value2[coef_tbl$p_value < 0.01] <- "< 0.01"
    coef_tbl$p_value <- p_value2
    
    out <- strsplit(as.character(fit$formula), "~")[[2]]
    
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
                
                form <- as.formula(paste(out, " ~ ", paste(covs[-i], collapse=" + "), 
                                         sep=""))
                
                fit2 <- geeglm(form,
                               id = id,
                               family = binomial,
                               data = df,
                               corstr = corstr)
                
                op <- anova(fit, fit2)[1,"P(>|Chi|)"]
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
                
                form <- as.formula(paste(out, " ~ ", paste(covs[-i], collapse=" + "), 
                                         sep=""))
                
                fit2 <- geeglm(form,
                               id = id,
                               family = binomial,
                               data = df,
                               corstr = corstr)
                
                op <- anova(fit, fit2)[1,"P(>|Chi|)"]
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
    
    estname <- "Estimate"
    if (family == "binomial") estname <- "Odds Ratio"
    if (family == "poisson")  estname <- "Rate Ratio"
    if (family == "gaussian") estname <- "Difference"
    
    tbl <- tbl[,c(ncol(tbl), 2:ncol(tbl)-1)]
    if (overallp == TRUE){
        names(tbl) <- c("Variable", estname, "95% CI", "Wald p-value", "Chisq p-value")
    }
    if (overallp == FALSE){
        names(tbl) <- c("Variable", estname, "95% CI", "p-value")
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