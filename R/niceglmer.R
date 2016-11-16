
#' Lisa's GLMER Regression Table Function 
#'
#' This function creates a nice looking regression table from a glmer object 
#' @param df Dataframe. 
#' @param fit GLMER binomial model object.
#' @param family Model family in quotes ("guassian", "binomial", "poisson"). 
#' @param covs Vector of covariate names in model.
#' @param intercept If TRUE the intercept will be included in the table. Default is FALSE.
#' @param ref If TRUE, the reference category gets its own line (left blank). Default is FALSE.
#' @param labels Covariate labels - default is NA (variable names are used).
#' @param blanks If TRUE, blank lines will be inserted separating covariates - default is FALSE.
#' @param overallp If TRUE, a likelihood ratio test pvalue (using drop1 Chisq tests) will be calculated for each variable. Default is TRUE.
#' @param est.dec Number of decimal places for estimates - default is 2.
#' @param ci.dec Number of decimal places for CI - default is 2.
#' @param pval.dec Number of decimal places for pvalues - default is 3.
#' @param estname Option to override estimate column name. Default is NA.
#' @param htmlTable Whether to use htmlTable package to display table (instead of xtable). Default = FALSE.
#' @param color Hex color to use for htmlTable output. Default = "#EEEEEE" (grey).
#' @param exp Option to exponentiate coefficients and CI's. Default is NA (exp for logistic and poisson families only).
#' @keywords nicetable mixed effects regression glmer lmer  
#' @importFrom xtable xtable
#' @export 
niceglmer <- function(df, 
                      fit, 
                      covs,
                      family = "gaussian",
                      exp = NA,
                      estname = NA,
                      intercept = FALSE,
                      ref = FALSE,
                      labels = NA, 
                      blanks = FALSE,
                      overallp = FALSE,
                      est.dec = 2,
                      ci.dec = 2,
                      pval.dec = 3,
                      htmlTable = FALSE,
                      color = "#EEEEEE",
                      printRMD = TRUE){
  
    # df = cdata
    # fit = x
    # covs = preds
    # intercept = FALSE
    # ref = FALSE
    # labels = NA
    # blanks = FALSE
    # overallp = FALSE
    # est.dec = 4
    # ci.dec = 4
    # pval.dec = 3
    # printRMD = TRUE
    # estname = NA
    # family = "binomial"
    # exp = NA
  
  ciformat <- paste("%.", ci.dec, "f", sep="")
  
  if (is.na(exp) & family %in% c("binomial", "poisson")){
    exp <- TRUE
  }
  if (is.na(exp) & family %in% c("gaussian")){
    exp <- FALSE
  }
  
  ### if no estimate name is specified pick a reasonable name for each situation
  if (is.na(estname) & family == "binomial" & exp == TRUE)  estname <- "OR"
  if (is.na(estname) & family == "poisson"  & exp == TRUE)  estname <- "RR"
  if (is.na(estname) & family == "gaussian" & exp == TRUE)  estname <- "Ratio"
  if (is.na(estname) & family == "binomial" & exp == FALSE) estname <- "Estimate"
  if (is.na(estname) & family == "poisson"  & exp == FALSE) estname <- "Estimate"
  if (is.na(estname) & family == "gaussian" & exp == FALSE) estname <- "Difference"
  
  if (exp){
    conf <- function(x){
      paste("[", 
            sprintf(ciformat, round(exp(x), ci.dec)[1]), ", ",
            sprintf(ciformat, round(exp(x), ci.dec)[2]) , "]", sep="")
    }
  }
  if (exp == FALSE){
    conf <- function(x){
      paste("[", 
            sprintf(ciformat, round(x, ci.dec)[1]), ", ",
            sprintf(ciformat, round(x, ci.dec)[2]) , "]", sep="")
    }
  }
  
  
  trim <- function(x) {
    gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
  }
  
  esformat <- paste("%.", est.dec, "f", sep="")
  
  orig <- data.frame(summary(fit)$coef)
  coef_tbl <- data.frame(summary(fit)$coef)
  if (exp)          coef_tbl$Estimate <- exp(summary(fit)$coef[,"Estimate"])
  if (exp == FALSE) coef_tbl$Estimate <-    (summary(fit)$coef[,"Estimate"])
  coef_tbl$Estimate <- sprintf(esformat, round(coef_tbl$Estimate, est.dec))
  
  names(coef_tbl)[grepl("Est", names(coef_tbl))] <- estname
  names(coef_tbl)[grepl("Pr", names(coef_tbl))] <- "p_value"
  
  cimat <- matrix(data = NA, ncol = 2, nrow = nrow(coef_tbl))
  cimat[,1] <- orig$Estimate - (qnorm(0.975)*orig$Std..Error)
  cimat[,2] <- orig$Estimate + (qnorm(0.975)*orig$Std..Error)
  
  cimat <- data.frame(cimat)
  rownames(cimat) <- rownames(coef_tbl)
  names(cimat) <- c("2.5 %", "97.5 %")
  
  coef_tbl$CI <- apply(cimat,1,conf)
  
  sformat <- paste("%.", pval.dec, "f", sep="")
  
  p_value2 <- sprintf(sformat, round(coef_tbl$p_value, pval.dec))
  if (pval.dec == 4) p_value2[coef_tbl$p_value < 0.0001] <- "&lt; 0.0001"
  if (pval.dec == 3) p_value2[coef_tbl$p_value < 0.001] <- "&lt; 0.001"
  if (pval.dec == 2) p_value2[coef_tbl$p_value < 0.01] <- "&lt; 0.01"
  if (pval.dec == 4) p_value2[coef_tbl$p_value > 0.9999] <- "&gt; 0.9999"
  if (pval.dec == 3) p_value2[coef_tbl$p_value > 0.999] <- "&gt; 0.999"
  if (pval.dec == 2) p_value2[coef_tbl$p_value > 0.99] <- "&gt; 0.99"
  coef_tbl$p_value <- p_value2
  
  # covs <- strsplit(as.character(fit$formula), "~")[[3]]
  # covs <- trim(unlist(strsplit(covs, "+", fixed=T)))
  
  coef_tbl <- coef_tbl[,c(estname, "CI", "p_value")]
  
  tbl <- NULL
  rgroup <- NULL
  
  rowlab <- NULL
  
  if (intercept == TRUE){
    tbl <- coef_tbl["(Intercept)",]
    if (overallp == TRUE) tbl$Overall_pvalue <- NA
    tbl$Variable <- "(Intercept)"
  }
  
  df <- data.frame(df)
  
  for (i in 1:length(covs)){
    
    if (class(df[,covs[i]]) == "numeric"){
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
    
    if (class(df[,covs[i]])  == "factor" | 
        class(df[,covs[i]]) == "character"){
      
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
    names(tbl) <- c("Variable", estname, "95% CI", "Wald p-value", "LR p-value")
  }
  if (overallp == FALSE){
    names(tbl) <- c("Variable", estname, "95% CI", "p-value")
  }
  if (printRMD == TRUE){
    
    if (overallp == TRUE){
      print(xtable(tbl, align="llccrr"), type='html', 
            include.rownames=F)
    }
    if (overallp == FALSE){
      print(xtable(tbl, align = "llccr"), type='html', 
            include.rownames=F)
    }
  }
  
  
  if (htmlTable == TRUE){
    
    final_html <- tbl
    
    ### stop htmlTable from treating everything as a factor
    for (i in 1:ncol(final_html)){
      final_html[,i] <- as.character(final_html[,i])
    }
    # ### remove blanks
    # final_html <- final_html[!is.na(final_html[,2]),]
    # ### get header rows
    # head <- rowlab
    # ### get non-header row
    # nohead <- rowlab == FALSE
    # ### indent non-header rows and remove *
    # final_html[nohead,"Variable"] <- paste("&nbsp; &nbsp; &nbsp;",
    #                                        substring(final_html[nohead,"Variable"], 3))
    # ### bold header rows
    # final_html[head,"Variable"] <- paste("<b>",
    #                                      final_html[head,"Variable"],
    #                                      "<b/>", sep="")
    head <- grepl("*", final_html[,"Variable"], fixed =T) == FALSE
    final_html[head,"Variable"] <- paste("<b>",
                                         final_html[head,"Variable"],
                                         "<b/>", sep="")
    final_html[,"Variable"] <- gsub("*", 
                                    "&nbsp; &nbsp; &nbsp;",
                                    final_html[,"Variable"],
                                    fixed = T)
    
    ### create htmlTable
    htmlver <- htmlTable(x = final_html[,2:ncol(final_html)],
                         rnames = final_html[,"Variable"],
                         rowlabel = paste("N =", summary(fit)$ngrps),
                         css.cell='border-collapse: collapse; padding: 4px;',
                         col.rgroup=rgroup)
    print(htmlver)
  }
  
  return(tbl)
  
}