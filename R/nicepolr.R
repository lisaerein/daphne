
#' Lisa's Propoportional Odds Logistic Regression Table Function
#'
#' This function creates a nice looking regression table for a polr {MASS} model.
#' Input either a polr object or outcome variable, vector of covariates, method, and type of analysis (bivariate or multiple regression).
#' The function returns a dataframe with regression coefficients.
#' By default a table is also printed via htmlTable - handy for R markdown html reports.
#' @param df Dataframe [REQUIRED].
#' @param fit polr model object [fit or method/covs/out are REQUIRED].
#' @param method Method name in quotes ("logistic", "probit", "loglog", "cloglog", "cauchit") [REQUIRED if fit = NA].
#' @param covs Vector of covariates to include in model [fit or method/covs/out are REQUIRED].
#' @param out Outcome for regression model [fit or method/covs/out are REQUIRED].
#' @param regtype Should the covariates be run separately ("uni") or together in a multiple regression model ("multi") [REQUIRED if no fit].
#' @param overallp If TRUE, a likelihood ratio test pvalue (using drop1 Chisq tests) will be calculated for each variable. Default is TRUE.
#' @param est.dec Number of decimal places for estimates. Default is 2.
#' @param ci.dec Number of decimal places for 95% confidence interval. Default is 2.
#' @param pval.dec Number of decimal places for pvalues. Default is 3.
#' @param estname Option to override default estimate column name. Default is NA.
#' @param exp Option to exponentiate coefficients and CI's. Default is NA (estimates are only exponentiated for binomial and poisson method models by default).
#' @param htmlTable Whether to use htmlTable package to display table (instead of xtable). Default is TRUE
#' @param color Hex color to use for htmlTable output. Default is "#EEEEEE" (grey).
#' @param printRMD Whether to print resulting table to Rmd via xtable. Default is FALSE
#' @keywords polr table ordered proportional odds logistic regression coefficients
#' @importFrom xtable xtable
#' @importFrom htmlTable htmlTable
#' @importFrom MASS polr
#' @export
nicepolr    <- function(df,
                        fit = NA,
                        method = "logistic",
                        covs = NA,
                        out = NA,
                        regtype = "multi",
                        exp = TRUE,
                        estname = "Estimate",
                        labels = NA,
                        overallp = FALSE,
                        est.dec = 2,
                        ci.dec = 2,
                        pval.dec = 3,
                        color = "#EEEEEE",
                        printRMD = FALSE,
                        htmlTable = TRUE){

  ### overallp feature does not work at this time!
  overallp = FALSE

  ### run separate models for univariate and 1 model for multivarite analyses
  ### if model fit is provided, make table as is
  if (!is.na(fit[1])) regtype <- "multi"

  if (regtype == "uni"  ) {
    nmods <- length(covs)
    covlist <- as.list(covs)
  }
  if (regtype == "multi") {
    nmods <- 1
    covlist <- list(covs)
  }

  tbl_uni <- NULL
  rgroup_uni <- NULL

  for (j in 1:nmods){

    ### if needed, run the model
    if (is.na(fit[1]) & !is.na(covlist[[j]][1]) & !is.na(out) & !is.na(method)){

      form <- as.formula(paste(out, "~", paste(covlist[[j]], collapse = " + ")))
      fit0 <- polr(formula = form, data = df, method = method, Hess = TRUE)
       y <- model.frame(fit0)
       fit <- polr(y, Hess = TRUE)

    }

    ### get method/covariates/types from glm object (fit) attributes
    nobs_fit <- nobs(fit)
    covs <- attr(terms(fit), "term.labels")
    out <- names(attr(terms(fit), "dataClasses")[1])
    nlevels_out <- nlevels(df[,out])

    class <- attr(terms(fit), "dataClasses")[-1]
    type <- rep(NA, length(covs))
    type[class %in% c("numeric", "integer")] <- 1
    type[class %in% c("factor", "character")] <- 2

    if (is.na(method)) method <- summary(fit)$method

    ########## formatting table of model coefficients ###############

    ciformat <- paste("%.", ci.dec, "f", sep="")

    ### if no estimate name is specified pick a reasonable name for each situation

    if (is.na(estname) & method == "logistic" & exp == TRUE )  estname <- "OR"
    if (is.na(estname) & method != "logistic" & exp == TRUE )  estname <- "exp(Estimate)"
    if (is.na(estname)                        & exp == FALSE)  estname <- "Estimate"

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

    coef_tbl <- data.frame(summary(fit)$coef)

    coef_tbl$p_value <- pnorm(abs(coef_tbl[, "t.value"]), lower.tail = FALSE) * 2

    if (exp)          coef_tbl$Estimate <- exp(summary(fit)$coef[,"Value"])
    if (exp == FALSE) coef_tbl$Estimate <-    (summary(fit)$coef[,"Value"])
    coef_tbl$Estimate <- sprintf(esformat, round(coef_tbl$Estimate, est.dec))

    names(coef_tbl)[grepl("Est", names(coef_tbl))] <- estname

    ### remove rows with ordered intercepts: last nlevels - 1 rows of table
    coef_tbl <- coef_tbl[1:(nrow(coef_tbl)-nlevels_out+1),]

    cimat <- data.frame(confint(fit))
    if (nrow(coef_tbl) == 1) coef_tbl$CI <- conf(t(cimat))
    if (nrow(coef_tbl) >  1) coef_tbl$CI <- apply(cimat,1,conf)

    sformat <- paste("%.", pval.dec, "f", sep="")

    p_value2 <- sprintf(sformat, round(coef_tbl$p_value, pval.dec))
    if (pval.dec == 4) p_value2[coef_tbl$p_value < 0.0001] <- "< 0.0001"
    if (pval.dec == 3) p_value2[coef_tbl$p_value < 0.001]  <- "< 0.001"
    if (pval.dec == 2) p_value2[coef_tbl$p_value < 0.01]   <- "< 0.01"
    if (pval.dec == 4) p_value2[coef_tbl$p_value > 0.9999] <- "> 0.9999"
    if (pval.dec == 3) p_value2[coef_tbl$p_value > 0.999]  <- "> 0.999"
    if (pval.dec == 2) p_value2[coef_tbl$p_value > 0.99]   <- "> 0.99"
    if (htmlTable){
      if (pval.dec == 4) p_value2[coef_tbl$p_value < 0.0001] <- "&lt; 0.0001"
      if (pval.dec == 3) p_value2[coef_tbl$p_value < 0.001]  <- "&lt; 0.001"
      if (pval.dec == 2) p_value2[coef_tbl$p_value < 0.01]   <- "&lt; 0.01"
      if (pval.dec == 4) p_value2[coef_tbl$p_value > 0.9999] <- "&gt; 0.9999"
      if (pval.dec == 3) p_value2[coef_tbl$p_value > 0.999]  <- "&gt; 0.999"
      if (pval.dec == 2) p_value2[coef_tbl$p_value > 0.99]   <- "&gt; 0.99"
    }
    coef_tbl$p_value <- p_value2

    coef_tbl <- coef_tbl[,c(estname, "CI", "p_value")]

    tbl <- NULL
    rgroup <- NULL

    rowlab <- NULL

    ### create variable for coefficient table line number
    line <- 0

    for (i in 1:length(covs)){

      if (type[i] == 1){
        ngroups <- 1
        line <- line + ngroups
        tmp <- coef_tbl[line,]
        if (overallp == TRUE) {
          op <- drop1(fit, test = "Chisq")[covs[i],"Pr(>Chi)"]
          op2 <- sprintf(sformat, round(op, pval.dec))
          if (pval.dec == 4) op2[op < 0.0001] <- "< 0.0001"
          if (pval.dec == 3) op2[op < 0.001]  <- "< 0.001"
          if (pval.dec == 2) op2[op < 0.01]   <- "< 0.01"
          if (pval.dec == 4) op2[op > 0.9999] <- "> 0.9999"
          if (pval.dec == 3) op2[op > 0.999]  <- "> 0.999"
          if (pval.dec == 2) op2[op > 0.99]   <- "> 0.99"
          if (htmlTable){
            if (pval.dec == 4) op2[op < 0.0001] <- "&lt; 0.0001"
            if (pval.dec == 3) op2[op < 0.001]  <- "&lt; 0.001"
            if (pval.dec == 2) op2[op < 0.01]   <- "&lt; 0.01"
            if (pval.dec == 4) op2[op > 0.9999] <- "&gt; 0.9999"
            if (pval.dec == 3) op2[op > 0.999]  <- "&gt; 0.999"
            if (pval.dec == 2) op2[op > 0.99]   <- "&gt; 0.99"
          }
          tmp$Overall_pvalue <- op2
        }

        if (is.na(labels[1])) tmp$Variable <- covs[i]
        if (!is.na(labels[1])) tmp$Variable <- labels[i]
        if (regtype == "uni") tmp$Variable <- paste(tmp$Variable, " (N = ", nobs_fit, ")", sep="")

        tbl <- rbind(tbl, tmp)

        rowlab <- c(rowlab, TRUE, rep(FALSE, (nrow(tmp)-1)))

      }

      if (type[i] == 2){

        df[,covs[i]] <- as.factor(df[,covs[i]])
        ngroups <- nlevels(df[,covs[i]])
        tmp <- coef_tbl[(line+1):(line+ngroups-1),]
        line <- line + ngroups - 1
        if (overallp == TRUE) tmp$Overall_pvalue <- NA
        title <- data.frame(tmp[1,])
        title[1,] <- NA
        if (overallp == TRUE) {
          op <- drop1(fit, test = "Chisq")[covs[i],"Pr(>Chi)"]
          op2 <- sprintf(sformat, round(op, pval.dec))

          if (pval.dec == 4) op2[op < 0.0001] <- "< 0.0001"
          if (pval.dec == 3) op2[op < 0.001]  <- "< 0.001"
          if (pval.dec == 2) op2[op < 0.01]   <- "< 0.01"
          if (pval.dec == 4) op2[op > 0.9999] <- "> 0.9999"
          if (pval.dec == 3) op2[op > 0.999]  <- "> 0.999"
          if (pval.dec == 2) op2[op > 0.99]   <- "> 0.99"
          if (htmlTable){
            if (pval.dec == 4) op2[op < 0.0001] <- "&lt; 0.0001"
            if (pval.dec == 3) op2[op < 0.001]  <- "&lt; 0.001"
            if (pval.dec == 2) op2[op < 0.01]   <- "&lt; 0.01"
            if (pval.dec == 4) op2[op > 0.9999] <- "&gt; 0.9999"
            if (pval.dec == 3) op2[op > 0.999]  <- "&gt; 0.999"
            if (pval.dec == 2) op2[op > 0.99]   <- "&gt; 0.99"
          }
          title$Overall_pvalue <- op2
        }
        if (is.na(labels[1])) title$Variable <- covs[i]
        if (!is.na(labels[1])) title$Variable <- labels[i]
        if (regtype == "uni") title$Variable <- paste(title$Variable, " (N = ", nobs_fit, ")", sep="")

        tmp$Variable <- paste("*", levels(df[,covs[i]])[2:nlevels(df[,covs[i]])], "vs.",
                              levels(df[,covs[i]])[1])

        tbl <- rbind(tbl, title, tmp)

        rowlab <- c(rowlab, TRUE, rep(FALSE, (nrow(tmp)-1)))

      }

      ### set colors for htmlTable striping
      if (i %% 2 == 0) rgroup <- c(rgroup, rep("none", ngroups))
      if (i %% 2 != 0) rgroup <- c(rgroup, rep(color, ngroups))
    }

    tbl <- tbl[,c(ncol(tbl), 2:ncol(tbl)-1)]

    tbl_uni <- rbind(tbl_uni, tbl)
    fit <- NA

    ### set colors for htmlTable striping
    if (j %% 2 == 0) rgroup_uni <- c(rgroup_uni, rep("none", nrow(tbl)))
    if (j %% 2 != 0) rgroup_uni <- c(rgroup_uni, rep(color, nrow(tbl)))

  }
  tbl <- tbl_uni

  if (overallp == TRUE){
    names(tbl) <- c("Variable", estname, "95% CI", "Wald p-value", "LR p-value")
  }
  if (overallp == FALSE){
    names(tbl) <- c("Variable", estname, "95% CI", "p-value")
  }
  if (printRMD){
    if (overallp == TRUE){
      print(xtable(tbl, align="llccrr"), type='html', include.rownames=F)
    }
    if (overallp == FALSE){
      print(xtable(tbl, align = "llccr"), type='html', include.rownames=F)
    }
  }


  if (htmlTable){
    final_html <- tbl

    ### stop htmlTable from treating everything as a factor
    for (i in 1:ncol(final_html)){
      final_html[,i] <- as.character(final_html[,i])
    }
    head <- grepl("*", final_html[,"Variable"], fixed =T) == FALSE
    final_html[head,"Variable"] <- paste("<b>",
                                         final_html[head,"Variable"],
                                         "<b/>", sep="")
    final_html[,"Variable"] <- gsub("*",
                                    "&nbsp; &nbsp; &nbsp;",
                                    final_html[,"Variable"],
                                    fixed = T)

    ### create htmlTable
    if (regtype == "uni"  ) rgroup <- rgroup_uni
    if (regtype == "uni"  ) rowlab <- ""
    if (regtype == "multi") rowlab <- paste("N =", nobs_fit)

    htmlver <- htmlTable(x = final_html[,2:ncol(final_html)],
                         rnames = final_html[,"Variable"],
                         rowlabel = rowlab,
                         css.cell='border-collapse: collapse; padding: 4px;',
                         col.rgroup=rgroup)
    print(htmlver)
  }

  return(tbl)

}
