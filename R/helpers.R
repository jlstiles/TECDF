library(boot)
############## somewhat general TMLE framework takes initial data, an estimation function and an update
############## function can be used for arbitrary TMLEs

#' @export
ci_gentmle <- function(gentmle_obj, level = 0.95) {
  
  n <- nrow(gentmle_obj$initdata$Q)
  n_ests <- length(gentmle_obj$tmleests)
  if (gentmle_obj$simultaneous.inference == TRUE){
    check = apply(gentmle_obj$Dstar, 2, FUN = function(IC) {
      uu = length(unique(IC))
      if (uu==1) return(0) else return(1)
    })
    if (any(check==0)) {
      z <- qnorm((1 + level)/2)
    } else {
      S = cor(gentmle_obj$Dstar)
      Z = rmvnorm(1000000, rep(0,ncol(gentmle_obj$Dstar)), S)
      Z_abs = apply(Z,1,FUN = function(x) max(abs(x)))
      z = quantile(Z_abs, level)
    }
  } else {
    z <- qnorm((1 + level)/2)
  }
  plyr::ldply(seq_len(n_ests), function(i) {
    est <- gentmle_obj$tmleests[i]
    se <- sqrt(gentmle_obj$ED2[i])/sqrt(n)
    lower <- est - z * se
    upper <- est + z * se
    data.frame(est = est, se = se, lower = lower, upper = upper)
  })
}



print.gentmle <- function(gentmle_obj) {
  cat(sprintf("TMLE ran for %d step(s)\n", gentmle_obj$steps))
  EDtext <- sprintf("E[%s]=%1.2e", names(gentmle_obj$ED), gentmle_obj$ED)
  cat(sprintf("The mean of the IC is %s\n", paste(EDtext, collapse = ", ")))
  
  cat("\n\n")
  print(ci_gentmle(gentmle_obj))
  
  cat("\n")
  
}

# truncation function for Q so logistic regression doesn't break on Y close to 0 or 1
#' @export
truncate <- function(x, lower = 0.01, upper = 1 - lower) {
  pmin(pmax(x, lower), upper)
}

# function to estimate logistic parametric submodel and get updated estimate logistic
# fluctuation
#' @export
logit_fluctuate <- function(tmledata, flucmod, truncate = 0) {
  suppressWarnings({
    fluc <- glm(flucmod, data = tmledata, family = "binomial")
  })
  list(eps = coef(fluc))
}

#' @export
ind_choose = function(CIs, incre) {
  right = CIs[,3]
  left = CIs[,2]
  if (incre) {
    ind = order(right)[1]
  } else {
    ind = order(left, decreasing = TRUE)[1]
  }
  return(ind)
}

#' @export
ci_form = function(est, SE, z_alpha) {
  df = as.data.frame(cbind(est, est - z_alpha*SE, est + z_alpha*SE))
  colnames(df) = c("ests", "left", "right")
  return(df)
}
