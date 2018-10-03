library(boot)
############## somewhat general TMLE framework takes initial data, an estimation function and an update
############## function can be used for arbitrary TMLEs


#' @export 
make_kernel = function(degree, R){
  if (is.null(degree)) {
    kern = function(x, R, veck) 1/(2*R)*as.numeric(-R <= x & R >= x) 
    kern_cdf = function(x, R, veck) (1/(2*R))*as.numeric(x > -R)*(pmin(x ,R)+R)
    veck = 1
  } else {
    kk = degree/2-2
    area_row = vapply(0:(kk+2), FUN = function(i) 2*R^(2*i+1)/(2*i+1), FUN.VALUE = 1)
    zero_row = vapply(0:(kk+2), FUN = function(i) R^(2*i), FUN.VALUE = 1)
    deriv_row = c(0,vapply(0:(kk+1), FUN = function(i) 2*(i + 1)*R^(2*i+1), FUN.VALUE = 1))
    if (kk>0) {
      orth_rows = lapply(seq(0,max((2*kk-2),0),2), FUN = function(r) {
        vapply(0:(kk+2), FUN = function(i) 2*R^(2*i+3+r)/(2*i+3+r), FUN.VALUE = 1)
      })
      orth_rows = do.call(rbind, orth_rows) 
      mm = rbind(area_row, zero_row, deriv_row, orth_rows)
    } else mm = rbind(area_row, zero_row, deriv_row)
    
    mm_inv = solve(mm)
    veck = mm_inv %*% c(1, rep(0,kk+2))
    kern = function(x, R, veck) {
      ll = lapply(1:length(veck), FUN = function(c) veck[c]*x^(2*c-2))
      w = Reduce("+", ll)*(x > -R & x < R)
      return(w)
    }
    
    kern_cdf = function(x, R, veck) {
      u = pmin(x, R)
      ll = lapply(1:length(veck), FUN = function(c) veck[c]*(u^(2*c-1) + R^(2*c-1))/(2*c-1))
      w = Reduce("+", ll)*as.numeric(x > -R)
      return(w)
    }
  }
  
  return(list(veck = veck, R = R,  kern = kern, kern_cdf = kern_cdf))
}

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


#' @export
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
