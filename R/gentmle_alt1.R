library(boot)
############## somewhat general TMLE framework takes initial data, an estimation function and an update
############## function can be used for arbitrary TMLEs

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


#' @title gentmle
#' @description General TMLE function that takes care of the bookkeeping of estimation and update steps.
#'
#' @param estimate_fun Function for estimation step
#' @param update_fun, Function for update step
#' @param max_iter, Maximum number of iteration steps
#' @param ..., Extra arguments that can be passed to update_fun and estimate_fun
#'
#' @export
gentmle_alt1 <- function(initdata, estimate_fun, update_fun, max_iter = 100, N=NULL,
                         t,h, k, kernel_cdf=NULL, ...) {
  
  # create the kernel according to specs
  R = k$range
  
  if (is.null(k$degree)) {
    kernel = list(kern = function(x, R, veck) .5*as.numeric(-1<=x&1>=x), 
    kern_cdf = function(x, R, veck) (1/(2*R))*as.numeric(x > -R)*(pmin(x ,R) + R))
    veck = 1
  } else {
    deg = k$degree
    mm = vapply(seq(1,deg,2), FUN = function(r) {
      vapply(seq(r,(r+deg+1),2), FUN = function(x) 2*R^x, FUN.VALUE = 1)/seq(r,(r+deg+1),2)
    }, FUN.VALUE = rep(1,(deg+3)/2))
    
    mm = cbind(mm, vapply(seq(0,deg+1,2), FUN = function(x) R^x, FUN.VALUE = 1)) 
    mm = t(mm)
    mm_inv = solve(mm)
    veck = mm_inv %*% c(1,rep(0, (deg+1)/2))
    kernel = list(kern = function(x, R, veck) {
      ll = lapply(1:length(veck), FUN = function(c) veck[c]*x^(2*c-2))
      w = Reduce("+", ll)*(x > -R & x < R)
      return(w)
    }, 
    kern_cdf = function(x, R, veck) {
      u = pmin(x, R)
      ll = lapply(1:length(veck), FUN = function(c) veck[c]*(u^(2*c-1) + R^(2*c-1))/(2*c-1))
      w = Reduce("+", ll)*as.numeric(x > -R)
      return(w)
    })
  }
  kernel$R = R
  kernel$veck = veck
  converge <- F
  n=length(initdata$Y)
  # cat(sprintf('bw: %f\n',bw))
  eststep <- estimate_fun(tmledata=initdata, t=t, h=h, kernel=kernel)

  
  initests <- eststep$ests
  if (is.null(N)) N=n
  order <- 1/N
  
  for (j in seq_len(max_iter)) {
    #
    #         if (any(apply(eststep$HAW,2,FUN = function(x) all(x==0))==TRUE))
    #         {
    #         ED <- sapply(eststep$Dstar, mean)
    #         break}
    updatestep <- update_fun(tmledata=eststep)
    eststep <- estimate_fun(tmledata=updatestep, t=t, h=h, kernel=kernel)
    
    
    ED <- apply(eststep$Dstar,2,mean)
    sigma <- apply(eststep$Dstar,2,sd)
    # cat(sprintf('ED_psi=%e ED_sigma=%e psi=%f sigma2=%f\n coef_h=%f coef_Cy=%f
    # coef_Cg=%f\n',ED[1],ED[2],eststep$ests[1],sigma2=eststep$ests[2],updatestep$coefs[1],updatestep$coefs[2],updatestep$coefs[3]))
    
    if (all(abs(ED) < sigma*order)) {
      converge <- T
      break
    }
  }
  
  ED2 <- apply(eststep$Dstar, 2, FUN = function(x) mean(x^2))
  result <- list(initdata = initdata, Q = eststep$Q, initests = initests, tmleests = eststep$ests,
                 steps = j, coefs = updatestep$coefs,Dstar = eststep$Dstar, ED = ED, ED2 = ED2)
  
  return(result)
}


#' @export
ci_gentmle <- function(gentmle_obj, level = 0.95) {

    n <- nrow(gentmle_obj$initdata)
    n_ests <- length(gentmle_obj$tmleests)
    ldply(seq_len(n_ests), function(i) {

        est <- gentmle_obj$tmleests[i]
        sd <- sqrt(gentmle_obj$ED2[i])/sqrt(n)
        z <- (1 + level)/2
        lower <- est - qnorm(z) * sd
        upper <- est + qnorm(z) * sd
        data.frame(parameter = names(est), est = est, sd = sd, lower = lower, upper = upper)
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
