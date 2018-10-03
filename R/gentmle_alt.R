library(boot)
############## somewhat general TMLE framework takes initial data, an estimation function and an update
############## function can be used for arbitrary TMLEs

#' @title gentmle
#' @description General TMLE function that takes care of the bookkeeping of estimation and update steps.
#'
#' @param estimate_fun Function for estimation step
#' @param update_fun, Function for update step
#' @param max_iter, Maximum number of iteration steps
#' @param N, controls order criterion, currently not used
#' @param simultaneous.inference, passed to the 
#' @param ..., Extra arguments that can be passed to update_fun and estimate_fun
#' @export
gentmle_alt <- function(initdata, estimate_fun, update_fun, max_iter = 100, N=NULL
                        ,simultaneous.inference = TRUE, kernel = NULL, blip = NULL, h = NULL, ...) {
  
  # create the kernel according to specs
  converge <- F
  n=length(initdata$Y)
  # cat(sprintf('bw: %f\n',bw))
  if (!is.null(kernel)) {
    eststep <- estimate_fun(tmledata=initdata, b=blip, h=h, kernel=kernel)
  } else {
    eststep <- estimate_fun(tmledata=initdata)
  }
  
  
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
    if (!is.null(kernel)) {
      eststep <- estimate_fun(tmledata=updatestep, b=blip, h=h, kernel=kernel)
    } else {
      eststep <- estimate_fun(tmledata=updatestep)
    }
    
    
    ED <- apply(eststep$Dstar,2,mean)
    sigma <- apply(eststep$Dstar,2,sd)
    # cat(sprintf('ED_psi=%e ED_sigma=%e psi=%f sigma2=%f\n coef_h=%f coef_Cy=%f
    # coef_Cg=%f\n',ED[1],ED[2],eststep$ests[1],sigma2=eststep$ests[2],updatestep$coefs[1],updatestep$coefs[2],updatestep$coefs[3]))
    
    if (all(abs(ED) <= sigma*order)) {
      converge <- T
      break
    }
  }
  
  ED2 <- apply(eststep$Dstar, 2, FUN = function(x) mean(x^2))
  result <- list(initdata = initdata, Q = eststep$Q, initests = initests, tmleests = eststep$ests,
                 steps = j, coefs = updatestep$coefs,Dstar = eststep$Dstar, ED = ED, ED2 = ED2,
                 simultaneous.inference = simultaneous.inference)
  
  return(result)
}

