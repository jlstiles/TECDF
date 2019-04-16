#' @title TECDF
#' @description computes kernel smoothed treatment effect or Treatment Effect CDF
#' @param initdata a list with elements Q, a data.frame with columns named QAW, Q1W, Q0W for initial 
#' predictions for the outcome, outcome under A=1 and under A=0, resp. A, a vector of binary 
#' treatment assignments and Y, the outcome and g1W, a vector of propensity scores. 
#' @param kernel, see make_kernel
#' @param blip, a vector of treatment effect value(s).
#' @param h, the bandwidth
#' @param max_iter, Maximum number of iteration steps
#' @param simultaneous.inference, do you want to compute simultaneous CI's (see ci_gentmle) 
#' @example /inst/quick_example.R  
#' @export
TECDF = function(initdata=initdata, kernel = kernel, TE, h, 
                   max_iter = 100, simultaneous.inference = TRUE)  {
  estimate_fun = TEdist_estimate
  update_fun = TEdist_update
    
    # create the kernel according to specs
    converge <- F
    n=length(initdata$Y)
    # cat(sprintf('bw: %f\n',bw))
    if (!is.null(kernel)) {
      eststep <- estimate_fun(tmledata=initdata, b = TE, h=h, kernel=kernel)
    } else {
      eststep <- estimate_fun(tmledata=initdata)
    }
    
    
    initests <- eststep$ests
    
    order <- 1/n
    
    for (j in seq_len(max_iter)) {
      #
      #         if (any(apply(eststep$HAW,2,FUN = function(x) all(x==0))==TRUE))
      #         {
      #         ED <- sapply(eststep$Dstar, mean)
      #         break}
      updatestep <- update_fun(tmledata=eststep)
      if (!is.null(kernel)) {
        eststep <- estimate_fun(tmledata=updatestep, b=TE, h=h, kernel=kernel)
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
    names(ED) = names(ED2) = names(eststep$ests) = paste0("F(", 1:length(ED),")")
    
    result <- list(initdata = initdata, Q = eststep$Q, initests = initests, tmleests = eststep$ests,
                   steps = j, coefs = updatestep$coefs,Dstar = eststep$Dstar, ED = ED, ED2 = ED2,
                   simultaneous.inference = simultaneous.inference)
    
    return(result)
  }


