
TECDF1 = function(initdata=initdata, kernel = kernel, TE, h, 
                   max_iter = 100, simultaneous.inference = TRUE)  {
  estimate_fun = TEdist_estimate
  update_fun = TEdist_update
  risk = c(NA)  
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
      Y_temp = eststep$Y
      Q_temp = eststep$Q$QAW
      loss = -Y_temp*log(Q_temp) - (1-Y_temp)*log(1-Q_temp)
      risk[j] = sum(loss)
      if (all(abs(ED) <= sigma*order)) {
        converge <- T
        break
      }
    }
    
    ED2 <- apply(eststep$Dstar, 2, FUN = function(x) mean(x^2))
    names(ED) = names(ED2) = names(eststep$ests) = paste0("F(", 1:length(ED),")")
    
    result <- list(initdata = initdata, Q = eststep$Q, initests = initests, tmleests = eststep$ests,
                   steps = j, coefs = updatestep$coefs,Dstar = eststep$Dstar, ED = ED, ED2 = ED2,
                   simultaneous.inference = simultaneous.inference, risk=risk)
    
    return(result)
  }


