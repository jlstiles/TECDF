

# not in use
gentmle_alt <- function(initdata, estimate_fun, update_fun, max_iter = 100,N=NULL, ...) {

    converge <- F
    n=nrow(initdata)
    # cat(sprintf('bw: %f\n',bw))
    eststep <- estimate_fun(initdata)

    initests <- eststep$ests
    if (is.null(N)) N=n
    order <- 1/N

    for (j in seq_len(max_iter)) {
      ED <- sapply(eststep$Dstar, mean)
      sig <- sapply(eststep$Dstar,sd)
      if (all(abs(ED) < sig*order)) {
        converge <- T
        updatestep = list()
        updatestep$coefs = NA
        break
      }
      updatestep <- update_fun(tmledata=eststep$tmledata,Dstar=eststep$Dstar)
      if (any(is.na(updatestep$tmledata$Qk))) {
        ED <- sapply(eststep$Dstar, mean)
        break}
      if (any(round(updatestep$tmledata$Qk,9)==1)|any(round(updatestep$tmledata$Qk,9)==0))
      {
        ED <- sapply(eststep$Dstar, mean)
        break}
      
      
      eststep <- estimate_fun(updatestep$tmledata,updatestep$coefs)

      ED <- sapply(eststep$Dstar, mean)
        # cat(sprintf('ED_psi=%e ED_sigma=%e psi=%f sigma2=%f\n coef_h=%f coef_Cy=%f
        # coef_Cg=%f\n',ED[1],ED[2],eststep$ests[1],sigma2=eststep$ests[2],updatestep$coefs[1],updatestep$coefs[2],updatestep$coefs[3]))

        if (all(abs(ED) < order)) {
            converge <- T
            break
        }



    }

    ED2 <- sapply(eststep$Dstar, function(x) mean(x^2))
    ED3 <- sapply(eststep$Dstar, function(x) mean(x^3))
    result <- list(initdata = initdata, tmledata = eststep$tmledata, initests = initests, tmleests = eststep$ests,
        steps = j, coefs = updatestep$coefs,Dstar = eststep$Dstar, ED = ED, ED2 = ED2, ED3 = ED3)

    class(result) <- "gentmle"

    return(result)
}


