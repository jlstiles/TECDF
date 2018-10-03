# not in use
gentmle_alt1 <- function(initdata, estimate_fun, update_fun, max_iter = 100, N=NULL,
                         t,h, k, simultaneous.inference = TRUE, ...) {
  
  # create the kernel according to specs
  R = k$range
  
  if (is.null(k$degree)) {
    kernel = list(kern = function(x, R, veck) 1/(2*R)*as.numeric(-R <= x & R >= x), 
    kern_cdf = function(x, R, veck) (1/(2*R))*as.numeric(x > -R)*(pmin(x ,R)+R))
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


