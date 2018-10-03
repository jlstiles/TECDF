
# Not for use unless simultaneously estimating params with different bandwidths
blipdist_estimate3 <- function(tmledata, b, h, kernel) {
  nn <- length(tmledata$Y)
  B = tmledata$Q[,"Q1W"]-tmledata$Q[,"Q0W"]
  
  if (length(h)==1) {
    int = vapply(b, FUN = function(x0) {
      w = with(kernel, kern_cdf((B - x0)/h, R=R, veck=veck))
      return(w)
    } ,FUN.VALUE=rep(1,nn))
    tmledata$ests = apply(int, 2, mean)
    
    tmledata$HAW = vapply(b, FUN = function(x) {
      (1/h)*with(kernel, kern((B-x)/h, R=R, veck=veck))*with(tmledata, A/g1W-(1-A)/(1-g1W))
    } ,FUN.VALUE= rep(1,nn))
    tmledata$H0W = vapply(b, FUN = function(x) {
      (1/h)*with(kernel, kern((B-x)/h, R=R, veck=veck))*(-1/(1-tmledata$g1W))
    }, FUN.VALUE= rep(1,nn))
    tmledata$H1W = vapply(b, FUN = function(x) {
      (1/h)*with(kernel, kern((B-x)/h, R=R, veck=veck))*(1/tmledata$g1W)
    }, FUN.VALUE= rep(1,nn))
    
    tmledata$Dstar <- vapply(1:length(b),FUN = function(x) with(tmledata, HAW[,x]*(Y-Q[,"QAW"])+
                                                                  int[,x]-ests[x]), FUN.VALUE = rep(1,nn))
  } else {
    int = vapply(h, FUN = function(x) {
      w = with(kernel, kern_cdf((B - b)/x, R=R, veck=veck))
      return(w)
    } ,FUN.VALUE=rep(1,nn))
    tmledata$ests = apply(int, 2, mean)
    
    tmledata$HAW = vapply(h, FUN = function(x) {
      (1/x)*with(kernel, kern((B-b)/x, R=R, veck=veck))*with(tmledata, A/g1W-(1-A)/(1-g1W))
    } ,FUN.VALUE= rep(1,nn))
    tmledata$H0W = vapply(h, FUN = function(x) {
      (1/x)*with(kernel, kern((B-b)/x, R=R, veck=veck))*(-1/(1-tmledata$g1W))
    }, FUN.VALUE= rep(1,nn))
    tmledata$H1W = vapply(h, FUN = function(x) {
      (1/x)*with(kernel, kern((B-b)/x, R=R, veck=veck))*(1/tmledata$g1W)
    }, FUN.VALUE= rep(1,nn))
    
    tmledata$Dstar <- vapply(1:length(h),FUN = function(x) with(tmledata, HAW[,x]*(Y-Q[,"QAW"])+
                                                                  int[,x]-ests[x]), FUN.VALUE = rep(1,nn))
  }
  return(tmledata)
}


# not using this one
blipdist_estimate <- function(tmledata, t, h, kernel, ...) {
  nn <- length(tmledata$Y)
  B = tmledata$Q[,"Q1W"]-tmledata$Q[,"Q0W"]
  int = vapply(t, FUN = function(x0) {
    # lower = blip-h
    k = function(x) (1/h)*kernel((x-x0)/h)
    D1 = vapply(B, FUN = function(blip) {
      upper = min(blip, (x0+h))
      lower = x0-h
      if (kernel(0)==kernel(.1)) i = .5*(blip>(x0-h))*(upper-lower)/h else{
        upper = (min(blip, (x0+h))-x0)/h
        i = .75*(-upper^3/3+min(blip, (x0+h))/h-1/3-(x0-h)/h)*(blip>(x0-h))}
      return(i)
    }, FUN.VALUE=1)
    return(D1)
  },FUN.VALUE=rep(1,nn))
  tmledata$ests = apply(int, 2, mean)
  
  A = tmledata$A
  g1W = tmledata$g1W
  
  tmledata$HAW = vapply(t, FUN = function(x) {
    (1/h)*(kernel((B-x)/h))*(A/g1W-(1-A)/(1-g1W))
  } ,FUN.VALUE= rep(1,nn))
  tmledata$H0W = vapply(t, FUN = function(x) {
    (1/h)*(kernel((B-x)/h))*(-1/(1-g1W))
  }, FUN.VALUE= rep(1,nn))
  tmledata$H1W = vapply(t, FUN = function(x) {
    (1/h)*(kernel((B-x)/h))/g1W
  }, FUN.VALUE= rep(1,nn))
  HAW = tmledata$HAW
  psi = tmledata$ests
  tmledata$Dstar <- apply(HAW,2,FUN = function(x) with(tmledata, x*(Y-Q[,"QAW"])))+
    int-matrix(rep(psi,nn),byrow=TRUE,nrow=nn)
  return(tmledata)
}

# not usin this one
blipdist_estimate1 <- function(tmledata, t, h, kernel, ...) {
  nn <- length(tmledata$Y)
  B = tmledata$Q[,"Q1W"]-tmledata$Q[,"Q0W"]
  int = vapply(t, FUN = function(x0) {
    # lower = blip-h
    k = function(x) (1/h)*kernel((x-x0)/h)
    D1 = vapply(B, FUN = function(blip) {
      i = integrate(f = kernel, lower = -Inf, upper = (blip-x0)/h)$value
    }, FUN.VALUE=1)
    return(D1)
  },FUN.VALUE=rep(1,nn))
  tmledata$ests = apply(int, 2, mean)
  
  A = tmledata$A
  g1W = tmledata$g1W
  
  tmledata$HAW = vapply(t, FUN = function(x) {
    (1/h)*(kernel((B-x)/h))*(A/g1W-(1-A)/(1-g1W))
  } ,FUN.VALUE= rep(1,nn))
  tmledata$H0W = vapply(t, FUN = function(x) {
    (1/h)*(kernel((B-x)/h))*(-1/(1-g1W))
  }, FUN.VALUE= rep(1,nn))
  tmledata$H1W = vapply(t, FUN = function(x) {
    (1/h)*(kernel((B-x)/h))/g1W
  }, FUN.VALUE= rep(1,nn))
  HAW = tmledata$HAW
  psi = tmledata$ests
  tmledata$Dstar <- apply(HAW,2,FUN = function(x) with(tmledata, x*(Y-Q[,"QAW"])))+
    int-matrix(rep(psi,nn),byrow=TRUE,nrow=nn)
  return(tmledata)
}
