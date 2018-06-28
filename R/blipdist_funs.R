#' @export
blipdist_update <- function(tmledata, Q.trunc = 1e-04) {
  eps_q <- 0
  # fluctuate Q
  deriv = colMeans(as.data.frame(tmledata$Dstar))
  norm = sqrt(sum(deriv^2))
  H = tmledata$HAW%*%deriv/norm
  H1 = as.vector(tmledata$H1W%*%deriv/norm)
  H0 = as.vector(tmledata$H0W%*%deriv/norm)
  # H = tmledata$HAW
  Qtrunc <- apply(tmledata$Q,2,truncate, Q.trunc)
  data <- data.frame(Y = tmledata$Y, H = H, offset =  qlogis(Qtrunc[,"QAW"]))
  data1 <- data0 <- data
  data1$H <- H1
  data1$offset <- qlogis(Qtrunc[,"Q1W"])
  data0$H <- H0
  data0$offset <- qlogis(Qtrunc[,"Q0W"])
  
  # qfluc <- glm(Y ~ 1, data = data, offset = offset, weights = H, family='binomial')
  if (all(H==0)) eps_q = 0 else {
    qfluc <- glm(Y ~ -1 + H, data = data, offset = offset, family='binomial')
    eps_q <- qfluc$coef
    QAW <- predict(qfluc, type = 'response')
    Q1W <- predict(qfluc, newdata = data1, type = 'response')
    Q0W <- predict(qfluc, newdata = data0, type = 'response')
    tmledata$Q = data.frame(QAW=QAW,Q0W=Q0W,Q1W=Q1W)
  }
  # qfluc <- glm(tmledata$Y ~ -1 + H + offset(qlogis(Qtrunc[,"QAW"])),family='binomial')

  # QAW <- as.vector(plogis(qlogis(Qtrunc[,"QAW"]) + H*eps_q))
  # Q1W <- as.vector(with(tmledata, plogis(qlogis(Qtrunc[,"Q1W"]) + eps_q*H1W%*%deriv/norm)))
  # Q0W <- as.vector(with(tmledata, plogis(qlogis(Qtrunc[,"Q0W"]) + eps_q*H0W%*%deriv/norm)))
  tmledata$coefs = c(eps_q)
  return(tmledata)
}

#' @export
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


#' @export
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


#' @export
blipdist_estimate2 <- function(tmledata, t, h, kernel) {
  nn <- length(tmledata$Y)
  B = tmledata$Q[,"Q1W"]-tmledata$Q[,"Q0W"]
  
  int = vapply(t, FUN = function(x0) {
      w = with(kernel, kern_cdf((B - x0)/h, R=R, veck=veck))
      return(w)
  } ,FUN.VALUE=rep(1,nn))
  tmledata$ests = apply(int, 2, mean)

  tmledata$HAW = vapply(t, FUN = function(x) {
    (1/h)*with(kernel, kern((B-x)/h, R=R, veck=veck))*with(tmledata, A/g1W-(1-A)/(1-g1W))
  } ,FUN.VALUE= rep(1,nn))
  tmledata$H0W = vapply(t, FUN = function(x) {
    (1/h)*with(kernel, kern((B-x)/h, R=R, veck=veck))*(-1/(1-tmledata$g1W))
  }, FUN.VALUE= rep(1,nn))
  tmledata$H1W = vapply(t, FUN = function(x) {
    (1/h)*with(kernel, kern((B-x)/h, R=R, veck=veck))*(1/tmledata$g1W)
  }, FUN.VALUE= rep(1,nn))
  
  tmledata$Dstar <- vapply(1:length(t),FUN = function(x) with(tmledata, HAW[,x]*(Y-Q[,"QAW"])+
    int[,x]-ests[x]), FUN.VALUE = rep(1,nn))
  return(tmledata)
}







