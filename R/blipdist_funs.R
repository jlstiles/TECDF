
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



blipdist_estimate <- function(tmledata, b, h, kernel) {
  nn <- length(tmledata$Y) 
  B = tmledata$Q[,"Q1W"]-tmledata$Q[,"Q0W"]
  
  int = vapply(b, FUN = function(x0) {
      w = 1-with(kernel, kern_cdf((B - x0)/h, R=R, veck=veck))
      return(w)
  } ,FUN.VALUE=rep(1,nn))
  tmledata$ests = apply(int, 2, mean)

  tmledata$HAW = vapply(b, FUN = function(x) {
    (-1/h)*with(kernel, kern((B-x)/h, R=R, veck=veck))*with(tmledata, A/g1W-(1-A)/(1-g1W))
  } ,FUN.VALUE= rep(1,nn))
  tmledata$H0W = vapply(b, FUN = function(x) {
    (-1/h)*with(kernel, kern((B-x)/h, R=R, veck=veck))*(-1/(1-tmledata$g1W))
  }, FUN.VALUE= rep(1,nn))
  tmledata$H1W = vapply(b, FUN = function(x) {
    (-1/h)*with(kernel, kern((B-x)/h, R=R, veck=veck))*(1/tmledata$g1W)
  }, FUN.VALUE= rep(1,nn))
  
  tmledata$Dstar <- vapply(1:length(b),FUN = function(x) with(tmledata, HAW[,x]*(Y-Q[,"QAW"])+
    int[,x]-ests[x]), FUN.VALUE = rep(1,nn))
  return(tmledata)
}







