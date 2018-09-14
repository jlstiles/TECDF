

#' @export
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

