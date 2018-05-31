#' @export 
get.zscore = function(Dstar, alpha) {
  
  if (ncol(Dstar)==1) return(1.96544) else{
    sigma = cor(Dstar)
    means = rep(0,ncol(sigma))
    zs = rmvnorm(n=2000000,mean= means, sigma=sigma, method= "chol")
    zabs = apply(zs,1,FUN = function(x) max(abs(x)))
    zscore = quantile(zabs, probs = 1-alpha)
    return(zscore)}
}