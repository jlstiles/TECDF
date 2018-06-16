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

#' @export 
gendata.blip=function(n, d, g0, Q0){
  var_names = paste0("W", 1:d)
  W_list = lapply(1:d, FUN = function(i) rnorm(n))
  names(W_list) = var_names
  pscores = do.call(g0, W_list)
  A=rbinom(n,1,pscores)
  AW_list = W_list
  AW_list$A = A
  Q = do.call(Q0, AW_list)
  Y=rbinom(n,1,Q)
  AW_list1 = AW_list0 = AW_list
  AW_list1$A = rep(1,n)
  AW_list0$A = rep(0,n)
  blip = do.call(Q0, AW_list1) - do.call(Q0, AW_list0)
  df = cbind(A = A, as.data.frame(W_list), Y = Y)
  return(list(df=df,blip=blip))
}

#' @export 
gentmledata = function(n, d, g0, Q0) {
  data = gendata.blip(n, d, g0, Q0)$df
  data0 = data
  data0$A = 0
  data1 = data
  data1$A = 1
  newX = rbind(data,data1,data0)
  
  fitQ = glm(Y~A*W1+A*W2+W3+W4,data=data, family = "binomial")
  Q0W = predict(fitQ,newdata = data0, type = 'response')
  Q1W = predict(fitQ, newdata = data1, type = 'response')
  QAW = predict(fitQ, newdata = data, type= 'response')
  
  Q = cbind(QAW, Q0W, Q1W)
  datag = data
  datag$Y = NULL
  fitg = glm(A~., data=datag, family='binomial')
  g1W = predict(fitg, type = 'response')
  
  tmledata = list(Q=Q,Y=data$Y,A=data$A,g1W=g1W)
  return(tmledata)
}
