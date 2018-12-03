
get.zscore = function(Dstar, alpha) {
  
  if (ncol(Dstar)==1) return(1.96544) else{
    sigma = cor(Dstar)
    means = rep(0,ncol(sigma))
    zs = rmvnorm(n=2000000,mean= means, sigma=sigma, method= "chol")
    zabs = apply(zs,1,FUN = function(x) max(abs(x)))
    zscore = quantile(zabs, probs = 1-alpha)
    return(zscore)}
}


gendata.blip=function(n, d, g0, Q0){
  # d=1
  # n=10
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


get.truth = function(t, h, kernel, d, g0, Q0) {
  # create the kernel according to specs
  B = gendata.blip(3e6, d, g0, Q0)$blip
  vapply(t, FUN = function(x0) {
    t(vapply(h, FUN = function(bw) {
      int = with(kernel, kern_cdf(x=(B - x0)/bw, R = R, veck = veck))
      truth_h = mean(int)
      truth = mean(B>=x0)
      return(c(truth, truth_h))
    }, FUN.VALUE = c(1,1)))
  }, FUN.VALUE = matrix(rep(1, 2*length(h)), ncol=2))
}

gentmledata = function(n, d, g0, Q0, V, formu = NULL) {
  # n=100
  # d=1
  # V=10
  data = gendata.blip(n, d, g0, Q0)$df
  data0 = data
  data0$A = 0
  data1 = data
  data1$A = 1
  datag = data
  datag$Y = NULL
  
  if (is.null(formu)) {
    covs = colnames(data)[!colnames(data) %in% c("A","Y")]
    formuQ = formula(paste0("Y ~ ", paste0("A*(", paste(covs, "", collapse = "+"), ")")))
    formug = formula("A ~.")
  } else  {
    formuQ = formu$Q
    formug = formu$g
  }
  
  if (V == 1) {
    fitQ = glm(formuQ,data=data, family = "binomial")
    Q0W = predict(fitQ,newdata = data0, type = 'response')
    Q1W = predict(fitQ, newdata = data1, type = 'response')
    QAW = predict(fitQ, newdata = data, type= 'response')
    
    Q = cbind(QAW, Q0W, Q1W)
    datag = data
    datag$Y = NULL
    fitg = glm(formug, data=datag, family='binomial')
    g1W = predict(fitg, type = 'response')
    
    tmledata = list(Q=Q,Y=data$Y,A=data$A,g1W=g1W)
  } else {
    folds = make_folds(n=n, V=V)
    fold_preds = lapply(folds, FUN = function(fold) {
      fitQ = glm(formuQ,data=data[fold$training_set,], family = "binomial")
      Q0W = predict(fitQ,newdata = data0[fold$validation_set,], type = 'response')
      Q1W = predict(fitQ, newdata = data1[fold$validation_set,], type = 'response')
      QAW = predict(fitQ, newdata = data[fold$validation_set,], type= 'response')
      
      Q = cbind(QAW, Q0W, Q1W)
      fitg = glm(formug, data=datag[fold$training_set,], family='binomial')
      g1W = predict(fitg, newdata = datag[fold$validation_set, ], type = 'response')
      
      return(list(Q=Q,Y=data$Y[fold$validation_set],A=data$A[fold$validation_set],g1W=g1W))
    })
    tmledata = list(Q = do.call(rbind, lapply(fold_preds, FUN = function(x) x$Q)),
                    Y = unlist(lapply(fold_preds, FUN = function(x) x$Y)),
                    A = unlist(lapply(fold_preds, FUN = function(x) x$A)),
                    g1W = unlist(lapply(fold_preds, FUN = function(x) x$g1W)))
  }
  return(tmledata)
}

