gendataTE = function (n, g0, Q0) {
  W1 = runif(n, -3, 3)
  W2 = rnorm(n)
  W3 = runif(n)
  W4 = rnorm(n)
  A = rbinom(n, 1, g0(W1, W2, W3, W4))
  Y = rbinom(n, 1, Q0(A, W1, W2, W3, W4))
  Q1W = Q0(rep(1,n), W1, W2, W3, W4)
  Q0W = Q0(rep(0,n), W1, W2, W3, W4)
  data.frame(A, W1, W2, W3, W4, Y, TE = Q1W-Q0W)
}
pop = gendataTE(1e6, g0 = g0_linear, Q0 = Q0_trig)
var(pop$TE)

lapply(seq(-1,1,.1), function(x) {
  mean(pop$TE<=x)
})

TEvals = seq(-1,1,.1)[8:15]

# make polynomial kernel of order 6. Note, you can only input even degrees for the kernel which 
# will be the highest degree
k=make_kernel(order=2,R=2)
k

simTECDF = function(n) {
  data = gendata(n, g0 = g0_linear, Q0 = Q0_trig)
  h = 1/n^.33
  data1 = data0 = data
  data1$A = 1
  data0$A = 0
  fitQ = glm(Y~A*(W1+W2+W3+W4), data = data, family = binomial)
  QAW = predict(fitQ, type = "response")
  Q1W = predict(fitQ, newdata = data1, type = "response")
  Q0W = predict(fitQ, newdata = data0, type = "response")
  fitg = glm(A~W1+W2+W3+W4, data = data, family = binomial)
  g1W = predict(fitg, type = 'response')
  
  initdata = list(Q = data.frame(QAW,Q1W,Q0W), g1W=g1W, A=data$A, Y = data$Y)
  est.info = TECDF(initdata = initdata, kernel = k, TE = TEvals, h = h, 
                   max_iter = 1000, simultaneous.inference = TRUE)
  
  
  plot(TEvals, est.info$tmleests)
  # tmle estimates
  est.info$tmleests
  # steps to convergence
  est.info$steps
  # mean of the influence curves
  est.info$ED
  
  return(list(initdata=initdata, est.info = est.info, 
              converge = est.info$converge))
}

res = lapply(1:100, function(x) simTECDF(500))
converge = which(!unlist(lapply(res, function(L) L$converge)))
converge[19]

mono = unlist(lapply(res, function(L) {
  risk = L$est.info$risk
  mono = !any(unlist(lapply(2:(L$est.info$steps-1), function(i) {
    risk[i]>risk[i-1]
  })))
  mono
}))

which(!mono)

bads = lapply(res[converge], function(L) {
  df = data.frame(step = 1:L$est.info$steps,risk = L$est.info$risk)
  p = ggplot(df, aes(x=step,y=risk))+geom_point()
  return(list(p = p, ED = L$est.info$ED))
  })
badbads = lapply(lapply(res[converge], function(L) L$est.info$ED), 
       function(L) any(abs(L)>.01))

converge[which(unlist(badbads))]
res[[78]]$initdata
h = 1/500^.33
badres = TECDF(initdata = res[[78]]$initdata, kernel = k, TE = TEvals, h = h, 
      max_iter = 10000, simultaneous.inference = TRUE)

badres$steps
plot(1:badres$steps, badres$risk)
# for simultaneous inference, default set to 5%
ci = ci_gentmle(est.info, level = 0.95)
ci
# number of se's used for simultaneous inference at type I error rate of 5%
(ci[1,4]-ci[1,1])/ci[1,2]
