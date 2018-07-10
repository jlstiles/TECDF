library(ggplot2)
library(cowplot)
library(parallel)
library(doSNOW)
library(cateSurvival)

sim_abbr = "lin"
# Define th DGP functions for SCM
g0 = function(W1) plogis(.2+.2*W1)
Q0 = function(A,W1) plogis(A + 2.5*A*W1 + W1)

# define kernel
k = list()
k$degree = NULL
k$range = 1

# getting two kinds of truths, plotting them
# make a grid of true parameters
true = gendata.blip(2000000,d = 1, g0, Q0)
hist(true$blip, breaks=200)
pscores = with(true$df, g0(W1))
hist(pscores, breaks = 200)
var(true$blip)
mean(true$blip)
m = min(true$blip)
M = max(true$blip)
m = -0.195
M = 0.319

# function start here
##
##
ii = 0
rr = 1
num_draws = 550
# for (rr in 5:5) {
  n=1000
  ii = ii + 1
  print(ii)
  t = seq(m, M, .01)
  bw = n^-.2
  
  truths = truth.get(t=t, h=bw, k, d=1, g0=g0, Q0=Q0)
  truth_h = truths$truth_h
  truth = truths$truth
  
  true_df = data.frame(truth = c(truth_h, truth), 
                       CATE = rep(t,2),
                       type = c(rep("smoothed", length(truth)),
                                rep("true", length(truth))))
  # gg_true = ggplot(true_df, aes(x=CATE, y=truth, color = type)) + geom_line()
  
  # setting up a simulation
  getres = function(n, t, h, k, truth, d = 1, g0, Q0, formu = NULL) {
    data=gentmledata_hal(n, d = 1, g0, Q0, V = 10, RCT = FALSE)
    reshal_simul = CATEsurv_plot(t = t, h = h, k = k, truth = truth, n = n, tmledata = data$tmledata)
    resglm_simul = CATEsurv_plot(t = t, h = h, k = k, truth = truth, n = n, tmledata = data$tmledata1)
    reshal = lapply(1:length(t), FUN = function(x) {
      res = CATEsurv_plot(t = t[x], h = h, k = k, truth = truth[x], n = n, tmledata = data$tmledata)
      return(res)
    })
    resglm = lapply(1:length(t), FUN = function(x) {
      res = CATEsurv_plot(t = t[x], h = h, k = k, truth = truth[x], n = n, tmledata = data$tmledata1)
      return(res)
    })
    risk = data$risk
    rownames(risk) = c("Q", "g")
    return(list(reshal_simul = reshal_simul$info, 
                resglm_simul = resglm_simul$info, 
                reshal = reshal,
                resglm = resglm,
                risk = risk, supnorm = data$supnorm
    ))
  }
  # tests = getres(n, t[blip], bw, k = k, 
  #        truth = truth_h[blip], d = 1, 
  #        g0 = g0, Q0 = Q0)
  
  if (n >= 1000) {
    cl = makeCluster(12, type = "SOCK")
  } else cl = makeCluster(detectCores(), type = "SOCK")
  
  registerDoSNOW(cl)
  blip = seq(6, 48, 6)
  B=num_draws
  
  allresults=foreach(i=1:B,
                     .packages=c("cateSurvival","mvtnorm","ggplot2", "cowplot", "hal9001")
                     ,.errorhandling='remove'
  )%dopar%
  {getres(n, t[blip], bw, k = k, 
          truth = truth_h[blip], d = 1, 
          g0 = g0, Q0 = Q0)}
  
  fname = paste0("unifCVhalglm_",sim_abbr,rr,"_",n,".RData")
  save(allresults, file = fname)
# }


