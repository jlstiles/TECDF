library(ggplot2)
library(cowplot)
library(parallel)
library(doSNOW)
library(cateSurvival)

# Define th DGP functions for SCM
g0 = function(W1) plogis(-.1-.5*sin(W1) - .4*(abs(W1)>1)*W1^2)
Q0 = function(A,W1) plogis(.3*A + 5*A*sin(W1)^2 - A*cos(W1))

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

# function start here
##
##
ii = 0
num_draws = 24
for (n in c(100, 250, 500, 1000, 2500, 5000)) {
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
  gg_true = ggplot(true_df, aes(x=CATE, y=truth, color = type)) + geom_line()

  # setting up a simulation
  getres = function(n, t, h, k, truth, d = 1, g0, Q0, formu = NULL) {
    data=gentmledata_hal(n, d = 1, g0, Q0, V = 10, RCT = FALSE)
    reshal = CATEsurv_plot(t = t, h = h, k = k, truth = truth, n = n, tmledata = data$tmledata)
    res = CATEsurv_plot(t = t, h = h, k = k, truth = truth, n = n, tmledata = data$tmledata1)
    return(list(reshal = reshal$info, res = res$info))
  }

  if (n >= 10000) {
    cl = makeCluster(8, type = "SOCK")
  } else cl = makeCluster(detectCores(), type = "SOCK")
  
  registerDoSNOW(cl)
  blip = seq(10,49,13)
  B=num_draws
  
  allresults=foreach(i=1:B,
                     .packages=c("cateSurvival","mvtnorm","ggplot2", "cowplot", "hal9001")
                     ,.errorhandling='remove'
  )%dopar%
  {getres(n, t[blip], bw, k = k, 
          truth = truth_h[blip], d = 1, 
          g0 = g0, Q0 = Q0)}
  
  allresults_hal = lapply(allresults, FUN = function(x) x$reshal)
  allresults_glm = lapply(allresults, FUN = function(x) x$res)
  
  res_hal = data.matrix(do.call(rbind, allresults_hal))
  res_hal = as.data.frame(res_hal)
  res_glm = data.matrix(do.call(rbind, allresults_glm))
  res_glm = as.data.frame(res_glm)
  
  B = length(allresults)
  L = length(blip)
  base_seq = seq(1,L*B,L)
  rows_res = unlist(lapply(1:length(blip), FUN = function(x) base_seq+x-1))
  res_hal = res_hal[rows_res,]
  res_glm = res_glm[rows_res,]
  
  cover_hal = unlist(lapply(1:B, FUN = function(x) {
    cover_hal = all(allresults_hal[[x]]$cover==1)
    return(cover_hal)
  }))
  
  cover_glm = unlist(lapply(1:B, FUN = function(x) {
    cover_glm = all(allresults_glm[[x]]$cover==1)
    return(cover_glm)
  }))
  
  coverage_hal = mean(cover_hal)
  coverage_glm = mean(cover_glm)
  
  assign(paste0("res", n, "unif_simulhal"), list(coverage = coverage_hal, 
                                              B = B, 
                                              h = bw, 
                                              res = res_hal,
                                              blip = t[blip],
                                              true_df = true_df,
                                              plot_true = gg_true))
  
  assign(paste0("res", n, "unif_simulglm"), list(coverage = coverage_glm, 
                                                 B = B, 
                                                 h = bw, 
                                                 res = res_glm,
                                                 blip = t[blip],
                                                 true_df = true_df,
                                                 plot_true = gg_true))

  if (n >= 10000) {
    cl = makeCluster(8, type = "SOCK")
  } else cl = makeCluster(detectCores(), type = "SOCK")
  
  registerDoSNOW(cl)
  blip = seq(10,49,13)
  B=num_draws
  
  for (b in blip) {
    # b = blip[1]
    allresults=foreach(i=1:B,
                       .packages=c("cateSurvival","mvtnorm","ggplot2", "cowplot", "hal9001")
                       ,.errorhandling='remove'
    )%dopar%
    {getres(n, t[b], bw, k = k, 
            truth = truth_h[b], d = 1, 
            g0 = g0, Q0 = Q0)}
    
    allresults_hal = lapply(allresults, FUN = function(x) x$reshal)
    allresults_glm = lapply(allresults, FUN = function(x) x$res)
    
    res_hal = data.matrix(do.call(rbind, allresults_hal))
    res_hal = as.data.frame(res_hal)
    res_glm = data.matrix(do.call(rbind, allresults_glm))
    res_glm = as.data.frame(res_glm)
    
    B = length(allresults)
    cover_hal = unlist(lapply(1:B, FUN = function(x) {
      cover = all(allresults_hal[[x]]$cover==1)
      return(cover)
    }))
    
    coverage_hal = mean(cover_hal)
    
    cover_glm = unlist(lapply(1:B, FUN = function(x) {
      cover = all(allresults_glm[[x]]$cover==1)
      return(cover)
    }))
    
    coverage_glm = mean(cover_glm)
    
    assign(paste0("res", n, "unifhal", b), list(coverage = coverage_hal, 
                                                   B = B, 
                                                   h = bw, 
                                                   res = res,
                                                   blip = t[b]))
    
    assign(paste0("res", n, "unifglm", b), list(coverage = coverage_hal, 
                                             B = B, 
                                             h = bw, 
                                             res = res,
                                             blip = t[b]))

  }
  
  if (n == 5000) {
    save(res100unifhal13, res100unifhal18, res100unifhal23, res100unifhal28, 
         res250unifhal13, res250unifhal18, res250unifhal23, res250unifhal28,
         res500unifhal13, res500unifhal18, res500unifhal23, res500unifhal28,
         res1000unifhal13, res1000unifhal18, res1000unifhal23, res1000unifhal28,
         res2500unifhal13, res2500unifhal18, res2500unifhal23, res2500unifhal28,
         res5000unifhal13, res5000unifhal18, res5000unifhal23, res5000unifhal28,
         res100unif_simulhal, res250unif_simulhal, res500unif_simulhal, res1000unif_simulhal,
         res2500unif_simulhal, res5000unif_simulhal,
         res250unifglm13, res250unifglm18, res250unifglm23, res250unifglm28,
         res500unifglm13, res500unifglm18, res500unifglm23, res500unifglm28,
         res1000unifglm13, res1000unifglm18, res1000unifglm23, res1000unifglm28,
         res2500unifglm13, res2500unifglm18, res2500unifglm23, res2500unifglm28,
         res5000unifglm13, res5000unifglm18, res5000unifglm23, res5000unifglm28,
         res100unif_simulglm, res250unif_simulglm, res500unif_simulglm, res1000unif_simulglm,
         res2500unif_simulglm, res5000unif_simulglm,
         g0, Q0, file = "kernel_sim_unifallCVhalglm.RData")
  } 
  
  if (n == 50000) {
    save(res10000unif13, res10000unif18, res10000unif23, res10000unif28,
         res25000unif13, res25000unif18, res25000unif23, res25000unif28,
         res50000unif13, res50000unif18, res50000unif23, res50000unif28,
         res10000unif_simul, res25000unif_simul, 
         res50000unif_simul, g0, Q0, file = "kernel_sim_unifallCVhalglm1.RData")
    }
}

