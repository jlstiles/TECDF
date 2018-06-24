library(ggplot2)
library(cowplot)
library(parallel)
library(doSNOW)
library(cateSurvival)

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
var(true$blip)
mean(true$blip)
m = min(true$blip)
M = max(true$blip)

# function start here
##
##
ii = 0
for (n in c(10000, 25000, 50000)) {
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
    tmledata=gentmledata(n, d = 1, g0, Q0, V = 10, formu = NULL)
    res = CATEsurv_plot(t = t, h = h, k = k, truth = truth, n = n, tmledata = tmledata)
    return(res$info)
  }

  if (n >= 10000) {
    cl = makeCluster(8, type = "SOCK")
  } else cl = makeCluster(24, type = "SOCK")
  
  registerDoSNOW(cl)
  blip = seq(13,28,5)
  B=5000
  
  allresults=foreach(i=1:B,
                     .packages=c("cateSurvival","mvtnorm","ggplot2", "cowplot")
                     ,.errorhandling='remove'
  )%dopar%
  {getres(n, t[blip], bw, k = k, 
          truth = truth_h[blip], d = 1, 
          g0 = g0, Q0 = Q0)}
  
  res = data.matrix(do.call(rbind, allresults))
  res = as.data.frame(res)
  
  B = length(allresults)
  L = length(blip)
  base_seq = seq(1,L*B,L)
  rows_res = unlist(lapply(1:length(blip), FUN = function(x) base_seq+x-1))
  res = res[rows_res,]
  
  cover = unlist(lapply(1:B, FUN = function(x) {
    cover = all(allresults[[x]]$cover==1)
    return(cover)
  }))
  
  coverage = mean(cover)
  
  assign(paste0("res", n, "unif_simul"), list(coverage = coverage, 
                                              B = B, 
                                              h = bw, 
                                              res = res,
                                              blip = t[blip],
                                              true_df = true_df,
                                              plot_true = gg_true))

  if (n >= 10000) {
    cl = makeCluster(8, type = "SOCK")
  } else cl = makeCluster(24, type = "SOCK")
  
  registerDoSNOW(cl)
  blip = seq(13,28,5)
  B=5000
  
  for (b in blip) {
    # b = blip[1]
    allresults=foreach(i=1:B,
                       .packages=c("cateSurvival","mvtnorm","ggplot2", "cowplot")
                       ,.errorhandling='remove'
    )%dopar%
    {getres(n, t[b], bw, k = k, 
            truth = truth_h[b], d = 1, 
            g0 = g0, Q0 = Q0)}
    
    res = data.matrix(do.call(rbind, allresults))
    res = as.data.frame(res)
    
    B = length(allresults)
    cover = unlist(lapply(1:B, FUN = function(x) {
      cover = all(allresults[[x]]$cover==1)
      return(cover)
    }))
    
    coverage = mean(cover)
    
    assign(paste0("res", n, "unif", b), list(coverage = coverage, 
                                                   B = B, 
                                                   h = bw, 
                                                   res = res,
                                                   blip = t[b]))

  }
  
  if (n == 5000) {
    save(res100unif13, res100unif18, res100unif23, res100unif28, 
         res250unif13, res250unif18, res250unif23, res250unif28,
         res500unif13, res500unif18, res500unif23, res500unif28,
         res1000unif13, res1000unif18, res1000unif23, res1000unif28,
         res2500unif13, res2500unif18, res2500unif23, res2500unif28,
         res5000unif13, res5000unif18, res5000unif23, res5000unif28,
         res100unif_simul, res250unif_simul, res500unif_simul, res1000unif_simul,
         res2500unif_simul, res5000unif_simul,
         g0, Q0, file = "kernel_sim_unifallCV.RData")
    rm(res100unif13, res100unif18, res100unif23, res100unif28, 
       res250unif13, res250unif18, res250unif23, res250unif28,
       res500unif13, res500unif18, res500unif23, res500unif28,
       res1000unif13, res1000unif18, res1000unif23, res1000unif28,
       res2500unif13, res2500unif18, res2500unif23, res2500unif28,
       res5000unif13, res5000unif18, res5000unif23, res5000unif28,
       res100unif_simul, res250unif_simul, res500unif_simul, res1000unif_simul,
       res2500unif_simul, res5000unif_simul)
  } 
  
  if (n == 50000) {
    save(res10000unif13, res10000unif18, res10000unif23, res10000unif28,
         res25000unif13, res25000unif18, res25000unif23, res25000unif28,
         res50000unif13, res50000unif18, res50000unif23, res50000unif28,
         res10000unif_simul, res25000unif_simul, 
         res50000unif_simul, g0, Q0, file = "kernel_sim_unifallCV1.RData")
    }
}




