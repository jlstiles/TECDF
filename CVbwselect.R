library(ggplot2)
library(cowplot)
library(parallel)
library(doSNOW)
library(cateSurvival)
library(Iso)
library(foreach)

# Define th DGP functions for SCM
g0 = function(W1) plogis(.2 + .2*W1)
Q0 = function(A,W1) plogis(A + 2.5*A*W1 + W1)

# define kernel

kernel4 = make_kernel(degree = 4, R = 2)
kernel6 = make_kernel(degree = 4, R = 2)
kernel8 = make_kernel(degree = 4, R = 2)
kernel10 = make_kernel(degree = 4, R = 2)
kernel12 = make_kernel(degree = 4, R = 2)

# getting two kinds of truths, plotting them
# make a grid of true parameters
true = gendata.blip(2000000,d = 1, g0, Q0)
hist(true$blip, breaks=200)
pscores = with(true$df, g0(W1))
hist(pscores, breaks = 200)
var(true$blip)
mean(true$blip>.4)
mean(true$df$Y)
m = min(true$blip)
M = max(true$blip)
m = -0.195
M = 0.319

for (n in c(1000, 2500, 5000, 10000, 25000, 50000)) {
  blips = seq(m, M, .01)
  bw = n^-.2
  step = round(bw/20, 3)
  bw_seq = seq(step, 20*step, step)
  big_ind = length(bw_seq)

  truths_h = lapply(bw_seq, FUN = function(h) truth.get(t = blips, h = h, k = list(deg = 4, range = 2),
                                                        d=1, g0, Q0)$truth_h)
  
  truthname = paste0("truths_h_", n, ".RData")
  save(truths_h, blips, g0, Q0, file = truthname)
  for (a in seq(6,48,6)) {
    blip = blips[a]
    truth = mean(true$blip> blip)
    B = 1000
    cl_size = ifelse(n > 10000, 12, 24)
    cl = makeCluster(cl_size, type = "SOCK")
    registerDoSNOW(cl)
    
    allresults=foreach(i=1:B,
                       .packages=c("cateSurvival","ggplot2","cowplot","mvtnorm","Iso")
                       ,.errorhandling='remove'
    )%dopar%
    {sim_bwselect(n=n, blip, truth = truth, truths_h = truths_h[[big_ind]][a], 
                  bw_seq = bw_seq, g0, Q0, kernel)}
    
    nname = (paste0("bwselect_", n, "_", a, ".RData"))
    save(allresults, bw_seq, blip, file = nname)
    stopCluster(cl)
  }
}