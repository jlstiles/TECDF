library(ggplot2)
library(cowplot)
library(parallel)
library(doSNOW)
library(cateSurvival)
library(Iso)
library(foreach)

getwd()

# Define th DGP functions for SCM
g0 = function(W1) plogis(.2 + .2*W1)
Q0 = function(A,W1) plogis(A + 2.5*A*W1 + W1)

# define kernel

kernel4 = make_kernel(degree = 4, R = 2)
kernel6 = make_kernel(degree = 6, R = 4)
kernel8 = make_kernel(degree = 8, R = 5)
kernel10 = make_kernel(degree = 10, R = 5.6)
kernel12 = make_kernel(degree = 12, R = 5.6)

kernel_list = list(kernel4, kernel6, kernel8, kernel10, kernel12)
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


kernels = list(c(4,2), c(6, 4), c(8, 5), c(10, 5.6), c(12, 5.6))
for (i in 1:5){
  kernel = kernel_list[[i]]
  for (n in c(1000, 2500, 5000, 10000, 25000, 50000)) {
    blips = seq(m, M, .01)
    # power = -1/(2*(kernels[[i]][1]) - 3)
    # bw = n^power
    bw = n^-.2
    step = round(bw/20, 3)
    bw_seq = seq(step, 20*step, step)
    big_ind = length(bw_seq)
    
    truths_h = lapply(bw_seq, FUN = function(h) truth.get(t = blips, h = h, 
                                                          k = list(deg = kernels[[i]][1], 
                                                                   range = kernels[[i]][2]),
                                                          d=1, g0, Q0)$truth_h)
    
    truthname = paste0("results_selector/truths_h_", n,"_kernel",i, ".RData")
    save(truths_h, blips, g0, Q0, kernel, file = truthname)
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
      {sim_bwselect1(n=n, blip, truth = truth, truths_h = truths_h[[big_ind]][a], 
                    bw_seq = bw_seq, g0, Q0, kernel)}
      
      nname = (paste0("results_selector/bwselect_", n, "_", a,"_kernel", i, ".RData"))
      save(allresults, bw_seq, blip, kernel, file = nname)
      stopCluster(cl)
    }
  }
}