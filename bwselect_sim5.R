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
kernel6 = make_kernel(degree = 6, R = 4)
kernel8 = make_kernel(degree = 8, R = 5)
kernel10 = make_kernel(degree = 10, R = 5.6)
kernel12 = make_kernel(degree = 12, R = 5.6)

kernel_list = list(kernel4, kernel6, kernel8, kernel10, kernel12)
# getting two kinds of truths, plotting them
# make a grid of true parameters
true = gendata.blip(2000000,d = 1, g0, Q0)
m = min(true$blip)
M = max(true$blip)
m = -0.195
M = 0.319
blips = seq(m, M, .01)

# up to here we have kept the simulation the same as CVbwselect2
# but here we will get simultaneous inference and proceed as before

for (j in 5) {
  # j=1
  kernel = kernel_list[[j]]
  for (size in c(25000,50000)) {
    # n=10000
    degree = length(kernel_list$veck)/2+1
    bw = size^-(1/(2*degree+3))
    step = round(bw/20, 3)
    bw_seq = seq(step, 20*step, step)

    for (a in seq(6,48,6)) {
      # a = 48
      blip = blips[a]
      truth = mean(true$blip> blip)
      B = 1000
      cl = makeCluster(4, type = "SOCK")
      registerDoSNOW(cl)
      
      allresults=foreach(i=1:B,
                         .packages=c("cateSurvival","ggplot2","cowplot","mvtnorm","Iso")
                         ,.errorhandling='remove'
      )%dopar%
      {
        info = sim_bwselect(size, blip, bw_seq, g0, Q0, kernel, zscore = NULL)
        return(list(ests = info$ests, SE = info$SE))
      }
      nname = (paste0("results_selector/bwselect_", size, "_", a,"_kernel", j, ".RData"))
      save(allresults, file = nname)
    }
  }
}

