library(ggplot2)
library(cowplot)
library(parallel)
library(doSNOW)
library(cateSurvival)
library(Iso)
library(foreach)

# j=3
# tnames = lapply(c(1000,2500,5000,10000,25000,50000), FUN = function(size) {
#   paste0("/Users/jlstiles/Dropbox/cateSurvival/results/kernel",j,"_bwselect_new/truths_h_",
#          size,"_kernel",j,".RData")
# })

# i=1
# load(tnames[[i]])

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
  kernel = kernel_list[[j]]
  for (n in c(1000,2500,5000,10000,25000,50000)) {
    bw = n^-.2
    step = round(bw/20, 3)
    bw_seq = seq(step, 20*step, step)
    r = length(bw_seq)
  
    for (a in seq(6,48,6)) {
      blip = blips[a]
      truth = mean(true$blip> blip)
      B = 1000
      cl_size = 12
      cl = makeCluster(cl_size, type = "SOCK")
      registerDoSNOW(cl)
      
      allresults=foreach(i=1:B,
                         .packages=c("cateSurvival","ggplot2","cowplot","mvtnorm","Iso")
                         ,.errorhandling='remove'
      )%dopar%
      {
        tmledata = gentmledata(n, d = 1, g0, Q0, V = 1, formu = NULL)
        tmle_info = gentmle_alt3(initdata=tmledata, estimate_fun = blipdist_estimate3,
                               update_fun = blipdist_update, max_iter = 1000, kernel = kernel,
                               simultaneous.inference = TRUE, blip = blip, h = bw_seq)
        ci = ci_gentmle(tmle_info)
        steps = tmle_info$steps
        return(list(ci = ci, steps = steps))
    }
    nname = (paste0("results_selector1/bwselect1_", n, "_", a,"_kernel", j, ".RData"))
    }
  }
}





tmle_info1 = gentmle_alt3(initdata=tmledata, estimate_fun = blipdist_estimate3,
                         update_fun = blipdist_update, max_iter = 1000, kernel = kernel,
                         simultaneous.inference = TRUE, blip = blip, h = bw_seq[1])

tmle_info1$tmleests

tmle_info$steps
tmle_info$tmleests
tmle_info$initests
abs(tmle_info$ED)-
sqrt(tmle_info$ED2)/n

truths = unlist(lapply(truths_h, FUN = function(x) x[seq(6,48,6)[a]]))
ci = ci_gentmle(tmle_info)
ci$lower <= truths & ci$upper >= truths 

