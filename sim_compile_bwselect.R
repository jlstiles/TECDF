library(ggplot2)
library(cowplot)
library(parallel)
library(doSNOW)
library(cateSurvival)
library(gridExtra)
# This was done to add the CIs conveniently to the results RData files
# for (j in 4) {
# 
#   fnames = lapply(c(1000,2500,5000,10000,25000,50000), FUN = function(size) {
#     paste0("/Users/jlstiles/Dropbox/cateSurvival/results/kernel",j,"_bwselect_new/bwselect_",
#            size,"_",seq(6,48,6),"_kernel",j,".RData")
#   })
# 
#   for (i in 1:6) {
#     for (a in 1:length(seq(6,48,6))){
#       load(fnames[[i]][[a]])
#       CIs = lapply(allresults, FUN = function(x) x$CIs_all[,1:3])
#       save(allresults, kernel, blip, bw_seq, CIs, file = fnames[[i]][[a]])
#     }
#   }
# }


info_allkerns = lapply(1:4, FUN = function(j) {
  # j=1
  fnames = lapply(c(1000,2500,5000,10000,25000,50000), FUN = function(size) {
    paste0("/Users/jlstiles/Dropbox/cateSurvival/results/kernel",j,"_bwselect_new/bwselect_",
           size,"_",seq(6,48,6),"_kernel",j,".RData")
  })
  
  tnames = lapply(c(1000,2500,5000,10000,25000,50000), FUN = function(size) {
    paste0("/Users/jlstiles/Dropbox/cateSurvival/results/kernel",j,"_bwselect_new/truths_h_",
           size,"_kernel",j,".RData")
  })
  
  info_allsizes = lapply(1:6, FUN = function(i) {
    # i=1
    load(tnames[[i]])
    true = gendata.blip(2000000,d = 1, g0, Q0)
    blip_select = blips[seq(6,48,6)]
    truth = unlist(lapply(1:8, FUN = function(a) mean(true$blip> blip_select[a])))
    info_allblips = lapply (1:8, FUN = function(a) {
      # a=7
      load(fnames[[i]][[a]])
      S0h = truths_h[[20]][seq(6,48,6)[a]]
      S0 = truth[a]
      allests = do.call(rbind, lapply(CIs, FUN = function(x) {
        x[,1]
      }))
      z_alpha = 1.96
      SE_true = apply(allests, 2, sd)
      
      # undebug(bwselect_m)
      # tester = bwselect_m(CIs[[222]], truth = list(S0 = truth[a], S0h = S0h),
      #            SE_true = SE_true, z_alpha = 1.96)
      # tester
      
      CI_info = lapply(CIs, FUN = function(CI) {
        bwselect_m(CI, truth = list(S0 = S0, S0h = S0h), 
                   SE_true = SE_true, z_alpha = z_alpha)
      }
      )
      
      CI_all = lapply(CI_info, FUN = function(x) x[2*1:7-1])
      ind_all = lapply(CI_info, FUN = function(x) unlist(x[2*1:6]))
      ind_avg = colMeans(do.call(rbind, ind_all))
      
      coverage_all = lapply(CI_all, FUN = function(x) {
        unlist(lapply(x, FUN = function(ci) {
          return(as.numeric(S0 >= ci[2] & S0 <= ci[3]))
        }))
      })
      coverage_all = do.call(rbind, coverage_all)
      coverage = colMeans(coverage_all)
      
      MSE = colMeans(do.call(rbind, lapply(CI_all, FUN = function(x) {
        unlist(lapply(x, FUN = function(ci) {
          return((ci[1] - S0)^2)
        }))
      })))
      
      SE = colMeans(do.call(rbind, lapply(CI_all, FUN = function(x) {
        unlist(lapply(x, FUN = function(ci) {
          return((ci[3] - ci[2])/(2*z_alpha))
        }))
      })))
      
      mono_avg = mean(unlist(lapply(CI_info, FUN = function(x) x$mono)))
      
      return(list(mono_avg = mono_avg, SE = SE, MSE = MSE, 
                  coverage = coverage, ind_avg = ind_avg))
      
    })
    return(info_allblips)
  })
  return(info_allsizes)
})

info_allkerns[[1]][[4]][[3]]
save(info_allkerns, file = "bwselect_allkerns.RData") 

cov_allkerns = lapply(info_allkerns, FUN = function(x) {
  lapply(x, FUN = function(y) {
    do.call(rbind, lapply(y, FUN = function(z) {
      z$coverage
    }))
  })
})
cov_allkerns

mono_allkerns = lapply(info_allkerns, FUN = function(x) {
  lapply(x, FUN = function(y) {
    do.call(rbind, lapply(y, FUN = function(z) {
      z$mono
    }))
  })
})
mono_allkerns

ind_allkerns = lapply(info_allkerns, FUN = function(x) {
  lapply(x, FUN = function(y) {
    do.call(rbind, lapply(y, FUN = function(z) {
      z$ind_avg
    }))
  })
})
ind_allkerns

SE_allkerns = lapply(info_allkerns, FUN = function(x) {
  lapply(x, FUN = function(y) {
    do.call(rbind, lapply(y, FUN = function(z) {
      z$SE
    }))
  })
})
SE_allkerns

