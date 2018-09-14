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

coverage_allkerns = lapply(1:4, FUN = function(j) {
  fnames = lapply(c(1000,2500,5000,10000,25000,50000), FUN = function(size) {
    paste0("/Users/jlstiles/Dropbox/cateSurvival/results/kernel",j,"_bwselect_new/bwselect_",
           size,"_",seq(6,48,6),"_kernel",j,".RData")
  })
  
  tnames = lapply(c(1000,2500,5000,10000,25000,50000), FUN = function(size) {
    paste0("/Users/jlstiles/Dropbox/cateSurvival/results/kernel",j,"_bwselect_new/truths_h_",
           size,"_kernel",j,".RData")
  })
  
  coverage_all = lapply(1:6, FUN = function(i) {
    load(tnames[[i]])
    true = gendata.blip(2000000,d = 1, g0, Q0)
    blip_select = blips[seq(6,48,6)]
    truth = unlist(lapply(1:8, FUN = function(a) mean(true$blip> blip_select[a])))
    lapply(1:8, FUN = function(a) {
      load(fnames[[i]][[a]])
      CIs_jl_info = lapply(CIs, FUN = function(x) bwselect_jl(x, 5))
      
      CIs_jl = as.data.frame(do.call(rbind, lapply(CIs_jl_info, FUN = function(ll) {
        c(ll$ind, ll$CI, ll$ind_plus, ll$CI_plus)
      })))
      colnames(CIs_jl) = c("ind_jl", "est_jl", "lower_jl", "upper_jl",
                           "ind_jlp", "est_jlp", "lower_jlp", "upper_jlp")
      
      
      CIs_m_info = lapply(CIs, FUN = function(x) bwselect_m(x, 5))
      
      CIs_m = as.data.frame(do.call(rbind, lapply(CIs_m_info, FUN = function(ll) {
        c(ll$ind, ll$CI, ll$ind_plus, ll$CI_plus)
      })))
      colnames(CIs_m) = c("ind_m", "est_m", "lower_m", "upper_m",
                          "ind_mp", "est_mp", "lower_mp", "upper_mp")
      
      CI_res = cbind(CIs_jl, CIs_m)[,c(3:4,7:8, 11:12, 15:16)]
      
      cov_jl = CI_res[,1] <= truth[a] & CI_res[,2] >= truth[a]
      cov_jlp = CI_res[,3] <= truth[a] & CI_res[,4] >= truth[a]
      cov_m = CI_res[,5] <= truth[a] & CI_res[,6] >= truth[a]
      cov_mp = CI_res[,7] <= truth[a] & CI_res[,8] >= truth[a]
      coverage = c(cov_jl = mean(cov_jl), 
                   cov_jlp = mean(cov_jlp), 
                   cov_m = mean(cov_m), 
                   cov_mp = mean(cov_mp))
      return(coverage)
    }
    )
  })
  
  coverage_all = lapply(coverage_all, FUN = function(ll) do.call(rbind, ll))
  return(coverage_all)
})

coverage_allkerns[[4]]
