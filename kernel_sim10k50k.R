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
for (n in c(10000, 50000)) {
  ii = ii + 1
  print(ii)
  h_vec = seq(.01, 0.4, by = .01)
  t = seq(m, M, .01)
  ideal = n^-.2
  bw = order(abs(h_vec - ideal))[1]
  
  truths = truth.get(t=t, h=h_vec[bw], k, d=1, g0=g0, Q0=Q0)
  truth_h = truths$truth_h
  truth = truths$truth
  
  true_df = data.frame(truth = c(truth_h, truth), 
                       CATE = rep(t,2),
                       type = c(rep("smoothed", length(truth)),
                                rep("true", length(truth))))
  gg_true = ggplot(true_df, aes(x=CATE, y=truth, color = type)) + geom_line()
  
  # setting up a simulation
  getres = function(n, t, h, k, truth, d = 1, g0, Q0, formu = NULL) {
    tmledata=gentmledata(n, d = 1, g0, Q0, formu = NULL)
    res = CATEsurv_plot(t = t, h = h, k = k, truth = truth, n = n, tmledata = tmledata)
    return(res$info)
  }
  
  cl = makeCluster(8, type = "SOCK")
  registerDoSNOW(cl)
  blip = seq(13,28,5)
  B=5000
  allresults=foreach(i=1:B,
                     .packages=c("cateSurvival","mvtnorm","ggplot2", "cowplot")
                     ,.errorhandling='remove'
  )%dopar%
  {getres(n, t[blip], h_vec[bw], k = k, 
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
  
  num = L-1
  plots = lapply(0:num, FUN = function(D) {
    res_temp = res[1:B+D*B,]
    S_t = res_temp$truth[1]
    inds = c(1,4)
    ests = c(unlist(lapply(inds, FUN = function(x) res_temp[,x])))
    types = c("TMLE", "Initial")
    type = c(unlist(lapply(types[1:2], FUN = function(x) rep(x,B))))
    
    inds = inds[order(types)]
    colors = c("red","blue", "green", "yellow", "orange") 
    
    plotdf = data.frame(ests = ests, type = type)
    ggover = ggplot(plotdf,aes(x=ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("CATE survival, t=", t[blip[D+1]],", bw=",h_vec[bw]))
    ggover = ggover+geom_vline(xintercept = S_t,color="black")+
      geom_vline(xintercept=mean(res_temp[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(res_temp[,inds[2]]),color = colors[2])
    plot = ggover
  })
 
  assign(paste0("res", n, "unif_simul"), list(plots = plots, 
                                              coverage = coverage, 
                                              B = B, 
                                              h = h_vec[bw], 
                                              res = res,
                                              blip = t[blip],
                                              true_df = true_df,
                                              plot_true = gg_true))

  cl = makeCluster(8, type = "SOCK")
  registerDoSNOW(cl)
  blip = seq(13,28,5)
  
  for (b in blip) {
    # b = blip[1]
    allresults=foreach(i=1:B,
                       .packages=c("cateSurvival","mvtnorm","ggplot2", "cowplot")
                       ,.errorhandling='remove'
    )%dopar%
    {getres(n, t[b], h_vec[bw], k = k, 
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
    
    
    S_t = res$truth[1]
    inds = c(1,4)
    ests = c(unlist(lapply(inds, FUN = function(x) res[,x])))
    types = c("TMLE", "Initial")
    type = c(unlist(lapply(types[1:2], FUN = function(x) rep(x,B))))
    
    inds = inds[order(types)]
    colors = c("red","blue", "green", "yellow", "orange") 
    
    plotdf = data.frame(ests = ests, type = type)
    ggover = ggplot(plotdf,aes(x=ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle(paste0("CATE survival, t=", t[b],", bw=",h_vec[bw]))
    ggover = ggover+geom_vline(xintercept = S_t,color="black")+
      geom_vline(xintercept=mean(res[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(res[,inds[2]]),color = colors[2])
    
    assign(paste0("res", n, "unif", b), list(plot = ggover, 
                                                   coverage = coverage, 
                                                   B = B, 
                                                   h = h_vec[bw], 
                                                   res = res,
                                                   blip = t[b]))

  }
}

save(res10000unif13, res10000unif18, res10000unif23, res10000unif28,
     res50000unif13, res50000unif18, res50000unif23, res50000unif28,
     res10000unif_simul, res50000unif_simul, g0, Q0, file = "kernel_sim_unif10k50k.RData")

