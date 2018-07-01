library(ggplot2)
library(cowplot)
library(parallel)
library(doSNOW)
library(cateSurvival)
library(gridExtra)
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

# set n, bw, blips
n=100
t = seq(m, M, .01)
bw = n^-.2
blip = seq(13,28,5)

# function to create plots for simultaneous pts
plots = plot_create(res10000unif_simul, 10000)
res10000unif_simul$plots = plots

plot_create = function(result, n) {
  B = result$B
  CATE = result$blip
  bw = result$h
  num = length(CATE) - 1
  res =result$res
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
      ggtitle("CATE 'survival' Sampling Dists", 
              subtitle = paste0("n = ", n, ", t = ",round(CATE[D+1],4),", bw = ", round(bw,4)))  
    ggover = ggover+geom_vline(xintercept = S_t,color="black")+
      geom_vline(xintercept=mean(res_temp[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(res_temp[,inds[2]]),color = colors[2])
    caption = paste0("5000 bootstrap samples were drawn \n",
                     "truth smoothed parameter at the black line \n",
                     "kernel smoothed P(CATE > t)")
    ggover=ggdraw(add_sub(ggover,caption,x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                          vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                          colour = "black", size = 12, angle = 0, lineheight = 0.9))
    plot = ggover
    return(plot)
  })
  return(plots)
}  


results = list(res10000unif13,
               res10000unif18,
               res10000unif23,
               res10000unif28)

plots = lapply(results, plot_create1, n = 10000)
plots
res10000unif13$plot = plots[[1]]
res10000unif18$plot = plots[[2]]
res10000unif23$plot = plots[[3]]
res10000unif28$plot = plots[[4]]


plot_create1 = function(result, n)  {
  bw = result$h
  CATE = result$blip
  res = result$res
  B = result$B
  coverage = result$coverage
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
    ggtitle("CATE 'survival' Sampling Dists", 
            subtitle = paste0("n = ", n, ", t = ",round(CATE,4),", bw = ", round(bw,4)))  
ggover = ggover+geom_vline(xintercept = S_t,color="black")+
    geom_vline(xintercept=mean(res[,inds[1]]),color = colors[1])+
    geom_vline(xintercept=mean(res[,inds[2]]),color = colors[2])
  caption = paste0("coverage of the true smoothed parameter is ", round(coverage,3),"\n",
                   "5000 bootstrap samples were drawn \n",
                   "truth smoothed parameter at the black line \n",
                   "kernel smoothed P(CATE > t)")
  ggover=ggdraw(add_sub(ggover,caption,x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                        vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                        colour = "black", size = 12, angle = 0, lineheight = 0.9))
  return(ggover)
}

save(res100unif13, res100unif18, res100unif23, res100unif28, 
     res250unif13, res250unif18, res250unif23, res250unif28,
     res500unif13, res500unif18, res500unif23, res500unif28,
     res1000unif13, res1000unif18, res1000unif23, res1000unif28,
     res2500unif13, res2500unif18, res2500unif23, res2500unif28,
     res5000unif13, res5000unif18, res5000unif23, res5000unif28,
     res10000unif23, res10000unif28,res25000unif13, res25000unif18, res25000unif23,
     res25000unif28, res50000unif13, res50000unif18, res50000unif23,
     res50000unif28,
     res100unif_simul, res250unif_simul, res500unif_simul, 
     res1000unif_simul, res2500unif_simul, res5000unif_simul, 
     res10000unif13,
     res10000unif18, res10000unif_simul, res25000unif_simul, res50000unif_simul,
     g0, Q0, file = "kernel_sim_unifallCV.RData")
 
plot25000unif_simul = arrangeGrob(res25000unif_simul$plots[[1]], res25000unif_simul$plots[[2]], 
                     res25000unif_simul$plots[[3]], res25000unif_simul$plots[[4]], 
                ncol=2,nrow=2,widths = c(3.5,3.5),heights = c(3,3))
plot25000unif_simul
ggsave("plot25000unif_simul.png", plot = plot25000unif_simul, device = NULL,    path = NULL, 
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)


plot10000unif = arrangeGrob(res10000unif13$plot, res10000unif18$plot, 
                                  res10000unif23$plot, res10000unif28$plot, 
                                  ncol=2,nrow=2,widths = c(3.5,3.5),heights = c(3,3))
plot10000unif
ggsave("plot10000unif.png", plot = plot10000unif, device = NULL,    path = NULL, 
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

