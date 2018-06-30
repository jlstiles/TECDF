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

# function to create plots for simultaneous pts

fnames = paste0("sim_unifCVhalglm_",1:3,"_1000_49.RData")
res_glm = lapply(seq_along(fnames), FUN = function(i) {
  load(fnames[[i]])
  assign(paste0("glm", i), glmstuff)
  })
res_hal = lapply(seq_along(fnames), FUN = function(i) {
  load(fnames[[i]])
  assign(paste0("hal", i), halstuff)
})


info = info_compile(res_hal, res_glm, n = 1000)
plots = info$plots
info$coverage_hal
info$coverage_glm
info$plots

save(info, g0, Q0, k, file = "sim_unifCVhalglm_1000_49.RData")

info_compile = function(res_hal, res_glm, n) {
  
  CATE = res_glm[[1]]$blip
  bw = res_glm[[1]]$h
  num = length(CATE) - 1
  
  res_halsep = lapply(0:num, FUN = function(D) {
  hals = lapply(res_hal, FUN = function(x) {
    B = x$B
    x$res[1:B+D*B,]}
    ) 
  res = do.call(rbind, hals)
  return(res)
  })
  
  res_glmsep = lapply(0:num, FUN = function(D) {
    glms = lapply(res_glm, FUN = function(x) {
      B = x$B
      x$res[1:B+D*B,]}
    ) 
    res = do.call(rbind, glms)
    return(res)
  })
  
  B = sum(unlist(lapply(res_hal, FUN = function(x) x$B)))
  cover_hal = sum(unlist(lapply(res_hal, FUN = function(x) x$coverage*x$B/B)))
  cover_glm = sum(unlist(lapply(res_glm, FUN = function(x) x$coverage*x$B/B)))
  
  plots = lapply(0:num, FUN = function(D) {
    res_temp = cbind(res_halsep[[D+1]], res_glmsep[[D+1]])
    S_t = res_temp[1,5]
    inds = c(1,4, 7, 10)
    ests = c(unlist(lapply(inds, FUN = function(x) res_temp[,x])))
    types = c("TMLE_hal", "Initial_hal", "TMLE_glm", "Initial_glm")
    type = c(unlist(lapply(types, FUN = function(x) rep(x,B))))
    
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
      geom_vline(xintercept=mean(res_temp[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(res_temp[,inds[3]]),color = colors[3])+
      geom_vline(xintercept=mean(res_temp[,inds[4]]),color = colors[4])
    
    if (length(CATE)==1) {
    caption = paste0("5000 bootstrap samples were drawn \n",
                     "coverage for tmle_glm = ", round(cover_glm, 3), "\n",
                     "coverage for tmle_hal = ", round(cover_hal, 3),"\n",
                     "truth smoothed parameter at the black line \n",
                     "kernel smoothed P(CATE > t)")
    } else {
      caption = paste0("5000 bootstrap samples were drawn \n",
                       "truth smoothed parameter at the black line \n",
                       "kernel smoothed P(CATE > t)")
    }
    ggover=ggdraw(add_sub(ggover,caption,x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                          vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                          colour = "black", size = 12, angle = 0, lineheight = 0.9))
    plot = ggover
    return(plot)
  })
  
  res_hal = do.call(rbind, lapply(res_hal, FUN = function(x) x$res))
  res_glm = do.call(rbind, lapply(res_glm, FUN = function(x) x$res))
  return(list(plots = plots, coverage_hal = cover_hal, coverage_glm = cover_glm,
              B = B, n = n, bw = bw, blip = CATE, res_hal = res_hal, res_glm = res_glm))
}  


plots = plot_create1(halstuff, glmstuff, n = 1000)
plots

plot1000unifsimul_halglm = arrangeGrob(info$plots[[1]], 
                                        info$plots[[2]], 
                                        info$plots[[3]], 
                                        info$plots[[4]], 
                                        ncol=2,nrow=2,widths = c(3.5,3.5),heights = c(3,3))
plot1000unifsimul_halglm
ggsave("plot1000unifsimul_halglm.png", plot = plot1000unifsimul_halglm, device = NULL,    path = NULL, 
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

