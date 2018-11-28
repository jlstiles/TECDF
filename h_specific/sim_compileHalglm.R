library(ggplot2)
library(cowplot)
library(parallel)
library(doSNOW)
library(cateSurvival)
library(gridExtra)

# for well-spec
sim_abbr = "lin"
# Define th DGP functions for SCM
g0 = function(W1) plogis(.2+.2*W1)
Q0 = function(A,W1) plogis(A + 2.5*A*W1 + W1)

m = -0.195
M = 0.319
t = seq(m, M, .01)
blip = seq(6, 48, 6)

# for hal vs glm 
sim_abbr = ""
# Define th DGP functions for SCM
g0 = function(W1) plogis(-.1-.5*sin(W1) - .4*(abs(W1)>1)*W1^2)
Q0 = function(A,W1) plogis(.3*A + 5*A*sin(W1)^2 - A*cos(W1))

m = -0.168
M = .495
t = seq(m, M, .01)
blip = seq(8,64,8)

# define kernel
kernel = make_kernel(degree = 4,R = 2)
n=2500
truth_h = get.truth(t = t[blip], h = n^-.2, kernel = kernel, d = 1, g0 = g0, Q0 = Q0)
truth_h 

# FXING wrong kernel being used for truth CVhalglm0_2500 and CVhalglm_1000
# allresults = lapply(allresults, FUN = function(x) {
#   # x = allresults[[1]]
#   x$reshal_simul[,5] = truth_h
#   x$resglm_simul[,5] = truth_h
#   for (a in 1:length(blip)) {
#     x$reshal[[a]]$info[5] =  truth_h[a]
#     x$resglm[[a]]$info[5] =  truth_h[a]
# }
#     return(x)
# })

# allresults = lapply(allresults, FUN = function(x) {
#   x$reshal_simul[,6] = (x$reshal_simul[,2] <=  x$reshal_simul[,5]) & 
#     (x$reshal_simul[,3] >=  x$reshal_simul[,5])
#   x$resglm_simul[,6] = (x$resglm_simul[,2] <=  x$resglm_simul[,5]) & 
#     (x$resglm_simul[,3] >=  x$resglm_simul[,5])
#   for (a in  1:length(blip)) {
#     tt = (x$reshal[[a]]$info[2] <=  x$reshal[[a]]$info[5]) & 
#       (x$reshal[[a]]$info[3] >=  x$reshal[[a]]$info[5])
#     x$reshal[[a]]$info[6] = tt
#     tt = (x$resglm[[a]]$info[2] <=  x$resglm[[a]]$info[5]) & 
#       (x$resglm[[a]]$info[3] >=  x$resglm[[a]]$info[5])
#     x$resglm[[a]]$info[6] = tt
#   }
#   return(x)
# })



# loading in RData files
allresults1 = allresults
allresults = append(allresults1, allresults)

# save(allresults, file = "results/CVhalglm0_2500.RData")
# for hal vs glm
results = cateSurvival::get_results(allresults, n=2500, L = 8, blips = t[blip])
plots = results$plots
cov_simul = c(results$cover_hal_simul, results$cover_glm_simul)
cov_ind = cbind(results$cover_hal, results$cover_glm)
cov = rbind(cov_simul, cov_ind)
colnames(cov) = c("TMLE_hal", "TMLE_glm")
rownames(cov)[2:9] = paste0("blip = ", t[blip])
stargazer::stargazer(cov)
cov

plots = lapply(plots, FUN = function(p) {
  capt = paste0("Truth = Pr(blip > t) at black line")
  p = ggdraw(add_sub(p,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                     vpadding = grid::unit(1, "lines"), fontfamily = "", 
                     fontface = "plain",colour = "black", size = 14, angle = 0, 
                     lineheight = 0.9))
})
plots
MSE = do.call(rbind,lapply(results$mse, FUN = function(x) x[c(3,5),3]))
MSE


ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                      vpadding = grid::unit(1, "lines"), fontfamily = "", 
                      fontface = "plain",colour = "black", size = 14, angle = 0, 
                      lineheight = 0.9))
# results1$plots[[8]]
unlist(results$cover_hal)
unlist(results$cover_glm)
# allresults[[100]]$risk
# allresults[[10]]$supnorm

riskg = lapply(allresults, FUN = function(x) x$risk[2,])
riskg = do.call(rbind, riskg)
colMeans(riskg)

risk = lapply(allresults, FUN = function(x) x$risk[1,])
risk = do.call(rbind, risk)
colMeans(risk)

# analysis of IC approx of variance
sd_hal_IC = lapply(1:8, FUN = function(d){
  unlist(lapply(1:length(allresults), FUN = function(x) {
    (allresults[[x]]$reshal[[d]]$info[3] - allresults[[x]]$reshal[[d]]$info[2])/(2*1.96)
  }))
})

sd_hal_IC = unlist(lapply(sd_hal_IC, mean))
sd_hal_true = vapply(1:8, FUN = function(x) sd(results$est[[x]][,3]), 1)
sd_hal_IC
sd_hal_true

sd_glm_IC = lapply(1:8, FUN = function(d){
  unlist(lapply(1:length(allresults), FUN = function(x) {
    (allresults[[x]]$resglm[[d]]$info[3] - allresults[[x]]$resglm[[d]]$info[2])/(2*1.96)
  }))
})

sd_glm_IC = unlist(lapply(sd_glm_IC, mean))
sd_glm_true = vapply(1:8, FUN = function(x) sd(results$est[[x]][,5]), 1)
sd_glm_IC
sd_glm_true

# for well-spec glm
results = cateSurvival::get_results_well(allresults, n=1000, L = 8)
results$plots
results$cover_well_simul
results$cover_well
results$mse

plot1000unifsimul_halglm = arrangeGrob(info$plots[[1]], 
                                        info$plots[[2]], 
                                        info$plots[[3]], 
                                        info$plots[[4]], 
                                        ncol=2,nrow=2,widths = c(3.5,3.5),heights = c(3,3))
plot1000unifsimul_halglm
ggsave("plot1000unifsimul_halglm.png", plot = plot1000unifsimul_halglm, device = NULL,    path = NULL, 
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

