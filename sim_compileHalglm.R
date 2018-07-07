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
k = list()
k$degree = NULL
k$range = 1

# loading in RData files
allresults1 = allresults
allresults = append(allresults1, allresults)

# for hal vs glm
results = cateSurvival::get_results(allresults, n=25000, L = 8)

# for well-spec glm
results = cateSurvival::get_results_well(allresults, n=25000, L = 8)
results$plots


plot1000unifsimul_halglm = arrangeGrob(info$plots[[1]], 
                                        info$plots[[2]], 
                                        info$plots[[3]], 
                                        info$plots[[4]], 
                                        ncol=2,nrow=2,widths = c(3.5,3.5),heights = c(3,3))
plot1000unifsimul_halglm
ggsave("plot1000unifsimul_halglm.png", plot = plot1000unifsimul_halglm, device = NULL,    path = NULL, 
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

