library(ggplot2)
library(cowplot)
library(parallel)
library(doSNOW)
library(cateSurvival)

detectCores()

cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)

B=1000
time = proc.time()
tester_fcn = function(x) sd(rnorm(1e6))
allresults=foreach(i=1:B) %dopar% tester_fcn(1)
time = proc.time() - time
time

cl = makeCluster(1, type = "SOCK")
registerDoSNOW(cl)

B=1000
time = proc.time()
tester_fcn = function(x) sd(rnorm(1e6))
allresults=foreach(i=1:B) %dopar% tester_fcn(1)
time = proc.time() - time
time