# library(reshape2)
# library(plyr)
# library(dplyr)
# library(SuperLearner)
library(ggplot2)
library(cowplot)
library(parallel)
library(doSNOW)

R = 5
degree = 4
kk = degree/2-2
area_row = vapply(0:(kk+2), FUN = function(i) 2*R^(2*i+1)/(2*i+1), FUN.VALUE = 1)
zero_row = vapply(0:(kk+2), FUN = function(i) R^(2*i), FUN.VALUE = 1)
deriv_row = c(0,vapply(0:(kk+1), FUN = function(i) 2*(i + 1)*R^(2*i+1), FUN.VALUE = 1))
if (kk>0) {
  orth_rows = lapply(seq(0,max((2*kk-2),0),2), FUN = function(r) {
    vapply(0:(kk+2), FUN = function(i) 2*R^(2*i+3+r)/(2*i+3+r), FUN.VALUE = 1)
  })
  orth_rows = do.call(rbind, orth_rows) 
  mm = rbind(area_row, zero_row, deriv_row, orth_rows)
} else mm = rbind(area_row, zero_row, deriv_row)

mm_inv = solve(mm)
veck = mm_inv %*% c(1, rep(0,kk+2))
veck

# The gaussian kernels are useless here, take a while to integrate and 
# then there's round off error so I choose explicit integrals, which is lightning fas
# and such kernels are formulaic to construct as per below--see tenth above

# unif_cdf = function(x) .5*as.numeric(x > -1)*(pmin(x ,1) + 1)
# Unif = function(x) .5*as.numeric(-1<=x&1>=x)
# 
# deg = 4
# R = 2
# mm = vapply(seq(1,deg,2), FUN = function(r) {
#   vapply(seq(r,(r+deg+1),2), FUN = function(x) 2*R^x, FUN.VALUE = 1)/seq(r,(r+deg+1),2)
# }, FUN.VALUE = rep(1,(deg+3)/2))
# 
# mm = cbind(mm, vapply(seq(0,deg+1,2), FUN = function(x) R^x, FUN.VALUE = 1))
# mm = t(mm)
# mm_inv = solve(mm)
# veck = mm_inv %*% c(1,rep(0, (deg + 1)/2))

kern = function(x, R, veck) {
  ll = lapply(1:length(veck), FUN = function(c) veck[c]*x^(2*c-2))
  w = Reduce("+", ll)*(x > -R & x < R)
  return(w)
}

kern_cdf = function(x, R, veck) {
  u = pmin(x, R)
  ll = lapply(1:length(veck), FUN = function(c) veck[c]*(u^(2*c-1) + R^(2*c-1))/(2*c-1))
  w = Reduce("+", ll)*as.numeric(x > -R)
  return(w)
}

# integrate to 1?
area = integrate(kern, lower = -R, upper = R, R = R, veck = veck, subdivisions = 10000)$value
area

# plot
s = seq(-R,R,.001)
y = kern(s, R=R, veck=veck)
plot = qplot(x=s,y=y)+geom_line()
plot

# check orthogonality (first val should be the area of 1)

test_fcn = as.data.frame(vapply(0:(degree-4), FUN = function(r) {
  test_fcn = function(x) (x^r)*kern(x, R=R, veck = veck)
  test_int = integrate(test_fcn, lower = -R, upper = R,subdivisions = 10000)
  return(c(test_int$abs.error, test_int$value))
}, FUN.VALUE = c(1,1)))
rownames(test_fcn) = c("abs_error", "integral")
colnames(test_fcn) = as.character(0:(degree - 4))
test_fcn

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

h_vec = seq(.01, 0.3, by = .01)
t = seq(m, M, .01)
n=1e4

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
gg_true

tmledata=gentmledata(n, d = 1, g0, Q0, formu = NULL)
blip = seq(10,25,5)
# undebug(CATEsurv_plot)
test = CATEsurv_plot(t = t[blip], h = h_vec[bw], k = k, 
                     truth = truth_h[blip], n = n, tmledata = tmledata)
test$info
test$plot
test$steps

# setting up a simulation
getres = function(n, t, h, k, truth, d = 1, g0, Q0, formu = NULL) {
 tmledata=gentmledata(n, d = 1, g0, Q0, formu = NULL)
 res = CATEsurv_plot(t = t, h = h, k = k, truth = truth, n = n, tmledata = tmledata)
 return(res$info)
 }

# undebug(getres)
# undebug(CATEcdf_plot)
# undebug(gentmle_alt1)
# undebug(blipdist_estimate2)
# debug(blipdist_update)
test = getres(n, t[blip], h_vec[bw], k = k, truth = truth_h[blip], d = 1, g0 = g0, Q0 = Q0)
test

B=1000

cl = makeCluster(detectCores(), type = "SOCK")
registerDoSNOW(cl)
blip = seq(13,28,5)
blip = c(5,30)
blip = 5
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

mean(cover)
head(res)


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
plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]
