library(plyr)
library(dplyr)
library(mvtnorm)
library(ggplot2)
library(foreach)
library(doSNOW)
library(SuperLearner)
library(parallel)
library(reshape)
library(reshape2)
library(cowplot)
# library(origami)
library(geepack)
library(mvtnorm)
library(gentmle)
library(orthopolynom)

# playing with integrating functions
f2 = function(x) (1/sqrt(2*pi))*exp(-x^2/2)
f4 = function(x) (1/sqrt(2*pi))*exp(-x^2/2)*(1.5-x^2/2)
f6 = function(x) (1/8)*(1/sqrt(2*pi))*exp(-x^2/2)*(15-10*x^2+x^4)
f8 = function(x) (1/48)*(1/sqrt(2*pi))*exp(-x^2/2)*(105-105*x^2+21*x^4-x^6)

kernels = list(f2,f4,f6,f8)

# verify integral 1
lapply(kernels, FUN =function(x) {
  integrate(x, lower = -100, upper = 100, subdivisions = 10000)$value
} )

# verify orthogonality
ff2 = function(x) x^3*f8(x)
int2 = integrate(ff2, lower = -10, upper = 10, subdivisions = 10000)$value
int2

ff6 = function(x) f8(x)*x^6
int6 = integrate(ff6, lower = -10, upper = 10, subdivisions = 10000)$value
int6

# view of kernel orthogonal to degree 6 polyns
x = seq(-5,5,.001)
y = f8(x)
plot(x,y)

chebs = chebyshev.u.polynomials(10, normalized=TRUE)
chebs
polyL =  legendre.polynomials(10,normalized=TRUE)
polyL


g6 = function(x) (x^6-(126/99)*x^4+(35/99)*x^2)*10395/128
g4 =function(x) (105/8)*((5/7)*x^2-x^4)
g2 = function(x) 1.5*x^2

integrate(g4, lower =-1, upper = 1, subdivisions = 1000)$value

x=seq(-1,1,.01)
plot(x,g4(x))

poly = chebyshev.c.polynomials(3, normalized=FALSE)
poly
cf2 = function(x) -1 + 0.5*x^2 
cf3 = function(x) -1.5*x + 0.5*x^3
ff = function(x)  cf2(x)*cf3(x)
chebyshev.c.inner.products(2)

U = function(x) .5*as.numeric(-1<=x&1>=x)
UU = function(x) (1-x^2)*as.numeric(-1<x&x<1)*.75

g0=function(W1,W2,W3,W4) {plogis(-.8*W1+.39*W2+.08*W3-.12*W4-.15)}
# g0=function(W1,W2,W3,W4) {plogis(-.28*W1+1*W2+.08*W3-.12*W4-1)}
# Q0=function(A,W1,W2,W3,W4) {plogis(3*A-1*W1*A+1.2*W2-1.9*W3+.4*W4*A*W2-cos(W4))}
# Q0=function(A,W1,W2,W3,W4) {plogis(.5*A-1*W1+1.2*W2-1.9*W3+.4*W4)}
Q0=function(A,W1,W2,W3,W4) {plogis(A-1*W1+1.2*W2-1.9*W3+.4*W4+.8*A*W2+0.7*A*W1)}
# Q0=function(A,W1,W2,W3,W4) {plogis(A-1*W1+1.2*W2-1.9*W3+.4*W4+.3*W2+.5*A*W1)}
gendata.blip=function(n){
  U1 = runif(n,0,1)
  W1= -(U1<.5)+(U1>=.5)
  # W1=rnorm(n)
  W2=rnorm(n)
  W3=rbinom(n,1,.3)
  W4=rnorm(n)
  A=rbinom(n,1,g0(W1,W2,W3,W4))
  Y=rbinom(n,1,Q0(A,W1,W2,W3,W4))
  blip = Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4)
  return(list(df=data.frame(A,W1,W2,W3,W4,Y),blip=blip))
}

gentmledata = function(n) {
  data = gendata.blip(n)$df
  data0 = data
  data0$A = 0
  data1 = data
  data1$A = 1
  newX = rbind(data,data1,data0)

  fitQ = glm(Y~A*W1+A*W2+W3+W4,data=data, family = "binomial")
  Q0W = predict(fitQ,newdata = data0, type = 'response')
  Q1W = predict(fitQ, newdata = data1, type = 'response')
  QAW = predict(fitQ, newdata = data, type= 'response')

  Q = cbind(QAW, Q0W, Q1W)
  datag = data
  datag$Y = NULL
  fitg = glm(A~., data=datag, family='binomial')
  g1W = predict(fitg, type = 'response')

  tmledata = list(Q=Q,Y=data$Y,A=data$A,g1W=g1W)
  return(tmledata)}

get.zscore = function(Dstar, alpha) {

  if (ncol(Dstar)==1) return(1.96544) else{
  sigma = cor(Dstar)
  means = rep(0,ncol(sigma))
  zs = rmvnorm(n=2000000,mean= means, sigma=sigma, method= "chol")
  zabs = apply(zs,1,FUN = function(x) max(abs(x)))
  zscore = quantile(zabs, probs = 1-alpha)
  return(zscore)}
}

# do an estimate
#
#
#
#
n=10000
tmledata=gentmledata(n)
kernel = U

# t=c(-.3,0,.3,.6)
h=.1
t=c(-.3,-.2,-.1,0,.1,.2,.3,.4,.5,.6,.7)
# t=.4
ff= gentmle::gentmle_alt1(initdata=tmledata, estimate_fun = blipdist_estimate,
             update_fun = blipdist_update,max_iter = 1000, t=t, h=h,
             kernel = U)

zscore = get.zscore(Dstar = ff$Dstar, .05)
zscore

pm = zscore*apply(ff$Dstar, 2, sd)*sqrt(n-1)/n
pm

whole.curve = ff$tmleests
whole.curve

initests = ff$initests

true = gendata.blip(1000000)
pscores = with(true$df, g0(W1,W2,W3,W4))
blips = true$blip

hist(pscores, breaks = 200)
hist(blips, breaks=200)
var(blips)

# get a sample S(t)
#
#
#
x0=.25
S_t = round(mean(blips>x0), 2)
S_t

with(dens, polygon(x=x, y= y, col="gray"))

dd <- with(dens,data.frame(x,y))
dd
plot2 = ggplot(data = dd, aes(x,y)) + geom_line()+
  geom_ribbon(data=subset(dd,x>x0),aes(ymax=y),ymin=0,
              fill="gray",colour=NA,alpha=0.5)+
  scale_x_continuous(breaks= c(-1,0,0.25,1), labels = 
                       c(-1,0,"0.25", 1))+
  ggtitle("CATE Distribution",subtitle = "S(0.25) = 0.29 = prob of treatment effect beyond 0.25") +
  labs(x="CATE", y = NULL)
plot2=ggdraw(add_sub(plot2,paste0("S(.25)=",S_t), x= .65, y = 4.5, hjust = 0, vjust = 0.5,
                     vpadding = grid::unit(1, "lines"), fontfamily = "", 
                     fontface = "plain",colour = "black", size = 12, angle = 0, 
                     lineheight = 0.9))
capt = "E[Y | (A, W1, W2, W3, W4)] = A-1*W1+1.2*W2-1.9*W3+.4*W4+.8*A*W2+.7*A*W1"
plot2=ggdraw(add_sub(plot2,capt, x= .05, y = 1, hjust = 0, vjust = 0.5,
                     vpadding = grid::unit(1, "lines"), fontfamily = "", 
                     fontface = "plain",colour = "black", size = 12, angle = 0, 
                     lineheight = 0.9))
plot2
ggsave("~/Dropbox/Quals/v0321/blipSurv.jpg", plot = plot2, 
       device = NULL, path = NULL, scale = 1, width = NA, height = NA, 
       units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

## Plot different bandwidth truths and the true param curve
#
#
#

truth_h = vapply(t, FUN = function(x0){
  truth = vapply(blips, FUN = function(x) {
    i = .5*(x>(x0-h))*(min(x,x0+h)-x0+h)/h
    # upper = (min(x, (x0+h))-x0)/h
    # lower = -1
    # i=(pnorm(upper)-pnorm(lower))*(x>x0-h)/(pnorm(1)-pnorm(-1))
    return(i)}, FUN.VALUE = 1)
  pt = mean(truth)
  return(pt)
}, FUN.VALUE=1)

truth = vapply(t, FUN = function(x) {
  mean(blips>x)
}, FUN.VALUE = 1)

bw = c(.1,.2,.3)
truths = vapply(bw, FUN = function(h) {
  vapply(t, FUN = function(x0){
  truth = vapply(blips, FUN = function(x) {
    i = .5*(x>(x0-h))*(min(x,x0+h)-x0+h)/h
    # upper = (min(x, (x0+h))-x0)/h
    # lower = -1
    # i=(pnorm(upper)-pnorm(lower))*(x>x0-h)/(pnorm(1)-pnorm(-1))
    return(i)}, FUN.VALUE = 1)
  pt = mean(truth)
  return(pt)
}, FUN.VALUE=1)
}, FUN.VALUE = rep(1, length(t)))


df = data.frame(t = rep(t,4),
                truth = c(truths[,1],truths[,2],truths[,3],truth),
                bandwidth = c(rep("h = 0.1",length(t)),
                         rep("h = 0.2",length(t)),
                         rep("h = 0.3",length(t)), 
                         rep("true", length(t))))

df
surv_params = ggplot(data = df, aes(x=t,y=truth, color = bandwidth)) + geom_line()+
  ggtitle(paste0("Different Parameters for the true ", "'survival'", " curve"),
          subtitle = "bandwidth = 0.1,0.2,0.3 Uniform[-1,1]")+labs(x="CATE", y = "S(CATE)")
surv_params

ggsave("~/Dropbox/Quals/v0321/surv_params.jpg", plot = surv_params,
       device = NULL, path = NULL, scale = 1, width = NA, height = NA,
       units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

# Plot an estimate
#
#
#
#
ff$initests
ff$steps
colMeans(ff$Dstar)
mean(blips)
var(blips)

df = data.frame(t = rep(t,4),
                left = rep(whole.curve-pm,4),
                right = rep(whole.curve + pm, 4),
                truth = c(whole.curve,initests, truths[,1],truth),
                type = c(rep("est", length(t)),rep("init est", length(t)),
                         rep("smoothed",length(t)), rep("true", length(t))))

df
surv_est = ggplot(data = df, aes(x=t,y=truth, color = type)) + geom_line()+
  geom_ribbon(data=df,aes(ymax=right,ymin=left),
              fill="gray",colour=NA,alpha=0.5)+labs(y = "S(CATE)", x = "CATE")+
  ggtitle(paste0("estimating CATE ", "'survival'", " curve"),
          subtitle = "bandwidth = 0.1, Uniform[-1,1] kernel, n = 10000")

capt = paste0("S(t) = prob CATE exceeds t.\n",
              "'smoothed' means the true parameter for bandwidth 0.1\n",
              "'true' is the true S(t) and 'est' is the estimate of smoothed parameter\n",
              "Shaded region represents simultaneous 95% CI to cover the whole curve\n",
              "both pscore and outcome model were well-specified")
surv_est=ggdraw(add_sub(surv_est,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                        vpadding = grid::unit(1, "lines"), fontfamily = "",
                        fontface = "plain",colour = "black", size = 12, angle = 0,
                        lineheight = 0.9))
surv_est

ggsave("~/Dropbox/Quals/v0321/surv_est.jpg", plot = surv_est,
       device = NULL, path = NULL, scale = 1, width = NA, height = NA,
       units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

# plot examples of blip dists with blip vars
#
#
#

library(Simulations)


Q0_sm = function (A, W1, W2, W3, W4) 
{
  plogis(0.14 * (2 * A + .5 * A * W1 + 2 * A * W3 * W4 + W2 * 
                   W1 + W3 * W4 + 5 * A * cos(W4)))
}

Qs = list(Q0_sm, Q0_2, Q0_1, Q0_trig)
plots = lapply(Qs, FUN = function(x) {
  data = gendata(1e6,g0_linear, x)
  blips = with(data, x(1,W1,W2,W3,W4) - x(0,W1,W2,W3,W4))
  v = var(blips)
  a = mean(blips)
  stdev = sqrt(v)
  plot = ggplot(data.frame(CATE = blips), aes(x=CATE))+
    geom_density(color = "#CCCCCC",fill = "#CCCCCC")+
    ggtitle("A CATE Distribution")
  capt = paste0("CATE variance = ",round(v,3),
                "\nCATE st dev = ",round(stdev,3),
                "\nCATE mean = ",round(a,3))
  plot=ggdraw(add_sub(plot,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                          vpadding = grid::unit(1, "lines"), fontfamily = "",
                          fontface = "plain",colour = "black", size = 12, angle = 0,
                          lineheight = 0.9))
  
  return(plot)
})


var_examples=marrangeGrob(list(plots[[1]],plots[[3]],
                     plots[[2]],plots[[4]]),
                ncol=2,nrow=2,widths = c(3.5,3.5),heights = c(1,1))
ggsave("~/Dropbox/Quals/v0321/var_examples.jpg", plot = var_examples, device = NULL,    path = NULL, 
       scale = 1, width = 8, height = 6, units = c("in"), dpi = 300, limitsize = TRUE)


# Plot the cv_advert
#
#
#
#

load("~/Dropbox/Quals/v0321/case2bSL2.RData")
load("~/Dropbox/Quals/v0321/case2bSL1.RData")
load("~/Dropbox/Quals/v0321/case2bCVSL2.RData")
colnames(res_list[[1]])

g0 = g0_linear
Q0 = Q0_trig
testdata=gendata(1000000, g0=g0, Q0 = Q0)
blip_true = with(testdata,Q0(1,W1,W2,W3,W4)-Q0(0,W1,W2,W3,W4))
propensity = with(testdata, g0(W1,W2,W3,W4))
ATE0 = mean(blip_true)
var0 = var(blip_true)

res_list = list(results_case2bCVSL2, results_case2bSL2, results_case2bSL1)
B = c(unlist(lapply(res_list, FUN = function(x) nrow(x))), nrow(results_case2bCVSL2))
B
type= c(rep("CV-TMLE SL2", B[1]), rep("TMLE SL2", B[2]), rep("TMLE SL1", B[3]), rep("LR", B[1]))
types = c("CV-TMLE SL2", "TMLE SL2", "TMLE SL1", "LR")


ests = c(res_list[[1]][,25],res_list[[2]][,25],res_list[[3]][,25],res_list[[1]][,40])
means = c(mean(res_list[[1]][,25]),mean(res_list[[2]][,25]),
          mean(res_list[[3]][,25]),mean(res_list[[1]][,40]))
means
means = means[order(types)]
colors = c("red","blue", "green","orange")

plotdf = data.frame(ests=ests, type=type)
# ticks = c(round(c(mean(res_LR[,1]), mean(res_LR[,7])), 3), seq(.05, .20, .05))
ggover = ggplot(plotdf,aes(x=ests, color = type, fill=type)) + 
  geom_density(alpha=.5)+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  labs(title = "ATE sampling distributions",
       subtitle = "TMLE is very robust for ATE")+
  theme(axis.title.x = element_blank(), plot.title = element_text(size = rel(1.5)))
ggover = ggover+geom_vline(xintercept = ATE0,color="black")+
  geom_vline(xintercept=means[1],color = colors[1])+
  geom_vline(xintercept=means[2],color = colors[2])+
  geom_vline(xintercept=means[3],color = colors[3])+
  geom_vline(xintercept=means[4],color = colors[4])

capt = paste0("Truth is at black vline.\n",
              "Logistic Regression suffers a little bit")

ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                      vpadding = grid::unit(1, "lines"), fontfamily = "", 
                      fontface = "plain",colour = "black", size = 14, angle = 0, 
                      lineheight = 0.9))

ggover

ggsave("~/Dropbox/Quals/v0321/ATEtrig.jpg", plot = ggover, device = NULL,    path = NULL, 
       scale = 1, width = 8, height = 6, units = c("in"), dpi = 300, limitsize = TRUE)


# plot the mixed results for mispec g and barQ
#
#
#

g0 = function (W1, W2) {
  plogis(.4*(-0.4 * W1*W2 + 0.63 * W2^2 -.66*cos(W1) - 0.25))
}

Q0 = function (A, W1, W2) {
  plogis(0.2 * W1*W2 + 0.1 * W2^2 - .8*A*(cos(W1) + .5*A*W1*W2^2) - 0.35)
}

Q0 = function (A, W1, W2) {
  plogis(0.1 * W1*W2 + 1.5*A*cos(W1) + 0.15*W1 - .4*W2*(abs(W2) > 1) -1*W2*(abs(W2 <=1)))
}

gendata.fcn = function (n, g0, Q0) 
{
  W1 = runif(n, -3, 3)
  W2 = rnorm(n)
  A = rbinom(n, 1, g0(W1, W2))
  Y = rbinom(n, 1, Q0(A, W1, W2))
  data.frame(A, W1, W2, Y)
}

pop = gendata.fcn(1e6, g0, Q0)
pscores = with(pop, g0(W1, W2))
pscore_plot = ggplot(data.frame(pscores=pscores), aes(x=pscores)) + 
  geom_density(fill="#FF6600")+
  ggtitle("Propensity Score Distribution")
capt = "g0(A|W) = expit(.4*(-0.4 * W1*W2 + 0.63 * W2^2 -.66*cos(W1) - 0.25))"
pscore_plot

pscore_plot=ggdraw(add_sub(pscore_plot,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                      vpadding = grid::unit(1, "lines"), fontfamily = "", 
                      fontface = "plain",colour = "black", size = 14, angle = 0, 
                      lineheight = 0.9))
pscore_plot

ggsave("~/Dropbox/Quals/v0321/pscore.jpg", plot = pscore_plot, device = NULL,    path = NULL, 
       scale = 1, width = 8, height = 4, units = c("in"), dpi = 300, limitsize = TRUE)


blips = with(pop, Q0(1, W1, W2) - Q0(0, W1, W2))
blips_plotGood = ggplot(data.frame(CATE=blips), aes(x=CATE)) + geom_density(fill="blue")+
  ggtitle("CATE Distribution, decent case for CATE variance")
capt = paste0("Q0(A,W) = expit(0.1 * W1*W2 + 1.5*A*cos(W1) + 0.15*W1 -\n",
              "0.4*W2*(abs(W2) > 1) -1*W2*(abs(W2) <=1)))")

blips_plotGood=ggdraw(add_sub(blips_plotGood,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                           vpadding = grid::unit(1, "lines"), fontfamily = "", 
                           fontface = "plain",colour = "black", size = 14, angle = 0, 
                           lineheight = 0.9))
blips_plotGood
ggsave("~/Dropbox/Quals/v0321/blips_caseGood.jpg", plot = blips_plotGood, device = NULL,    path = NULL, 
       scale = 1, width = 8, height = 4, units = c("in"), dpi = 300, limitsize = TRUE)


Q0 = function (A, W1, W2) {
  plogis(0.2 * W1*W2 + 0.1 * W2^2 - .8*A*(cos(W1) + .5*A*W1*W2^2) - 0.35)
}
blips = with(pop, Q0(1, W1, W2) - Q0(0, W1, W2))
blips_plotBad = ggplot(data.frame(CATE=blips), aes(x=CATE)) + geom_density(fill="blue")+
  ggtitle("CATE Distribution, bad case for CATE variance")
capt = "Q0(A,W) = expit(0.2 * W1*W2 + 0.1 * W2^2 - .8*A*(cos(W1) + .5*A*W1*W2^2) - 0.35)"

blips_plotBad=ggdraw(add_sub(blips_plotBad,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                              vpadding = grid::unit(1, "lines"), fontfamily = "", 
                              fontface = "plain",colour = "black", size = 14, angle = 0, 
                              lineheight = 0.9))
blips_plotBad
ggsave("~/Dropbox/Quals/v0321/blips_caseBad.jpg", plot = blips_plotBad, device = NULL,    path = NULL, 
       scale = 1, width = 8, height = 4, units = c("in"), dpi = 300, limitsize = TRUE)


var0 = var(blips)
ATE0 = mean(blips)
setwd("~/Dropbox/Quals/v0321")


load("case_d24shit.RData")
results = data.matrix(data.frame(do.call(rbind, ALL)))

load("case_d23shit.RData")
results1 = data.matrix(data.frame(do.call(rbind, ALL)))

load("case_d22shit.RData")
results2 = data.matrix(data.frame(do.call(rbind, ALL)))

load("case_d21shit.RData")
results3 = data.matrix(data.frame(do.call(rbind, ALL)))

load("case_d2shit.RData")
results4 = data.matrix(data.frame(do.call(rbind, ALL)))


results = rbind(results, results1, results2, results3, results4)
nrow(results)

colnames(results)
varinds = c(1, 45, 54)
ateinds = varinds+3
varAll = c(varinds, 7, 57)
ateAll = c(ateinds, 8, 58)

cov.check(results, ATE0, ateinds)
cov.check(results, var0, varinds)

mean(results[,6] - results[,5])
mean(results[,50] - results[,49])

B=nrow(results)
param = var0
inds = c(varinds, 7)
# inds = c(ateinds, 8)
type= c(rep("CV-TMLE SL 1step",B), rep("CV-TMLE LR 1step",B), rep("Delta LR 1step",B), 
        rep("SL initial",B))
types = c("CV-TMLE SL 1step", "CV-TMLE LR 1step", "Delta LR 1step", "SL initial")
ests = unlist(lapply(inds, FUN = function(x) results[,x]))

inds = inds[order(types)]
colors = c("red","blue", "yellow", "green")

plotdf = data.frame(ests=ests, type=type)

ggover = ggplot(plotdf,aes(x=ests, color = type, fill=type)) + 
  geom_density(alpha=.5)+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  theme(axis.title.x = element_blank())+
  ylim(0,100)+
  ggtitle("CATE variance sampling distributions")
ggover = ggover+geom_vline(xintercept = param,color="black")+
  geom_vline(xintercept=mean(results[,inds[1]]),color = colors[1])+
  geom_vline(xintercept=mean(results[,inds[2]]),color = colors[2])+
  geom_vline(xintercept=mean(results[,inds[3]]),color = colors[3])+
  geom_vline(xintercept=mean(results[,inds[4]]),color = colors[4])

capt = paste0("Truth is at black vline.\n",
              "CV-TMLE SL covers at 30.6%, others castastrophic")
ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                      vpadding = grid::unit(1, "lines"), fontfamily = "", 
                      fontface = "plain",colour = "black", size = 10, angle = 0, 
                      lineheight = 0.9))
ggover

ggsave("~/Dropbox/Quals/v0321/misspecBad.jpg", plot = ggover, device = NULL,    path = NULL, 
       scale = 1, width = 8, height = 6, units = c("in"), dpi = 300, limitsize = TRUE)


SLcoef = colMeans(results[,grep("coef", colnames(results))])
SLrisk = colMeans(results[,grep("risk", colnames(results))])[1:(length(grep("risk", colnames(results)))-2)]
# rownames(SLres)=NULL
SLres = cbind(SLcoef, SLrisk)
stargazer(SLres, summary = FALSE)
SLres

results[1:20,grep("risk", colnames(results))]
results[1:20,grep("coef", colnames(results))]
nrow(results)


# Do a depiction of the Pn^* gotten from Pn
#
#
#
#

x = seq(-.2,.05,by=.001)
fx = function(x) {
  z = rep(-999, length(x))
  z[x>= 0] = (x[x>=0])^2
  z[x< 0] = .4*(x[x< 0])^2
  return(z)
}
y = fx(x)

plotdf = data.frame(x=x, y=y)

ggover = ggplot(plotdf,aes(x=x, y=y)) + geom_line() + geom_point(x=-.1, y = .004)+
  theme_bw(base_size = 18)+
  geom_point(x=0,y=0, color = "red") + labs(y="Loss at P_{n,e}", x = "e")+
  scale_x_continuous(breaks= c(-.1,0), labels = c(expression(epsilon=0),"e*"))+
  scale_y_continuous(breaks = NULL)+
  ggtitle("One-step TMLE Mapping Visual")

capt = paste0("The empirical loss is minimized at P_{n,e*} = P_n^*\n",
              "The desired TMLE mapping from initial estimate P_n^0 (P_{n,0})\n", 
              "to P_n^* such that P_n D*(P_n^*) = 0.  For the one-step TMLE\n",
              "the derivative of the empirical loss wrt e is the same as the\n",
              "empirical mean of the efficient influence curve approximation.")
ggover=ggdraw(add_sub(ggover,capt, x= 0, y = 0.3, hjust = 0, vjust = 0.5,
                      vpadding = grid::unit(1, "lines"), fontfamily = "", 
                      fontface = "plain",colour = "black", size = 18, angle = 0, 
                      lineheight = 0.9))
ggover

ggsave("~/Dropbox/Quals/v0321/onestepVisual.jpg", plot = ggover, device = NULL,    path = NULL, 
       scale = 1, width = 8, height = 4, units = c("in"), dpi = 300, limitsize = TRUE)


# visual of noise
#
#
#
#
library(Simulations)
W = seq(-5,5,.01)
Q1W = plogis(.2*1 + W)
n = 100
Z = rnorm(n)
mm = n^(-.3)*Z
ss = abs(n^(-.3)*Z)
bias = rnorm(length(mm),mm,ss)
Qnoise100 = plogis(.2*1 + Z + 2*bias)
d100 = data.frame(Z=Z[1:100],Qnoise=Qnoise100[1:100])

n = 5000
Z = rnorm(n)
mm = n^(-.3)*Z
ss = abs(n^(-.3)*Z)
bias = rnorm(length(mm),mm,ss)
Qnoise1000 = plogis(.2*1 + Z + 2*bias)
d1000 = data.frame(Z=Z[1:100],Qnoise=Qnoise1000[1:100])

colors = c("orange","blue") 
df = rbind(d1000,d100)
df$size = c(rep("n=5000",100), rep("n=100",100))
head(df)
plotdf = data.frame(W=W, Q1W = Q1W)
noiseplot = ggplot(plotdf, aes(x=W, y= Q1W)) + geom_line() + 
  geom_line(data = df, mapping = aes(x=Z, y=Qnoise, color = size))+
  scale_color_manual(values=colors)
noiseplot

ggsave("~/Dropbox/Quals/v0321/noiseplot.jpg", plot = noiseplot, device = NULL,    path = NULL, 
       scale = 1, width = 5, height = 3, units = c("in"), dpi = 300, limitsize = TRUE)








getres = function(t=.5,h=.1,kernel=U,true=.2313561, TG=.240668,n){
  tmledata=gentmledata(n)
  n = nrow(tmledata$Q)
  res = gentmle_alt1(initdata=tmledata, estimate_fun = blipdist_estimate,
               update_fun = blipdist_update,max_iter = 1000, t=t, h=h,
               kernel = kernel)
  est = res$tmleests
  inits = res$initests
  steps = res$steps
  sigma = sqrt((n-1)/n)*sd(res$Dstar)
  left = est-1.96*sigma/sqrt(n)
  right = est+ 1.96*sigma/sqrt(n)
  cover = left<=true&right>=true
  result= c(est=est,left=left,right=right,init=inits,steps=steps,cover=cover)
  return(result)
}
B=100

cl = makeCluster(4, type = "SOCK")
registerDoSNOW(cl)

allresults=foreach(i=1:B,
                   .packages=c("plyr","gentmle"))%dopar%{getres(t=.4,h=.15,kernel=U,true=truth,
                                                      n=4000)}

data.matrix(do.call(rbind, allresults))
cover=unlist(lapply(allresults,FUN=function(x) x[[1]]$cover))
mean(cover)


ests = unlist(lapply(allresults,FUN=function(x) x[[1]]$est))
mean(ests)
cover=unlist(lapply(allresults,FUN=function(x) x[[1]]$cover))
mean(cover)
right=unlist(lapply(allresults105,FUN=function(x) x[[1]]$right))
mean(right)
left=unlist(lapply(allresults105,FUN=function(x) x[[1]]$left))
mean(left)
inits=unlist(lapply(allresults105,FUN=function(x) x[[1]]$init))
mean(inits)
truth.05[6]
hist(ests,breaks=100)
steps = unlist(lapply(allresults,FUN=function(x) x[[1]]$steps))
steps

hist(estsUU,breaks=100)

allresults[[2]][[1]]$cover

ff$Q[1:10,]

t=c(-.3,-.2,-.1,0,.1,.2,.3,.4,.5,.6,.7)
h=.05
t=.4
true = gendata.blip(1000000)
hist(true$blip, breaks=200)
blip_df = data.frame(blip=true$blip)
blipdistexample = ggplot(blip_df, aes(x=blip))+geom_density(fill="blue")+
  ggtitle("example of CATE distribution")
ggsave("~/Dropbox/blipCumDist/blipdistexample.png",
       plot = blipdistexample, device = NULL, path = NULL, scale = 1, width = 6,
       height = 4, units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

truth = vapply(t, FUN = function(x0){
  truth = vapply(true$blip, FUN = function(x) {
  i = .5*(x>(x0-h))*(min(x,x0+h)-x0+h)/h
    # upper = (min(x, (x0+h))-x0)/h
    # lower = -1
    # i=(pnorm(upper)-pnorm(lower))*(x>x0-h)/(pnorm(1)-pnorm(-1))
  return(i)}, FUN.VALUE = 1)
  pt = mean(truth)
  return(pt)
  }, FUN.VALUE=1)

truth.05 =truth
target = vapply(t, FUN = function(t) mean(true$blip>t), FUN.VALUE=1)
target

truth1 = truth
truth5 =truth
truth2 =truth
target = vapply(t, FUN = function(t) mean(true$blip>t), FUN.VALUE=1)
target

truth.05f2 = truth
truth.1f2 =truth
truth.15f2 =truth
truth.2f2 =truth

len = length(t)-1
len =11
plotdf = data.frame(t=rep(t[2:12],5), blip.surv = c(whole.curve, truth[2:12],
                                           target[2:12],lower,upper),
                    type = c(rep("ests",len),rep("true.smooth",len),
                             rep("true.target",len),rep("lower.est",len),
                             rep("upper.est",len)))

plotf2 = ggplot(plotdf, aes(x=t,y=blip.surv,color=type))+geom_line()

plotf2 <- plotf2+ labs(title="Estimating whole blip survival curve")
plotf2 = ggdraw(add_sub(plotf2,"does not do well at tails of blip probs.\nWell-specified linear models were estimated", x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                         vpadding = grid::unit(1, "lines"), fontfamily = "", fontface = "plain",
                         colour = "black", size = 8, angle = 0, lineheight = 0.9))
plotf2
ggsave("C:/Users/jlstiles/Dropbox/blipCumDist/wholecurve.png", plot = plotf2, device = NULL, path = NULL, scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

zscoref2 = get.zscore(ff$Dstar,alpha=.05)
zscoref2
sds = apply(ff$Dstar,2,FUN = function(x) sqrt((999)/1000)*sd(x))
ses = sds/sqrt(1000)
ses
lower = whole.curve-zscoref2*ses
upper = whole.curve + zscoref2*ses
lower
upper


hist2 = hist1
h.2 = hist(hist2,breaks=100)
hist3 = hist2
h.3 = hist(hist3,breaks=100)
truth.3
hist05 = ests
truth.05

bw = as.character(c(rep(.05,1000),rep(.2,1000),rep(.3,1000)))
histdf = data.frame(ests = c(hist05,hist2,hist3),bandwidth = bw)
gghist = ggplot(histdf, aes(x=ests,color=bandwidth))+geom_density()
gghist
ggsave("C:/Users/jlstiles/Dropbox/blipCumDist/hists.png", plot = gghist, device = NULL, path = NULL, scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

# ggover

result = gentmle_alt1(initdata=tmledata, estimate_fun = blipdist_estimate,
                      update_fun = blipdist_update,max_iter = 100, t=t, h=h,
                      kernel = kernel)


onestep.estimator <- function(tmledata, t, bandwidth, kernel, ...) {

  nn <- length(tmledata$Y)
  B = tmledata$Q[,"Q1W"]-tmledata$Q[,"Q0W"]
  int = vapply(t, FUN = function(blip) {
    k = function(x) (1/bandwidth)*kernel((x-blip)/bandwidth)
    D1 = vapply(B, FUN = function(x) {
      i = .5*(x>(t-h))*(min(x,t+h)-t+h)/h
      return(i)
      # integrate(k, lower = -1, upper = x, subdivisions = 1000)$value
    }, FUN.VALUE=1)
    return(D1)
  },FUN.VALUE=rep(1,nn))
  est = apply(int, 2, mean)
  k = function(x,tt) (1/bandwidth)*kernel((x-tt)/bandwidth)
  A = tmledata$A
  g1W = tmledata$g1W
  HAW = -vapply(t, FUN = function(x) {
    (B<(x+h)&B>(x-h))*(k(B,x)-1)*(A/g1W-(1-A)/(1-g1W))
    }, FUN.VALUE= rep(1,nn))
  psi = tmledata$ests
  Dstar <- apply(HAW,2,FUN = function(x) with(tmledata, x*(Y-Q[,"QAW"])))+
    int-matrix(rep(est,nn),byrow=TRUE,nrow=nn)
  est = est+mean(Dstar)
  return(est)
}


getres1 = function(t=.5,h=.1,kernel=U,true=.2313561, TG=.240668){
  tmledata=gentmledata(1000)
  res = gentmle_alt1(initdata=tmledata, estimate_fun = blipdist_estimate1,
                     update_fun = blipdist_update,max_iter = 1000, t=t, h=h,
                     kernel = kernel)
  est = res$tmleests
  inits = res$initests
  steps = res$steps
  sigma = sqrt((n-1)/n)*sd(res$Dstar)
  left = est-1.96*sigma/sqrt(n)
  right = est+ 1.96*sigma/sqrt(n)
  cover = left<=true&right>=true
  result= list(est=est,left=left,right=right,init=inits,steps=steps,cover=cover)
  return(list(result))
}

