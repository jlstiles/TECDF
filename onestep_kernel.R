library(reshape2)
library(plyr)
library(dplyr)
library(mvtnorm)
library(ggplot2)
library(foreach)
library(doSNOW)
library(SuperLearner)
library(parallel)
library(reshape)
library(plyr)
library(cowplot)
# library(origami)
library(geepack)
library(mvtnorm)
# playing with integrating functions
f2 = function(x) (1/sqrt(2*pi))*exp(-x^2/2)
f4 = function(x) (1/sqrt(2*pi))*exp(-x^2/2)*(1.5-x^2/2)
f6 = function(x) (1/8)*(1/sqrt(2*pi))*exp(-x^2/2)*(15-10*x^2+x^4)
f8 = function(x) (1/48)*(1/sqrt(2*pi))*exp(-x^2/2)*(105-105*x^2+21*x^4-x^6)

s2 = integrate(f2, lower = -100, upper = 100, subdivisions = 10000)$value
s4 = integrate(f4, lower = -100, upper = 100, subdivisions = 10000)$value
s6 = integrate(f6, lower = -100, upper = 100, subdivisions = 10000)$value
s8 = integrate(f8, lower = -100, upper = 100, subdivisions = 10000)$value


f2 = function(x) (1/sqrt(2*pi))*exp(-x^2/2)
f4 = function(x) (1/sqrt(2*pi))*exp(-x^2/2)*(1.5-x^2/2)*
f6 = function(x) (1/8)*(1/sqrt(2*pi))*exp(-x^2/2)*(15-10*x^2+x^4)
f8 = function(x) (1/48)*(1/sqrt(2*pi))*exp(-x^2/2)*(105-105*x^2+21*x^4-x^6)

ff2 = function(x) x^3*f8(x)
int2 = integrate(ff2, lower = -10, upper = 10, subdivisions = 10000)$value
int2

ff6 = function(x) f8(x)*x^6
int6 = integrate(ff6, lower = -10, upper = 10, subdivisions = 10000)$value
int6

x = seq(-5,5,.001)
y = f8(x)
plot(x,y)

ttt = function(x) f8(x)*x^6

f1 = function(x) (2-2*(1-x))
f2 = function(x) 3-8*(1-x)+4*(1-x)^2

library(orthopolynom)
chebs = chebyshev.u.polynomials(10, normalized=TRUE)
chebs
polyL =  legendre.polynomials(10,normalized=TRUE)
polyL

g6 = function(x) (x^6-(126/99)*x^4+(35/99)*x^2)*10395/128
g4 =function(x) (105/8)*((5/7)*x^2-x^4)
ff = function(x) UU(x)*x^2
g2 = function(x) 1.5*x^2
integrate(ff, lower =-1, upper = 1, subdivisions = 1000)$value

x=seq(-1,1,.01)
plot(x,g4(x))

poly = chebyshev.c.polynomials(3, normalized=FALSE)

chebyshev.c.inner.products(2)

ff = function(x) poly[[4]]
ff(222)

x=seq(-1,1,.001)
plot(x,f2(x))

dnorm(.2)/(pnorm(1)-pnorm(-1))
f2(.2)

U = function(x) .5*as.numeric(-1<=x&1>=x)

UU = function(x) (1-x^2)*as.numeric(-1<x&x<1)*.75

UU(-1)




g0=function(W1,W2,W3,W4) {plogis(-.8*W1+.39*W2+.08*W3-.12*W4-.15)}
# g0=function(W1,W2,W3,W4) {plogis(-.28*W1+1*W2+.08*W3-.12*W4-1)}
# Q0=function(A,W1,W2,W3,W4) {plogis(3*A-1*W1*A+1.2*W2-1.9*W3+.4*W4*A*W2-cos(W4))}
# Q0=function(A,W1,W2,W3,W4) {plogis(.5*A-1*W1+1.2*W2-1.9*W3+.4*W4)}
Q0=function(A,W1,W2,W3,W4) {plogis(A-1*W1+1.2*W2-1.9*W3+.4*W4+A*W2+2*A*W1)}
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

# undebug(blipdist_update)
# undebug(blipdist_estimate)
# undebug(gentmle_alt1)

get.zscore = function(Dstar, alpha) {

  if (ncol(Dstar)==1) return(1.96544) else{
  sigma = cor(Dstar)
  means = rep(0,ncol(sigma))
  zs = rmvnorm(n=2000000,mean= means, sigma=sigma, method= "chol")
  zabs = apply(zs,1,FUN = function(x) max(abs(x)))
  zscore = quantile(zabs, probs = 1-alpha)
  return(zscore)}
}

n=10000
tmledata=gentmledata(n)
kernel = UU

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

true = gendata.blip(1000000)
hist(true$blip, breaks=200)
truth_h = vapply(t, FUN = function(x0){
  truth = vapply(true$blip, FUN = function(x) {
    i = .5*(x>(x0-h))*(min(x,x0+h)-x0+h)/h
    # upper = (min(x, (x0+h))-x0)/h
    # lower = -1
    # i=(pnorm(upper)-pnorm(lower))*(x>x0-h)/(pnorm(1)-pnorm(-1))
    return(i)}, FUN.VALUE = 1)
  pt = mean(truth)
  return(pt)
}, FUN.VALUE=1)

truth = vapply(t, FUN = function(x) {
  mean(true$blip>x)
}, FUN.VALUE = 1)


df = data.frame(t = rep(t,3),
                left = rep(whole.curve-pm,3),
                right = rep(whole.curve + pm, 3),
                truth = c(whole.curve,truth_h,truth),
                type = c(rep("est", length(t)),rep("smoothed",length(t)), rep("true", length(t))))

df
surv_est = ggplot(data = df, aes(x=t,y=truth, color = type)) + geom_line()+
  geom_ribbon(data=df,aes(ymax=right,ymin=left),
              fill="gray",colour=NA,alpha=0.5)+labs(y = "S(t)")+
  ggtitle(paste0("estimating CATE ", "'survival'", " curve"),
          subtitle = "bandwidth = 0.1, Uniform[-1,1] kernel, n = 10000")

capt = paste0("S(t) = prob CATE exceeds t.\n",
              "'smoothed' means the true parameter for bandwidth 0.1\n",
              "'true' is the true S(t) and 'est' is the estimate of smoothed parameter\n",
              "Shaded region represents simultaneous 95% CI to cover the whole curve\n",
              "both pscore and outcome model were well-specified")
surv_est=ggdraw(add_sub(surv_est,capt, x= 0, y = 0.8, hjust = 0, vjust = 0.5,
                     vpadding = grid::unit(1, "lines"), fontfamily = "",
                     fontface = "plain",colour = "black", size = 12, angle = 0,
                     lineheight = 0.9))
surv_est

ggsave("~/Dropbox/Quals/v0321/surv_est.jpg", plot = surv_est,
       device = NULL, path = NULL, scale = 1, width = NA, height = NA,
       units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE)

ff$initests
ff$steps
colMeans(ff$Dstar)

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
blipdistexample = ggplot(blip_df, aes(x=blip))+geom_density(fill="blue")+ggtitle("example of blip distribution")
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

