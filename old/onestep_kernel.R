# library(reshape2)
# library(plyr)
# library(dplyr)
# library(SuperLearner)
library(ggplot2)
library(cowplot)
library(parallel)
# playing with integrating functions
kernel_list = list(gaus2 = function(x) (1/sqrt(2*pi))*exp(-x^2/2),
                   gaus4 = function(x) (1/sqrt(2*pi))*exp(-x^2/2)*(1.5-x^2/2),
                   gaus6 = function(x) (1/8)*(1/sqrt(2*pi))*exp(-x^2/2)*(15-10*x^2+x^4),
                   gaus8 = function(x) (1/48)*(1/sqrt(2*pi))*exp(-x^2/2)*(105-105*x^2+21*x^4-x^6),
                   Unif = function(x) .5*as.numeric(-1<=x&1>=x),
                   fourth = function(x) (0.7031250 - 0.5859375*x^2 + 0.1025391*x^4)*as.numeric(-2<=x&x<=2),
                   sixth = function(x) (0.683593750 - 0.531684028*x^2 + 
                                          0.106336806*x^4 - 0.006188915*x^6)*as.numeric(-3<=x&x<=3),
                   eighth = function(x) (0.897216797 - 1.196289063*x^2 + 0.438639323*x^4 - 
                                           0.060341917*x^6 + 0.002793607*x^8)*as.numeric(-3<=x&x<=3),
                   tenth = function(x) (veck[1] + veck[2]*x^2 + veck[3]*x^4 + 
                                          veck[4]*x^6 + veck[5]*x^8 + veck[6]*x^10)
                   *as.numeric(-3.5<=x&x<=3.5)
                   
                     
)

# The gaussian kernels are useless here, take a while to integrate and 
# then there's round off error so I choose explicit integrals, which is lightning fas
# and such kernels are formulaic to construct as per below--see tenth above

deg = 9
R = 3.5
mm = vapply(seq(1,deg,2), FUN = function(r) {
  vapply(seq(r,(r+deg+1),2), FUN = function(x) 2*R^x, FUN.VALUE = 1)/seq(r,(r+deg+1),2)
}, FUN.VALUE = rep(1,(deg+3)/2))

mm = cbind(mm, vapply(seq(0,deg+1,2), FUN = function(x) R^x, FUN.VALUE = 1)) 
mm = t(mm)
mm_inv = solve(mm)
veck = mm_inv %*% c(1,rep(0, (deg + 1)/2))
veck

unif_cdf = function(x) .5*as.numeric(x > -1)*(pmin(x ,1) + 1)


F0 = list(unif_cdf = unif_cdf, fourth_cdf = fourth_cdf, 
          sixth_cdf = sixth_cdf, eighth_cdf = eighth_cdf, tenth_cdf) 


# integrate to 1?
areas = lapply(kernel_list, FUN = function(x) {
  integrate(x, lower = -Inf, upper = Inf, subdivisions = 10000)$value
} )

areas

# generate_plots
test = integrate(kernel_list[[9]], -3.5, 3.5, subdivisions = 10000)
test

plot_kernels = lapply(kernel_list, FUN = function(x) {
  s = seq(-3.5,3.5,.001)
  y = x(s)
  plot = qplot(x=s,y=y)+geom_line()
  return(plot)
})
plot_kernels[[9]]

# check orthogonality

test_fcn = function(x) (x^2)*kernel_list$tenth(x)
int2 = integrate(test_fcn, lower = -3.5, upper = 3.5,subdivisions = 10000)
int2$abs.error
int2$value

xx = seq(-10,10,length.out = 100)
yy = test_fcn(xx)
plot(xx,yy)

ff6 = function(x) kernel_list[[4]](x)*x^4
int6 = integrate(ff6, lower = -10, upper = 10, subdivisions = 10000)
int6$abs.error
int6$value

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
  return(tmledata)
}

# make a grid of true parameters
true = gendata.blip(1000000)
hist(true$blip, breaks=200)
h_vec = seq(.01, 0.3, by = .01)
t=c(-.3,-.2,-.1,0,.1,.2,.3,.4,.5,.6,.7)

truth = vapply(t, FUN = function(x) mean(true$blip>x), FUN.VALUE = 1)


p.05 = ggplot(plot_df, aes(x = t, y = h0.05, color = names)) + geom_line()
p.10 = ggplot(plot_df, aes(x = t, y = h0.1, color = names)) + geom_line()
p.15 = ggplot(plot_df, aes(x = t, y = h0.15, color = names)) + geom_line()
p.20 = ggplot(plot_df, aes(x = t, y = h0.2, color = names)) + geom_line()
p.25 = ggplot(plot_df, aes(x = t, y = h0.25, color = names)) + geom_line()
p.30 = ggplot(plot_df, aes(x = t, y = h0.3, color = names)) + geom_line()

p.05
p.10
p.15
p.20
p.25
p.30

# note: k is a list with elts, degree and range for polyn kernel and kernel cdf construction
ff= gentmle_alt1(initdata=tmledata, estimate_fun = blipdist_estimate2,
                 update_fun = blipdist_update,max_iter = 1000, t=t, h=h[4],
                 k = k)

ff$steps
ff$initests
ff$tmleests

n=1000
tmledata=gentmledata(n)

k = list()
k$degree = NULL
k$range = 1
h = h_vec[15]
h
res = kernel_plot(t = t, h = h, k = k, truth = truth, n = n, tmledata = tmledata)
res$steps
res$info
res$plot

# debug(gentmle_alt1)
# debug(blipdist_estimate2)
# undebug(kernel_plot)

kernel_plot = function(t, h, k,truth, n, tmledata = NULL) {
  if (is.null(tmledata)) tmledata=gentmledata(n)

  ff= gentmle_alt1(initdata=tmledata, estimate_fun = blipdist_estimate2,
                   update_fun = blipdist_update,max_iter = 1000, t=t, h=h,
                   k = k)
  
  zscore = get.zscore(Dstar = ff$Dstar, .05)
  # zscore
  
  pm = zscore*apply(ff$Dstar, 2, sd)*sqrt(n-1)/n
  
  whole.curve = ff$tmleests
  init.whole.curve = ff$initests
  
  df = data.frame(t = rep(t,3),
                  left = rep(whole.curve-pm,3),
                  right = rep(whole.curve + pm, 3),
                  true = c(whole.curve,truth, init.whole.curve),
                  type = c(rep("est", length(t)), 
                           rep("true", length(t)),
                           rep("init", length(t))
                  )
  )
  
  # df
  if (is.null(k$deg)) {
    kern_name = paste0("unif[" , -k$range, ", ", k$range, "]")
  } else {
    kern_name = paste0("deg " , k$degree, ", range = ", k$range)
  }
  surv_est = ggplot(data = df, aes(x=t,y=true, color = type)) + geom_line()+
    scale_x_continuous(breaks = t)+
    geom_ribbon(data=df,aes(ymax=right,ymin=left),
                fill="gray",colour=NA,alpha=0.5)+labs(y = "S(t)")+
    ggtitle(paste0("estimating CATE ", "'survival'", " curve"),
            subtitle = paste0("bandwidth = ", h, " n = ",n, " kernel: ", kern_name))
  
  capt = "S(t) = prob CATE exceeds t."
  surv_est=ggdraw(add_sub(surv_est,capt, x= 0, y = 0.8, hjust = 0, vjust = 0.5,
                          vpadding = grid::unit(1, "lines"), fontfamily = "",
                          fontface = "plain",colour = "black", size = 12, angle = 0,
                          lineheight = 0.9))
  left = df$left[1:length(t)] 
  right = df$right[1:length(t)]
  cover = truth <= right & truth >= left
  info = data.frame(est = whole.curve, left = left, right =right, 
                    truth = truth,
                    cover = cover)
  
  
  return(list(info = info, steps = ff$steps, plot = surv_est))
}

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

