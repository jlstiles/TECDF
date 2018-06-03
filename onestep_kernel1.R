# library(reshape2)
# library(plyr)
# library(dplyr)
# library(SuperLearner)
library(ggplot2)
library(cowplot)
library(parallel)

# The gaussian kernels are useless here, take a while to integrate and 
# then there's round off error so I choose explicit integrals, which is lightning fas
# and such kernels are formulaic to construct as per below--see tenth above

unif_cdf = function(x) .5*as.numeric(x > -1)*(pmin(x ,1) + 1)
Unif = function(x) .5*as.numeric(-1<=x&1>=x)

deg = 11
R = 5
mm = vapply(seq(1,deg,2), FUN = function(r) {
  vapply(seq(r,(r+deg+1),2), FUN = function(x) 2*R^x, FUN.VALUE = 1)/seq(r,(r+deg+1),2)
}, FUN.VALUE = rep(1,(deg+3)/2))

mm = cbind(mm, vapply(seq(0,deg+1,2), FUN = function(x) R^x, FUN.VALUE = 1)) 
mm = t(mm)
mm_inv = solve(mm)
veck = mm_inv %*% c(1,rep(0, (deg + 1)/2))

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

# check orthogonality

test_fcn = function(x) (x^10)*kern(x, R=R, veck = veck)
test_int = integrate(test_fcn, lower = -R, upper = R,subdivisions = 10000)
test_int$abs.error
test_int$value

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

n=1000
tmledata=gentmledata(n)

k = list()
k$degree = 11
k$range = 5
h = h_vec[20]
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

