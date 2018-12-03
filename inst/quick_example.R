data("data_example")

head(data_example$Q)
data_example$Y[1:6]
data_example$A[1:6]
data_example$g1W[1:6]

t = seq(m, M, .01)
blip = t[seq(12, 48, 6)]

# make polynomial kernel of order 9
k=make_kernel(degree=12,R=5)

est.info = blipCDF(initdata = data_example, kernel = k, blip = blip, h = .1, 
                   max_iter = 1000, simultaneous.inference = TRUE)

plot(blip, est.info$tmleests)
# tmle estimates
est.info$tmleests
# steps to convergence
est.info$steps
# mean of the influence curves
est.info$ED

# for simultaneous inference, default set to 5%
ci = ci_gentmle(est.info, level = 0.95)
ci
# number of se's used for simultaneous inference at type I error rate of 5%
(ci[1,4]-ci[1,1])/ci[1,2]
