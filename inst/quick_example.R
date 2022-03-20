data("data_example")

head(data_example$Q)
data_example$Y[1:6]
data_example$A[1:6]
data_example$g1W[1:6]

TE = seq(-0.08,0.3, .06)

# make polynomial kernel of order 6. Note, you can only input odd degrees for the kernel which 
# will be the highest degree for which any polynomial of less than or equal to this
# degree will be orthogonal to the kernel
k=make_kernel(order=6,R=5)

est.info = TECDF(initdata = data_example, kernel = k, TE = TE, h = .1, 
                   max_iter = 1000, simultaneous.inference = TRUE)

plot(TE, est.info$tmleests)
# tmle estimates
est.info$tmleests
# steps to convergence
est.info$steps
# mean of the influence curves
est.info$ED

plot(1:415, est.info$risk)
# for simultaneous inference, default set to 5%
ci = ci_gentmle(est.info, level = 0.95)
ci
# number of se's used for simultaneous inference at type I error rate of 5%
(ci[1,4]-ci[1,1])/ci[1,2]
