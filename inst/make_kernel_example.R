# The order of the kernel is the degree of the highest order zero moment. Kernels will have support between
# -R and R.  Order must be an odd number because all kernels aere are rthogonal to odd polynomials 
# since they are symmetric
order = 7
R = 2
k=TECDF:::make_kernel(order,R)

# check it is a kernel
area = with(k, integrate(kern, lower = -R, upper = R, R = R, veck = veck, subdivisions = 10000)$value)
area

# plot
s = seq(-R,R,.001)
y = with(k, kern(s, R=R, veck=veck))
plot = plot(s,y)
plot

# check orthogonality to a polynomial of degree less or equal than the order 
test_fcn = as.data.frame(vapply(0:(order+1), FUN = function(r) {
  test_fcn = function(x) (x^r)*with(k, kern(x, R=R, veck = veck))
  test_int = integrate(test_fcn, lower = -R, upper = R,subdivisions = 10000)
  return(c(test_int$abs.error, test_int$value))
}, FUN.VALUE = c(1,1)))
rownames(test_fcn) = c("abs_error", "integral")
colnames(test_fcn) = as.character(0:(order))
# We see the integral of the kernel times an 8th degree polynomial is non trivial
test_fcn

