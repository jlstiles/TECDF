# making polynomial kernel degree 12, which is an order 10 kernel with support -5 to 5 and 
# smooth at -5 and 5. 
order = 10
R = 5
k=blipCDF:::make_kernel(order,R)

# check it is a kernel
area = with(k, integrate(kern, lower = -R, upper = R, R = R, veck = veck, subdivisions = 10000)$value)
area

# plot
s = seq(-R,R,.001)
y = with(k, kern(s, R=R, veck=veck))
plot = plot(s,y)
plot

# check orthogonality to a polynomial less than the order 
test_fcn = function(x) (x^9+x^8)*with(k, kern(x, R=R, veck = veck))
test_int = integrate(test_fcn, lower = -R, upper = R,subdivisions = 10000)
test_int

