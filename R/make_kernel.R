#' @title make_kernel
#' @description constructs polynomial kernels up to order 10
#' @param order, the degree of the 1st non-zero moment, an even number since all these kernels
#' are orthogonal to odd polynomials.  If NULL then a uniform kernel is constructed
#' @param R, support is -R to R and kernel is smooth at the boundary
#' @return  a list containing coefficients of the even polynomial kernel, veck, Range, R,  
#' and functions for the kernel and its cdf, kern and kern_cdf.  
#' @example /inst/make_kernel_example.R 
#' @export
make_kernel = function(order, R){
  if (order/2 != floor(order/2)) stop("order must be even")

  if (is.null(order)) {
    kern = function(x, R, veck) 1/(2*R)*as.numeric(-R <= x & R >= x) 
    kern_cdf = function(x, R, veck) (1/(2*R))*as.numeric(x > -R)*(pmin(x ,R)+R)
    veck = 1
  } else {
    kk = (order+2)/2-2
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
  }
  
  return(list(veck = veck, R = R,  kern = kern, kern_cdf = kern_cdf))
}