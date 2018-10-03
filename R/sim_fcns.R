#' @export
sim_bwselect = function(n, blip, bw_seq, g0, Q0, kernel, zscore = NULL)
{
  tmledata = gentmledata(n, d = 1, g0, Q0, V=1, formu = NULL)
  r = length(bw_seq)
  # undebug(gentmle_alt3)
  select_info = lapply(bw_seq, FUN = function(h) {
    info = gentmle_alt(initdata=tmledata, estimate_fun = blipdist_estimate,
                        update_fun = blipdist_update, max_iter = 1000, kernel = kernel,
                        simultaneous.inference = FALSE, blip = blip, h = h)
    return(return(list(est = info$tmleests, IC = info$Dstar, SE = info$ED2)))
  })
  
  ests = vapply(1:r, FUN = function(x) {
    select_info[[x]]$est
  }, FUN.VALUE = 1)
  
  SE = vapply(1:r, FUN = function(x) select_info[[x]]$SE/sqrt(n), FUN.VALUE = 1)
  
  IC = vapply(1:r, FUN = function(x) select_info[[x]]$IC, FUN.VALUE = rep(1,n))
  
  if (is.null(zscore)) return(list(ests = ests, SE = SE)) else {
    zscore = get.zscore(IC, .05) 
    ci = cbind(ests = ests, left = ests - zscore*SE, right = ests + zscore*SE)
    return(list(ci = ci, zscore = zscore))
  }
}

  

  
