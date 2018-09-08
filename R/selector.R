#' @export
sim_bwselect1 = function(n, blip, truth, truths_h, bw_seq, g0, Q0, kernel)
{
  tmledata = gentmledata(n, d = 1, g0, Q0, V=10, formu = NULL)
  big_ind = length(bw_seq)
  # undebug(gentmle_alt3)
  test_bwselect = lapply(bw_seq, FUN = function(h) {
    gentmle_alt3(initdata=tmledata, estimate_fun = blipdist_estimate2,
                 update_fun = blipdist_update, max_iter = 1000, kernel = kernel,
                 simultaneous.inference = FALSE, blip = blip, h = h)
  })

  est_vec = unlist(lapply(test_bwselect, FUN = function(x) x$tmleests))
  est_derivs = vapply(2:length(bw_seq), FUN = function(i) {
    (est_vec[i]-est_vec[i-1])/.01
  }, FUN.VALUE = 1)
  # est_derivs
  
  SE_vec = unlist(lapply(1:length(bw_seq), FUN = function(i) {
    x = test_bwselect[[i]]
    sqrt((x$ED2 - x$ED^2)/n)
  }))
  
  SE_derivs = vapply(2:length(SE_vec), FUN = function(i) {
    (SE_vec[i]-SE_vec[i-1])/.01
  }, FUN.VALUE = 1)
  # SE_derivs
  
  # start at a sensible spot, then depending on increasing or decr psi in h
  # we travel until no longer helpful
  for(i in 6:big_ind) {
    ans = all(sign(est_vec[i:(i-4)] - est_vec[(i-1):(i-5)])>=0) |
      all(sign(est_vec[i:(i-4)] - est_vec[(i-1):(i-5)])<=0)
    if (ans) {
      start_h = i-5 
      break
    } else start_h = big_ind
  }

  decre = est_vec[start_h+5] < est_vec[start_h+4]
  for (i in start_h:(big_ind-1)){
    if (decre) {
      end_h = i + 1
      if (!(est_vec[i] >= est_vec[i+1])) {
        end_h = i
        break
      }
    } else {
      end_h = i + 1
      if (!(est_vec[i] <= est_vec[i+1])) {
        end_h = i
        break
      }
    }
  }

  if (decre) ind = order((est_vec - 1.96*SE_vec)[start_h:end_h], decreasing = TRUE)[1] else {
    ind = order((est_vec + 1.96*SE_vec)[start_h:end_h])[1]
  }
  
  
  ind = ifelse((start_h:end_h)[ind] == big_ind, big_ind, (start_h:end_h)[ind]+1)
  
  CI_lower = (est_vec - 1.96*SE_vec)[1:big_ind]
  CI_upper = (est_vec + 1.96*SE_vec)[1:big_ind]
  
  CIs = cbind(est_vec[1:big_ind], CI_lower, CI_upper)
  cover = truth >= CIs[,2] & truth <= CIs[,3]
  
  # do linear regresssion and check positivity of coeffcoefficient
  x = CIs[,1]
  y = bw_seq[1:big_ind]
  incre = coef(lm(y~x))[2] >= 0
  
  ests = CIs[,1]
  if (incre) {
    estsI = pava(y=ests, decreasing = FALSE)
  } else {
    estsI = pava(y=ests, decreasing = TRUE)
  }
  
  SE_vecI = pava(y=SE_vec, decreasing = TRUE)
  
  if (decre) ind_alt = order((est_vec - 1.96*SE_vecI)[start_h:end_h], decreasing = TRUE)[1] else {
    ind_alt = order((est_vec + 1.96*SE_vecI)[start_h:end_h])[1]
  }
  ind_alt = (start_h:end_h)[ind_alt]
  
  CI_lowerI = estsI - 1.96*SE_vecI
  CI_upperI = estsI + 1.96*SE_vecI
  CIsI = data.frame(ests = estsI, lower = CI_lowerI, upper = CI_upperI)
  
  CI_lowerI_alt = ests - 1.96*SE_vecI
  CI_upperI_alt = ests + 1.96*SE_vecI
  CIsI_alt = data.frame(ests = ests, lower = CI_lowerI_alt, upper = CI_upperI_alt)
  
  coverI_alt = truth >= CIsI_alt[,2] & truth <= CIsI_alt[,3]
  
  if (incre) {
    m = order(CIsI[,3])[1]
    indI = max(which(CIsI[,3]==CIsI[m,3]))
    indI = ifelse(indI == big_ind, indI, indI + 1)
  } else {
    m = order(CIsI[,2], decreasing = TRUE)[1]
    indI = max(which(CIsI[,2]==CIsI[m,2]))
    indI = ifelse(indI == big_ind, indI, indI + 1)
  }
  
  CIs_all = cbind(CIs, CIsI, CIsI_alt)
  return(list(CIs_all = CIs_all
              ,CIs = rbind(CIs[ind, ], CIsI_alt[ind_alt,], CIs[indI, ],  CIsI_alt[indI,], CIs[big_ind,])
              ,cover = c(cover[ind], coverI_alt[ind_alt], cover[indI], coverI_alt[indI], cover[big_ind])
              ,truths = c(truth = truth, truth_h = truths_h)
              ,incre = incre
              ,inds = c(ind, ind_alt, indI) 
              ,range = c(start_h, end_h)))
}

#' @export
sim_bwselect = function(n, blip, truth, truths_h, bw_seq, g0, Q0, kernel)
{
  tmledata = gentmledata(n, d = 1, g0, Q0, V=10, formu = NULL)
  big_ind = length(bw_seq)
  # undebug(gentmle_alt3)
  test_bwselect = lapply(bw_seq, FUN = function(h) {
    gentmle_alt3(initdata=tmledata, estimate_fun = blipdist_estimate2,
                 update_fun = blipdist_update, max_iter = 1000, kernel = kernel,
                 simultaneous.inference = FALSE, blip = blip, h = h)
  })
  
  est_vec = unlist(lapply(test_bwselect, FUN = function(x) x$tmleests))
  est_derivs = vapply(2:length(bw_seq), FUN = function(i) {
    (est_vec[i]-est_vec[i-1])/.01
  }, FUN.VALUE = 1)
  # est_derivs
  
  SE_vec = unlist(lapply(1:length(bw_seq), FUN = function(i) {
    x = test_bwselect[[i]]
    sqrt((x$ED2 - x$ED^2)/n)
  }))
  
  SE_derivs = vapply(2:length(SE_vec), FUN = function(i) {
    (SE_vec[i]-SE_vec[i-1])/.01
  }, FUN.VALUE = 1)
  # SE_derivs
  
  # start at a sensible spot, then depending on increasing or decr psi in h
  # we travel until no longer helpful
  for(i in 6:big_ind) {
    ans = all(sign(est_vec[i:(i-4)] - est_vec[(i-1):(i-5)])>=0) |
      all(sign(est_vec[i:(i-4)] - est_vec[(i-1):(i-5)])<=0)
    if (ans) {
      start_h = i-5 
      break
    } else start_h = big_ind
  }
  
  decre = est_vec[start_h+5] < est_vec[start_h+4]
  for (i in start_h:(big_ind-1)){
    if (decre) {
      end_h = i + 1
      if (!(est_vec[i] >= est_vec[i+1])) {
        end_h = i
        break
      }
    } else {
      end_h = i + 1
      if (!(est_vec[i] <= est_vec[i+1])) {
        end_h = i
        break
      }
    }
  }
  
  if (decre) ind = order((est_vec - 1.96*SE_vec)[start_h:end_h], decreasing = TRUE)[1] else {
    ind = order((est_vec + 1.96*SE_vec)[start_h:end_h])[1]
  }
  
  ind = ifelse((start_h:end_h)[ind] == big_ind, big_ind, (start_h:end_h)[ind]+1)
  
  CI_lower = (est_vec - 1.96*SE_vec)[1:big_ind]
  CI_upper = (est_vec + 1.96*SE_vec)[1:big_ind]
  
  CIs = cbind(est_vec[1:big_ind], CI_lower, CI_upper)
  cover = truth >= CIs[,2] & truth <= CIs[,3]
  
  
  # do linear regresssion and check positivity of coeffcoefficient
  x = CIs[,1]
  y = bw_seq[1:big_ind]
  incre = coef(lm(y~x))[2] >= 0
  
  ests = CIs[,1]
  if (incre) {
    estsI = pava(y=ests, decreasing = FALSE)
  } else {
    estsI = pava(y=ests, decreasing = TRUE)
  }
  
  SE_vecI = pava(y=SE_vec, decreasing = TRUE)
  CI_lowerI = estsI - 1.96*SE_vecI
  CI_upperI = estsI + 1.96*SE_vecI
  CIsI = data.frame(ests = estsI, lower = CI_lowerI, upper = CI_upperI)
  if (incre) {
    m = order(CIsI[,3])[1]
    indI = max(which(CIsI[,3]==CIsI[m,3]))
    indI = ifelse(indI == big_ind, indI, indI + 1)
  } else {
    m = order(CIsI[,2], decreasing = TRUE)[1]
    indI = max(which(CIsI[,2]==CIsI[m,2]))
    indI = ifelse(indI == big_ind, indI, indI + 1)
  }
  coverI = truth >= CIs[indI,2] & truth <= CIs[indI,3]
  
  CIs_all = cbind(CIs, CIsI)
  return(list(CIs_all = CIs_all
              ,CIs = rbind(CIs[ind, ], CIs[indI, ], CIs[big_ind,])
              ,cover = c(cover[ind], coverI, cover[big_ind])
              ,truths = c(truth = truth, truth_h = truths_h)
              ,incre = incre
              ,inds = c(ind, indI) 
              ,range = c(start_h, end_h)))
}
