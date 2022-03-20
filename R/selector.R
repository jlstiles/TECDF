#' @export
bwselect_jl = function(ests, SEs, len, plus = TRUE, z_alpha = 1.96)
{
  # form the CI's' according to z_alpha
  CIs = ci_form(ests, SEs, z_alpha)
  # start at a sensible spot, then depending on increasing or decr psi in h
  # we travel until no longer helpful
  
  r = nrow(CIs)
  for(i in len:r) {
    ans = all(sign(CIs[i:(i - len+2), 1] - CIs[(i - 1):(i - len+1), 1]) >= 0) |
      all(sign(CIs[i:(i - len+2), 1] - CIs[(i-1):(i - len+1), 1]) <= 0)
    if (ans) {
      start_h = i-len+1 
      decre = CIs[start_h+1, 1] < CIs[start_h, 1] 
      no_mono = FALSE
      end_h = i
      j=i
      continue = TRUE
      while (continue &  j<r) {
        if (decre) continue = CIs[(j+1), 1] - CIs[j, 1] < 0 else {
          continue = CIs[(j+1), 1] - CIs[j, 1] >= 0
        }
        if (continue) end_h = j+1
        j=j+1
      }
      break
    } else no_mono = TRUE
  }
  
  if (no_mono) {
    ind = r
    if (plus) {
      ind_plus = r
      CI_plus = CIs[ind_plus,]
    }
  } else {
    if (decre) ind = order(CIs[start_h:end_h,2], decreasing = TRUE)[1] else {
      ind = order(CIs[start_h:end_h,3])[1]
    }
    ind = (start_h:end_h)[ind]
    if (plus) {
      SE_vec = ((CIs[,3] - CIs[,2])/(2*z_alpha))[start_h:end_h]
      SE_vecI = pava(y=SE_vec, decreasing = TRUE)
      
      if (decre) ind_temp = order((CIs[start_h:end_h,1] - z_alpha*SE_vecI), decreasing = TRUE)[1] else {
        ind_temp = order((CIs[start_h:end_h,1] + z_alpha*SE_vecI))[1]
      }
      ind_plus = (start_h:end_h)[ind_temp]
      if (SE_vecI[ind_temp] >= SE_vec[ind_temp]) {
        CI_plus = c(CIs[ind_plus,1], 
                    CIs[ind_plus,1]-z_alpha*SE_vecI[ind_temp],
                    CIs[ind_plus,1]+z_alpha*SE_vecI[ind_temp])
      } else {
        CI_plus = CIs[ind_plus,]
      }
    }
  }
  if (plus) {
    return(list(CI = CIs[ind,], ind = ind, CI_plus = CI_plus, ind_plus = ind_plus))
  } else return(list(CI = CIs[ind,], ind = ind))
}


bwselect_m = function(ests, SEs, truth = NULL, SE_true = NULL, z_alpha = 1.96)
{
  
  # setting some overall parameters
  r = length(ests)
  x = 1:r
  CIs = ci_form(ests, SEs, z_alpha)
  incre = coef(lm(ests~x))[2] >= 0
  
  # compute SE's
  SEsI = pava(y=SEs, decreasing = TRUE)
  
  w = 1/SEsI/(sum(1/SEsI))
  
  #  computing various estimates
  estsI_wt = pava(y=ests, decreasing = !incre, w=w)
  
  # CIs with isotonic SE and true ests
  CIs_var = ci_form(ests, SEsI, z_alpha)
  # CIs with isotonic ests and istonics SEs
  CIsI_var = ci_form(estsI_wt, SEsI, z_alpha)
  
  # if we know the true params, use the info to determine incre
  if (!is.null(truth)) {
    incre0 = truth$S0 <= truth$S0h
    if (incre == incre0) {
      estsI0_wt =  estsI_wt
    } else estsI0_wt = pava(y=ests, decreasing = !incre0, w=w)
    
    # CIs with isotonic ests from true incre indication and isotonic SEs      
    CIsI0_var = ci_form(estsI0_wt, SEsI, z_alpha)
    
    # choose the index based on those CIs
    ind_I0_var = ind_choose(CIsI0_var, incre0)
    # bump it up 1
    ind_I0_var = ifelse(ind_I0_var==r, r, ind_I0_var+1)
    
    # if isotonic SE is bigger than estimated SE then choose that CI, otherwise choose the estimated CI
    if (SEsI[ind_I0_var] >= SEs[ind_I0_var]) {
      CI_I0_var = CIs_var[ind_I0_var,]
    } else {
      CI_I0_var = CIs[ind_I0_var,]
    }
    
  }
  
  # If we know the true SE use it
  if (!is.null(SE_true)) {
    w0 = 1/SE_true/(sum(1/SE_true))
    # istonize the ests using true vars as weights and the estimated incre
    estsI_wt0 = pava(y=ests, decreasing = !incre, w=w0)
    # isotonize using the true incr and true vars as wts
    if (incre == incre0) {
      estsI0_wt0 = estsI_wt0
    } else estsI0_wt0 = pava(y=ests, decreasing = !incre0, w=w0)
    
    # The below two CIs might be the same in incre = incre0
    CIsI_var0 = ci_form(estsI_wt0, SE_true, z_alpha)
    CIsI0_var0 = ci_form(estsI0_wt0, SE_true, z_alpha)
    
    # CIs with ests using true var
    CIs_var0 = ci_form(ests, SE_true, z_alpha)
    
    ind_I_var0 = ind_choose(CIsI_var0, incre)
    ind_I0_var0 = ind_choose(CIsI0_var0, incre0)
    
    CI_I_var0 = CIsI_var0[ind_I_var0,]
    CI_I0_var0 = CIsI0_var0[ind_I0_var0,]
    
    CI_jl_var0_info = bwselect_jl(CIs_var0[,1], SE_true, len = 5, plus = FALSE, z_alpha = z_alpha)
    ind_jl_var0 = ifelse(CI_jl_var0_info$ind==r, r, CI_jl_var0_info$ind-1)
    CI_jl_var0 = CIs_var0[ind_jl_var0,]
  }
  
  # get jl method estimates
  CI_jl_info = bwselect_jl(ests = ests, SEs = SEs, len = 5, plus = TRUE, z_alpha = z_alpha)
  CI_jl = CI_jl_info$CI
  CI_jl_var = CI_jl_info$CI_plus
  ind_jl = CI_jl_info$ind
  ind_jl_var = CI_jl_info$ind_plus
  
  ind_I_var = ind_choose(CIsI_var, incre)
  ind_I_var = ifelse(ind_I_var==r, r, ind_I_var+1)
  
  if (SEsI[ind_I_var] >= SEs[ind_I_var]) {
    CI_I_var = CIs_var[ind_I_var,]
  } else {
    CI_I_var = CIs[ind_I_var,]
  }
  
  CIr = CIs[r,]
  if (!is.null(SE_true) & !is.null(truth)) {
    return(list(CI_jl = CI_jl, ind_jl = ind_jl, 
                CI_jl_var = CI_jl_var, ind_jl_var = ind_jl_var,
                CI_jl_var0 = CI_jl_var0, ind_jl_var0 = ind_jl_var0,
                CI_I_var = CI_I_var, ind_I_var = ind_I_var,
                CI_I0_var = CI_I_var, ind_I0_var = ind_I0_var,
                CI_I_var0 = CI_I_var0, ind_I_var0 = ind_I_var0,
                CI_I0_var0 = CI_I0_var0, ind_I0_var0 = ind_I0_var0,
                CIr = CIr,
                mono = incre==incre0))
  } else {
    return(list(CI_jl = CI_jl, ind_jl = ind_jl, 
                CI_jl_var = CI_jl_var, ind_jl_var = ind_jl_var,
                CI_I_var = CI_I_var, ind_I_var = ind_I_var,
                CIr = CIr,
                mono = incre==incre0))
  }
}


