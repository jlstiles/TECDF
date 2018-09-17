#' @export
bwselect_jl = function(CIs, len, plus = TRUE, z_alpha = 1.96)
{
  # start at a sensible spot, then depending on increasing or decr psi in h
  # we travel until no longer helpful
  r = nrow(CIs)
  for(i in (len+1):r) {
    ans = all(sign(CIs[i:(i + 1 - len),1] - CIs[(i-1):(i - len),1]) >= 0) |
      all(sign(CIs[i:(i-4),1] - CIs[(i-1):(i-5),1]) <= 0)
    if (ans) {
      start_h = i-5 
      break
    } else start_h = r
  }
  
  if (start_h != r) {
  decre = CIs[start_h+5,1] < CIs[start_h+4,1]
  for (i in start_h:(r-1)){
    if (decre) {
      end_h = i + 1
      if (!(CIs[i,1] >= CIs[(i+1),1])) {
        end_h = i
        break
      }
    } else {
      end_h = i + 1
      if (!(CIs[i,1] <= CIs[(i+1),1])) {
        end_h = i
        break
      }
    }
  }
  
  if (decre) ind = order(CIs[start_h:end_h,2], decreasing = TRUE)[1] else {
    ind = order(CIs[start_h:end_h,3])[1]
  }
  ind = ifelse((start_h:end_h)[ind] == r, r, (start_h:end_h)[ind]+1)
  
  if (plus) {
    SE_vec = (CIs[,3] - CIs[,2])/(2*z_alpha)
    SE_vecI = pava(y=SE_vec, decreasing = TRUE)
    
    if (decre) ind_plus = order((CIs[,1] - z_alpha*SE_vecI)[start_h:end_h], decreasing = TRUE)[1] else {
      ind_plus = order((CIs[,1] + z_alpha*SE_vecI)[start_h:end_h])[1]
    }
    
    ind_plus = (start_h:end_h)[ind_plus]
    ind_plus = ifelse(ind_plus==r, r, ind_plus+1)
    if (SE_vecI[ind_plus] >= SE_vec[ind_plus]) {
      CI_plus = c(CIs[ind_plus,1], 
                  CIs[ind_plus,1]-z_alpha*SE_vecI[ind_plus],
                  CIs[ind_plus,1]+z_alpha*SE_vecI[ind_plus])
    } else {
      CI_plus = CIs[ind_plus,]
    }
    
    return(list(CI = CIs[ind,], ind = ind, CI_plus = CI_plus, ind_plus = ind_plus))
  } else return(list(CI = CIs[ind,], ind = ind))
  } else {
    if (plus) {
      return(list(CI = CIs[r,], ind = r, CI_plus = CIs[r,], ind_plus = r))
    } else return(list(CI = CIs[r,], ind = r))
  } 
}

#' @export
bwselect_m = function(CIs, truth = NULL, SE_true = NULL, z_alpha = 1.96)
{
  
  # setting some overall parameters
  r = nrow(CIs)
  y = 1:r
  ests = CIs[,1]
  incre = coef(lm(y~ests))[2] >= 0
  
  # compute SE's
  SEs = (CIs[,3] - CIs[,2])/(2*z_alpha)
  SEsI = pava(y=SEs, decreasing = TRUE)
  
  w = 1/SEsI/(sum(1/SEsI))
  
  #  computing various estimates
  estsI_wt = pava(y=ests, decreasing = !incre, w=w)
  CIs_var = ci_form(ests, SEsI, z_alpha)
  CIsI_var = ci_form(estsI_wt, SEsI, z_alpha)
  
  if (!is.null(truth)) {
    incre0 = truth$S0 <= truth$S0h
    if (incre == incre0) {
      estsI0_wt =  estsI_wt
      } else estsI0_wt = pava(y=ests, decreasing = !incre0, w=w)
      
    CIsI0_var = ci_form(estsI0_wt, SEsI, z_alpha)
    
    ind_I0_var = ind_choose(CIsI0_var, incre0)
    ind_I0_var = ifelse(ind_I0_var==r, r, ind_I0_var+1)
    
    if (SEsI[ind_I0_var] >= SEs[ind_I0_var]) {
      CI_I0_var = CIs_var[ind_I0_var,]
    } else {
      CI_I0_var = CIs[ind_I0_var,]
    }
    
  }
  
  if (!is.null(SE_true)) {
    w0 = 1/SE_true/(sum(1/SE_true))
    estsI_wt0 = pava(y=ests, decreasing = !incre, w=w0)
    if (incre == incre0) {
      estsI0_wt0 = estsI_wt0
      } else estsI0_wt0 = pava(y=ests, decreasing = !incre0, w=w0)
    
    CIsI_var0 = ci_form(estsI_wt0, SE_true, z_alpha)
    CIsI0_var0 = ci_form(estsI0_wt0, SE_true, z_alpha)
    CIs_var0 = ci_form(ests, SE_true, z_alpha)
    
    ind_I_var0 = ind_choose(CIsI_var0, incre)
    ind_I0_var0 = ind_choose(CIsI0_var0, incre0)
    
    CI_I_var0 = CIsI_var0[ind_I_var0,]
    CI_I0_var0 = CIsI0_var0[ind_I0_var0,]
    
    CI_jl_var0_info = bwselect_jl(CIs_var0, len = 5, plus = FALSE, z_alpha = z_alpha)
    ind_jl_var0 = ifelse(CI_jl_var0_info$ind==r, r, CI_jl_var0_info$ind-1)
    CI_jl_var0 = CIs_var0[ind_jl_var0,]
  }
  
  # get jl method estimates
  CI_jl_info = bwselect_jl(CIs, len = 5, plus = TRUE, z_alpha = z_alpha)
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


