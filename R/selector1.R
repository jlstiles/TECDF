#' @export
bwselect_jl = function(CIs, len, plus = TRUE)
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
    SE_vec = (CIs[,3] - CIs[,2])/(2*1.96)
    SE_vecI = pava(y=SE_vec, decreasing = TRUE)
    
    if (decre) ind_plus = order((CIs[,1] - 1.96*SE_vecI)[start_h:end_h], decreasing = TRUE)[1] else {
      ind_plus = order((CIs[,1] + 1.96*SE_vecI)[start_h:end_h])[1]
    }
    
    ind_plus = (start_h:end_h)[ind_plus]
    ind_plus = ifelse(ind_plus==r, r, ind_plus+1)
    if (SE_vecI[ind_plus] >= SE_vec[ind_plus]) {
      CI_plus = c(CIs[ind_plus,1], 
                  CIs[ind_plus,1]-1.96*SE_vecI[ind_plus],
                  CIs[ind_plus,1]+1.96*SE_vecI[ind_plus])
    } else {
      CI_plus = CIs[ind_plus,]
    }
    
    return(list(CI = CIs[ind,], ind = ind, CI_plus = CI_plus, ind_plus = ind_plus))
  } else return(list(CI = CIs[ind,], ind = ind))
}

#' @export
bwselect_m = function(CIs, plus = TRUE)
{
  # start at a sensible spot, then depending on increasing or decr psi in h
  # we travel until no longer helpful
  r = nrow(CIs)
  ests = CIs[,1]
  y = 1:r
  incre = coef(lm(y~ests))[2] >= 0

  SE_vec = (CIs[,3] - CIs[,2])/(2*1.96)
  SE_vecI = pava(y=SE_vec, decreasing = TRUE)
  
  w = (1/SE_vecI)/(sum(1/SE_vecI))
  
  if (incre) {
    estsI = pava(y=ests, decreasing = FALSE)
    estsI_plus = pava(y=ests, decreasing = FALSE,w=w)
  } else {
    estsI = pava(y=ests, decreasing = TRUE)
    estsI_plus = pava(y=ests, decreasing = TRUE,w=w)
  }
  
  right = estsI + 1.96*SE_vecI
  left = estsI - 1.96*SE_vecI
  
  right_plus = estsI_plus + 1.96*SE_vecI
  left_plus = estsI_plus - 1.96*SE_vecI
  
  if (incre) {
    m = order(right)[1]
    ind = max(which(right==right[m]))
  } else {
    m = order(left, decreasing = TRUE)[1]
    ind = max(which(left==left[m]))
  }
  
  if (incre) {
    m_plus = order(right_plus)[1]
    ind_plus = max(which(right_plus==right_plus[m]))
  } else {
    m_plus = order(left_plus, decreasing = TRUE)[1]
    ind_plus = max(which(left_plus==left_plus[m]))
  }
  
  ind = ifelse(ind==r, r, ind+1)
  ind_plus = ifelse(ind_plus==r, r, ind_plus+1)
  if (SE_vec[ind] >= SE_vecI[ind]) CI = CIs[ind,] else {
    CI = c(ests[ind], ests[ind] - 1.96*SE_vecI[ind], ests[ind] + 1.96*SE_vecI[ind])
  }
  
  if (SE_vec[ind_plus] >= SE_vecI[ind_plus]) CI_plus = CIs[ind_plus,] else {
    CI_plus = c(ests[ind_plus], ests[ind_plus] - 1.96*SE_vecI[ind_plus], 
           ests[ind_plus] + 1.96*SE_vecI[ind_plus])
  }
  
  return(list(CI = CI, ind = ind, CI_plus = CI_plus, ind_plus = ind_plus))
}
