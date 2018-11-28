#' @export 
get.zscore = function(Dstar, alpha) {
  
  if (ncol(Dstar)==1) return(1.96544) else{
    sigma = cor(Dstar)
    means = rep(0,ncol(sigma))
    zs = rmvnorm(n=2000000,mean= means, sigma=sigma, method= "chol")
    zabs = apply(zs,1,FUN = function(x) max(abs(x)))
    zscore = quantile(zabs, probs = 1-alpha)
    return(zscore)}
}

#' @export 
gendata.blip=function(n, d, g0, Q0){
  # d=1
  # n=10
  var_names = paste0("W", 1:d)
  W_list = lapply(1:d, FUN = function(i) rnorm(n))
  names(W_list) = var_names
  pscores = do.call(g0, W_list)
  A=rbinom(n,1,pscores)
  AW_list = W_list
  AW_list$A = A
  Q = do.call(Q0, AW_list)
  Y=rbinom(n,1,Q)
  AW_list1 = AW_list0 = AW_list
  AW_list1$A = rep(1,n)
  AW_list0$A = rep(0,n)
  blip = do.call(Q0, AW_list1) - do.call(Q0, AW_list0)
  df = cbind(A = A, as.data.frame(W_list), Y = Y)
  return(list(df=df,blip=blip))
}

#' @export 
get.truth = function(t, h, kernel, d, g0, Q0) {
  # create the kernel according to specs
  B = gendata.blip(3e6, d, g0, Q0)$blip
  vapply(t, FUN = function(x0) {
    t(vapply(h, FUN = function(bw) {
      int = with(kernel, kern_cdf(x=(B - x0)/bw, R = R, veck = veck))
      truth_h = mean(int)
      truth = mean(B>=x0)
      return(c(truth, truth_h))
    }, FUN.VALUE = c(1,1)))
  }, FUN.VALUE = matrix(rep(1, 2*length(h)), ncol=2))
}



#' @export 
CATEsurv_plot = function(t, h, k,truth, n, tmledata = NULL, fun, ...) {
  if (is.null(tmledata)) tmledata=gentmledata(n=n, d=d, g0=g0, Q0=Q0, formu=formu)
  if (length(t)>1) simultaneous.inference = TRUE else simultaneous.inference = FALSE
  
  ff= fun(initdata=tmledata, estimate_fun = blipdist_estimate2,
                   update_fun = blipdist_update, max_iter = 1000, t=t, h=h,
                   k = k, simultaneous.inference = simultaneous.inference)
  

  ci = ci_gentmle(ff, level = 0.95)
  whole.curve = ff$tmleests
  init.whole.curve = ff$initests
  
  if (length(t)>1) {
    df = data.frame(t = rep(t,3),
                    left = rep(ci[,3],3),
                    right = rep(ci[,4],3),
                    true = c(whole.curve,truth, init.whole.curve),
                    type = c(rep("est", length(t)), 
                             rep("true", length(t)),
                             rep("init", length(t))
                    )
    )
    
    # df
    if (is.null(k$deg)) {
      kern_name = paste0("unif[" , -k$range, ", ", k$range, "]")
    } else {
      kern_name = paste0("deg " , k$degree, ", range = ", k$range)
    }
    # surv_est = ggplot2::ggplot(data = df, aes(x=t,y=true, color = type)) + geom_line()+
    #   scale_x_continuous(breaks = t)+
    #   geom_ribbon(data=df,aes(ymax=right,ymin=left),
    #               fill="gray",colour=NA,alpha=0.5)+labs(y = "S(t)")+
    #   ggtitle(paste0("estimating CATE ", "'survival'", " curve"),
    #           subtitle = paste0("bandwidth = ", h, " n = ",n, " kernel: ", kern_name))
    # 
    # capt = "S(t) = prob CATE exceeds t."
    # surv_est=cowplot::ggdraw(add_sub(surv_est,capt, x= 0, y = 0.8, hjust = 0, vjust = 0.5,
    #                         vpadding = grid::unit(1, "lines"), fontfamily = "",
    #                         fontface = "plain",colour = "black", size = 12, angle = 0,
    #                         lineheight = 0.9))
    left = df$left[1:length(t)] 
    right = df$right[1:length(t)]
    cover = truth <= right & truth >= left
    info = data.frame(est = whole.curve, left = left, right =right, 
                      init = init.whole.curve,
                      truth = truth,
                      cover = cover)
    out = list(info = info, steps = ff$steps)
  } else {
    df = data.frame(t = rep(t,3),
                    left = rep(ci[,3],3),
                    right = rep(ci[,4],3),
                    true = c(whole.curve,truth, init.whole.curve),
                    type = c(rep("est", length(t)), 
                             rep("true", length(t)),
                             rep("init", length(t))
                    )
    )
    
    left = df$left[1:length(t)] 
    right = df$right[1:length(t)]
    cover = truth <= right & truth >= left
    info = data.frame(est = whole.curve, left = left, right =right, 
                      init = init.whole.curve,
                      truth = truth,
                      cover = cover)
    out = list(info = info, steps = ff$steps)
  }
  return(out)
}

#' @export 
gentmledata = function(n, d, g0, Q0, V, formu = NULL) {
  # n=100
  # d=1
  # V=10
  data = gendata.blip(n, d, g0, Q0)$df
  data0 = data
  data0$A = 0
  data1 = data
  data1$A = 1
  datag = data
  datag$Y = NULL
  
  if (is.null(formu)) {
    covs = colnames(data)[!colnames(data) %in% c("A","Y")]
    formuQ = formula(paste0("Y ~ ", paste0("A*(", paste(covs, "", collapse = "+"), ")")))
    formug = formula("A ~.")
  } else  {
    formuQ = formu$Q
    formug = formu$g
  }
  
  if (V == 1) {
    fitQ = glm(formuQ,data=data, family = "binomial")
    Q0W = predict(fitQ,newdata = data0, type = 'response')
    Q1W = predict(fitQ, newdata = data1, type = 'response')
    QAW = predict(fitQ, newdata = data, type= 'response')
    
    Q = cbind(QAW, Q0W, Q1W)
    datag = data
    datag$Y = NULL
    fitg = glm(formug, data=datag, family='binomial')
    g1W = predict(fitg, type = 'response')
    
    tmledata = list(Q=Q,Y=data$Y,A=data$A,g1W=g1W)
  } else {
    folds = make_folds(n=n, V=V)
    fold_preds = lapply(folds, FUN = function(fold) {
      fitQ = glm(formuQ,data=data[fold$training_set,], family = "binomial")
      Q0W = predict(fitQ,newdata = data0[fold$validation_set,], type = 'response')
      Q1W = predict(fitQ, newdata = data1[fold$validation_set,], type = 'response')
      QAW = predict(fitQ, newdata = data[fold$validation_set,], type= 'response')
      
      Q = cbind(QAW, Q0W, Q1W)
      fitg = glm(formug, data=datag[fold$training_set,], family='binomial')
      g1W = predict(fitg, newdata = datag[fold$validation_set, ], type = 'response')
      
      return(list(Q=Q,Y=data$Y[fold$validation_set],A=data$A[fold$validation_set],g1W=g1W))
    })
    tmledata = list(Q = do.call(rbind, lapply(fold_preds, FUN = function(x) x$Q)),
                    Y = unlist(lapply(fold_preds, FUN = function(x) x$Y)),
                    A = unlist(lapply(fold_preds, FUN = function(x) x$A)),
                    g1W = unlist(lapply(fold_preds, FUN = function(x) x$g1W)))
  }
  return(tmledata)
}

#' @export 
gentmledata_hal = function(n, d, g0, Q0, V, RCT = FALSE, formu = NULL) {
  # prepare inputs for hal
  info = gendata.blip(n, d, g0, Q0)
  data = info$df
  blip = info$blip
  
  X = data
  X$Y = NULL
  Y = data$Y
  X0 = X
  X0$A = 0
  X1 = X
  X1$A = 1
  
  # define the formula for glm
  if (is.null(formu)) {
    covs = colnames(data)[!colnames(data) %in% c("A","Y")]
    formuQ = formula(paste0("Y ~ ", paste0("A*(", paste(covs, "", collapse = "+"), ")")))
    formug = formula("A ~.")
  } else  {
    formuQ = formu$Q
    formug = formu$g
  }
  
  # 1 as in fitQ1 as opposed to fitQ means glm instead of hal
  
  if (V == 1) {
    # fitQ = glm(formuQ,data=data, family = "binomial")
    fitQ = fit_hal(X = X, Y = Y, degrees = NULL, fit_type = "glmnet",
                   n_folds = 10, use_min = TRUE, family = "binomial",
                   return_lasso = FALSE, yolo = TRUE)
    Q0W = predict(fitQ, new_data = X0, type = 'response')
    Q1W = predict(fitQ, new_data = X1, type = 'response')
    QAW = predict(fitQ, new_data = X, type = 'response')
    
    Q = cbind(QAW, Q0W, Q1W)
    
    fitQ1 = glm(formula = formuQ, data = data, family = 'binomial')
    Q0W = predict(fitQ1, newdata = X0, type = 'response')
    Q1W = predict(fitQ1, newdata = X1, type = 'response')
    QAW = predict(fitQ1, type = 'response')
    
    Q1 = cbind(QAW, Q0W, Q1W)
    
    if (RCT) {
      datag = X
      fitg = glm(formug, data=datag, family='binomial')
      g1W = predict(fitg, type = 'response')
      g1W1 = g1W
    } else {
      Xg = X
      A = Xg$A
      Xg$A = NULL
      fitg = fit_hal(X = Xg, Y = A, degrees = NULL, fit_type = "glmnet",
                     n_folds = 10, use_min = TRUE, family = "binomial",
                     return_lasso = FALSE, yolo = TRUE)
      g1W = predict(fitg, new_data = Xg, type = 'response')
      
      fitg1 = glm(formula = formug, data = X, family = 'binomial')
      g1W1 = predict(fitg1, type = 'response')
    }
    tmledata = list(Q=Q,Y=data$Y,A=data$A,g1W=g1W)
    tmledata1 = list(Q=Q1,Y=data$Y,A=data$A,g1W=g1W1)
  } else {
    folds = make_folds(n=n, V=V)
    fold_preds = lapply(folds, FUN = function(fold) {
      fitQ = fit_hal(X = X[fold$training_set,], Y = Y[fold$training_set], 
                     degrees = NULL, fit_type = "glmnet",
                     n_folds = 10, use_min = TRUE, family = "binomial",
                     return_lasso = FALSE, yolo = TRUE)
      Q0W = predict(fitQ,new_data = X0[fold$validation_set,], type = 'response')
      Q1W = predict(fitQ, new_data = X1[fold$validation_set,], type = 'response')
      QAW = predict(fitQ, new_data = X[fold$validation_set,], type= 'response')
      
      Q = cbind(QAW, Q0W, Q1W)
      
      fitQ1 = glm(formula = formuQ, data = data[fold$training_set,], family = 'binomial')
      Q0W = predict(fitQ1, newdata = X0[fold$validation_set,], type = 'response')
      Q1W = predict(fitQ1, newdata = X1[fold$validation_set,], type = 'response')
      QAW = predict(fitQ1, newdata = X[fold$validation_set,], type = 'response')
      
      Q1 = cbind(QAW, Q0W, Q1W)
      
      if (RCT) {
        datag = X[fold$training_set,]
        fitg = glm(formug, data=datag, family='binomial')
        g1W = predict(fitg, type = 'response')
        g1W1 = g1W
      } else {
        Xg = X[fold$training_set,]
        A = Xg$A
        Xg$A = NULL
        newX = X
        newX$A = NULL
        fitg = fit_hal(X = Xg, Y = A, degrees = NULL, fit_type = "glmnet",
                       n_folds = 10, use_min = TRUE, family = "binomial",
                       return_lasso = FALSE, yolo = TRUE)
        g1W = predict(fitg, new_data = newX[fold$validation_set,], type = 'response')
        
        fitg1 = glm(formula = formug, data = X[fold$training_set,], family = 'binomial')
        g1W1 = predict(fitg1, newdata = X[fold$validation_set, ], type = 'response')
      }
      
      return(list(Q=Q,Y=data$Y[fold$validation_set],A=data$A[fold$validation_set],g1W=g1W,
                  Q1=Q1,g1W1=g1W1, inds = fold$validation_set))
    })
    tmledata = list(Q = do.call(rbind, lapply(fold_preds, FUN = function(x) x$Q)),
                    Y = unlist(lapply(fold_preds, FUN = function(x) x$Y)),
                    A = unlist(lapply(fold_preds, FUN = function(x) x$A)),
                    g1W = unlist(lapply(fold_preds, FUN = function(x) x$g1W)))
    
    tmledata1 = list(Q = do.call(rbind, lapply(fold_preds, FUN = function(x) x$Q1)),
                     Y = unlist(lapply(fold_preds, FUN = function(x) x$Y)),
                     A = unlist(lapply(fold_preds, FUN = function(x) x$A)),
                     g1W = unlist(lapply(fold_preds, FUN = function(x) x$g1W1)))
    risk_hal = mean(with(tmledata, -Y*log(Q[,"QAW"])-(1-Y)*log(1-Q[,"QAW"])))
    risk_glm = mean(with(tmledata1, -Y*log(Q[,"QAW"])-(1-Y)*log(1-Q[,"QAW"])))
    risk_halg = mean(with(tmledata, -A*log(g1W)-(1-A)*log(1-g1W)))
    risk_glmg = mean(with(tmledata1, -A*log(g1W)-(1-A)*log(1-g1W)))
    
    supnorm_hal = with(tmledata, max(abs(Q[,"Q1W"] - Q[,"Q0W"] - blip)))
    supnorm_glm = with(tmledata1, max(abs(Q[,"Q1W"] - Q[,"Q0W"] - blip)))
    
    risk = c(risk_hal = risk_hal, risk_glm = risk_glm)
    riskg = c(risk_hal = risk_halg, risk_glm = risk_glmg)
    risk = rbind(risk, riskg)
    
    supnorm = c(sup_hal = supnorm_hal, sup_glm = supnorm_glm)
  }
  return(list(tmledata = tmledata, tmledata1 = tmledata1, risk = risk, supnorm = supnorm))
}


#' @export 
get_results_survtoCDF = function(allresults, n, L,blips) {
  cover_hal_simul = mean(unlist(lapply(allresults, FUN = function(x) {
    all(x$reshal_simul[,6]==TRUE)
  })))
  
  cover_glm_simul = mean(unlist(lapply(allresults, FUN = function(x) {
    all(x$resglm_simul[,6]==TRUE)
  })))
  
  cover_hal_simul
  cover_glm_simul
  
  cover_hal = lapply(1:L, FUN = function(nn) {
    cover = mean(unlist(lapply(allresults, FUN = function(x) {
      x$reshal[[nn]]$info[,6]==TRUE
    })))
    return(cover)
  })
  
  cover_glm = lapply(1:L, FUN = function(nn) {
    cover = mean(unlist(lapply(allresults, FUN = function(x) {
      x$resglm[[nn]]$info[,6]==TRUE
    })))
    return(cover)
  })
  
  # cover_hal
  # cover_glm
  
  
  # MSE
  # allresults[[10]]$reshal_simul[1,"truth"]
  
  est_hal_simul = lapply(1:L, FUN = function(b){
    unlist(lapply(allresults, FUN = function(x) {
      1-x$reshal_simul[b,1]
    }))
  })
  
  est_glm_simul = lapply(1:L, FUN = function(b){
    unlist(lapply(allresults, FUN = function(x) {
      1-x$resglm_simul[b,1]
    }))
  })
  
  est_simul = lapply(1:L, FUN = function(x) {
    xx = cbind(est_hal_simul[[x]], est_glm_simul[[x]])
    colnames(xx) = c("tmle_hal_simul", "tmle_glm_simul")
    return(xx)
  })
  
  truths = vapply(1:L, FUN = function(b) 1-allresults[[1]]$reshal_simul[b,"truth"], FUN.VALUE = 1)
  
  est_hal = lapply(1:L, FUN = function(b){
    xx = lapply(c(1,4) , FUN = function(i) {
      unlist(lapply(allresults, FUN = function(x) {
        1-x$reshal[[b]]$info[1,i]
      }))
    })
    xx = do.call(cbind, xx)
    colnames(xx) = c("tmle_hal", "init_hal")
    return(xx)
  })
  
  est_glm = lapply(1:L, FUN = function(b){
    xx = lapply(c(1,4) , FUN = function(i) {
      unlist(lapply(allresults, FUN = function(x) {
        1-x$resglm[[b]]$info[1,i]
      }))
    })
    xx = do.call(cbind, xx)
    colnames(xx) = c("tmle_glm", "init_glm")
    return(xx)
  })
  
  est = lapply(1:L, FUN = function(x) cbind(est_simul[[x]], est_hal[[x]], est_glm[[x]]))
  # head(est[[8]])
  
  mse = lapply(1:L, FUN = function(b) {
    t(apply(est[[b]], 2, FUN = function(x){
      variance = var(x)
      bias = mean(x-truths[b])
      mse = variance + bias^2                 
      return(c(variance = variance, bias = bias, mse = mse))
    }))
  })
  
  mse
  # head(est[[1]])
  
  plots = lapply(1:L, FUN = function(x) {
    # x=1
    res_temp = est[[x]]
    S_t = truths[x]
    B = nrow(res_temp)
    bw = n^-.2
    inds = 1:6
    ests = c(unlist(lapply(inds, FUN = function(x) res_temp[,x])))
    types = c("TMLE_hal_simul", "TMLE_glm_simul", "TMLE_hal", "Initial_hal", "TMLE_glm", "Initial_glm")
    type = c(unlist(lapply(types, FUN = function(x) rep(x,B))))
    
    inds = inds[order(types)]
    colors = c("red","blue", "green", "yellow", "orange","violet") 
    
    plotdf = data.frame(ests = ests, type = type)
    ggover = ggplot(plotdf,aes(x=ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle("CDF of blip Sampling Dists", 
              subtitle = paste0("n = ", n, ", t = ",blips[x],", bw = ", round(bw,4)))  
    ggover = ggover+geom_vline(xintercept = S_t,color="black")+
      geom_vline(xintercept=mean(res_temp[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(res_temp[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(res_temp[,inds[3]]),color = colors[3])+
      geom_vline(xintercept=mean(res_temp[,inds[4]]),color = colors[4])+
      geom_vline(xintercept=mean(res_temp[,inds[5]]),color = colors[5])+
      geom_vline(xintercept=mean(res_temp[,inds[6]]),color = colors[6])
  })
  
  out = list(cover_hal_simul = cover_hal_simul, cover_glm_simul = cover_glm_simul,
             cover_hal = cover_hal, cover_glm = cover_glm, mse = mse, est = est, plots = plots)
  return(out)
}

#' @export 
get_results_well_survtoCDF = function(allresults, n, L, suffix = "well", blips) {
  cover_well_simul = mean(unlist(lapply(allresults, FUN = function(x) {
    all(x$resglm_simul[,6]==TRUE)
  })))
  
  cover_well = lapply(1:L, FUN = function(nn) {
    cover = mean(unlist(lapply(allresults, FUN = function(x) {
      x$resglm[[nn]]$info[,6]==TRUE
    })))
    return(cover)
  })
  
  
  # cover_hal
  # cover_glm
  
  
  # MSE
  # allresults[[10]]$reshal_simul[1,"truth"]
  
  est_well_simul = lapply(1:L, FUN = function(b){
    unlist(lapply(allresults, FUN = function(x) {
      1-x$resglm_simul[b,1]
    }))
  })

  
  truths = 1-vapply(1:L, FUN = function(b) allresults[[1]]$resglm_simul[b,"truth"], FUN.VALUE = 1)

  est_well = lapply(1:L, FUN = function(b){
    xx = lapply(c(1,4) , FUN = function(i) {
      unlist(lapply(allresults, FUN = function(x) {
        1-x$resglm[[b]]$info[1,i]
      }))
    })
    xx = do.call(cbind, xx)
    colnames(xx) = c("tmle_well", "init_well")
    return(xx)
  })
  
  est = lapply(1:L, FUN = function(x) {
    ggg = cbind(est_well_simul[[x]], est_well[[x]])
    colnames(ggg)[1] = 'tmle_well_simul'
    return(ggg)
    })
  # head(est[[8]])
  
  mse = lapply(1:L, FUN = function(b) {
    t(apply(est[[b]], 2, FUN = function(x){
      variance = var(x)
      bias = mean(x-truths[b])
      mse = variance + bias^2                 
      return(c(variance = variance, bias = bias, mse = mse))
    }))
  })
  
  mse
  # head(est[[1]])
  
  plots = lapply(1:L, FUN = function(x) {
    # x=1
    res_temp = est[[x]]
    S_t = truths[x]
    B = nrow(res_temp)
    bw = n^-.2
    inds = 1:3
    ests = c(unlist(lapply(inds, FUN = function(x) res_temp[,x])))
    types = c("TMLE_glm_simul", "TMLE_glm", "Initial_glm")
    type = c(unlist(lapply(types, FUN = function(x) rep(x,B))))
    
    inds = inds[order(types)]
    colors = c("red","blue", "green") 
    
    plotdf = data.frame(ests = ests, type = type)
    ggover = ggplot(plotdf,aes(x=ests, color = type, fill=type)) + 
      geom_density(alpha=.5)+
      scale_fill_manual(values=colors)+
      scale_color_manual(values=colors)+
      theme(axis.title.x = element_blank())+
      ggtitle("CDF of blip Sampling Dists", 
              subtitle = paste0("n = ", n, ", t = ",blips[x],", bw = ", round(bw,4)))  
    ggover = ggover+geom_vline(xintercept = S_t,color="black")+
      geom_vline(xintercept=mean(res_temp[,inds[1]]),color = colors[1])+
      geom_vline(xintercept=mean(res_temp[,inds[2]]),color = colors[2])+
      geom_vline(xintercept=mean(res_temp[,inds[3]]),color = colors[3])
  })
  
  out = list(cover_well_simul = cover_well_simul, cover_well = cover_well, mse = mse, est = est, plots = plots)
  return(out)
}


#' @export
ind_choose = function(CIs, incre) {
  right = CIs[,3]
  left = CIs[,2]
  if (incre) {
    ind = order(right)[1]
  } else {
    ind = order(left, decreasing = TRUE)[1]
  }
  return(ind)
}

#' @export
ci_form = function(est, SE, z_alpha) {
  df = as.data.frame(cbind(est, est - z_alpha*SE, est + z_alpha*SE))
  colnames(df) = c("ests", "left", "right")
  return(df)
}