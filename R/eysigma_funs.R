############## TMLE targeting Ey1 and Sigma^2(Ey1) simultaneously
#' @export
eysigma_update <- function(tmledata, Q.trunc = 0.001, g.trunc = 0.001, ...) {
    eps_q <- 0
    # logit_fluctuate Q
    tmledata$Qktrunc <- with(tmledata, truncate(Qk, Q.trunc))
    qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + CyA + offset(qlogis(Qktrunc)), subset)
    eps_q <- qfluc$eps
    tmledata$Qk <- with(tmledata, plogis(qlogis(Qktrunc) + H1 * eps_q[1] + Cy1 * eps_q[2]))
    
    # logit_fluctuate g
    tmledata$gktrunc <- with(tmledata, truncate(gk, g.trunc))
    gfluc <- logit_fluctuate(tmledata, A ~ -1 + Ca + offset(qlogis(gktrunc)))
    eps_g <- gfluc$eps
    tmledata$gk <- with(tmledata, plogis(qlogis(gktrunc) + Ca * eps_g))
    
    list(tmledata = tmledata, coefs = c(eps_q, eps_g))
    
}

#' @export
eysigma_estimate <- function(tmledata, ...) {
    
    psi <- mean(tmledata$Qk)
    
    tmledata$H1 <- with(tmledata, (1/gk))
    tmledata$HA <- with(tmledata, (A * H1))
    
    tmledata$Cy1 <- with(tmledata, (1/gk) * ((1 - 2 * Qk)/gk + 2 * (Qk - psi)))
    tmledata$CyA <- with(tmledata, A * Cy1)
    
    tmledata$Ca <- with(tmledata, Qk * (1 - Qk)/(gk^2))
    
    # influence curves
    Dstar_psi <- with(tmledata, HA * (Y - Qk) + Qk - psi)
    
    varterm <- with(tmledata, Qk * (1 - Qk)/gk + (Qk - psi)^2)
    sigma2 <- mean(varterm)
    D_Qw <- varterm - sigma2
    D_Qbar <- with(tmledata, CyA * (Y - Qk))
    D_gbar <- with(tmledata, -Ca * (A - gk))
    Dstar_sigma <- D_Qw + D_Qbar + D_gbar
    
    
    
    list(tmledata = tmledata, ests = c(psi = psi, sigma2 = sigma2), Dstar = list(Dstar_psi = Dstar_psi, 
        Dstar_sigma = Dstar_sigma))
} 
