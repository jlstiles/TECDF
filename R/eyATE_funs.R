   
  
############## TMLE targeting EYATE
#' @export
eyATE_update <- function(tmledata, Q.trunc = 1e-04, ...) {
    subset <- with(tmledata, which(0 < Qk & Qk < 1))
    eps_q <- 0
    
    #fluctuate Q
    tmledata$Qktrunc <- with(tmledata, truncate(Qk, Q.trunc))
    qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + offset(qlogis(Qktrunc)))
    eps_q <- qfluc$eps
    tmledata$Qk <- with(tmledata, plogis(qlogis(Qktrunc) + HA * eps_q))
    tmledata$Q1k <- with(tmledata, plogis(qlogis(Q1k) + H1 * eps_q))
    tmledata$Q0k <- with(tmledata, plogis(qlogis(Q0k) + H0 * eps_q))
    # tmledata$Qk=qfluc$update
    
    
    list(tmledata = tmledata, coefs = c(eps_q))
    
}

#' @export
eyATE_estimate <- function(tmledata, ...) {
    
    psi <- mean(tmledata$Q1k - tmledata$Q0k)
    
    tmledata$H1 <- with(tmledata, (1/gk))
    tmledata$H0 <- with(tmledata, -1/(1 - gk))
    tmledata$HA <- with(tmledata, (A * H1 + (1 - A) * H0))
    
    # influence curves
    Dstar_psi <- with(tmledata, HA * (Y - Qk) + Q1k - Q0k - psi)
    
    list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))
    
} 
