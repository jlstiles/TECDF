
eysigmaATE_update <- function(tmledata, Q.trunc = 1e-04) {
    # fix points where Q is already 0 or 1 - perfect prediction
    subset <- with(tmledata, which(0 < Qk & Qk < 1))
    eps_q <- 0
    # fluctuate Q
    tmledata$Qktrunc <- with(tmledata, truncate(Qk, Q.trunc))
    qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + CYA + offset(qlogis(Qktrunc)))
    eps_q <- qfluc$eps
    tmledata$Qk <- with(tmledata, plogis(qlogis(Qktrunc) + HA * eps_q[1] + CYA * eps_q[2]))
    tmledata$Q1k <- with(tmledata, plogis(qlogis(Q1k) + H1 * eps_q[1] + CY1 * eps_q[2]))
    tmledata$Q0k <- with(tmledata, plogis(qlogis(Q0k) + H0 * eps_q[1] + CY0 * eps_q[2]))
    
    
    list(tmledata = tmledata, coefs = c(eps_q))
    
}

eysigmaATE_estimate <- function(tmledata, ...) {
    nn=length(tmledata$Q0k)
    psi <- mean(tmledata$Q1k - tmledata$Q0k)
    sigma <- var(tmledata$Q1k - tmledata$Q0k)*(nn-1)/nn
    tmledata$CYA <- with(tmledata, 2 * (Q1k - Q0k - psi) * (A/gk - (1 - A)/(1 - gk)))
    tmledata$CY1 <- with(tmledata, 2 * (Q1k - Q0k - psi)/gk)
    tmledata$CY0 <- with(tmledata, -2 * (Q1k - Q0k - psi)/(1 - gk))
    tmledata$HA <- with(tmledata, (A/gk - (1 - A)/(1 - gk)))
    tmledata$H1 <- with(tmledata, 1/gk)
    tmledata$H0 <- with(tmledata, -1/(1 - gk))
    
    # influence curves
    Dstar_psi <- with(tmledata, HA * (Y - Qk) + Q1k - Q0k - psi)
    Dstar_sigma <- with(tmledata, CYA * (Y - Qk) + (Q1k - Q0k - psi)^2 - sigma)
    
    list(tmledata = tmledata, ests = c(psi = psi, sigma = sigma), Dstar = list(Dstar_psi = Dstar_psi, 
        Dstar_sigma = Dstar_sigma))
} 
