
############## TMLE targeting EY1
#' @export
ey1_update <- function(tmledata, Q.trunc = 0.001, ...) {
    subset <- with(tmledata, which(0 < Qk & Qk < 1 & A == 1))
    eps_q <- 0
    
    # fluctuate Q
    tmledata$Qktrunc <- with(tmledata, truncate(Qk, Q.trunc))
    qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + offset(qlogis(Qktrunc)))
    eps_q <- qfluc$eps
    tmledata$Qk <- with(tmledata, plogis(qlogis(Qktrunc) + H1 * eps_q))
    # tmledata$Qk=qfluc$update
    
    
    list(tmledata = tmledata, coefs = c(eps_q))
    
}

#' @export
ey1_estimate <- function(tmledata, ...) {
    
    psi <- mean(tmledata$Qk)
    
    tmledata$H1 <- with(tmledata, (1/gk))
    tmledata$HA <- with(tmledata, (A * H1))
    
    # influence curves
    Dstar_psi <- with(tmledata, HA * (Y - Qk) + Qk - psi)
    
    list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))
    
}

############## TMLE targeting EY1 (more canonical than above, but should be equivalent)

#' @export
ey1_update2 <- function(tmledata, Q.trunc = 0.001, ...) {
    eps_q <- 0
    # fluctuate Q
    tmledata$QAktrunc <- with(tmledata, truncate(QAk, Q.trunc))
    tmledata$Q1ktrunc <- with(tmledata, truncate(Q1k, Q.trunc))
    # tmledata$Q0ktrunc=with(tmledata,truncate(Q0k, Q.trunc))
    qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + offset(qlogis(QAktrunc)), subset)
    eps_q <- qfluc$eps
    tmledata$Q1k <- with(tmledata, plogis(qlogis(Q1k) + H1 * eps_q))
    tmledata$QAk <- with(tmledata, plogis(qlogis(QAk) + HA * eps_q))
    # tmledata$Qk=qfluc$update
    
    
    list(tmledata = tmledata, coefs = c(eps_q))
    
}

#' @export
ey1_estimate2 <- function(tmledata, ...) {
    
    psi <- mean(tmledata$Q1k)
    
    tmledata$H1 <- with(tmledata, (1/gk))
    tmledata$HA <- with(tmledata, (A * H1))
    
    # influence curves
    Dstar_psi <- with(tmledata, HA * (Y - QAk) + Q1k - psi)
    
    list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))
    
} 
