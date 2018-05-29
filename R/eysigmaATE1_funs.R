#' @export
eysigmaATE1_update <- function(tmledata, Q.trunc = 1e-04,Dstar) {
    # fix points where Q is already 0 or 1 - perfect prediction
    subset <- with(tmledata, which(0 < Qk & Qk < 1))
    eps_q <- 0
    # fluctuate Q
    HA = cbind(tmledata$HAW,tmledata$CYAW)
    H1 = cbind(tmledata$H1W,tmledata$CY1W)
    H0 = cbind(tmledata$H0W,tmledata$CY0W)
    deriv = colMeans(as.data.frame(Dstar))
    norm = sqrt(sum(deriv^2))
    H = HA%*%deriv/norm
    tmledata$Qktrunc <- with(tmledata, truncate(Qk, Q.trunc))
    qfluc <- logit_fluctuate(tmledata, Y ~ -1 + H + offset(qlogis(Qktrunc)))
    eps_q <- qfluc$eps
    tmledata$Qk <- with(tmledata, plogis(qlogis(Qktrunc) + H*eps_q))
    tmledata$Q1k <- as.vector(with(tmledata, plogis(qlogis(Q1k) + eps_q*H1%*%deriv/norm)))
    tmledata$Q0k <- as.vector(with(tmledata, plogis(qlogis(Q0k) + eps_q*H0%*%deriv/norm)))


    list(tmledata = tmledata, coefs = c(eps_q))

}


#' @export
eysigmaATE1_estimate <- function(tmledata, ...) {
    n=nrow(tmledata)
    psi <- mean(tmledata$Q1k - tmledata$Q0k)
    sigma <- var(tmledata$Q1k - tmledata$Q0k)*(n-1)/n
    tmledata$CYAW <- with(tmledata, 2 * (Q1k - Q0k - psi) * (A/gk - (1 - A)/(1 - gk)))
    tmledata$CY1W <- with(tmledata, 2 * (Q1k - Q0k - psi)/gk)
    tmledata$CY0W <- with(tmledata, -2 * (Q1k - Q0k - psi)/(1 - gk))
    tmledata$HAW <- with(tmledata, (A/gk - (1 - A)/(1 - gk)))
    tmledata$H1W <- with(tmledata, 1/gk)
    tmledata$H0W <- with(tmledata, -1/(1 - gk))

    # influence curves
    Dstar_psi <- with(tmledata, HAW * (Y - Qk) + Q1k - Q0k - psi)
    Dstar_sigma <- with(tmledata, CYAW * (Y - Qk) + (Q1k - Q0k - psi)^2 - sigma)
    list(tmledata = tmledata, ests = c(psi = psi, sigma = sigma), Dstar = list(Dstar_psi = Dstar_psi,
        Dstar_sigma = Dstar_sigma))
}
