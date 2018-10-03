
# stochastic intervention update for the mean
#' @export
eysi_update <- function(tmledata, Q.trunc = 0.001, ...) {

    subset <- with(tmledata, which(Qk > 0 & Qk < 1))
    eps_q <- 0

    # fluctuate Q
    tmledata$Qktrunc <- with(tmledata, truncate(Qk, Q.trunc))
    qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + offset(qlogis(Qktrunc)))
    # qfluc <- logit_fluctuate(tmledata, Y ~ -1 + HA + offset(Qk))
    eps_q <- qfluc$eps
    tmledata$Q_a <- with(tmledata, plogis(qlogis(Q_a) + H * eps_q))
    # tmledata$Qk=qfluc$update

    list(tmledata = tmledata, coefs = c(eps_q))

}


# Stochastic intervention estimate for the mean
#' @export
eysi_estimate <- function(tmledata, ...) {

    # assign probs under stochastic intervention
    n <- nrow(tmledata)
    A_vals <- vals_from_factor(tmledata$A)
    tmledata$pAstar <- sapply(A_vals, function(A_val) gstar(rep(A_val, n), tmledata$pA))

    # compute the parameter under stochastic intervention for fixed g
    psi <- mean(rowSums(tmledata$Q_a * tmledata$pAstar))

    # which A popped up in reality
    ind <- factor_to_indicators(tmledata$A)

    # clever coordinate for each treatment
    tmledata$H <- with(tmledata, pAstar/pA)

    # clever coordinate for the fitted QAW
    tmledata$HA <- rowSums(ind * tmledata$H)

    # create QAW--our predictions over the treatments that occurred
    tmledata$Qk <- rowSums(ind * tmledata$Q_a)

    # influence curves
    tmledata$empirical <- with(tmledata, rowSums(pAstar * Q_a))
    Dstar_psi <- with(tmledata, HA * (Y - Qk) + empirical - psi)

    list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))

}

