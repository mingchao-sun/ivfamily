# Binary treatment and binary instrument: LATE, compliance shares and complier characteristics.
# Under independence, the exclusion restriction and monotonicity (Imbens and Angrist 1994) the Wald ratio is the average
# effect for compliers. The shares and group means use the raw sample (no controls).

iv_late_shares <- function(d, z) {
    n1 <- sum(z == 1); n0 <- sum(z == 0)
    p1 <- mean(d[z == 1]); p0 <- mean(d[z == 0])
    pa <- p0; pn <- 1 - p1; pc <- p1 - p0
    se_a <- sqrt(p0 * (1 - p0) / n0); se_n <- sqrt(p1 * (1 - p1) / n1)
    se_c <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    data.frame(group = c("Always-takers", "Never-takers", "Compliers"), share = c(pa, pn, pc), se = c(se_a, se_n, se_c),
               stringsAsFactors = FALSE)
}

# Means of covariates in the three groups: always-takers (treated although the instrument is off), never-takers
# (untreated although the instrument is on) and compliers from the overall mean.
iv_late_profile <- function(X, d, z) {
    sh <- iv_late_shares(d, z)
    pa <- sh$share[1]; pn <- sh$share[2]; pc <- sh$share[3]
    out <- lapply(colnames(X), function(v) {
        x <- X[, v]
        m_all <- mean(x)
        m_at <- mean(x[d == 1 & z == 0]); m_nt <- mean(x[d == 0 & z == 1])
        m_c <- (m_all - pa * m_at - pn * m_nt) / pc
        data.frame(variable = v, all = m_all, always = m_at, never = m_nt, compliers = m_c, ratio = m_c / m_all, stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}

iv_check_binary <- function(x, what) {
    u <- sort(unique(x))
    if (!all(u %in% c(0, 1)) || length(u) < 2)
        stop(sprintf("%s must be a binary variable coded 0/1 with both values present.", what), call. = FALSE)
}
