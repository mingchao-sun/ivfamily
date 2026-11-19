# Sample data of the marginal treatment effects analysis: schooling choice with essential heterogeneity (generalised Roy model).
# Run from the module root:  Rscript dev/make_data_mte.R
set.seed(20261128)
n <- 6000
score <- round(rnorm(n), 3); female <- rbinom(n, 1, 0.5)
distance <- round(rexp(n, 1 / 6), 2); tuition <- round(runif(n, 1, 4), 2)             # instruments: distance to a college and tuition
v <- rnorm(n)                                                                            # unobserved resistance to college: U_D = Phi(V)
p <- pnorm(1.5 + 0.5 * score - 0.10 * distance - 0.35 * tuition + 0.1 * female)
college <- as.integer(p >= pnorm(v))
u0 <- 0.4 * v + rnorm(n, 0, 0.5)
gain <- -0.30 * v + rnorm(n, 0, 0.25)                                                    # returns fall with resistance: those eager to go gain most
y0 <- 2 + 0.3 * score + 0.1 * female + u0
y1 <- y0 + 0.45 + 0.10 * score + gain
logwage <- round(ifelse(college == 1, y1, y0), 3)
write.csv(data.frame(id = seq_len(n), logwage, college = ifelse(college == 1, "yes", "no"), score, female, distance, tuition), "data/sim_mte.csv", row.names = FALSE)
cat("share in college", round(mean(college), 3), " range of the propensity score", round(range(p), 3), "\n")
cat("true ATE (at the mean score, 0.45 + 0.10 mean(score))", round(0.45 + 0.1 * mean(score), 3), " true MTE(u) = 0.45 + 0.10 score - 0.30 qnorm(u)\n")
