# Shift-share IV: the engine against the R package ShiftShareSE (ADH China-shock data and the sample data), then the jamovi
# analysis headless on the sample data and edge cases. Run from the module root:  Rscript tests/shiftshare.R
suppressMessages({library(jmvcore); library(R6); library(ggplot2)})
for (f in c("R/00jmv.R", "R/common.R", "R/engine_ss.R", "R/guide.R", "R/shiftshare.h.R", "R/shiftshare.b.R")) source(f)
nfail <- 0
chk <- function(name, got, want, tol) {
    d <- max(abs(got - want), na.rm = TRUE); ok <- is.finite(d) && d < tol
    cat(sprintf("%-60s %s  max |diff| = %.1e\n", name, if (ok) "ok  " else "FAIL", d)); if (!ok) nfail <<- nfail + 1
}
if (requireNamespace("ShiftShareSE", quietly = TRUE)) {
    suppressMessages(library(ShiftShareSE))
    r <- ADH$reg; W <- ADH$W
    ctl <- c("t2", "l_shind_manuf_cbp", "l_sh_popedu_c", "l_sh_popfborn", "l_sh_empl_f", "l_sh_routine33", "l_task_outsource")
    Z <- cbind(1, as.matrix(r[, ctl]))
    ref <- ivreg_ss(d_sh_empl_mfg ~ t2 + l_shind_manuf_cbp + l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f + l_sh_routine33 + l_task_outsource | shock,
                    X = IV, data = r, W = W, weights = weights, region_cvar = statefip, method = "all")
    me <- iv_ss_ivreg(r$d_sh_empl_mfg, r$shock, r$IV, W, Z, w = r$weights, region_cvar = r$statefip)
    chk("ADH weighted, controls: coefficient", me$beta, ref$beta, 1e-10)
    chk("ADH weighted: five standard errors", unname(me$se), unname(ref$se), 1e-10)
    chk("ADH weighted: p-values", unname(me$p), unname(ref$p), 1e-10)
    chk("ADH weighted: confidence limits (incl. AKM0)", c(me$lo, me$hi), c(ref$ci.l, ref$ci.r), 1e-9)
    sc <- substr(sprintf("%04d", ADH$sic), 1, 2)
    Z2 <- cbind(1, as.matrix(r[, c("t2", "l_shind_manuf_cbp")]))
    ref2 <- ivreg_ss(d_sh_empl_mfg ~ t2 + l_shind_manuf_cbp | shock, X = IV, data = r, W = W, weights = weights, region_cvar = statefip, sector_cvar = sc, method = "all")
    me2 <- iv_ss_ivreg(r$d_sh_empl_mfg, r$shock, r$IV, W, Z2, w = r$weights, region_cvar = r$statefip, sector_cvar = sc)
    chk("ADH with sector clusters: standard errors", unname(me2$se), unname(ref2$se), 1e-10)
    ref3 <- ivreg_ss(d_sh_empl_mfg ~ 1 | shock, X = IV, data = r, W = W, region_cvar = statefip, method = "all")
    me3 <- iv_ss_ivreg(r$d_sh_empl_mfg, r$shock, r$IV, W, matrix(1, nrow(r), 1), region_cvar = r$statefip)
    chk("ADH unweighted, no controls: standard errors", unname(me3$se), unname(ref3$se), 1e-10)
    rr <- reg_ss(shock ~ t2 + l_shind_manuf_cbp, X = IV, data = r, W = W, weights = weights, region_cvar = statefip, method = "all")
    mr <- iv_ss_reg(r$shock, r$IV, W, Z2, w = r$weights, region_cvar = r$statefip)
    chk("ADH first stage (reg_ss): coefficient and standard errors", c(mr$beta, unname(mr$se)), c(rr$beta, unname(rr$se)), 1e-10)
    # sample data through the long-form preparation
    d <- read.csv("data/sim_shiftshare.csv")
    pp <- iv_ss_prepare(d, "region", "sector", "share", "shock", "y", "d", "x1", "pop", "state", "sector_group")
    ref4 <- ivreg_ss(y ~ x1 | d, X = z, data = data.frame(y = pp$y, d = pp$d, z = pp$z, x1 = pp$Z[, "x1"], w = pp$w, st = pp$region_cvar),
                     W = pp$W, weights = w, region_cvar = st, sector_cvar = pp$sector_cvar, method = "all")
    me4 <- iv_ss_ivreg(pp$y, pp$d, pp$z, pp$W, pp$Z, pp$w, region_cvar = pp$region_cvar, sector_cvar = pp$sector_cvar)
    chk("sample data (long form): coefficient and standard errors", c(me4$beta, unname(me4$se)), c(ref4$beta, unname(ref4$se)), 1e-10)
} else cat("ShiftShareSE not installed: engine checks skipped\n")

d <- read.csv("data/sim_shiftshare.csv")
run <- function(label, data, expect_error = FALSE, ...) {
    a <- list(...)
    r <- tryCatch({ res <- do.call(shiftshare, c(list(data = data), a)); invisible(capture.output(print(res))); "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)))
    good <- if (expect_error) grepl("^ERROR", r) else r == "ok"
    cat(sprintf("%-46s %s\n", label, if (good) substr(r, 1, 90) else paste("UNEXPECTED:", r))); if (!good) nfail <<- nfail + 1
}
base <- list(outcome = "y", endog = "d", region = "region", sector = "sector", share = "share", shock = "shock")
opt <- list(controls = NULL, weights = NULL, regionCluster = NULL, sectorCluster = NULL)
call <- function(label, data, expect_error = FALSE, ...) { a <- modifyList(c(base, opt), list(...), keep.null = TRUE); do.call(run, c(list(label, data, expect_error), a)) }
call("minimal", d)
call("all options", d, controls = "x1", weights = "pop", regionCluster = "state", sectorCluster = "sector_group", showGuide = TRUE, showHelp = TRUE)
call("sum of shares control, h0 = 0.5, 90%", d, controls = "x1", sumShares = TRUE, h0 = 0.5, level = 90)
call("no stages table", d, showStages = FALSE)
d2 <- d; d2$y[d2$region == 3] <- NA
call("region with a missing outcome", d2)
d3 <- d; d3$y[1] <- d3$y[1] + 1
call("region-level variable varies in a region", d3, TRUE)
d4 <- d; d4$shock[1] <- d4$shock[1] + 1
call("shock varies in a sector", d4, TRUE)
call("duplicate region-sector rows", rbind(d, d[1:3, ]), TRUE)
d5 <- d; d5$share[1] <- -0.1
call("negative share", d5, TRUE)
call("missing share variable", d, share = NULL)
d6 <- d[d$region <= 5, ]
call("too few regions", d6, TRUE)
cat(if (nfail == 0) "\nAll shift-share checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
