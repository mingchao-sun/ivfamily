# jamovi analysis class: shift-share (Bartik) instrumental variables with exposure-robust standard errors.

shiftshareClass <- R6::R6Class(
    "shiftshareClass",
    inherit = shiftshareBase,
    private = list(

        .run = function() {
            o <- self$options
            iv_show_text(self$results, o, "shiftshare")
            req <- c("outcome", "endog", "region", "sector", "share", "shock")
            if (!all(vapply(req, function(k) length(o[[k]]) > 0, TRUE))) {
                self$results$notes$setContent(paste0("<p>Select the outcome, the endogenous regressor, the region, the sector, the exposure share and the shock. ",
                                                     "The data are in long form: one row per region and sector, with region-level variables repeated on every row of a region ",
                                                     "and the shock repeated on every row of a sector.</p>"))
                for (nm in c("info", "iv", "stages")) self$results[[nm]]$setVisible(FALSE)
                return()
            }
            alpha <- 1 - o$level / 100
            pp <- iv_guard(iv_ss_prepare(self$data, o$region, o$sector, o$share, o$shock, o$outcome, o$endog, o$controls, o$weights,
                                         o$regionCluster, o$sectorCluster, isTRUE(o$sumShares)))
            res <- iv_guard(iv_ss_ivreg(pp$y, pp$d, pp$z, pp$W, pp$Z, pp$w, beta0 = o$h0, alpha = alpha, region_cvar = pp$region_cvar, sector_cvar = pp$sector_cvar))
            notes <- character(0)

            # ---- sample and instrument
            it <- self$results$info
            it$deleteRows()
            ess <- iv_ss_effective_sectors(pp$W, pp$w)
            ssum <- pp$share_sum
            rows <- list(c("Regions", format(pp$L, big.mark = ",")), c("Sectors", format(pp$K, big.mark = ",")),
                         c("Effective number of sectors", iv_fmt(ess, 1)),
                         c("Sum of shares per region: mean (min to max)", sprintf("%s (%s to %s)", iv_fmt(mean(ssum), 3), iv_fmt(min(ssum), 3), iv_fmt(max(ssum), 3))))
            if (!is.null(pp$region_cvar)) rows[[length(rows) + 1]] <- c("Region clusters", format(length(unique(pp$region_cvar)), big.mark = ","))
            if (!is.null(pp$sector_cvar)) rows[[length(rows) + 1]] <- c("Sector clusters", format(length(unique(pp$sector_cvar)), big.mark = ","))
            for (i in seq_along(rows)) it$addRow(rowKey = i, values = list(item = rows[[i]][1], value = rows[[i]][2]))
            if (pp$n_dropped_rows > 0) notes <- c(notes, sprintf("%d row(s) with missing values were dropped%s.", pp$n_dropped_rows,
                                                                  if (pp$n_dropped_regions > 0) sprintf(", including %d region(s) with a missing region-level value", pp$n_dropped_regions) else ""))
            if (stats::sd(ssum) > 1e-8 && !pp$sum_shares_added)
                notes <- c(notes, "The shares do not add up to the same number in every region. The instrument is then also driven by the sum of the shares; tick 'Control for the sum of the shares' to hold it fixed.")
            if (ess < 20) notes <- c(notes, sprintf("The effective number of sectors is only %s: the exposure-robust standard errors rely on many sectors and can be unreliable with so few.", iv_fmt(ess, 1)))

            # ---- IV table
            tb <- self$results$iv
            tb$deleteRows()
            keep <- !is.na(res$se)
            ks <- which(keep)
            fin <- function(x) if (is.finite(x)) x else NaN
            for (i in seq_along(ks)) {
                k <- ks[i]
                tb$addRow(rowKey = i, values = list(method = names(res$se)[k], estimate = res$beta, se = fin(res$se[[k]]), p = res$p[[k]], lo = fin(res$lo[[k]]), hi = fin(res$hi[[k]])))
            }
            tb$setNote("a", paste0("Two-stage least squares with the shift-share instrument (the sum over sectors of share times shock) and the controls. ",
                                   "AKM standard errors allow regions with similar exposure shares to have correlated errors (Adao, Kolesar and Morales, 2019); AKM0 imposes the null value of the coefficient ",
                                   "and has better coverage when the instrument is weak. Only AKM and AKM0 should be used for inference; the others treat regions as independent."))
            if (!is.finite(res$se[["AKM0"]]) || !all(is.finite(c(res$lo[5], res$hi[5]))))
                notes <- c(notes, "The AKM0 confidence set is unbounded or empty, which signals a weak instrument at this level.")

            # ---- first stage and reduced form
            st <- self$results$stages
            st$deleteRows()
            fs <- iv_guard(iv_ss_reg(pp$d, pp$z, pp$W, pp$Z, pp$w, 0, alpha, pp$region_cvar, pp$sector_cvar))
            rf <- iv_guard(iv_ss_reg(pp$y, pp$z, pp$W, pp$Z, pp$w, 0, alpha, pp$region_cvar, pp$sector_cvar))
            add <- function(i, lab, r) st$addRow(rowKey = i, values = list(item = lab, coef = r$beta, seEhw = r$se[["Robust (EHW)"]], seAkm = r$se[["AKM"]],
                                                                          tAkm = r$beta / r$se[["AKM"]], f = (r$beta / r$se[["AKM"]])^2, r2 = r$partial_r2))
            add(1, "First stage: instrument on the regressor", fs)
            add(2, "Reduced form: instrument on the outcome", rf)
            st$setNote("f", "F (AKM) is the squared t statistic with the AKM standard error: an exposure-robust F for the strength of the instrument. The IV coefficient is the ratio of the two coefficients.")

            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        }
    )
)
