# jamovi analysis class: binary treatment, binary instrument; LATE, compliance shares and complier characteristics.

ivlateClass <- R6::R6Class(
    "ivlateClass",
    inherit = ivlateBase,
    private = list(

        .run = function() {
            opts <- self$options
            iv_show_text(self$results, opts, "ivlate")
            if (is.null(opts$outcome) || is.null(opts$treatment) || is.null(opts$instrument)) {
                self$results$notes$setContent("<p>Select an outcome, a binary treatment and a binary instrument.</p>")
                return()
            }
            if (opts$vce == "cluster" && is.null(opts$cluster))
                jmvcore::reject("Select a cluster variable for cluster-robust standard errors.")
            data <- self$data
            notes <- character(0)
            # two-level factors become 0/1 with the second level (alphabetical or as ordered in jamovi) as 1
            for (v in c(opts$treatment, opts$instrument)) {
                x <- data[[v]]
                if (is.factor(x) || is.character(x) || is.logical(x)) {
                    lv <- levels(factor(x))
                    if (length(lv) != 2) jmvcore::reject(sprintf("%s must have exactly two values.", v))
                    data[[v]] <- as.numeric(factor(x, levels = lv)) - 1
                    notes <- c(notes, sprintf("%s: '%s' is coded 0 and '%s' is coded 1.", v, lv[1], lv[2]))
                }
            }
            pp <- iv_guard(iv_prepare(data, opts$outcome, opts$treatment, opts$instrument, opts$controls, NULL, opts$cluster, NULL, extra = opts$profile))
            iv_guard(iv_check_binary(pp$X1[, 1], opts$treatment))
            iv_guard(iv_check_binary(pp$Z1[, 1], opts$instrument))
            res <- iv_guard(iv_fit(pp, "2sls", vce = opts$vce, small = opts$small))
            rf <- iv_guard(iv_reduced_form(pp, res))
            fs <- iv_guard(iv_first_stage(pp, res, opts$treatment))

            crit <- if (is.finite(res$df)) stats::qt(1 - (1 - opts$level / 100) / 2, res$df) else stats::qnorm(1 - (1 - opts$level / 100) / 2)
            pv <- function(est, se) if (is.finite(res$df)) 2 * stats::pt(-abs(est / se), res$df) else 2 * stats::pnorm(-abs(est / se))
            ct <- iv_coef_table(res, opts$level); ct <- ct[ct$term == opts$treatment, ]
            et <- self$results$effects
            et$addRow(rowKey = 1, values = list(item = "Reduced form: effect of the instrument on the outcome", estimate = rf$estimate[1], se = rf$se[1],
                                                p = rf$p[1], lo = rf$estimate[1] - crit * rf$se[1], hi = rf$estimate[1] + crit * rf$se[1]))
            et$addRow(rowKey = 2, values = list(item = "First stage: effect of the instrument on the treatment", estimate = fs$coef[1], se = fs$se[1],
                                                p = pv(fs$coef[1], fs$se[1]), lo = fs$coef[1] - crit * fs$se[1], hi = fs$coef[1] + crit * fs$se[1]))
            et$addRow(rowKey = 3, values = list(item = "LATE: effect of the treatment for compliers (Wald / 2SLS)", estimate = ct$estimate, se = ct$se,
                                                p = ct$p, lo = ct$lo, hi = ct$hi))
            et$setNote("w", "LATE = reduced form / first stage, the 2SLS estimate with a binary treatment and a binary instrument. It describes compliers, the people whose treatment is moved by the instrument.")

            d <- pp$X1[, 1]; z <- pp$Z1[, 1]
            sh <- iv_late_shares(d, z)
            st <- self$results$shares
            for (i in seq_len(nrow(sh))) st$addRow(rowKey = i, values = list(group = sh$group[i], share = sh$share[i], se = sh$se[i]))
            st$setNote("m", "Shares use the raw sample (no controls) and assume monotonicity: the instrument never pushes anyone out of treatment.")
            if (sh$share[3] <= 0) notes <- c(notes, "The estimated share of compliers is not positive: the instrument has no first stage or monotonicity fails.")
            else if (min(sh$share) < 0) notes <- c(notes, "A group share is negative, which contradicts monotonicity.")
            if (sh$share[3] < 0.05 && sh$share[3] > 0) notes <- c(notes, "Fewer than 5% compliers: the instrument is weak and the LATE is imprecise.")

            if (length(opts$profile) > 0 && !is.null(pp$extra)) {
                pr <- iv_late_profile(pp$extra, d, z)
                tb <- self$results$profile
                for (i in seq_len(nrow(pr))) tb$addRow(rowKey = i, values = as.list(pr[i, c("variable", "all", "always", "never", "compliers", "ratio")]))
                tb$setNote("c", "Complier means are derived from the overall mean and the always- and never-taker means; they have no standard errors here.")
            }
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        }
    )
)
