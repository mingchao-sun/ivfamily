# jamovi analysis class: first-stage and weak-instrument diagnostics.

ivfirstClass <- R6::R6Class(
    "ivfirstClass",
    inherit = ivfirstBase,
    private = list(

        .run = function() {
            opts <- self$options
            iv_show_text(self$results, opts, "ivfirst")
            if (!iv_ready(opts)) {
                self$results$notes$setContent("<p>Select an outcome, at least one endogenous regressor and at least as many excluded instruments.</p>")
                return()
            }
            st <- iv_setup(opts, self$data)
            pp <- st$pp; res <- st$res
            nm <- res$endog
            nz <- length(pp$instr)

            jt <- self$results$joint
            rowi <- 0
            add <- function(stat, value, df = "", p = NaN) { rowi <<- rowi + 1; jt$addRow(rowKey = rowi, values = list(stat = stat, value = value, df = df, p = p)) }
            cd <- iv_guard(iv_cragg_donald(res))
            add("Cragg-Donald Wald F (homoskedastic)", cd)
            if (opts$vce != "classical") {
                kp <- iv_guard(iv_kleibergen_paap(pp, res))
                add("Kleibergen-Paap rk Wald F", kp$F)
                add("Kleibergen-Paap rk Wald chi-square (H0: underidentified)", kp$chi2, as.character(kp$df), stats::pchisq(kp$chi2, kp$df, lower.tail = FALSE))
            }
            if (length(nm) == 1) {
                ef <- iv_guard(iv_effective_f(pp, res))
                add("Montiel Olea-Pflueger effective F (robust)", ef$F)
                for (k in names(ef$crit))
                    add(sprintf("  critical value, bias threshold %s (5%% level)", k), ef$crit[[k]])
            }
            jt$setNote("rule", paste0("With a single endogenous regressor the rule of thumb is F above 10; the effective F and the Stock-Yogo and Montiel Olea-Pflueger ",
                                      "critical values give sharper guidance. The effective F is compared with the simplified Montiel Olea-Pflueger critical values, which hold for 2SLS and LIML ",
                                      "and are conservative: reject weak instruments for a given bias threshold when F exceeds the value. Kleibergen-Paap statistics allow robust and cluster errors."))

            et <- self$results$each
            for (j in seq_along(nm)) {
                fs <- iv_guard(iv_first_stage(pp, res, nm[j]))
                sw <- if (length(nm) > 1) iv_guard(iv_sw_f(pp, res, j)) else NULL
                et$addRow(rowKey = j, values = list(term = nm[j], F = fs$F, df1 = fs$df1, df2 = if (is.finite(fs$df2)) fs$df2 else NaN, p = fs$p,
                                                   r2 = fs$partial_r2, swF = if (is.null(sw)) NaN else sw$F, swDf = if (is.null(sw)) NA_integer_ else sw$df1))
            }
            if (length(nm) > 1)
                et$setNote("sw", "Conditional F: first-stage F for each regressor after removing the other endogenous regressors (Sanderson and Windmeijer 2016); it tests whether the instruments identify that regressor separately.")

            if (isTRUE(opts$showReduced)) {
                rf <- iv_guard(iv_reduced_form(pp, res))
                rt <- self$results$reduced
                for (i in seq_along(rf$term))
                    rt$addRow(rowKey = i, values = list(term = rf$term[i], estimate = rf$estimate[i], se = rf$se[i], p = rf$p[i]))
                rt$setNote("f", sprintf("Joint F for the excluded instruments in the reduced form: %s (p = %s).", iv_fmt(rf$F, 3), iv_fmt(rf$pF, 4)))
            }
            notes <- character(0)
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            if (nz == length(nm)) notes <- c(notes, "The model is exactly identified: there are as many instruments as endogenous regressors.")
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        }
    )
)
