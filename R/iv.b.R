# jamovi analysis class: IV regression (2SLS, LIML, Fuller, k-class, two-step GMM).

ivClass <- R6::R6Class(
    "ivClass",
    inherit = ivBase,
    private = list(

        .run = function() {
            opts <- self$options
            iv_show_text(self$results, opts, "iv")
            if (is.null(opts$outcome) || length(opts$endog) == 0 || length(opts$instruments) == 0) {
                self$results$notes$setContent(
                    "<p>Select an outcome, at least one endogenous regressor and at least as many excluded instruments.</p>")
                return()
            }
            pp <- iv_guard(iv_prepare(self$data, opts$outcome, opts$endog, opts$instruments, opts$controls,
                                      opts$absorb, opts$cluster, opts$weights))
            if (opts$vce %in% c("cluster", "cluster2") && length(opts$cluster) < (if (opts$vce == "cluster2") 2 else 1))
                jmvcore::reject(if (opts$vce == "cluster2") "Two-way clustering needs two cluster variables."
                                else "Select a cluster variable for cluster-robust standard errors.")
            if (opts$method == "gmm" && opts$vce == "classical")
                jmvcore::reject("Two-step GMM needs robust or cluster-robust standard errors; change Standard errors.")

            res <- iv_guard(iv_fit(pp, opts$method, alpha = opts$alpha, kvalue = opts$kvalue, vce = opts$vce, small = opts$small))
            ct <- iv_coef_table(res, opts$level)
            if (!isTRUE(opts$showControls)) ct <- ct[ct$type == "Endogenous", , drop = FALSE]
            if (nrow(ct) > 0) {
                tb <- self$results$coef
                for (i in seq_len(nrow(ct)))
                    tb$addRow(rowKey = i, values = as.list(ct[i, c("term", "type", "estimate", "se", "stat", "p", "lo", "hi")]))
                tb$setNote("vce", paste0(iv_vce_label(res), if (is.finite(res$df)) sprintf(" Test statistics are t with %s degrees of freedom.", iv_fmt(res$df, 0))
                                                              else " Test statistics are normal."))
            }

            if (isTRUE(opts$showOls)) {
                ols <- iv_guard(iv_fit(pp, "ols", vce = opts$vce, small = opts$small))
                tb <- self$results$ols
                for (j in seq_along(res$endog)) {
                    e <- res$endog[j]
                    tb$addRow(rowKey = j, values = list(term = e, ols = ols$beta[[e]], olsSe = ols$se[[e]], iv = res$beta[[e]],
                                                         ivSe = res$se[[e]], diff = res$beta[[e]] - ols$beta[[e]]))
                }
            }

            fF <- self$results$firstF
            fc <- self$results$first
            nr <- 0
            for (j in seq_along(res$endog)) {
                fs <- iv_guard(iv_first_stage(pp, res, res$endog[j]))
                fF$addRow(rowKey = j, values = list(term = res$endog[j], F = fs$F, df1 = fs$df1, df2 = if (is.finite(fs$df2)) fs$df2 else NaN,
                                                    p = fs$p, r2 = fs$partial_r2))
                for (q in seq_along(fs$coef)) {
                    nr <- nr + 1
                    fc$addRow(rowKey = nr, values = list(term = res$endog[j], instrument = colnames(res$mm$Z1)[q], coef = fs$coef[[q]], se = fs$se[[q]]))
                }
            }
            fF$setNote("f", paste0("F tests that the excluded instruments are jointly zero in the first stage, using the chosen standard errors ",
                                   "(with robust or cluster-robust errors this is the Wald F, not the Kleibergen-Paap statistic). ",
                                   "A common rule of thumb is F above 10 for a single endogenous regressor."))

            sm <- self$results$summary
            items <- list(
                c("Estimator", iv_method_label(opts$method, res$k)),
                c("Observations", format(res$N, big.mark = ",")),
                c("Dropped for missing values", format(pp$n_drop, big.mark = ",")),
                c("Endogenous regressors / excluded instruments", sprintf("%d / %d", length(pp$endog), length(pp$instr))))
            if (!is.null(res$G)) items[[length(items) + 1]] <- c("Clusters", format(res$G, big.mark = ","))
            if (length(pp$fes)) items[[length(items) + 1]] <- c("Absorbed fixed effects", paste(names(pp$fes), collapse = ", "))
            if (!is.null(res$J))
                items[[length(items) + 1]] <- c(sprintf("%s (df = %d)", res$J$label, res$J$df),
                                                 sprintf("%s, p = %s", iv_fmt(res$J$stat, 3), iv_fmt(res$J$p, 3)))
            else if (length(pp$instr) > length(pp$endog))
                items[[length(items) + 1]] <- c("Overidentification test", "not available with this estimator or too few clusters")
            for (i in seq_along(items)) sm$addRow(rowKey = i, values = list(item = items[[i]][1], value = items[[i]][2]))

            notes <- character(0)
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            if (!is.null(res$G) && res$G < 30)
                notes <- c(notes, sprintf("Only %d clusters: cluster-robust standard errors can be unreliable with fewer than about 30 clusters.", res$G))
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        }
    )
)
