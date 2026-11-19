# jamovi analysis class: overidentification, orthogonality (C) and endogeneity tests.

ivoverClass <- R6::R6Class(
    "ivoverClass",
    inherit = ivoverBase,
    private = list(

        .run = function() {
            opts <- self$options
            iv_show_text(self$results, opts, "ivover")
            # A variable can sit in one box only, so the instruments to test and the regressors to test are in addition to the
            # other boxes: they are used as instruments / endogenous regressors too.
            o <- list(outcome = opts$outcome, endog = unique(c(opts$endog, opts$endogTest)), instruments = unique(c(opts$instruments, opts$suspects)),
                      controls = opts$controls, absorb = opts$absorb, cluster = opts$cluster, weights = opts$weights, vce = opts$vce, small = opts$small)
            if (!iv_ready(o)) {
                self$results$notes$setContent("<p>Select an outcome, at least one endogenous regressor and at least as many excluded instruments (instruments to test count as instruments).</p>")
                return()
            }
            st <- iv_setup(o, self$data)
            pp <- st$pp; res <- st$res
            nm <- res$endog; ins <- pp$instr
            lab <- function(v) paste(v, collapse = ", ")

            ot <- self$results$overid
            if (!is.null(res$J)) {
                ot$addRow(rowKey = 1, values = list(stat = res$J$label, value = res$J$stat, df = res$J$df, p = res$J$p))
                ot$setNote("h0", "H0: all excluded instruments are valid (uncorrelated with the error). Rejection means at least one instrument is invalid or effects differ across instruments.")
            } else if (length(ins) == length(nm)) {
                ot$setNote("exact", "The model is exactly identified: validity of the instruments cannot be tested.")
            } else {
                ot$setNote("na", "Not available: too few clusters for the cluster-robust weight matrix.")
            }

            if (length(opts$suspects) > 0) {
                ct <- self$results$ctest
                r <- iv_guard(iv_orthog_test(pp, res, opts$suspects))
                ct$addRow(rowKey = 1, values = list(label = lab(opts$suspects), C = r$C, df = r$df, p = r$p))
                ct$setNote("c", "C = difference between the Hansen J of the full instrument set and of the set without the tested instruments. H0: the tested instruments are valid, given that the others are.")
            }

            et <- self$results$endog
            tested <- if (length(opts$endogTest) > 0) opts$endogTest else nm
            r <- iv_guard(iv_endog_test(pp, res, tested))
            et$addRow(rowKey = 1, values = list(label = lab(tested), C = r$C, df = r$df, p = r$p))
            et$setNote("c", "H0: the tested regressors can be treated as exogenous (OLS is consistent). Rejection supports instrumenting them. Valid only if the remaining instruments are valid.")

            if (isTRUE(opts$showEach)) {
                tb <- self$results$each
                k <- 0
                for (z in ins) {
                    rr <- tryCatch(iv_orthog_test(pp, res, z), error = function(e) NULL)
                    k <- k + 1
                    if (is.null(rr)) tb$addRow(rowKey = k, values = list(label = z, C = NaN, df = 1L, p = NaN))
                    else tb$addRow(rowKey = k, values = list(label = z, C = rr$C, df = rr$df, p = rr$p))
                }
                tb$setNote("n", "Not available (blank) when dropping the instrument would leave the model underidentified or, with cluster errors, when clusters are too few.")
            }

            notes <- character(0)
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        }
    )
)
