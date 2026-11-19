# jamovi analysis class: side-by-side comparison of OLS, 2SLS, LIML, Fuller and two-step GMM.

ivcmpClass <- R6::R6Class(
    "ivcmpClass",
    inherit = ivcmpBase,
    private = list(

        .run = function() {
            opts <- self$options
            iv_show_text(self$results, opts, "ivcmp")
            if (!iv_ready(opts)) {
                self$results$notes$setContent("<p>Select an outcome, at least one endogenous regressor and at least as many excluded instruments.</p>")
                return()
            }
            pp <- iv_guard(iv_prepare(self$data, opts$outcome, opts$endog, opts$instruments, opts$controls, opts$absorb, opts$cluster, opts$weights))
            if (opts$vce %in% c("cluster", "cluster2") && length(opts$cluster) < (if (opts$vce == "cluster2") 2 else 1))
                jmvcore::reject(if (opts$vce == "cluster2") "Two-way clustering needs two cluster variables." else "Select a cluster variable for cluster-robust standard errors.")

            list_m <- list(c("ols", "OLS", isTRUE(opts$showOls)), c("2sls", "2SLS", TRUE), c("liml", "LIML", isTRUE(opts$showLiml)),
                           c("fuller", sprintf("Fuller (alpha = %s)", iv_fmt(opts$alpha, 2)), isTRUE(opts$showFuller)), c("gmm", "Two-step GMM", isTRUE(opts$showGmm)))
            tb <- self$results$est
            rows <- list(); notes <- character(0); k <- 0
            for (m in list_m) {
                if (m[3] != "TRUE") next
                if (m[1] == "gmm" && opts$vce == "classical") { notes <- c(notes, "Two-step GMM was skipped: it needs robust or cluster-robust standard errors."); next }
                private$.checkpoint()
                res <- tryCatch(iv_fit(pp, m[1], alpha = opts$alpha, vce = opts$vce, small = opts$small), error = function(e) conditionMessage(e))
                if (is.character(res)) { notes <- c(notes, sprintf("%s could not be estimated: %s", m[2], res)); next }
                ct <- iv_coef_table(res, opts$level)
                ct <- ct[ct$type == "Endogenous", , drop = FALSE]
                jtxt <- if (!is.null(res$J)) sprintf("%s %s (p = %s)", res$J$label, iv_fmt(res$J$stat, 2), iv_fmt(res$J$p, 3)) else ""
                for (i in seq_len(nrow(ct))) {
                    k <- k + 1
                    tb$addRow(rowKey = k, values = list(method = m[2], term = ct$term[i], estimate = ct$estimate[i], se = ct$se[i], lo = ct$lo[i], hi = ct$hi[i],
                                                         j = if (i == 1) jtxt else ""))
                    rows[[k]] <- data.frame(method = m[2], term = ct$term[i], estimate = ct$estimate[i], lo = ct$lo[i], hi = ct$hi[i], stringsAsFactors = FALSE)
                }
            }
            tb$setNote("n", paste0(iv_vce_label(list(vce = opts$vce, small = opts$small)),
                                   " OLS treats the endogenous regressors as exogenous. LIML and Fuller are less biased than 2SLS with many or weak instruments. ",
                                   "Large differences between the IV estimators point to weak or invalid instruments."))
            if (isTRUE(opts$plot) && length(rows)) self$results$cmpPlot$setState(do.call(rbind, rows))
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        },

        .plotCmp = function(image, ggtheme, theme, ...) {
            df <- image$state
            if (is.null(df)) return(FALSE)
            df$method <- factor(df$method, levels = rev(unique(df$method)))
            g <- ggplot2::ggplot(df, ggplot2::aes(y = method, x = estimate, xmin = lo, xmax = hi)) +
                ggplot2::geom_pointrange() +
                ggplot2::labs(x = "Estimate and confidence interval", y = NULL) + ggtheme
            if (length(unique(df$term)) > 1) g <- g + ggplot2::facet_wrap(~ term, scales = "free_x")
            print(g)
            TRUE
        }
    )
)
