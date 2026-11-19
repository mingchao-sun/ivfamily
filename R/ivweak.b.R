# jamovi analysis class: weak-instrument-robust tests (AR, K, CLR) and confidence sets for one endogenous regressor.

ivweakClass <- R6::R6Class(
    "ivweakClass",
    inherit = ivweakBase,
    private = list(

        .run = function() {
            opts <- self$options
            iv_show_text(self$results, opts, "ivweak")
            if (!iv_ready(opts)) {
                self$results$notes$setContent("<p>Select an outcome, one endogenous regressor and at least one excluded instrument.</p>")
                return()
            }
            if (length(opts$endog) != 1)
                jmvcore::reject("The weak-instrument-robust tests are available for one endogenous regressor.")
            st <- iv_setup(opts, self$data)
            pp <- st$pp; res <- st$res
            a <- iv_guard(iv_ar_setup(pp, res))
            b0 <- opts$nullValue
            tt <- iv_guard(iv_weak_tests(a, b0))

            tb <- self$results$tests
            tb$addRow(rowKey = 1, values = list(test = "Anderson-Rubin (AR)", stat = tt$AR, df = tt$df, p = tt$p_ar))
            tb$addRow(rowKey = 2, values = list(test = "Kleibergen (K)", stat = tt$K, df = 1L, p = tt$p_k))
            tb$addRow(rowKey = 3, values = list(test = "Conditional likelihood ratio (CLR)", stat = tt$CLR, df = tt$df, p = tt$p_clr))
            tb$setNote("h0", sprintf("H0: the coefficient of %s equals %s. AR is chi-square with %d df; K with 1 df; the CLR p-value is conditional on the rank statistic (rk = %s).",
                                     res$endog, iv_fmt(b0, 4), tt$df, iv_fmt(tt$rk, 3)))

            alpha <- 1 - opts$level / 100
            bh <- res$beta[[res$endog]]; se <- res$se[[res$endog]]
            grid <- seq(bh - opts$gridWidth * se, bh + opts$gridWidth * se, length.out = 1201)
            ct <- self$results$csets
            ci <- iv_coef_table(res, opts$level)
            ci <- ci[ci$term == res$endog, ]
            ct$addRow(rowKey = 1, values = list(test = "Wald (2SLS, not robust to weak instruments)", set = sprintf("[%s, %s]", iv_fmt(ci$lo, 3), iv_fmt(ci$hi, 3))))
            sets <- list(AR = "p_ar", K = "p_k", CLR = "p_clr")
            labs <- c(AR = "Anderson-Rubin (AR)", K = "Kleibergen (K)", CLR = "Conditional likelihood ratio (CLR)")
            k <- 1
            for (nm in names(sets)) {
                private$.checkpoint()
                cs <- iv_guard(iv_cset(function(b) iv_weak_tests(a, b)[[sets[[nm]]]], grid, alpha))
                k <- k + 1
                ct$addRow(rowKey = k, values = list(test = labs[[nm]], set = iv_cset_text(cs)))
            }
            ct$setNote("g", sprintf("Sets of coefficients not rejected at the %s%% level, found by inverting the tests on a grid of %s to %s (2SLS estimate +/- %s standard errors). ",
                                    iv_fmt(opts$level, 1), iv_fmt(grid[1], 2), iv_fmt(grid[length(grid)], 2), iv_fmt(opts$gridWidth, 0)))
            ct$setNote("inf", "A bound shown as -Inf or Inf means the set reaches the edge of the grid and may be unbounded; widen the grid to check. A set can consist of several intervals.")

            if (isTRUE(opts$plot)) {
                pg <- seq(bh - min(opts$gridWidth, 8) * se, bh + min(opts$gridWidth, 8) * se, length.out = 241)
                pv <- lapply(pg, function(b) { r <- iv_weak_tests(a, b); c(r$p_ar, r$p_k, r$p_clr) })
                pv <- do.call(rbind, pv)
                self$results$pPlot$setState(list(b = pg, p = pv, alpha = alpha, bh = bh, null = b0))
            }
            notes <- character(0)
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            fs <- iv_first_stage(pp, res, res$endog)
            if (fs$F > 20) notes <- c(notes, "The instruments look strong (first-stage F above 20); the robust sets should be close to the conventional interval.")
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        },

        .plotP = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            df <- data.frame(b = rep(st$b, 3), p = as.vector(st$p),
                             test = factor(rep(c("AR", "K", "CLR"), each = length(st$b)), levels = c("AR", "K", "CLR")))
            g <- ggplot2::ggplot(df, ggplot2::aes(x = b, y = p, linetype = test, colour = test)) +
                ggplot2::geom_line(linewidth = 0.8) +
                ggplot2::geom_hline(yintercept = st$alpha, linetype = "dashed") +
                ggplot2::geom_vline(xintercept = st$bh, linetype = "dotted") +
                ggplot2::labs(x = paste("Hypothesised coefficient of", self$options$endog), y = "p-value", colour = NULL, linetype = NULL) +
                ggplot2::coord_cartesian(ylim = c(0, 1)) + ggtheme
            print(g)
            TRUE
        }
    )
)
