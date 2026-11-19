# jamovi analysis class: sensitivity of the IV estimate to violations of the exclusion restriction.

ivplausClass <- R6::R6Class(
    "ivplausClass",
    inherit = ivplausBase,
    private = list(

        .run = function() {
            opts <- self$options
            iv_show_text(self$results, opts, "ivplaus")
            if (!iv_ready(opts)) {
                self$results$notes$setContent("<p>Select an outcome, at least one endogenous regressor and at least as many excluded instruments.</p>")
                return()
            }
            st <- iv_setup(opts, self$data)
            pp <- st$pp; res <- st$res
            lev <- opts$level
            tip <- iv_guard(iv_plaus_tipping_exact(pp, res, lev))
            idx <- match(res$endog, names(res$beta))
            z <- stats::qnorm(1 - (1 - lev / 100) / 2)
            A <- iv_plaus_A(res)

            tt <- self$results$tip
            for (j in seq_along(tip)) tt$addRow(rowKey = j, values = list(term = res$endog[j], estimate = res$beta[[idx[j]]], tip = tip[j]))
            tt$setNote("t", paste0("All excluded instruments are allowed a direct effect on the outcome between -delta and delta (in outcome units per unit of the instrument). ",
                                   "At the tipping delta the union of confidence intervals just reaches zero; a zero tipping value means the ordinary interval already contains zero. ",
                                   "A direct effect of this size, relative to what you consider plausible, tells you how fragile the result is."))

            deltas <- list()
            for (j in seq_along(tip)) {
                scale <- max(tip[j], z * res$se[[idx[j]]] / sum(abs(A[j, ])))
                dj <- scale * c(0, 0.25, 0.5, 0.75, 1, 1.5)
                if (opts$delta > 0) dj <- sort(unique(c(dj, opts$delta)))
                deltas[[j]] <- dj
            }
            ut <- self$results$uci
            k <- 0
            for (j in seq_along(tip)) {
                for (dl in deltas[[j]]) {
                    private$.checkpoint()
                    u <- iv_guard(iv_plaus_uci_exact(pp, res, -dl, dl, lev))[j, ]
                    k <- k + 1
                    ut$addRow(rowKey = k, values = list(term = res$endog[j], delta = dl, lo = u$lo, hi = u$hi, zero = if (u$lo <= 0 && u$hi >= 0) "yes" else "no"))
                }
            }
            ut$setNote("u", "Conley, Hansen and Rossi (2012): the union, over all direct effects in the box, of the confidence intervals of the 2SLS coefficient; 2SLS and its standard errors are recomputed at every corner of the box. delta = 0 is the conventional interval (normal critical values).")

            if (opts$ltzSd > 0) {
                lt <- iv_guard(iv_plaus_ltz(res, opts$ltzMean, opts$ltzSd^2, lev))
                tb <- self$results$ltz
                for (j in seq_len(nrow(lt))) tb$addRow(rowKey = j, values = as.list(lt[j, c("term", "estimate", "se", "lo", "hi")]))
                tb$setNote("l", "Direct effects are assumed to be normal with the given mean and standard deviation, the same for every instrument and independent; the estimate shifts by the implied bias and the variance grows.")
            }

            if (isTRUE(opts$plot)) {
                sc <- max(deltas[[1]]) * 1.2
                dg <- seq(0, sc, length.out = 25)
                b <- t(vapply(dg, function(dl) { u <- iv_plaus_uci_exact(pp, res, -dl, dl, lev)[1, ]; c(u$lo, u$hi) }, c(0, 0)))
                self$results$plausPlot$setState(list(delta = dg, lo = b[, 1], hi = b[, 2], est = res$beta[[idx[1]]], tip = tip[1], term = res$endog[1]))
            }
            notes <- character(0)
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            if (length(pp$instr) == length(pp$endog))
                notes <- c(notes, "The model is exactly identified; the exclusion restriction cannot be tested, which is when this sensitivity analysis matters most.")
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        },

        .plotPlaus = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            df <- data.frame(delta = st$delta, lo = st$lo, hi = st$hi)
            g <- ggplot2::ggplot(df, ggplot2::aes(x = delta)) +
                ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.25) +
                ggplot2::geom_line(ggplot2::aes(y = lo)) + ggplot2::geom_line(ggplot2::aes(y = hi)) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
                ggplot2::geom_hline(yintercept = st$est, linetype = "dotted") +
                ggplot2::labs(x = "delta: largest direct effect of an instrument on the outcome", y = paste("Union interval for", st$term)) + ggtheme
            if (st$tip > 0) g <- g + ggplot2::geom_vline(xintercept = st$tip, linetype = "dotted")
            print(g)
            TRUE
        }
    )
)
