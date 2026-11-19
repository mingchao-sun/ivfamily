# jamovi analysis class: marginal treatment effects by local instrumental variables (polynomial).

ivmteClass <- R6::R6Class(
    "ivmteClass",
    inherit = ivmteBase,
    private = list(

        .run = function() {
            o <- self$options
            iv_show_text(self$results, o, "ivmte")
            if (is.null(o$outcome) || is.null(o$treat) || length(o$instruments) == 0) {
                self$results$notes$setContent("<p>Select the outcome, a binary treatment and at least one instrument to begin. Covariates are optional.</p>")
                for (nm in c("support", "pars", "tests", "curve")) self$results[[nm]]$setVisible(FALSE)
                return()
            }
            if (anyDuplicated(c(o$outcome, o$treat, o$instruments, o$covs, o$cluster))) jmvcore::reject("A variable can sit in one box only.")
            dat <- self$data
            tv <- dat[[o$treat]]
            d <- if (is.factor(tv)) { lv <- levels(droplevels(tv)); if (length(lv) != 2) jmvcore::reject("A nominal treatment must have exactly two levels (the second level is treated)."); as.numeric(droplevels(tv) == lv[2]) } else as.numeric(tv)
            if (!all(d %in% c(0, 1, NA))) jmvcore::reject("The treatment must be binary (0/1).")
            y <- dat[[o$outcome]]; if (is.factor(y)) jmvcore::reject("The outcome must be numeric (set its type to Continuous in the Data tab).")
            mm <- function(vars) {
                df <- dat[, vars, drop = FALSE]
                for (v in vars) if (is.factor(df[[v]])) df[[v]] <- droplevels(df[[v]]) else df[[v]] <- as.numeric(df[[v]])
                f <- stats::reformulate(paste0("`", vars, "`"))
                M <- stats::model.matrix(f, stats::model.frame(f, df, na.action = stats::na.pass))[, -1, drop = FALSE]
                M[, apply(M, 2, function(z) stats::sd(z, na.rm = TRUE)) > 0, drop = FALSE]
            }
            Z <- mm(o$instruments); X <- if (length(o$covs)) mm(o$covs) else NULL
            cl <- if (!is.null(o$cluster)) dat[[o$cluster]] else NULL
            keep <- stats::complete.cases(y, d, Z, if (!is.null(X)) X, if (!is.null(cl)) cl)
            n_drop <- length(y) - sum(keep)
            y <- as.numeric(y)[keep]; d <- d[keep]; Z <- Z[keep, , drop = FALSE]; if (!is.null(X)) X <- X[keep, , drop = FALSE]; if (!is.null(cl)) cl <- cl[keep]
            if (min(sum(d == 1), sum(d == 0)) < 30) jmvcore::reject("Each treatment group needs at least 30 units.")
            private$.checkpoint()
            res <- iv_guard(mte_boot(y, d, X, Z, o$degree, o$link, cl, o$nboot, o$level / 100, o$seed, function() private$.checkpoint()))
            est <- res$est
            sb <- self$results$support
            for (g in list(list("All", rep(TRUE, length(d))), list("Treated", d == 1), list("Untreated", d == 0)))
                sb$addRow(rowKey = g[[1]], values = list(group = g[[1]], n = sum(g[[2]]), mean = mean(est$p[g[[2]]]), min = min(est$p[g[[2]]]), max = max(est$p[g[[2]]])))
            sb$setNote("a", sprintf("The marginal treatment effect is identified only where the propensity score varies; between the lowest and the highest score above the curve is an extrapolation of the polynomial. %s%% of the units have a score below 0.01 or above 0.99.", iv_fmt(100 * est$support$share_outside, 1)))
            pb <- self$results$pars
            labs <- c(ATE = "Average treatment effect (ATE)", ATT = "Effect on the treated (ATT)", ATUT = "Effect on the untreated (ATUT)", `2SLS` = "2SLS estimate (a weighted average of marginal effects)")
            for (i in seq_len(nrow(res$parameters))) {
                pr <- res$parameters[i, ]
                pb$addRow(rowKey = i, values = list(what = labs[[rownames(res$parameters)[i]]], est = pr[["estimate"]], se = pr[["se"]], lo = pr[["lo"]], hi = pr[["hi"]], p = pr[["p"]]))
            }
            pb$setNote("a", sprintf("Standard errors and intervals from %d bootstrap samples that repeat the first stage and the second stage. The parameters average the marginal effect over the resistance to treatment u = 0.01, ..., 0.99, with the weights of Heckman and Vytlacil as in Stata mtefe.", res$nboot))
            tb <- self$results$tests
            tu <- res$test_unobserved; to <- res$test_observed
            tb$addRow(rowKey = 1, values = list(what = "Essential heterogeneity: does the effect vary with the unobserved resistance u?", chisq = tu[["chisq"]], df = tu[["df"]], p = tu[["p"]]))
            if (to[["df"]] > 0) tb$addRow(rowKey = 2, values = list(what = "Observed heterogeneity: does the effect vary with the covariates?", chisq = to[["chisq"]], df = to[["df"]], p = to[["p"]]))
            tb$setNote("a", "Wald tests with the bootstrap covariance matrix. If the first is not rejected the marginal effect is flat in u: the effect does not depend on how reluctant people are, and ATE, ATT, ATUT and the 2SLS estimate coincide.")
            cb <- self$results$curve
            sh <- which(abs(res$curve[, "u"] * 20 - round(res$curve[, "u"] * 20)) < 1e-9)
            for (j in seq_along(sh)) { r <- res$curve[sh[j], ]; cb$addRow(rowKey = j, values = list(u = r[["u"]], mte = r[["mte"]], se = r[["se"]], lo = r[["lo"]], hi = r[["hi"]])) }
            cb$setNote("a", "The effect at the mean of the covariates for people who are just indifferent between treatment and no treatment when their resistance to treatment is u: small u means those who take the treatment even when it is unattractive to most.")
            notes <- sprintf("%d units (%d treated). Polynomial of degree %d in u (the outcome equation has the propensity score up to the power %d); %s first stage.", length(y), sum(d == 1), o$degree, o$degree + 1, if (o$link == "probit") "probit" else "logit")
            if (n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", n_drop))
            if (est$support$share_outside > 0.1) notes <- c(notes, "Many scores are close to 0 or 1: the instruments are very strong or the model separates the groups, which makes the extremes of the curve rely on few units.")
            self$results$notes$setContent(paste0("<p>", paste(notes, collapse = " "), "</p>"))
            if (isTRUE(o$plot)) {
                self$results$mtePlot$setState(list(curve = res$curve, rng = est$support$all, ate = res$parameters["ATE", "estimate"]))
                self$results$psPlot$setState(list(p = est$p, d = d))
            }
        },

        .plotMte = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            df <- as.data.frame(st$curve)
            ins <- df$u >= st$rng[1] & df$u <= st$rng[2]
            g <- ggplot2::ggplot(df, ggplot2::aes(x = u, y = mte)) + ggtheme +
                ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "#6b9de8", alpha = 0.18) +
                ggplot2::geom_line(colour = "#6b9de8", linewidth = 0.9, linetype = "dashed") +
                ggplot2::geom_line(data = df[ins, ], colour = "#6b9de8", linewidth = 1.2) +
                ggplot2::geom_hline(yintercept = st$ate, colour = "#e6ac40", linewidth = 0.8) + ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = "dotted") +
                ggplot2::labs(x = "Resistance to treatment u", y = "Marginal treatment effect", caption = "Solid: within the range of the propensity score; dashed: extrapolated. Orange: ATE")
            print(g)
            TRUE
        },

        .plotPs = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            df <- data.frame(p = st$p, group = factor(ifelse(st$d == 1, "Treated", "Untreated"), levels = c("Untreated", "Treated")))
            g <- ggplot2::ggplot(df, ggplot2::aes(x = p, fill = group)) + ggtheme +
                ggplot2::geom_histogram(position = "identity", alpha = 0.55, bins = 30, colour = NA) +
                ggplot2::scale_fill_manual(values = c(Untreated = "#0072B2", Treated = "#D55E00"), name = NULL) + ggplot2::xlim(0, 1) +
                ggplot2::labs(x = "Propensity score P(treated | instruments, covariates)", y = "Count", caption = "The marginal effect is identified where both groups are present")
            print(g)
            TRUE
        }
    )
)
