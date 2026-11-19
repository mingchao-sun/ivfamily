# Helpers shared by the IV analyses. Nothing here depends on jamovi except iv_guard() and iv_show_text().

iv_fmt <- function(x, d = 3) formatC(x, format = "f", digits = d, drop0trailing = TRUE)

# Turn engine errors into readable jamovi messages.
iv_guard <- function(expr) tryCatch(expr, error = function(e) {
    m <- conditionMessage(e)
    if (grepl("singular|rank", m, ignore.case = TRUE) && !grepl("rank deficient|collinear|identify", m))
        m <- paste0("The model could not be estimated; the variables may be perfectly collinear or have no variation. (", m, ")")
    jmvcore::reject(m)
})

iv_method_label <- function(method, k) {
    switch(method, ols = "OLS", `2sls` = "2SLS", liml = sprintf("LIML (k = %s)", iv_fmt(k, 4)),
           fuller = sprintf("Fuller (k = %s)", iv_fmt(k, 4)), kclass = sprintf("k-class (k = %s)", iv_fmt(k, 4)),
           gmm = "Two-step GMM")
}

iv_vce_label <- function(res) {
    base <- switch(res$vce, classical = "Classical standard errors.", robust = "Heteroskedasticity-robust (HC1) standard errors.",
                   cluster = "Cluster-robust standard errors.", cluster2 = "Two-way cluster-robust standard errors.")
    if (isTRUE(res$small)) paste(base, "Small-sample corrections applied.") else base
}

# Guide and method notes shown at the top of the analysis when ticked in the Output section.
iv_show_text <- function(results, opts, key) {
    if (isTRUE(opts$showGuide)) results$guide$setContent(iv_guide_html(key))
    if (isTRUE(opts$showHelp)) results$helpText$setContent(iv_notes_html(key))
}

# Shared set-up of the diagnostic analyses: validates the selections, builds the sample and fits 2SLS with the chosen errors.
iv_setup <- function(opts, data) {
    pp <- iv_guard(iv_prepare(data, opts$outcome, opts$endog, opts$instruments, opts$controls, opts$absorb, opts$cluster, opts$weights))
    if (opts$vce %in% c("cluster", "cluster2") && length(opts$cluster) < (if (opts$vce == "cluster2") 2 else 1))
        jmvcore::reject(if (opts$vce == "cluster2") "Two-way clustering needs two cluster variables."
                        else "Select a cluster variable for cluster-robust standard errors.")
    res <- iv_guard(iv_fit(pp, "2sls", vce = opts$vce, small = opts$small))
    list(pp = pp, res = res)
}

iv_ready <- function(opts) !is.null(opts$outcome) && length(opts$endog) > 0 && length(opts$instruments) > 0
