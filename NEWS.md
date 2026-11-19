# ivfamily 1.0.0

* First public release (jamovi library submission). Contains everything listed below; the earlier 0.x entries are development history.

# ivfamily 0.3.0

* New: **Marginal Treatment Effects** (menu IV > Estimation): the effect of a binary treatment for people at each margin of choosing it (Heckman and Vytlacil, 2005), by local instrumental variables with a polynomial in the propensity score as in the Stata command `mtefe`: marginal effect curve with bootstrap bands, ATE, ATT and ATUT, the 2SLS estimate for comparison, tests of essential and observed heterogeneity, propensity score overlap plot. Recovers a known marginal effect in simulations (`tests/mte.R`).

# ivfamily 0.2.0

- **New analysis: Shift-Share (Bartik) IV** (menu IV > Estimation). IV with an instrument built from sector shocks and regional exposure shares, in long form (one row per region and sector). Homoskedastic, robust and region-cluster standard errors next to the exposure-robust AKM and AKM0 inference of Adao, Kolesar and Morales (2019), optional sector clusters, weights, controls and a control for the sum of the shares; first stage and reduced form with AKM standard errors and the effective number of sectors. Base-R engine (`R/engine_ss.R`); results equal the R package ShiftShareSE (ADH China-shock data and the sample data, `tests/shiftshare.R`). Sample data: Simulated shift-share design.

# ivfamily 0.1.0

First version. Seven analyses in the menu group IV:

* **IV Regression**: OLS comparison, 2SLS, LIML, Fuller, k-class, two-step GMM; classical, robust, cluster-robust and two-way cluster-robust standard errors; absorbed fixed effects (one or two sets); analytic weights; first-stage F; Sargan and Hansen J.
* **Compare IV Estimators**, **Binary Instrument and LATE**.
* **First Stage and Weak Instruments**: conditional F, Cragg-Donald, Kleibergen-Paap rk, effective F with simplified critical values, reduced form.
* **Weak-Instrument-Robust Inference**: AR, K and CLR tests and confidence sets for one endogenous regressor.
* **Overidentification and Endogeneity**: Hansen J, C tests, endogeneity test.
* **Plausibly Exogenous Instruments**: union of confidence intervals, tipping point, local-to-zero approximation.

Every analysis offers *Method notes* and *Data and setup guide* (Output section). Estimates agree with Stata `ivreg2`, `ivreghdfe`, `weakiv`, `weakivtest` and `plausexog` and with R `AER`, `sandwich` and `fixest` on the checks listed in README.md.
