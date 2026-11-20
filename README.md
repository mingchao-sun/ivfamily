# ivfamily — instrumental variables in jamovi

A jamovi module for instrumental-variable (IV) analysis: estimators, first-stage and weak-instrument diagnostics, weak-instrument-robust inference,
overidentification tests, plausibly exogenous instruments, LATE, shift-share designs and marginal treatment effects. The estimation engine is written in
base R (no dependency on `AER`, `ivreg` or `fixest`) and is checked against R packages and Stata (`ivreg2`, `ivreghdfe`, `weakiv`, `weakivtest`,
`plausexog`, `mtefe`). Chinese summary: [README.zh.md](README.zh.md). Japanese summary: [README.ja.md](README.ja.md).

## Analyses (jamovi menu: IV)

| Analysis | What it does |
|---|---|
| **IV Regression** | OLS comparison, 2SLS, LIML, Fuller, k-class and two-step GMM; classical, robust (HC1), cluster-robust and two-way cluster-robust standard errors; one or two sets of absorbed fixed effects; analytic weights; first-stage F and Sargan / Hansen J |
| **Compare IV Estimators** | OLS, 2SLS, LIML, Fuller and GMM side by side, with a forest plot |
| **Binary Instrument and LATE** | Binary treatment and binary instrument: reduced form, first stage, LATE, shares of always-takers, never-takers and compliers, covariate means by group (Abadie kappa) |
| **First Stage and Weak Instruments** | First-stage F and partial R2, Sanderson-Windmeijer conditional F, Cragg-Donald F, Kleibergen-Paap rk Wald F, Montiel Olea-Pflueger effective F with simplified critical values, reduced form |
| **Weak-Instrument-Robust Inference** | Anderson-Rubin, Kleibergen K and conditional likelihood ratio tests and confidence sets (one endogenous regressor), with a p-value plot |
| **Overidentification and Endogeneity** | Sargan / Hansen J, difference-in-Hansen (C) tests for chosen instruments or each instrument in turn, endogeneity test |
| **Plausibly Exogenous Instruments** | Conley, Hansen and Rossi (2012): union of confidence intervals, tipping point, local-to-zero approximation |
| **Shift-Share (Bartik) IV** | Instrument built from sector shocks and regional exposure shares; AKM and AKM0 exposure-robust inference (Adao, Kolesar and Morales, 2019), first stage, effective number of sectors |
| **Marginal Treatment Effects** | Local IV with a polynomial in the propensity score (Heckman and Vytlacil, 2005; the specification of Stata `mtefe`): marginal effect curve, ATE, ATT, ATUT, the 2SLS estimate for comparison, bootstrap inference, tests of essential and observed heterogeneity, propensity score overlap plot |

Every analysis has *Method notes* and *Data and setup guide* in its Output section.

## Install
- **jamovi library**: in jamovi open the modules button (top right), *jamovi library*, find *Instrumental Variables Family* and install (once the module is listed).
- **Sideload**: download the `.jmo` file and use *Install from file* in the same menu. A `.jmo` works only on the operating system, processor type and jamovi version it was built for.
- **From source** (needs R and the jamovi compiler `jmvtools`; adjust the path to your jamovi folder):
```r
jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")
```

## Use from R (Rj editor, scripts)
Every analysis is also an R function. Pass `NULL` for the arguments you do not use, because they have no default values. For example, after installing the module,
```r
ivfamily::iv(data = d, outcome = "lq", endog = "lp", instruments = c("tdiff", "tax_r"), controls = "li", absorb = NULL, cluster = NULL, weights = NULL, vce = "robust")
```

## Data format
One row per observation. A numeric outcome, one or more endogenous regressors (numeric; binary 0/1 for the LATE and marginal-treatment-effect analyses),
instruments, optional controls, cluster and fixed-effect variables and weights. Rows with missing values are dropped and reported. The shift-share analysis
takes long-form data with one row per region and sector.

## Example data
Shown in jamovi under **Open > Data Library > Instrumental Variables Family**; each entry says which analysis to use and which variable goes where.
Files in `data/` (`dev/make_data.R` rebuilds the AER-based files):

| File | Source | Used for |
|---|---|---|
| `cigarettes_sw.csv` | `AER::CigarettesSW`, 1995 | IV Regression, Plausibly Exogenous Instruments |
| `college_distance.csv` | `AER::CollegeDistance` | IV Regression |
| `fertility.csv` | `AER::Fertility` (Angrist and Evans 1998) | Binary Instrument and LATE |
| `sim_two_endogenous.csv` | simulated | First Stage, Overidentification |
| `sim_weak.csv` | simulated | Weak-Instrument-Robust Inference, Compare IV Estimators |
| `sim_shiftshare.csv` | simulated: 500 regions, 40 sectors | Shift-Share IV |
| `sim_mte.csv` | simulated schooling choice (generator `dev/make_data_mte.R`) | Marginal Treatment Effects |

## References
Papers and R packages behind each analysis appear in the References section of the jamovi output (defined in `jamovi/00refs.yaml`).

## Validation
`tests/reference.R` checks the engine against `AER`, `sandwich` and `fixest`; `tests/shiftshare.R` and `tests/mte.R` check the shift-share and marginal-treatment
designs (the latter recovers a known marginal effect in a generalised Roy model: parameters, curve and coverage of the bootstrap interval). The Stata comparisons are in
`../stata_validation/` (`compare_iv.R`, `compare_weak.R`, `compare_effective.R`, `compare_ar.R`, `compare_plaus_late.R`, `compare_mtefe_log.R`).
Robust LIML standard errors differ from `ivreg2` by about 0.07%; everything else agrees to at least five digits; the marginal treatment effect curve and its
ATE, ATT and ATUT agree with `mtefe` to 2e-7. `tests/smoke.R` runs every analysis headless (dozens of scenarios, including bad input).

## Known limitations
- Weights are analytic (as Stata's `aweight`): rescaling them changes nothing, replicating rows is not equivalent.
- Weak-instrument-robust inference and the effective F are for one endogenous regressor. With several, use the conditional F and the Kleibergen-Paap statistic in *First Stage and Weak Instruments*.
- Not included: tF adjustments, TSLS- and LIML-specific effective-F critical values, examiner designs, more than two sets of fixed effects.
- Marginal treatment effects need an instrument that moves the propensity score over a wide range; the curve is estimated (and shown solid) only where the propensity score has support.

## Notes for developers
- Layout: `R/engine_*.R` is plain R with no jamovi dependency; `R/*.b.R` and `jamovi/*.yaml` are the jamovi front end.
- Tests: `Rscript tests/smoke.R`, `tests/reference.R`, `tests/shiftshare.R`, `tests/mte.R`; `tests/ui_checklist.md` lists what to check by hand in jamovi. Run them in jamovi's own R with the bundled libraries as well.
- Ribbon icon: `icon/analysis-iv.svg` (74 x 44): instrument, regressor and outcome with a dashed unobserved confounder. jmvtools 28.17 cannot attach an icon to a module, so `dev/icon_preview.ps1` adds one CSS rule for the menu group to jamovi's client stylesheet (administrator PowerShell, restart jamovi afterwards).

## How to cite
Sun, M. *ivfamily: Instrumental Variables Family for jamovi* (jamovi module). https://github.com/mingchao-sun/ivfamily. Please also cite the papers behind the analysis you use (listed in the output).

## License
GPL-3, as stated in `DESCRIPTION`; third-party code and data are listed in `NOTICE.md`. Bug reports and suggestions: https://github.com/mingchao-sun/ivfamily/issues
