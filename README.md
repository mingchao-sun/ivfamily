# ivfamily: Instrumental Variables Family for jamovi

A jamovi module for instrumental-variable (IV) analysis, in the menu group **IV**. The estimation engine is written in base R
(no dependency on `AER`, `ivreg` or `fixest`) and is checked against R packages and Stata (`ivreg2`, `ivreghdfe`, `weakiv`,
`weakivtest`, `plausexog`).

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

Every analysis has *Method notes* and *Data and setup guide* in its Output section. Five example data sets are in the Data Library.

## Example data (`data/`)

| File | Source | Used for |
|---|---|---|
| `cigarettes_sw.csv` | `AER::CigarettesSW`, 1995 | IV Regression, Plausibly Exogenous Instruments |
| `college_distance.csv` | `AER::CollegeDistance` | IV Regression |
| `fertility.csv` | `AER::Fertility` (Angrist and Evans 1998) | Binary Instrument and LATE |
| `sim_two_endogenous.csv` | simulated (`dev/` and `stata_validation/`) | First Stage, Overidentification |
| `sim_weak.csv` | simulated | Weak-Instrument-Robust Inference, Compare IV Estimators |

`dev/make_data.R` rebuilds the AER-based files.

## Design notes and limits

- Weights are analytic (as Stata's `aweight`): rescaling them changes nothing, replicating rows is not equivalent.
- Weak-instrument-robust inference and the effective F are for one endogenous regressor. With several, use the conditional F
  and the Kleibergen-Paap statistic in *First Stage and Weak Instruments*.
- Not included: tF adjustments, TSLS- and LIML-specific effective-F critical values, shift-share and examiner designs,
  marginal treatment effects, more than two sets of fixed effects.
- Robust LIML standard errors differ from `ivreg2` by about 0.07%; everything else agrees to at least five digits
  (see `PLAN.md` at the project root and `tests/`).

## Build and test

```r
jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")
```

`Rscript tests/smoke.R` runs every analysis headless (about 50 scenarios, including bad input); `Rscript tests/reference.R`
checks the engine against `AER`, `sandwich` and `fixest`. Run both in jamovi's own R with the bundled libraries too
(see `PLAN.md`). `tests/ui_checklist.md` lists what to look at by hand in jamovi. The Stata comparisons are in
`../stata_validation/` (`compare_iv.R`, `compare_weak.R`, `compare_effective.R`, `compare_ar.R`, `compare_plaus_late.R`).

## Ribbon icon

`icon/analysis-iv.svg` (74 x 44, same size and colours as the other families): instrument, regressor and outcome as three circles with
arrows, the effect arrow in orange, and a dashed unobserved confounder acting on the regressor and the outcome.
jmvtools 28.17 cannot attach an icon to a module, so `dev/icon_preview.ps1` adds one CSS rule for the menu group IV to jamovi's client
stylesheet (run it in an administrator PowerShell, see the comments at the top; restart jamovi afterwards). `dev/icon_preview.html`
shows the icon at several sizes next to the other families.

## Marginal Treatment Effects (0.3.0)

`R/engine_mte.R` estimates marginal treatment effects (Heckman and Vytlacil, 2005) by local instrumental variables with a polynomial in the propensity score, in the specification of the Stata command `mtefe` (Andresen, 2018): marginal effect curve, ATE, ATT, ATUT, the 2SLS estimate for comparison, bootstrap inference, and Wald tests of essential and observed heterogeneity. `tests/mte.R` recovers a known marginal effect in a generalised Roy model (parameters, curve, coverage of the bootstrap interval); `stata_validation/18_mtefe.do` and `compare_mtefe.R` compare with `mtefe`. Example data: `data/sim_mte.csv` (generator `dev/make_data_mte.R`).
