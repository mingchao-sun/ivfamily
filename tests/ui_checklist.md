# ivfamily: checklist for testing in jamovi

Menu: IV (analysis group). Install: `jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`.
Open sample data from File > Open > Data Library > Instrumental Variables Family.

## 1. IV Regression, cigarette demand
Dataset: Cigarette demand (Stock and Watson). Outcome `lq`; Endogenous `lp`; Excluded instruments `tdiff`, `tax_r`; Controls `li`; Standard errors: Robust (HC1).
Expect: lp -1.277, s.e. 0.250; OLS table present; first-stage F 209.7; Hansen J 0.335 (df 1, p 0.563).
Then: switch Estimator to LIML, Fuller, k-class (k = 0.5), Two-step GMM (needs robust); tick Method notes and Data and setup guide.

## 2. IV Regression, college distance
Dataset: College distance. Outcome `wage`; Endogenous `education`; Instruments `distance`, `tuition`; Controls `gender`, `ethnicity` (set both to Nominal) and `unemp`; Robust.
Expect: education 0.915 (s.e. 0.163), OLS 0.0055, Hansen J about 110 (rejects).
Then: add a cluster variable (any nominal or ID variable) with Cluster-robust; add an absorbed fixed effect.

## 3. First Stage and Weak Instruments
Dataset: Simulated two endogenous regressors. Outcome `y`; Endogenous `x1`, `x2`; Instruments `z1`-`z4`; Controls `w`; Robust.
Expect: Kleibergen-Paap rk Wald F 74.68, first-stage F x1 95.33 and x2 82.25, conditional F 116.34 and 100.76.
With only `x1` as endogenous: effective F 96.03 and four critical values (29.1, 16.9, 10.4, 8.0).

## 4. Overidentification and Endogeneity
Same dataset and variables, but move `z3` out of Excluded instruments (keep `z1`, `z2`, `z4`) into "Instruments to test": a variable can sit in one box only, and the tested instruments are used as instruments too. Robust.
Expect: Hansen J 6.104 (df 2); C test for z3 5.637 (df 1, p 0.018); endogeneity test (x1, x2) 231.5 (df 2).

## 5. Weak-Instrument-Robust Inference
Dataset: Simulated weak instruments. Outcome `y`; Endogenous `x`; Instruments `z1`, `z2`, `z3`; Controls `w`; Standard errors: Robust. Hypothesised coefficient 0.
Expect: AR 6.990 (p 0.072), K 3.900 (p 0.048), CLR 4.887 (p 0.044); sets: AR [-0.171, 1.484], K [0.019, 1.383] U [2.259, Inf], CLR [0.056, 1.382]; the plot shows three p-value curves and a dashed line at 0.05.
Then: set Hypothesised coefficient to 0.5 (CLR p 0.215); try Classical and Cluster (cluster variable `g`); select two endogenous regressors to see the message.

## 6. Compare IV Estimators
Dataset: Simulated weak instruments. Outcome `y`; Endogenous `x`; Instruments `z1`-`z3`; Controls `w`; Robust. Tick Fuller and Two-step GMM.
Expect: a table of OLS, 2SLS, LIML, Fuller and GMM estimates (true coefficient 0.5; 2SLS about 0.91), a forest plot, overidentification text for 2SLS and GMM. With Classical errors GMM is skipped with a note.

## 7. Binary Instrument and LATE
Dataset: Family size and labour supply. Outcome `weeks`; Binary treatment `morekids`; Binary instrument `samesex`; Covariates to compare `age`, `boy1`, `afam`, `hispanic`, `other`; Robust.
Expect: reduced form -0.426, first stage 0.0675, LATE -6.314 (s.e. 1.275); shares 34.6% always-takers, 58.6% never-takers, 6.8% compliers; complier boy1 mean 0.397 against 0.514 overall.

## 8. Plausibly Exogenous Instruments
Dataset: Cigarette demand. Outcome `lq`; Endogenous `lp`; Instruments `tdiff`, `tax_r`; Controls `li`; Robust (large-sample, untick small-sample corrections to match Stata). delta 0.05; LTZ standard deviation 0.02, mean 0.01.
Expect: tipping delta 0.0074; union interval at delta 0.05 [-7.00, 4.76]; LTZ interval [-5.72, 1.09]; the plot shows a widening band with a dashed zero line.

## 9. Shift-Share (Bartik) IV
Dataset: Simulated shift-share design (long form: one row per region and sector, 4,000 rows). Outcome `y`; Endogenous regressor `d`; Region `region`; Sector (shock) `sector`; Exposure share `share`; Shock `shock`; Controls `x1`; Weights `pop`; Region cluster `state`; Sector cluster `sector_group`.
Set `region`, `sector`, `state` and `sector_group` to ID or Nominal; the others Continuous.
Expect: coefficient 0.513; standard errors 0.124 (homoskedastic), 0.145 (robust), 0.104 (state clusters), AKM 0.188, AKM0 0.392 with confidence set 0.046 to 1.583; first stage 0.938 (AKM s.e. 0.065, F 207); info table: 500 regions, 40 sectors, effective number of sectors 39.2, shares add up to 1.
Then: remove `sector_group` -> AKM 0.246 and AKM0 0.304 (AKM0 set 0.022 to 1.215); remove `state` -> the region-cluster row disappears; tick "Control for the sum of the shares" (no change here, the shares add up to 1); untick "First stage and reduced form" (the table disappears); null value 0.5; confidence level 90.
Messages to try: drop `share` (a hint asks for the six required variables); put `y` in Sector cluster (message: varies within a sector or a region); delete the rows of one region from the middle (still works); delete all but five regions (message: fewer than 10 regions).

## Things to look at in the interface
- Collapsed groups open and close; options that do not apply are disabled (Fuller alpha, k).
- Messages for wrong input are in plain words (no instruments, cluster without a variable, GMM with classical errors).
- Tables fill in as variables are dropped into the boxes and clear when a variable is removed.
- Number formats, p-values (< .001), long variable names, nominal controls that expand into dummies.
- The p-value plot in Weak-Instrument-Robust Inference: three curves, legend, axis labels, theme and dark mode.

## 10. Marginal Treatment Effects (IV > Estimation)
Dataset: Simulated schooling choice (marginal treatment effects). Outcome `logwage`; Binary treatment `college` (yes or no); Instruments `distance`, `tuition`; Covariates `score`, `female`; Degree of the polynomial 1; Probit; 300 bootstrap samples (takes several seconds).
Expect: propensity score table: all 6000 (mean 0.540), treated 3241 (0.660), untreated 2759 (0.398), lowest about 0 and highest 0.995; ATE 0.455 (s.e. 0.041, interval 0.376 to 0.535), ATT 0.610 (0.068), ATUT 0.274 (0.091), 2SLS 0.498 (0.044); tests: essential heterogeneity chi-square 4.14 (p 0.042), observed heterogeneity p 0.001; the curve table has 19 rows from u = 0.05 (0.762) to 0.95 (0.149) with 0.455 at 0.50. The first plot shows a falling curve with a shaded band, solid where the propensity score has support and dashed beyond, and an orange line at the ATE; the second plot shows two overlapping histograms.
Then: degree 2 (ATE 0.440, test p 0.12), degree 3, 4; logit; remove the covariates; add a cluster variable; untick the plot. Method notes and guide on/off.
Messages: a continuous treatment -> "must be binary"; a variable in two boxes; an instrument that is random noise -> either a message (second stage not identified, most bootstrap samples failing) or very wide intervals and large standard errors.
