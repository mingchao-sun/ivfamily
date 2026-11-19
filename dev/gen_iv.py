# Generates the .a/.r/.u.yaml files of ivfirst and ivover. Run from the module root: python dev/gen_iv.py
import os

BASE_OPTS = ["outcome", "endog", "instruments", "controls", "absorb", "cluster", "weights"]
CLEAR = ["outcome", "endog", "instruments", "controls", "absorb", "cluster", "weights", "vce", "small"]


def opt_vars():
    return """    - name: data
      type: Data

    - name: outcome
      title: Outcome
      type: Variable
      suggested:
          - continuous
      permitted:
          - numeric

    - name: endog
      title: Endogenous regressors
      type: Variables
      permitted:
          - numeric

    - name: instruments
      title: Excluded instruments
      type: Variables
      permitted:
          - numeric

    - name: controls
      title: Exogenous controls
      type: Variables
      suggested:
          - continuous
          - nominal
      permitted:
          - numeric
          - factor

    - name: absorb
      title: Absorbed fixed effects
      type: Variables
      permitted:
          - factor
          - id
          - numeric

    - name: cluster
      title: Cluster variables
      type: Variables
      permitted:
          - factor
          - id
          - numeric

    - name: weights
      title: Weights
      type: Variable
      permitted:
          - numeric

    - name: vce
      title: Standard errors
      type: List
      options:
          - name: classical
            title: Classical
          - name: robust
            title: Robust (HC1)
          - name: cluster
            title: Cluster-robust
          - name: cluster2
            title: Two-way cluster-robust
      default: classical

    - name: small
      title: Small-sample corrections (t and F distributions)
      type: Bool
      default: true

"""


def opt_tail(extra_bools):
    s = ""
    for n, t, d in extra_bools:
        s += f"    - name: {n}\n      title: {t}\n      type: Bool\n      default: {d}\n\n"
    s += """    - name: showHelp
      title: Method notes
      type: Bool
      default: false

    - name: showGuide
      title: Data and setup guide
      type: Bool
      default: false
"""
    return s


def head_a(name, title, sub, mtitle):
    return f"""---
name: {name}
title: {title}
jas: '1.2'
version: '0.1.0'
menuGroup: IV
menuSubgroup: {sub}
menuTitle: {mtitle}

options:
"""


def r_head(name, title):
    return f"""---
name: {name}
title: {title}
jrs: '1.1'

items:
    - name: guide
      title: Data and setup guide
      type: Html
      visible: (showGuide)

    - name: helpText
      title: Method notes
      type: Html
      visible: (showHelp)

    - name: notes
      title: Notes
      type: Html
      visible: true

"""


def table(name, title, cols, visible="true", extra_clear=()):
    s = f"    - name: {name}\n      title: '{title}'\n      type: Table\n      visible: {visible}\n      clearWith:\n"
    for c in list(CLEAR) + list(extra_clear):
        s += f"          - {c}\n"
    s += "      columns:\n"
    for n, t, ty in cols:
        s += f"          - name: {n}\n            title: '{t}'\n            type: {ty}\n"
        if n in ("p",):
            s += "            format: zto,pvalue\n"
    return s + "\n"


def u_vars(label_extra=""):
    s = """title: TITLE
name: NAME
jus: '3.0'
stage: 0
compilerMode: tame
children:
    - type: VariableSupplier
      persistentItems: false
      stretchFactor: 1
      children:
"""
    boxes = [("Outcome", "outcome", 1), ("Endogenous regressors", "endog", None), ("Excluded instruments", "instruments", None),
             ("Exogenous controls", "controls", None), ("Absorbed fixed effects (max 2)", "absorb", 2),
             ("Cluster variables (max 2)", "cluster", 2), ("Weights", "weights", 1)]
    for lab, nm, mx in boxes:
        s += f"          - type: TargetLayoutBox\n            label: {lab}\n            children:\n                - type: VariablesListBox\n                  name: {nm}\n                  isTarget: true\n                  height: small\n"
        if mx:
            s += f"                  maxItemCount: {mx}\n"
    return s


def u_tail(extra_boxes, extra_checks):
    s = ""
    s += """    - type: CollapseBox
      label: Standard errors
      collapsed: true
      children:
          - type: ComboBox
            name: vce
          - type: CheckBox
            name: small
"""
    if extra_checks:
        s += "    - type: CollapseBox\n      label: Tables\n      collapsed: true\n      children:\n"
        for c in extra_checks:
            s += f"          - type: CheckBox\n            name: {c}\n"
    s += """    - type: CollapseBox
      label: Output
      collapsed: true
      children:
          - type: CheckBox
            name: showHelp
            label: Method notes
          - type: CheckBox
            name: showGuide
            label: Data and setup guide
"""
    return s


def write(path, text):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


os.makedirs("jamovi", exist_ok=True)

# ---------------------------------------------------------------- ivfirst
write("jamovi/ivfirst.a.yaml",
      head_a("ivfirst", "First Stage and Weak Instruments", "Diagnostics", "First Stage and Weak Instruments") +
      opt_vars() + opt_tail([("showReduced", "Reduced form of the outcome on the excluded instruments", "true")]) + "...\n")
write("jamovi/ivfirst.r.yaml",
      r_head("ivfirst", "First Stage and Weak Instruments") +
      table("joint", "Joint statistics for all endogenous regressors",
            [("stat", "Statistic", "text"), ("value", "Value", "number"), ("df", "df", "text"), ("p", "p", "number")]) +
      table("each", "First stage, one endogenous regressor at a time",
            [("term", "Endogenous regressor", "text"), ("F", "F", "number"), ("df1", "df1", "integer"),
             ("df2", "df2", "number"), ("p", "p", "number"), ("r2", "Partial R2", "number"),
             ("swF", "Conditional F (Sanderson-Windmeijer)", "number"), ("swDf", "Conditional df1", "integer")]) +
      table("reduced", "Reduced form: outcome on the excluded instruments",
            [("term", "Instrument", "text"), ("estimate", "Estimate", "number"), ("se", "Std. error", "number"),
             ("p", "p", "number")], visible="(showReduced)") + "...\n")
write("jamovi/ivfirst.u.yaml",
      u_vars().replace("TITLE", "First Stage and Weak Instruments").replace("NAME", "ivfirst") + u_tail([], ["showReduced"]))

# ---------------------------------------------------------------- ivover
write("jamovi/ivover.a.yaml",
      head_a("ivover", "Overidentification and Endogeneity", "Diagnostics", "Overidentification and Endogeneity") +
      opt_vars() + """    - name: suspects
      title: Instruments to test (also used as instruments)
      type: Variables
      permitted:
          - numeric

    - name: endogTest
      title: Regressors to test (also endogenous; empty = all)
      type: Variables
      permitted:
          - numeric

""" + opt_tail([("showEach", "Drop each instrument in turn", "true")]) + "...\n")
write("jamovi/ivover.r.yaml",
      r_head("ivover", "Overidentification and Endogeneity") +
      table("overid", "Overidentification test (Sargan / Hansen J)",
            [("stat", "Statistic", "text"), ("value", "Value", "number"), ("df", "df", "integer"), ("p", "p", "number")],
            extra_clear=()) +
      table("ctest", "Orthogonality of the selected instruments (C test)",
            [("label", "Instruments tested", "text"), ("C", "C", "number"), ("df", "df", "integer"), ("p", "p", "number")],
            extra_clear=("suspects",)) +
      table("endog", "Endogeneity test (C test)",
            [("label", "Regressors tested", "text"), ("C", "C", "number"), ("df", "df", "integer"), ("p", "p", "number")],
            extra_clear=("endogTest",)) +
      table("each", "C test dropping each excluded instrument in turn",
            [("label", "Instrument dropped", "text"), ("C", "C", "number"), ("df", "df", "integer"), ("p", "p", "number")],
            visible="(showEach)") + "...\n")
u = u_vars().replace("TITLE", "Overidentification and Endogeneity").replace("NAME", "ivover")
u += """          - type: TargetLayoutBox
            label: Instruments to test (also used as instruments)
            children:
                - type: VariablesListBox
                  name: suspects
                  isTarget: true
          - type: TargetLayoutBox
            label: Regressors to test (also endogenous; empty = all)
            children:
                - type: VariablesListBox
                  name: endogTest
                  isTarget: true
"""
write("jamovi/ivover.u.yaml", u + u_tail([], ["showEach"]))

# ---------------------------------------------------------------- ivweak
write("jamovi/ivweak.a.yaml",
      head_a("ivweak", "Weak-Instrument-Robust Inference", "Diagnostics", "Weak-Instrument-Robust Inference") +
      opt_vars() + """    - name: nullValue
      title: Hypothesised coefficient (H0)
      type: Number
      default: 0

    - name: level
      title: Confidence level (%)
      type: Number
      default: 95
      min: 50
      max: 99.9

    - name: gridWidth
      title: Grid half-width (in 2SLS standard errors)
      type: Number
      default: 30
      min: 1

""" + opt_tail([("plot", "P-value plot", "true")]) + "...\n")
write("jamovi/ivweak.r.yaml",
      r_head("ivweak", "Weak-Instrument-Robust Inference") +
      table("tests", "Tests of the hypothesised coefficient",
            [("test", "Test", "text"), ("stat", "Statistic", "number"), ("df", "df", "integer"), ("p", "p", "number")],
            extra_clear=("nullValue",)) +
      table("csets", "Confidence sets",
            [("test", "Method", "text"), ("set", "Set", "text")], extra_clear=("level", "gridWidth")) +
      """    - name: pPlot
      title: P-values of the tests against the coefficient
      type: Image
      width: 640
      height: 420
      renderFun: .plotP
      visible: (plot)
      clearWith:
          - outcome
          - endog
          - instruments
          - controls
          - absorb
          - cluster
          - weights
          - vce
          - small
          - level
          - gridWidth

...\n""")
write("jamovi/ivweak.u.yaml",
      u_vars().replace("TITLE", "Weak-Instrument-Robust Inference").replace("NAME", "ivweak") + """    - type: CollapseBox
      label: Tests
      collapsed: true
      children:
          - type: TextBox
            name: nullValue
            label: Hypothesised coefficient (H0)
            format: number
          - type: TextBox
            name: level
            label: Confidence level (%)
            format: number
          - type: TextBox
            name: gridWidth
            label: Grid half-width (2SLS standard errors)
            format: number
""" + u_tail([], ["plot"]))

# ---------------------------------------------------------------- ivcmp
write("jamovi/ivcmp.a.yaml",
      head_a("ivcmp", "Compare IV Estimators", "Estimation", "Compare IV Estimators") +
      opt_vars() + """    - name: alpha
      title: Fuller alpha
      type: Number
      default: 1
      min: 0

    - name: level
      title: Confidence level (%)
      type: Number
      default: 95
      min: 50
      max: 99.9

""" + opt_tail([("showOls", "OLS", "true"), ("showLiml", "LIML", "true"), ("showFuller", "Fuller", "true"),
                ("showGmm", "Two-step GMM (robust or cluster errors only)", "false"), ("plot", "Forest plot", "true")]) + "...\n")
write("jamovi/ivcmp.r.yaml",
      r_head("ivcmp", "Compare IV Estimators") +
      table("est", "Estimates of the endogenous regressors",
            [("method", "Estimator", "text"), ("term", "Variable", "text"), ("estimate", "Estimate", "number"), ("se", "Std. error", "number"),
             ("lo", "CI lower", "number"), ("hi", "CI upper", "number"), ("j", "Overid. test", "text")],
            extra_clear=("alpha", "level", "showOls", "showLiml", "showFuller", "showGmm")) +
      """    - name: cmpPlot
      title: Estimates with confidence intervals
      type: Image
      width: 640
      height: 380
      renderFun: .plotCmp
      visible: (plot)
      clearWith:
          - outcome
          - endog
          - instruments
          - controls
          - absorb
          - cluster
          - weights
          - vce
          - small
          - alpha
          - level
          - showOls
          - showLiml
          - showFuller
          - showGmm

...\n""")
write("jamovi/ivcmp.u.yaml",
      u_vars().replace("TITLE", "Compare IV Estimators").replace("NAME", "ivcmp") + """    - type: CollapseBox
      label: Estimators
      collapsed: true
      children:
          - type: CheckBox
            name: showOls
          - type: CheckBox
            name: showLiml
          - type: CheckBox
            name: showFuller
          - type: TextBox
            name: alpha
            label: Fuller alpha
            format: number
            enable: (showFuller)
          - type: CheckBox
            name: showGmm
          - type: TextBox
            name: level
            label: Confidence level (%)
            format: number
""" + u_tail([], ["plot"]))

# ---------------------------------------------------------------- ivplaus
write("jamovi/ivplaus.a.yaml",
      head_a("ivplaus", "Plausibly Exogenous Instruments", "Diagnostics", "Plausibly Exogenous Instruments") +
      opt_vars() + """    - name: delta
      title: Direct effect of each instrument lies in [-delta, delta]
      type: Number
      default: 0
      min: 0

    - name: ltzMean
      title: LTZ prior mean of the direct effects
      type: Number
      default: 0

    - name: ltzSd
      title: LTZ prior standard deviation (0 = no LTZ)
      type: Number
      default: 0
      min: 0

    - name: level
      title: Confidence level (%)
      type: Number
      default: 95
      min: 50
      max: 99.9

""" + opt_tail([("plot", "Plot of the bounds against delta", "true")]) + "...\n")
write("jamovi/ivplaus.r.yaml",
      r_head("ivplaus", "Plausibly Exogenous Instruments") +
      table("tip", "Tipping point: smallest delta whose union interval contains zero",
            [("term", "Endogenous regressor", "text"), ("estimate", "2SLS estimate", "number"), ("tip", "Tipping delta", "number")],
            extra_clear=("level",)) +
      table("uci", "Union of confidence intervals for direct effects in [-delta, delta]",
            [("term", "Endogenous regressor", "text"), ("delta", "delta", "number"), ("lo", "Lower", "number"), ("hi", "Upper", "number"),
             ("zero", "Contains zero", "text")], extra_clear=("delta", "level")) +
      table("ltz", "Local-to-zero approximation",
            [("term", "Endogenous regressor", "text"), ("estimate", "Adjusted estimate", "number"), ("se", "Std. error", "number"),
             ("lo", "CI lower", "number"), ("hi", "CI upper", "number")], visible="(ltzSd)", extra_clear=("ltzMean", "ltzSd", "level")) +
      """    - name: plausPlot
      title: Bounds of the union interval against delta
      type: Image
      width: 640
      height: 400
      renderFun: .plotPlaus
      visible: (plot)
      clearWith:
          - outcome
          - endog
          - instruments
          - controls
          - absorb
          - cluster
          - weights
          - vce
          - small
          - level

...\n""")
write("jamovi/ivplaus.u.yaml",
      u_vars().replace("TITLE", "Plausibly Exogenous Instruments").replace("NAME", "ivplaus") + """    - type: CollapseBox
      label: Direct effects
      collapsed: true
      children:
          - type: TextBox
            name: delta
            label: delta (direct effect in [-delta, delta])
            format: number
          - type: TextBox
            name: ltzMean
            label: LTZ prior mean
            format: number
          - type: TextBox
            name: ltzSd
            label: LTZ prior standard deviation (0 = off)
            format: number
          - type: TextBox
            name: level
            label: Confidence level (%)
            format: number
""" + u_tail([], ["plot"]))
print("written")
