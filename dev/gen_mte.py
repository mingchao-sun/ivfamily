# Generates ivmte.a/.r/.u.yaml (Marginal Treatment Effects). Run from the module root:  python dev/gen_mte.py
NL = chr(10)


def q(s):
    return "'" + s.replace("'", "''") + "'"


def var(name, title, multi=False, perm=("numeric",), sug=None):
    s = f"    - name: {name}{NL}      title: {q(title)}{NL}      type: {'Variables' if multi else 'Variable'}{NL}"
    if sug:
        s += f"      suggested:{NL}" + "".join(f"          - {x}{NL}" for x in sug)
    s += f"      permitted:{NL}" + "".join(f"          - {x}{NL}" for x in perm) + NL
    return s


def num(name, title, default, mn=None, mx=None, typ="Number"):
    s = f"    - name: {name}{NL}      title: {q(title)}{NL}      type: {typ}{NL}      default: {default}{NL}"
    if mn is not None:
        s += f"      min: {mn}{NL}"
    if mx is not None:
        s += f"      max: {mx}{NL}"
    return s + NL


def lst(name, title, opts, default):
    s = f"    - name: {name}{NL}      title: {q(title)}{NL}      type: List{NL}      options:{NL}"
    for n, t in opts:
        s += f"          - name: {n}{NL}            title: {q(t)}{NL}"
    return s + f"      default: {default}{NL}{NL}"


def boo(name, title, default="false"):
    return f"    - name: {name}{NL}      title: {q(title)}{NL}      type: Bool{NL}      default: {default}{NL}{NL}"


OPTS = "outcome, treat, instruments, covs, cluster, degree, link, nboot, level, seed"

a = (f"---{NL}name: ivmte{NL}title: 'Marginal Treatment Effects'{NL}jas: '1.2'{NL}version: '1.0.0'{NL}menuGroup: IV{NL}menuSubgroup: Estimation{NL}"
     f"menuTitle: 'Marginal Treatment Effects'{NL}menuSubtitle: 'Heckman and Vytlacil: effects by resistance to treatment, ATE, ATT, ATUT'{NL}{NL}options:{NL}    - name: data{NL}      type: Data{NL}{NL}")
a += var("outcome", "Outcome") + var("treat", "Binary treatment (0/1)", False, ("numeric", "factor", "id"), ("nominal", "continuous"))
a += var("instruments", "Instruments (excluded from the outcome)", True, ("numeric", "factor"), ("continuous", "nominal"))
a += var("covs", "Covariates (optional)", True, ("numeric", "factor"), ("continuous", "nominal")) + var("cluster", "Cluster (optional)", False, ("numeric", "factor", "id"), ("nominal",))
a += num("degree", "Degree of the polynomial in the marginal treatment effect", 2, 1, 4, "Integer")
a += lst("link", "First stage (propensity score)", [("probit", "Probit"), ("logit", "Logit")], "probit")
a += num("nboot", "Bootstrap samples", 300, 50, 5000, "Integer") + num("level", "Confidence level (%)", 95, 50, 99.9) + num("seed", "Random seed", 1234, 1, 999999, "Integer")
a += boo("plot", "Marginal treatment effect plot", "true") + boo("showHelp", "Method notes") + boo("showGuide", "Data and setup guide") + "..." + NL


def col(name, title, typ, extra=""):
    return f"        - name: {name}{NL}          title: {q(title)}{NL}          type: {typ}{NL}" + (extra + NL if extra else "")


def tbl(name, title, cols, visible=None):
    s = f"    - name: {name}{NL}      title: {q(title)}{NL}      type: Table{NL}      rows: 0{NL}"
    if visible:
        s += f"      visible: {visible}{NL}"
    return s + f"      clearWith: [{OPTS}]{NL}      columns:{NL}" + "".join(cols) + NL


PC = "          format: zto,pvalue"
r = (f"---{NL}name: ivmte{NL}title: 'Marginal Treatment Effects'{NL}jrs: '1.1'{NL}{NL}items:{NL}"
     f"    - name: guide{NL}      title: Data and setup guide{NL}      type: Html{NL}      visible: (showGuide){NL}{NL}"
     f"    - name: helpText{NL}      title: Method notes{NL}      type: Html{NL}      visible: (showHelp){NL}{NL}"
     f"    - name: notes{NL}      title: Notes{NL}      type: Html{NL}      visible: true{NL}{NL}")
r += tbl("support", "Propensity score and common support", [col("group", "", "text"), col("n", "N", "integer"), col("mean", "Mean", "number"), col("min", "Lowest", "number"), col("max", "Highest", "number")])
r += tbl("pars", "Treatment effect parameters", [col("what", "", "text"), col("est", "Estimate", "number"), col("se", "Std. error", "number"), col("lo", "CI lower", "number"), col("hi", "CI upper", "number"), col("p", "p", "number", PC)])
r += tbl("tests", "Tests of heterogeneity", [col("what", "", "text"), col("chisq", "Chi-square", "number"), col("df", "df", "integer"), col("p", "p", "number", PC)])
r += tbl("curve", "Marginal treatment effect by resistance to treatment", [col("u", "Resistance u", "number"), col("mte", "MTE", "number"), col("se", "Std. error", "number"), col("lo", "CI lower", "number"), col("hi", "CI upper", "number")])
r += (f"    - name: mtePlot{NL}      title: 'Marginal treatment effect'{NL}      type: Image{NL}      width: 560{NL}      height: 380{NL}"
      f"      renderFun: .plotMte{NL}      visible: (plot){NL}      clearWith: [{OPTS}]{NL}{NL}")
r += (f"    - name: psPlot{NL}      title: 'Propensity score by treatment'{NL}      type: Image{NL}      width: 560{NL}      height: 340{NL}"
      f"      renderFun: .plotPs{NL}      visible: (plot){NL}      clearWith: [{OPTS}]{NL}{NL}")
r += f"refs:{NL}    - heckmanvytlacil2005{NL}    - andresen2018{NL}...{NL}"


def box(label, name, single=True):
    s = f"          - type: TargetLayoutBox{NL}            label: {q(label)}{NL}            children:{NL}                - type: VariablesListBox{NL}                  name: {name}{NL}                  isTarget: true{NL}                  height: small{NL}"
    if single:
        s += f"                  maxItemCount: 1{NL}"
    return s


def tb(name, label):
    return f"          - type: TextBox{NL}            name: {name}{NL}            label: {q(label)}{NL}            format: number{NL}"


def cb(name, label=None):
    return f"          - type: CheckBox{NL}            name: {name}{NL}" + (f"            label: {q(label)}{NL}" if label else "")


def cmb(name):
    return f"          - type: ComboBox{NL}            name: {name}{NL}"


def collapse(label, body, collapsed=True):
    return f"    - type: CollapseBox{NL}      label: {q(label)}{NL}      collapsed: {str(collapsed).lower()}{NL}      children:{NL}" + body


u = (f"title: 'Marginal Treatment Effects'{NL}name: ivmte{NL}jus: '3.0'{NL}stage: 0{NL}compilerMode: tame{NL}children:{NL}"
     f"    - type: VariableSupplier{NL}      persistentItems: false{NL}      stretchFactor: 1{NL}      children:{NL}")
u += box("Outcome", "outcome") + box("Binary treatment (0/1)", "treat") + box("Instruments (excluded from the outcome)", "instruments", False) + box("Covariates (optional)", "covs", False) + box("Cluster (optional)", "cluster")
u += collapse("Model", tb("degree", "Degree of the polynomial in the marginal treatment effect") + cmb("link"), False)
u += collapse("Inference", tb("nboot", "Bootstrap samples") + tb("level", "Confidence level (%)") + tb("seed", "Random seed") + cb("plot"))
u += collapse("Output", cb("showHelp", "Method notes") + cb("showGuide", "Data and setup guide"))

for ext, t in (("a", a), ("r", r), ("u", u)):
    with open(f"jamovi/ivmte.{ext}.yaml", "w", encoding="utf-8", newline="") as f:
        f.write(t)
print("written")
