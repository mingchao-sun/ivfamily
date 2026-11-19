# ivfamily 设计文档：工具变量（IV / 2SLS）家族

状态：设计草案（2026-11-15），尚未开工。模块名暂定 `ivfamily`，菜单组 `IV`，界面与帮助文档用英文。
沿用 rdfamily、scmfamily 的模板：一个家族一个模块，每个分析配 Method notes、Data and setup guide、示例数据、NEWS、NOTICE。

## 1. 目标

让 jamovi 用户在一个模块里完成 IV 研究的标准流程：估计，然后检验工具变量的强度和有效性，最后在工具变量弱或外生性存疑时给出稳健推断。
这个流程里每一步都有成熟的 R 包和 Stata 命令，所以验证可以从一开始就带着对照做。

## 2. 分析清单

### 第一批（v1.0，7 个）

| 名称 | 菜单标题 | 内容 | 主要输出 |
|---|---|---|---|
| `iv` | IV Regression | 2SLS、LIML、Fuller(α)、一般 k-class、两步 GMM；内生变量可多个；外生控制变量；可吸收一组或两组固定效应；加权 | 第二阶段系数表（含置信区间）、OLS 对照列、第一阶段系数表、样本量与聚类数 |
| `ivfirst` | First Stage and Weak Instruments | 每个内生变量的第一阶段、偏 R²、各种 F 统计量 | 经典 F、稳健 F（Kleibergen-Paap rk Wald F）、Montiel Olea-Pflueger 有效 F、Sanderson-Windmeijer 条件 F（多内生变量）、Stock-Yogo 临界值（仅同方差）、Cragg-Donald |
| `ivweak` | Weak-Instrument-Robust Inference | 不依赖工具变量强度的推断 | Anderson-Rubin、条件似然比 CLR、Kleibergen K（LM）检验及其置信集（网格反转）；单内生变量时给 tF 调整后的置信区间 |
| `ivover` | Overidentification and Endogeneity | 工具变量外生性与内生性检验 | Sargan、Hansen J、C（差分 Sargan，可指定工具变量子集）、Durbin-Wu-Hausman（稳健版）、对每个工具变量的单独 J 对比 |
| `ivplaus` | Plausibly Exogenous | 排他性约束的敏感性（Conley, Hansen, Rossi 2012） | 允许工具变量直接效应 γ 取值于区间时的置信区间：union of CIs 和 LTZ；γ 从 0 到临界值的曲线 |
| `ivlate` | Binary Instrument and LATE | 二值工具变量、二值内生变量的 LATE 解释 | 简化式（ITT）、第一阶段（遵从者比例）、Wald 比、遵从者/始终接受者/从不接受者比例、协变量下的 Abadie κ 加权、单调性提示 |
| `ivcmp` | Compare IV Estimators | 同一设定下并排比较 OLS、2SLS、LIML、Fuller、GMM、简化式 | 估计量表、森林图、OLS 与 IV 差异 |

### 第二批（v1.x）
- **Shift-share（Bartik）IV**：Borusyak-Hull-Jaravel、Adão-Kolesár-Morales 标准误。
- **多工具变量 / 多重弱**：JIVE、split-sample、Hansen-Hausman-Newey 多工具变量稳健标准误。
- **Examiner / judge 设计**：留一法工具变量。
- **MTE**（边际处理效应）。

## 3. 选项设计（`iv` 为例，其余共享）

| 区域 | 选项 |
|---|---|
| 变量 | Outcome、Endogenous regressors（可多个）、Instruments（可多个）、Exogenous controls（可多个，含因子变量）、Weights、Cluster（一个或两个变量）、Absorb fixed effects（一个或两个变量） |
| 估计 | Estimator（2SLS / LIML / Fuller / k-class / GMM 两步）、Fuller α、k 值、Small-sample correction（是否用 t 与 F 分布） |
| 方差 | Classical / HC1 / HC3 / Cluster / Two-way cluster；GMM 的权重矩阵 |
| 输出 | OLS 对照、第一阶段表、简化式表、置信水平 |
| Output（末尾折叠） | Method notes、Data and setup guide |

要求：每个选项组默认折叠；依赖选项在未勾选时禁用；报错信息指出具体变量和原因（秩不足、工具变量少于内生变量、聚类数过少、共线性）。

## 4. 实现方案

**引擎自写，不直接包装 AER/ivreg。** 理由有三。
1. 本项目已经吃过版本差异的亏：jamovi 随模块打包的包版本和本机不同（rdhte 0.1.0 对 0.2.0，rdrobust 3.0.0 对 4.1.1），并因此出过错误。自写引擎只依赖基础 R 的矩阵运算，行为稳定。
2. LIML、Fuller、k-class、有效 F、KP rk、AR/CLR 都需要拿到中间矩阵，现成包的返回值往往不全。
3. scmfamily 里 sdid、augsynth 也是这样自己移植的，验证方法已经成熟。

结构：
```
ivfamily/
├── DESCRIPTION  LICENSE  NEWS.md  NOTICE.md  README.md
├── jamovi/  0000.yaml  00refs.yaml  <name>.{a,r,u}.yaml
├── R/  engine_iv.R       # 2SLS / k-class / GMM，partialling-out，vcov（classical, HC, cluster, two-way）
│       engine_weak.R     # F 统计量族、AR、CLR、K、tF
│       engine_tests.R    # Sargan/Hansen、C、DWH
│       engine_plaus.R    # Conley-Hansen-Rossi
│       common.R  guide.R  method_notes.R  <name>.b.R
├── data/  示例数据
├── tests/  smoke.R  edge.R  reference.R  plots.R
└── dev/  gen_*.py（批量生成 yaml）、make_data.R、make_help.R
```

固定效应：先写 within 变换（交替投影，一组和两组），不依赖 fixest。高维固定效应（三组以上、百万级观测）不在 v1 范围。

依赖：jmvcore、R6、ggplot2；`Matrix`、`stats` 是 jamovi 自带。只有 `ivfirst` 的 Montiel Olea-Pflueger 有效 F 和 `ivweak` 的 CLR 需要特殊数值例程，自己实现（CLR 的 p 值用 Moreira 的条件分布数值积分）；若实现风险过高，可退而依赖 `ivmodel`，并列为可选。

## 5. 示例数据

| 数据 | 来源 | 用途 |
|---|---|---|
| 香烟需求 | `AER::CigarettesSW`（96 行） | 经典 2SLS，销售税为工具变量，课本复现 |
| 二胎与劳动供给 | `AER::Fertility`（254,654 行） | 二值内生变量、二值工具变量（头两胎同性别），`ivlate` |
| 大学距离 | `AER::CollegeDistance`（4,739 行） | 多控制变量、聚类、稳健 SE |
| 模拟：弱工具变量 | `dev/make_data.R` | 演示 AR/CLR 与 2SLS 的差别（一阶段 F 约 3） |
| 模拟：多工具变量，一个无效 | `dev/make_data.R` | 演示 Hansen J、C 检验、`ivplaus` |

发布前需要确认这些数据的许可证允许随模块再分发（AER 为 GPL，可再分发）。变量名以实际包内为准，建库时核对。

## 6. 验证计划（设计阶段就确定，不留到事后）

对照对象分三层：R 包、Stata 命令、已发表的结果。

| 分析 | R 参考 | Stata 参考 | 已发表结果 |
|---|---|---|---|
| iv：2SLS、稳健/聚类 SE | `AER::ivreg`、`fixest::feols` | `ivreg2`、`ivreghdfe`（均已装） | Stock-Watson 教材中香烟需求的 2SLS 系数 |
| iv：LIML、Fuller、GMM | `ivmodel`（k-class） | `ivreg2 ..., liml / fuller(1) / gmm2s` | — |
| iv：吸收固定效应 | `fixest::feols` | `ivreghdfe` | — |
| ivfirst：KP rk F | — | `ivreg2` 输出的 Kleibergen-Paap rk Wald F | — |
| ivfirst：有效 F | `weakivtest`（需确认可装） | `weakivtest`（SSC，需安装） | Andrews-Stock-Sun (2019) 示例 |
| ivfirst：条件 F | — | `ivreg2, first` 的 SW 统计量 | — |
| ivweak：AR、CLR、K | `ivmodel` | `weakiv`（SSC，需安装） | — |
| ivweak：tF | — | `tf`（SSC，需安装） | Lee et al. (2022) 示例 |
| ivover：Sargan、Hansen、C、DWH | `AER::summary(ivreg, diagnostics = TRUE)` | `ivreg2, endog()`、`overid` | — |
| ivplaus | — | `plausexog`（已装） | Conley et al. (2012) 示例 |
| ivlate | — | `ivreg2`、`kappalate`（SSC，需安装） | Angrist-Evans (1998) 的 Wald 估计 |

验证规则（从本项目吸取的教训）：
1. 所有测试在 **jamovi 自带 R 与随模块打包的库**下再跑一遍，不只跑本机 R。
2. 用界面实测每个分析一遍，按指南里的设定点完。
3. 每个对照脚本固定数据、固定随机数种子，容差写进脚本，结果写入 `PLAN.md`。
4. Stata 脚本里不要在注释中写 `/*`；批处理前先 `capture file close _all`。

## 7. 已知难点与风险

| 风险 | 说明 | 对策 |
|---|---|---|
| KP rk 统计量 | 实现复杂，聚类和异方差下的定义因软件而异 | 先对照 `ivreg2` 的输出，数值不一致时以 Stata 为准并写明差异 |
| 有效 F | Olea-Pflueger 的临界值要数值模拟，可能慢 | 提供解析近似选项，默认用 Nagar 近似 |
| CLR | p 值数值积分，要小心边界 | 与 `weakiv` 对照；失败时降级为 AR 并提示 |
| 多内生变量的弱工具检验 | 方法不统一 | v1 只实现 SW 条件 F，其余明确写进“局限” |
| 一与多聚类 | 聚类数少时 SE 不可靠 | 报告聚类数，少于 30 时给警告 |
| 解释责任 | 用户常把 2SLS 当作 ATE | Method notes 明确写 LATE 的解释条件和排他性假设 |

## 8. 里程碑

| 阶段 | 内容 | 交付物 |
|---|---|---|
| I0 | 脚手架：DESCRIPTION、0000.yaml、`gen_*.py`、示例数据 | 可安装的空模块 |
| I1 | `engine_iv.R` 与 `iv` 分析，对照 ivreg2/ivreghdfe/AER | iv 通过 R 与 Stata 对照 |
| I2 | `ivfirst`、`ivover` | 通过对照 |
| I3 | `ivweak`（AR、CLR、K、tF） | 通过 `weakiv`、`tf` 对照 |
| I4 | `ivplaus`、`ivlate`、`ivcmp` | 通过对照 |
| I5 | 帮助文档、界面实测、打包库下重跑测试 | v1.0 候选版 |

## 9. 需要你确认的决策

1. **自写引擎还是包装 AER/ivreg？**（推荐自写，理由见第 4 节。）
2. **固定效应范围**：v1 只做一组和两组 within，不做高维？（推荐如此。）
3. **v1 是否包含 `ivplaus` 和 `ivlate`？**（推荐包含，两者在应用文献里很常见，Stata 对照也现成。）
4. **Shift-share 放 v1 还是 v1.x？**（推荐 v1.x，它的推断部分要单独验证。）
5. **模块名**：`ivfamily` 还是别的？

## 10. 需要在本机安装的包（确认后再装）

- Stata（SSC）：`weakiv`、`weakivtest`、`tf`、`kappalate`；`ivreg2`、`ivreghdfe`、`ranktest`、`plausexog` 已装。
- R（CRAN）：`ivmodel`（k-class 与 AR/CLR 对照）；`weakivtest`（需先确认包是否存在于 CRAN）。`AER`、`fixest`、`sandwich` 已装。
