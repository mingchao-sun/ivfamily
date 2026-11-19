# ivfamily — 工具变量（jamovi 模块）

工具变量分析：估计量、第一阶段与弱工具变量诊断、弱工具稳健推断、过度识别检验、可能外生的工具变量、LATE、shift-share 设计和边际处理效应。 英文完整说明见 [README.md](README.md)。

## 分析（jamovi 菜单：IV）

| 分析 | 内容 |
|---|---|
| **IV 回归** | OLS、2SLS、LIML、Fuller、k 类和两步 GMM，多种标准误，吸收固定效应 |
| **比较 IV 估计量** | 各估计量并列与森林图 |
| **二值工具与 LATE** | 依从者比例与 LATE |
| **第一阶段与弱工具变量** | F、条件 F、Cragg-Donald、Kleibergen-Paap、有效 F |
| **弱工具稳健推断** | Anderson-Rubin、Kleibergen K、条件似然比检验与置信集 |
| **过度识别与内生性** | Sargan/Hansen J、差分 Hansen、内生性检验 |
| **可能外生的工具变量** | Conley、Hansen 和 Rossi (2012) |
| **Shift-Share (Bartik) IV** | AKM 暴露稳健推断 |
| **边际处理效应** | 局部 IV 与倾向得分多项式（与 Stata `mtefe` 同规格）：MTE 曲线、ATE、ATT、ATUT |

每个分析的输出部分都有"方法说明"和"数据与设置指南"。

引擎为基础 R 编写，已与 R 包（`AER`、`sandwich`、`fixest`）和 Stata（`ivreg2`、`ivreghdfe`、`weakiv`、`plausexog`、`mtefe`）对照，除稳健 LIML 标准误相差约 0.07% 外，均一致到至少五位数字。示例数据见 `data/`。

## 安装
- **jamovi 模块库**：在 jamovi 右上角打开模块按钮，进入 *jamovi library*，搜索并安装（模块上架后）。
- **侧载**：下载 `.jmo` 文件，在同一菜单选 *Install from file*。`.jmo` 只适用于构建它的操作系统、处理器类型和 jamovi 版本。
- **从源码**：需要 R 和 `jmvtools`：`jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`（路径按本机 jamovi 修改）。

## 引用与许可
引用：Sun, M. *ivfamily* (jamovi module). https://github.com/mingchao-sun/ivfamily 。请同时引用所用方法的原始文献（输出中列出）。

许可：见 `DESCRIPTION`（GPL）；第三方代码与数据见 `NOTICE`。问题与建议：https://github.com/mingchao-sun/ivfamily/issues
