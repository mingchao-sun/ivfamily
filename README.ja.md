# ivfamily — 操作変数（jamovi モジュール）

操作変数分析：推定量、第 1 段階と弱操作変数の診断、弱操作変数に頑健な推測、過剰識別検定、もっともらしい外生性、LATE、shift-share、限界処置効果。 詳しい英語の説明は [README.md](README.md)、中国語の要約は [README.zh.md](README.zh.md) を参照してください。

## 分析（jamovi メニュー：IV）

| 分析 | 内容 |
|---|---|
| **IV 回帰** | OLS、2SLS、LIML、Fuller、k クラス、2 段階 GMM。各種標準誤差、固定効果の吸収 |
| **IV 推定量の比較** | 各推定量を並べてフォレスト図 |
| **2 値操作変数と LATE** | コンプライアーの割合と LATE |
| **第 1 段階と弱操作変数** | F、条件付き F、Cragg-Donald、Kleibergen-Paap、有効 F |
| **弱操作変数に頑健な推測** | Anderson-Rubin、Kleibergen K、条件付き尤度比検定と信頼集合 |
| **過剰識別と内生性** | Sargan/Hansen J、差の Hansen、内生性検定 |
| **もっともらしい外生の操作変数** | Conley、Hansen、Rossi (2012) |
| **Shift-Share (Bartik) IV** | AKM 露出に頑健な推測 |
| **限界処置効果** | 傾向スコアの多項式による局所 IV（Stata `mtefe` と同仕様）：MTE 曲線、ATE、ATT、ATUT |

各分析の出力には「方法の説明」と「データと設定ガイド」があります。

エンジンは base R で書かれ、R パッケージ（`AER`、`sandwich`、`fixest`）と Stata（`ivreg2`、`ivreghdfe`、`weakiv`、`plausexog`、`mtefe`）と照合済みです。頑健 LIML の標準誤差が `ivreg2` と約 0.07% 異なる以外は、少なくとも 5 桁一致します。サンプルデータは `data/` を参照。

## インストール
- **jamovi ライブラリ**：jamovi 右上のモジュールボタンから *jamovi library* を開き、検索してインストール（モジュール公開後）。
- **サイドロード**：`.jmo` ファイルをダウンロードし、同じメニューの *Install from file* から導入します。`.jmo` はビルドした OS、プロセッサ、jamovi のバージョンでのみ動作します。
- **ソースから**：R と `jmvtools` が必要です：`jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`（パスはお使いの jamovi に合わせてください）。

## 引用とライセンス
引用：Sun, M. *ivfamily* (jamovi module). https://github.com/mingchao-sun/ivfamily 。使用した手法の原論文（出力に記載）も併せて引用してください。

ライセンス：`DESCRIPTION` を参照（GPL）。第三者のコードとデータは `NOTICE` に記載。不具合や提案：https://github.com/mingchao-sun/ivfamily/issues
