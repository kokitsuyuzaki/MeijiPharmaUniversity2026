# =============================================================================
# 一般化線形モデル（GLM）ハンズオン
# =============================================================================
#
# GLMはfamilyを変えるだけで、線形回帰・ロジスティック回帰・ポアソン回帰を
# 統一的な枠組みで扱える。本資料では以下の3つのデータを用いる:
#
#   1. ToothGrowth（線形回帰, family = gaussian）
#      - ビタミンCの用量・投与法と歯の成長の関係
#   2. MASS::birthwt（ロジスティック回帰, family = binomial）
#      - 低出生体重児のリスク因子
#   3. esoph（ポアソン回帰, family = poisson）
#      - 食道癌の症例数とアルコール・タバコ摂取量
#
# いずれも以下の共通パターンで分析する:
#   データ確認 → glm() → summary() → 可視化
# =============================================================================

# パッケージの読み込み
library("MASS")

# =============================================================================
# 1. 線形回帰（family = gaussian）: ToothGrowth
# =============================================================================
#
# ビタミンCの投与がモルモットの歯の成長に与える影響を調べた実験データ
# - len:  歯の長さ（連続値、応答変数）
# - supp: 投与法（OJ=オレンジジュース、VC=アスコルビン酸）
# - dose: 用量（0.5, 1.0, 2.0 mg/day）
#
# 連続値の応答変数 → 正規分布を仮定 → family = gaussian（線形回帰と等価）
# =============================================================================

# データの確認
data(ToothGrowth)
str(ToothGrowth)
head(ToothGrowth)
summary(ToothGrowth)

# 生データの可視化
# 用量が増えるほど歯が長くなる傾向、投与法による差も見られるか確認
options(repr.plot.width=10, repr.plot.height=5)
par(mfrow=c(1,2))

# 箱ひげ図: 用量ごとの歯の長さ
boxplot(len ~ dose, data=ToothGrowth,
    col="lightblue", main="用量と歯の長さ",
    xlab="ビタミンC用量 (mg/day)", ylab="歯の長さ")

# 箱ひげ図: 投与法ごとの歯の長さ
boxplot(len ~ supp, data=ToothGrowth,
    col=c("orange", "lightgreen"), main="投与法と歯の長さ",
    xlab="投与法", ylab="歯の長さ")

par(mfrow=c(1,1))

# モデルの構築
# glm(応答変数 ~ 説明変数, family = 分布族, data = データ)
# gaussian: 正規分布（恒等リンク関数）→ 通常の線形回帰と同じ
fit_gauss <- glm(len ~ dose + supp, family = gaussian, data = ToothGrowth)

# 結果の確認
# summary()の読み方:
#   Coefficients:
#     Estimate:  回帰係数の推定値
#       (Intercept): 切片（dose=0, supp=OJのときの予測値）
#       dose:        用量が1 mg/day増えたときの歯の長さの変化量
#       suppVC:      OJに対するVCの効果の差（OJが基準）
#     Std. Error:  標準誤差（推定値の不確かさ）
#     t value:     t統計量（Estimate / Std. Error）
#     Pr(>|t|):    p値（帰無仮説「係数=0」の検定）
#   Null deviance:    帰無モデル（切片のみ）の逸脱度
#   Residual deviance: 構築したモデルの逸脱度（小さいほど良い）
#   AIC: 赤池情報量基準（モデル比較に使用、小さいほど良い）
summary(fit_gauss)

# 回帰係数の信頼区間
# 95%信頼区間にゼロが含まれなければ、その変数は有意に効果がある
confint(fit_gauss)

# 予測値の可視化
# モデルの予測値と実測値を重ねて表示する
options(repr.plot.width=8, repr.plot.height=6)
dose_seq <- seq(0.5, 2.0, length.out=100)

# 投与法ごとに予測値を計算
pred_OJ <- predict(fit_gauss,
    newdata=data.frame(dose=dose_seq, supp="OJ"), type="response")
pred_VC <- predict(fit_gauss,
    newdata=data.frame(dose=dose_seq, supp="VC"), type="response")

# 実測値のプロット
plot(len ~ dose, data=ToothGrowth,
    col=ifelse(ToothGrowth$supp == "OJ", "orange", "green3"),
    pch=16, cex=1.2,
    main="線形回帰: ビタミンCと歯の成長",
    xlab="ビタミンC用量 (mg/day)", ylab="歯の長さ")

# 予測値の回帰直線
lines(dose_seq, pred_OJ, col="orange", lwd=2)
lines(dose_seq, pred_VC, col="green3", lwd=2)
legend("bottomright", legend=c("OJ（実測）", "VC（実測）",
    "OJ（予測）", "VC（予測）"),
    col=c("orange", "green3", "orange", "green3"),
    pch=c(16, 16, NA, NA), lty=c(NA, NA, 1, 1), lwd=2)

# =============================================================================
# 2. ロジスティック回帰（family = binomial）: birthwt
# =============================================================================
#
# 低出生体重児（2500g未満）のリスク因子を調べたデータ（Baystate Medical Center）
# - low:   低出生体重児かどうか（0=正常, 1=低体重、応答変数）
# - age:   母親の年齢
# - lwt:   最終月経時の母親の体重（ポンド）
# - smoke: 妊娠中の喫煙（0=なし, 1=あり）
# - ht:    高血圧の既往（0=なし, 1=あり）
#
# 二値の応答変数 → ベルヌーイ分布を仮定 → family = binomial（ロジスティック回帰）
# リンク関数はロジット: log(p/(1-p)) = β0 + β1*x1 + ...
# =============================================================================

# データの確認
data(birthwt)
str(birthwt)
head(birthwt)

# 因子型への変換（カテゴリカル変数を明示する）
birthwt$smoke <- factor(birthwt$smoke, labels=c("非喫煙", "喫煙"))
birthwt$ht    <- factor(birthwt$ht, labels=c("なし", "あり"))

summary(birthwt[, c("low", "age", "lwt", "smoke", "ht")])

# 生データの可視化
options(repr.plot.width=10, repr.plot.height=5)
par(mfrow=c(1,2))

# 喫煙と低出生体重の関係（モザイクプロット）
mosaicplot(table(birthwt$smoke, birthwt$low),
    main="喫煙と低出生体重",
    xlab="喫煙", ylab="低出生体重 (0=正常, 1=低体重)",
    col=c("lightblue", "salmon"))

# 母親の年齢と低出生体重の関係
boxplot(age ~ low, data=birthwt,
    col=c("lightblue", "salmon"),
    main="母親の年齢と低出生体重",
    xlab="低出生体重 (0=正常, 1=低体重)", ylab="母親の年齢")

par(mfrow=c(1,1))

# モデルの構築
# binomial: 二項分布（ロジットリンク関数）→ ロジスティック回帰
fit_binom <- glm(low ~ age + lwt + smoke + ht,
    family = binomial, data = birthwt)

# 結果の確認
# summary()の読み方（ロジスティック回帰特有の点）:
#   Coefficients:
#     Estimate: ログオッズ比（対数オッズの変化量）
#       (Intercept): 全説明変数=0（基準カテゴリ）のときの対数オッズ
#       age:     年齢が1歳上がるごとの対数オッズの変化
#       lwt:     体重が1ポンド増えるごとの対数オッズの変化
#       smoke喫煙: 非喫煙に対する喫煙の対数オッズの変化
#       htあり:  高血圧なしに対するありの対数オッズの変化
#     z value:  z統計量（正規近似による検定）
#     Pr(>|z|): p値
#
#   解釈のコツ: Estimateをexp()するとオッズ比になる
#     例: smoke喫煙のEstimate=0.5 → exp(0.5)=1.65
#         → 喫煙者は非喫煙者に比べてリスクが1.65倍
summary(fit_binom)

# オッズ比と信頼区間
# exp(係数) = オッズ比: 1より大→リスク増加、1より小→リスク減少
# 95%信頼区間に1が含まれなければ有意
cbind(
    OddsRatio = exp(coef(fit_binom)),
    exp(confint(fit_binom))
)

# 予測確率の可視化
# 母親の体重を変化させたときの低出生体重の予測確率
options(repr.plot.width=10, repr.plot.height=5)
par(mfrow=c(1,2))

lwt_seq <- seq(min(birthwt$lwt), max(birthwt$lwt), length.out=100)

# 喫煙者 vs 非喫煙者（年齢=平均、高血圧=なし）
pred_nonsmoker <- predict(fit_binom,
    newdata=data.frame(age=mean(birthwt$age), lwt=lwt_seq,
        smoke="非喫煙", ht="なし"),
    type="response")
pred_smoker <- predict(fit_binom,
    newdata=data.frame(age=mean(birthwt$age), lwt=lwt_seq,
        smoke="喫煙", ht="なし"),
    type="response")

plot(lwt_seq, pred_nonsmoker, type="l", col="blue", lwd=2, ylim=c(0, 1),
    main="予測確率: 喫煙の影響",
    xlab="母親の体重 (lb)", ylab="低出生体重の確率")
lines(lwt_seq, pred_smoker, col="red", lwd=2)
abline(h=0.5, lty=2, col="gray")
legend("topright", legend=c("非喫煙", "喫煙"),
    col=c("blue", "red"), lwd=2)

# 高血圧あり vs なし（年齢=平均、非喫煙）
pred_noht <- predict(fit_binom,
    newdata=data.frame(age=mean(birthwt$age), lwt=lwt_seq,
        smoke="非喫煙", ht="なし"),
    type="response")
pred_ht <- predict(fit_binom,
    newdata=data.frame(age=mean(birthwt$age), lwt=lwt_seq,
        smoke="非喫煙", ht="あり"),
    type="response")

plot(lwt_seq, pred_noht, type="l", col="blue", lwd=2, ylim=c(0, 1),
    main="予測確率: 高血圧の影響",
    xlab="母親の体重 (lb)", ylab="低出生体重の確率")
lines(lwt_seq, pred_ht, col="red", lwd=2)
abline(h=0.5, lty=2, col="gray")
legend("topright", legend=c("高血圧なし", "高血圧あり"),
    col=c("blue", "red"), lwd=2)

par(mfrow=c(1,1))

# =============================================================================
# 3. ポアソン回帰（family = poisson）: esoph
# =============================================================================
#
# フランスのイル=エ=ヴィレーヌ県における食道癌の症例対照研究データ
# - ncases:    食道癌の症例数（カウントデータ、応答変数）
# - ncontrols: 対照群の人数
# - agegp:     年齢層（25-34, 35-44, 45-54, 55-64, 65-74, 75+）
# - alcgp:     アルコール摂取量（0-39, 40-79, 80-119, 120+ g/day）
# - tobgp:     タバコ摂取量（0-9, 10-19, 20-29, 30+ g/day）
#
# カウントデータ → ポアソン分布を仮定 → family = poisson
# リンク関数は対数: log(μ) = β0 + β1*x1 + ...
# =============================================================================

# データの確認
data(esoph)
str(esoph)
head(esoph, 10)
summary(esoph)

# 生データの可視化
# 各要因ごとの症例数の分布を確認する
options(repr.plot.width=12, repr.plot.height=5)
par(mfrow=c(1,3))

# 年齢層ごとの症例数
boxplot(ncases ~ agegp, data=esoph, col="lightblue",
    main="年齢層と症例数",
    xlab="年齢層", ylab="症例数")

# アルコール摂取量ごとの症例数
boxplot(ncases ~ alcgp, data=esoph, col="lightyellow",
    main="アルコール摂取量と症例数",
    xlab="アルコール摂取量 (g/day)", ylab="症例数")

# タバコ摂取量ごとの症例数
boxplot(ncases ~ tobgp, data=esoph, col="lightgreen",
    main="タバコ摂取量と症例数",
    xlab="タバコ摂取量 (g/day)", ylab="症例数")

par(mfrow=c(1,1))

# モデルの構築
# poisson: ポアソン分布（対数リンク関数）→ ポアソン回帰
# 順序因子を数値に変換して使用（線形な用量反応関係を仮定）
esoph$alcgp_num <- as.numeric(esoph$alcgp)
esoph$tobgp_num <- as.numeric(esoph$tobgp)
esoph$agegp_num <- as.numeric(esoph$agegp)

fit_pois <- glm(ncases ~ agegp_num + alcgp_num + tobgp_num,
    family = poisson, data = esoph)

# 結果の確認
# summary()の読み方（ポアソン回帰特有の点）:
#   Coefficients:
#     Estimate: 対数リスク比（対数スケールでの変化量）
#       (Intercept): 全説明変数=0のときの対数期待値
#       agegp_num:   年齢層が1段階上がるごとの対数期待値の変化
#       alcgp_num:   アルコール摂取が1段階増えるごとの対数期待値の変化
#       tobgp_num:   タバコ摂取が1段階増えるごとの対数期待値の変化
#
#   解釈のコツ: Estimateをexp()するとリスク比（rate ratio）になる
#     例: alcgp_numのEstimate=0.5 → exp(0.5)=1.65
#         → アルコール摂取が1段階増えると期待症例数が1.65倍
#
#   過分散の確認:
#     Residual deviance / 残差自由度 ≈ 1 なら問題なし
#     1より大幅に大きい場合は過分散 → 負の二項分布などを検討
summary(fit_pois)

# 過分散の確認
# Residual deviance / df が1より大幅に大きければ過分散の可能性
deviance(fit_pois) / fit_pois$df.residual

# リスク比と信頼区間
# exp(係数) = リスク比: 1より大→リスク増加、1より小→リスク減少
cbind(
    RateRatio = exp(coef(fit_pois)),
    exp(confint(fit_pois))
)

# 予測値の可視化
options(repr.plot.width=10, repr.plot.height=5)
par(mfrow=c(1,2))

# アルコール摂取量ごとの予測症例数（年齢層・タバコ=中央値で固定）
alcgp_seq <- seq(1, 4, length.out=100)
pred_alc <- predict(fit_pois,
    newdata=data.frame(
        agegp_num=median(esoph$agegp_num),
        alcgp_num=alcgp_seq,
        tobgp_num=median(esoph$tobgp_num)),
    type="response")

plot(alcgp_seq, pred_alc, type="l", col="darkred", lwd=2,
    main="予測症例数: アルコールの影響",
    xlab="アルコール摂取量（カテゴリ）", ylab="予測症例数",
    xaxt="n")
axis(1, at=1:4, labels=levels(esoph$alcgp))

# 年齢層ごとの予測症例数（アルコール・タバコ=中央値で固定）
agegp_seq <- seq(1, 6, length.out=100)
pred_age <- predict(fit_pois,
    newdata=data.frame(
        agegp_num=agegp_seq,
        alcgp_num=median(esoph$alcgp_num),
        tobgp_num=median(esoph$tobgp_num)),
    type="response")

plot(agegp_seq, pred_age, type="l", col="darkblue", lwd=2,
    main="予測症例数: 年齢の影響",
    xlab="年齢層（カテゴリ）", ylab="予測症例数",
    xaxt="n")
axis(1, at=1:6, labels=levels(esoph$agegp))

par(mfrow=c(1,1))

# =============================================================================
# 3つのモデルの比較まとめ
# =============================================================================
#
# | モデル           | family    | リンク関数 | 応答変数   | 係数の解釈       |
# |------------------|-----------|------------|------------|------------------|
# | 線形回帰         | gaussian  | 恒等       | 連続値     | そのまま差分     |
# | ロジスティック回帰| binomial  | ロジット   | 二値(0/1)  | exp()→オッズ比   |
# | ポアソン回帰     | poisson   | 対数       | カウント   | exp()→リスク比   |
#
# 共通点:
#   - glm() の family を変えるだけで、同じ枠組みで異なる種類のデータを扱える
#   - summary() の構造は同じ（Estimate, Std. Error, 検定統計量, p値）
#   - confint() で信頼区間、predict(type="response") で予測値を得る
#
# 注意点:
#   - gaussian以外では、Estimateはリンク関数のスケールで出力される
#     → 元のスケールに戻すにはexp()（対数・ロジットリンクの場合）が必要
#   - ポアソン回帰では過分散に注意（Residual deviance / df >> 1）
#   - ロジスティック回帰ではオッズ比の信頼区間に1が含まれるか確認
# =============================================================================

# セッション情報
sessionInfo()
