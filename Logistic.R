# =============================================================================
# 二項分布（ロジスティック回帰）のアニメーション
# 回帰曲線上を二項分布がスライドする様子を可視化
# =============================================================================

library(ggplot2)
library(grid)

# --- データ取得 ---
# 久保(2012)「データ解析のための統計モデリング入門」第3章のデータ
# x: 体サイズ（連続値）、y: 種子数（カウントデータ）
tmpzip <- tempfile(fileext = ".zip")
download.file(
    "https://kuboweb.github.io/-kubo/stat/iwanamibook/kubobook_2012.zip",
    tmpzip, quiet = TRUE)
unzip(tmpzip, "kubobook_2012/poisson/d0.RData", exdir = tempdir())
load(file.path(tempdir(), "kubobook_2012/poisson/d0.RData"))
unlink(tmpzip)
N_trial <- max(d0$y)  # 試行数（二項分布のN）：データの最大値を使用

# --- モデルの当てはめ ---
fit <- glm(cbind(y, N_trial - y) ~ x,
    family = binomial(link = "logit"), data = d0)
a <- coef(fit)[["(Intercept)"]]
b <- coef(fit)[["x"]]

# --- テーマ設定 ---
theme_anim <- theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size = 16))

# --- フレーム生成 ---
xrange <- seq(0, 2, by = 0.05)
tmpdir <- file.path(tempdir(), "frames_logistic")
dir.create(tmpdir, showWarnings = FALSE)

for (j in seq_along(xrange)) {
    xj <- xrange[j]
    prob <- plogis(a + b * xj)

    # 二項分布のPMF
    kvec <- 0:N_trial
    pmf <- dbinom(kvec, N_trial, prob)
    pt_df <- data.frame(x = xj, y = prob * N_trial)
    bar_df <- data.frame(
        k = kvec, pmf = pmf,
        xmin = xj - pmf * 3, xmax = xj,
        ymin = kvec - 0.3, ymax = kvec + 0.3
    )

    # 左パネル: 散布図 + ロジスティック曲線（期待値スケール）+ 横向き棒グラフ
    p1 <- ggplot() +
        geom_rect(data = bar_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#F8766D", alpha = 0.8) +
        stat_function(
            fun = function(x) N_trial * plogis(a + b * x),
            color = "green3", linewidth = 1.5) +
        geom_point(data = pt_df, aes(x, y),
            size = 5, color = "green3") +
        geom_point(data = d0, aes(x, y), size = 3) +
        coord_cartesian(xlim = c(-1.5, 2.5), ylim = c(-2, 10)) +
        labs(x = "x", y = "y") + theme_anim

    # 右パネル: 二項分布のPMF（縦向き棒グラフ）
    p2 <- ggplot(bar_df, aes(x = k, y = pmf)) +
        geom_col(fill = "#F8766D", alpha = 0.8, width = 0.6) +
        coord_cartesian(xlim = c(-1, N_trial + 1), ylim = c(0, 1)) +
        labs(x = "y", y = "Probability") +
        annotate("text", x = N_trial * 0.8, y = 0.9,
            label = sprintf(
                "y ~ Bin(%d, q)\nlogit(q) = a + bx", N_trial),
            size = 4.5) +
        annotate("text", x = N_trial * 0.8, y = 0.72,
            label = "\u2193", size = 6) +
        annotate("text", x = N_trial * 0.8, y = 0.58,
            label = sprintf(
                "y ~ Bin(%d, q)\nlogit(q) = %.3f + %.3f x",
                N_trial, a, b),
            size = 4.5) +
        theme_anim

    png(file.path(tmpdir, sprintf("frame_%03d.png", j)),
        width = 960, height = 480)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1, 2)))
    print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
    dev.off()
}

# --- GIF生成 ---
frames <- file.path(tmpdir, sprintf("frame_%03d.png", seq_along(xrange)))
system(paste("convert -delay 10 -loop 0",
    paste(frames, collapse = " "), "Logistic.gif"))

# --- クリーンアップ ---
unlink(tmpdir, recursive = TRUE)
cat("Logistic.gif を生成しました\n")
