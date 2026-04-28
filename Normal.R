# =============================================================================
# 正規分布（線形回帰）のアニメーション
# 回帰直線上を正規分布がスライドする様子を可視化
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

# --- モデルの当てはめ ---
fit <- glm(y ~ x, family = gaussian, data = d0)
a <- coef(fit)[["(Intercept)"]]
b <- coef(fit)[["x"]]
s <- sqrt(summary(fit)$dispersion)

# --- テーマ設定 ---
theme_anim <- theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size = 16))

# --- フレーム生成 ---
xrange <- seq(0, 2, by = 0.05)
tmpdir <- file.path(tempdir(), "frames_normal")
dir.create(tmpdir, showWarnings = FALSE)

for (j in seq_along(xrange)) {
    xj <- xrange[j]
    mu <- a + b * xj

    # 正規分布（左パネル用: 90度回転して表示）
    yseq <- seq(mu - 4 * s, mu + 4 * s, length.out = 200)
    dens <- dnorm(yseq, mean = mu, sd = s)
    poly_df <- data.frame(x = xj - dens * 3, y = yseq)
    pt_df <- data.frame(x = xj, y = mu)

    # 左パネル: 散布図 + 回帰直線 + 分布
    p1 <- ggplot() +
        geom_polygon(data = poly_df, aes(x, y),
            fill = "#F8766D", alpha = 0.8) +
        geom_abline(intercept = a, slope = b,
            color = "green3", linewidth = 1.5) +
        geom_point(data = pt_df, aes(x, y),
            size = 5, color = "green3") +
        geom_point(data = d0, aes(x, y), size = 3) +
        coord_cartesian(xlim = c(-1.5, 2.5), ylim = c(-4, 12)) +
        labs(x = "x", y = "y") + theme_anim

    # 右パネル: 正規分布のPDF
    dist_df <- data.frame(y = yseq, dens = dens)
    dmax <- dnorm(0, 0, s)

    p2 <- ggplot(dist_df, aes(y, dens)) +
        geom_polygon(fill = "#F8766D", alpha = 0.8) +
        coord_cartesian(xlim = c(-4, 12), ylim = c(0, dmax * 1.8)) +
        labs(x = "y", y = "Density") +
        annotate("text", x = 9, y = dmax * 1.6,
            label = "y ~ N(\u03bc, \u03c3\u00b2)\n\u03bc = a + bx",
            size = 4.5) +
        annotate("text", x = 9, y = dmax * 1.3,
            label = "\u2193", size = 6) +
        annotate("text", x = 9, y = dmax * 1.0,
            label = sprintf(
                "y ~ N(\u03bc, %.2f\u00b2)\n\u03bc = %.3f + %.3f x",
                s, a, b),
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
    paste(frames, collapse = " "), "Normal.gif"))

# --- クリーンアップ ---
unlink(tmpdir, recursive = TRUE)
cat("Normal.gif を生成しました\n")
