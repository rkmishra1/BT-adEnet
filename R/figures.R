## Figures for bayes-tadenet.tex.  Rscript R/figures.R
##
## Palette: Okabe-Ito subset, colourblind-safe, contrast >= 3:1 on white.  Series
## are also distinguished by line type and plotting symbol, so identity never
## depends on colour alone.  Annotations are computed from the data.
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
out <- file.path(here, "..")

BLUE <- "#0072B2"
ORANGE <- "#D55E00"
GREEN <- "#009E73"
PURPLE <- "#CC79A7"
GRID <- "#E6E6E6"
AXIS <- "#4D4D4D"
logax <- function(side, at, zero_first = FALSE) {
  lab <- ifelse(at < 1000, format(at), paste0("10^", round(log10(at))))
  if (zero_first) lab[1] <- "0"
  axis(side, at = at, labels = parse(text = lab))
}
pan <- function() {
  par(mar = c(3.6, 3.8, 2.0, 0.8), mgp = c(2.3, 0.6, 0), tcl = -0.25,
      cex.axis = 0.85, cex.lab = 0.95, cex.main = 0.95, font.main = 1,
      col.axis = AXIS, col.lab = "black", fg = AXIS, las = 1)
}
grid_h <- function() grid(NA, NULL, col = GRID, lty = 1, lwd = 0.6)
have <- function(f) file.exists(file.path(out, f))

## ============================================ bounded sensitivity (E2)
if (have("e2-sensitivity.rds")) {
  e2 <- readRDS(file.path(out, "e2-sensitivity.rds"))$agg
  tk <- e2[e2$loss == "Tukey", ]
  tk <- tk[order(tk$t), ]
  gs <- e2[e2$loss == "Gaussian", ]
  gs <- gs[order(gs$t), ]
  xs <- function(v) pmax(v, 0.5)
  ser <- function(d, v, col, pch, lty) {
    lines(xs(d$t), d[[v]], col = col, lwd = 2, lty = lty)
    points(xs(d$t), d[[v]], col = col, pch = pch, bg = "white", cex = 0.9)
  }
  tick <- c(0.5, 1, 10, 100, 1000, 1e5)
  pdf(file.path(out, "fig-sensitivity.pdf"), width = 9.4, height = 2.9, pointsize = 10)
  par(mfrow = c(1, 3))
  pan()
  ## the squared-error displacement spans five orders of magnitude, so this panel
  ## is logarithmic in both axes and starts at t = 1, where the displacement is > 0
  tk1 <- tk[tk$t > 0, ]
  gs1 <- gs[gs$t > 0, ]
  plot(tk1$t, tk1$shift, type = "n", log = "xy", ylim = range(c(tk1$shift, gs1$shift)) * c(0.5, 2),
       xaxt = "n", xlab = "displacement t of a single observation",
       ylab = expression(group("||", Delta*hat(beta), "||")[2]), main = "(a) posterior mean")
  logax(1, c(1, 10, 100, 1000, 1e5))
  grid_h()
  ser(tk1, "shift", BLUE, 21, 1)
  ser(gs1, "shift", ORANGE, 24, 2)
  legend("topleft", c("Tukey", "squared error"), col = c(BLUE, ORANGE),
         pch = c(21, 24), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85, pt.bg = "white")
  plot(xs(tk$t), tk$sigma, type = "n", log = "xy", xaxt = "n",
       ylim = range(c(tk$sigma, gs$sigma)) * c(0.8, 1.6),
       xlab = "displacement t", ylab = expression(hat(sigma)), main = "(b) scale estimate")
  logax(1, tick, zero_first = TRUE)
  grid_h()
  ser(tk, "sigma", BLUE, 21, 1)
  ser(gs, "sigma", ORANGE, 24, 2)
  plot(xs(tk$t), tk$size, type = "n", log = "x", xaxt = "n",
       ylim = c(-0.4, max(c(tk$size, gs$size)) + 1),
       xlab = "displacement t", ylab = "variables selected", main = "(c) selected model size")
  logax(1, tick, zero_first = TRUE)
  grid_h()
  ser(tk, "size", BLUE, 21, 1)
  ser(gs, "size", ORANGE, 24, 2)
  dev.off()
}

## ======================================= learning rate: coverage, width, risk
if (have("calib2.rds")) {
  cb <- readRDS(file.path(out, "calib2.rds"))
  a1 <- cb$agg_sweep[order(cb$agg_sweep$w), ]
  wcal <- cb$agg_meth$w[cb$agg_meth$method == "BT-AdEnet"]
  w_risk <- a1$w[which.min(a1$excess)]
  vline <- function(v, col, lty) if (length(v) && is.finite(v)) abline(v = v, col = col, lty = lty)
  pdf(file.path(out, "fig-learning-rate.pdf"), width = 9.4, height = 2.9, pointsize = 10)
  par(mfrow = c(1, 3))
  pan()
  plot(a1$w, a1$cov, type = "n", log = "x",
       ylim = range(c(0.7, 1, a1$cov + 2 * a1$cov_se, a1$cov - 2 * a1$cov_se), na.rm = TRUE),
       xlab = "learning rate w", ylab = "marginal coverage", main = "(a) calibration")
  grid_h()
  abline(h = 0.95, col = AXIS, lty = 3)
  segments(a1$w, a1$cov - 2 * a1$cov_se, a1$w, a1$cov + 2 * a1$cov_se, col = BLUE)
  lines(a1$w, a1$cov, col = BLUE, lwd = 2)
  points(a1$w, a1$cov, col = BLUE, pch = 21, bg = "white")
  vline(cb$cross, GREEN, 2)
  vline(wcal, ORANGE, 4)
  legend("bottomleft", c("coverage, 2 s.e.", sprintf("crosses nominal (%.2f)", cb$cross),
                         sprintf("calibrated w (%.2f)", wcal)),
         col = c(BLUE, GREEN, ORANGE), lty = c(1, 2, 4), lwd = 2, pch = c(21, NA, NA),
         bty = "n", cex = 0.75, pt.bg = "white")
  rel <- a1$width / a1$width[a1$w == 1]
  plot(a1$w, rel, type = "n", log = "xy", xlab = "learning rate w",
       ylab = "interval width (relative to w = 1)", main = "(b) width scaling")
  grid_h()
  lines(a1$w, 1 / sqrt(a1$w), col = ORANGE, lwd = 5)
  lines(a1$w, rel, col = BLUE, lwd = 1.6)
  points(a1$w, rel, col = BLUE, pch = 21, bg = "white", cex = 0.9)
  legend("topright", c(expression(paste("predicted  ", w^{-1/2})), "observed"),
         col = c(ORANGE, BLUE), lty = 1, lwd = c(5, 1.6), pch = c(NA, 21),
         bty = "n", cex = 0.8, pt.bg = "white")
  plot(a1$w, a1$excess, type = "n", log = "x", xlab = "learning rate w",
       ylim = range(c(a1$excess + 2 * a1$excess_se, a1$excess - 2 * a1$excess_se), na.rm = TRUE),
       ylab = "excess Tukey risk", main = "(c) excess risk")
  grid_h()
  segments(a1$w, a1$excess - 2 * a1$excess_se, a1$w, a1$excess + 2 * a1$excess_se, col = BLUE)
  lines(a1$w, a1$excess, col = BLUE, lwd = 2)
  points(a1$w, a1$excess, col = BLUE, pch = 21, bg = "white")
  vline(w_risk, GREEN, 2)
  vline(cb$cross, ORANGE, 4)
  legend("topright", c(sprintf("risk-optimal (%.2f)", w_risk),
                       sprintf("coverage-optimal (%.2f)", cb$cross)),
         col = c(GREEN, ORANGE), lty = c(2, 4), lwd = 2, bty = "n", cex = 0.8)
  dev.off()
}

## ================================== leakage against p and against separation
if (have("ratio.rds")) {
  r <- readRDS(file.path(out, "ratio.rds"))
  r$Q <- (r$p - 6) / (r$prior_leak * r$w * r$lam1) / r$wmaxA
  r$cfg <- ifelse(r$floor == "p", "pfloor", paste0("iota", r$iota))
  md <- aggregate(cbind(leak, Q) ~ cfg + p, r, stats::median)
  sty <- list(iota0 = list(ORANGE, 24, 2, "uniform weights"),
              iota1 = list(BLUE, 21, 1, "adaptive, iota = 1"),
              iota2 = list(PURPLE, 23, 3, "adaptive, iota = 2"),
              pfloor = list(GREEN, 22, 1, "iota = 1, p-scaled floor"))
  gl <- function(d, s, xv) {
    d <- d[order(d$p), ]
    lines(d[[xv]], d$leak, col = s[[1]], lwd = 2, lty = s[[3]])
    points(d[[xv]], d$leak, col = s[[1]], pch = s[[2]], bg = "white", cex = 0.9)
  }
  pdf(file.path(out, "fig-leakage.pdf"), width = 6.6, height = 3.0, pointsize = 10)
  par(mfrow = c(1, 2))
  pan()
  plot(md$p, md$leak, type = "n", log = "xy", xlab = "dimension p",
       ylab = "posterior leakage", main = "(a) against dimension")
  grid_h()
  for (k in names(sty)) gl(md[md$cfg == k, ], sty[[k]], "p")
  legend("topleft", sapply(sty, `[[`, 4), col = sapply(sty, `[[`, 1), pch = sapply(sty, `[[`, 2),
         lty = sapply(sty, `[[`, 3), lwd = 2, bty = "n", cex = 0.7, pt.bg = "white")
  plot(md$Q, md$leak, type = "n", log = "xy", xlab = "separation ratio Q",
       ylab = "posterior leakage", main = "(b) against separation")
  grid_h()
  for (k in names(sty)) gl(md[md$cfg == k, ], sty[[k]], "Q")
  dev.off()
}

## ================================ real data, per-split error ratios
if (have("realdata.rds")) {
  rd <- readRDS(file.path(out, "realdata.rds"))$res
  ds <- c("Boston housing", "Ames housing", "Glass", "Riboflavin")
  short <- c("Boston", "Ames", "Glass", "Riboflavin")
  cols <- c(GREEN, BLUE, ORANGE, PURPLE)
  ratio_for <- function(metric, other) {
    w <- reshape(rd[, c("dataset", "arm", "rep", metric)], idvar = c("dataset", "rep"),
                 timevar = "arm", direction = "wide")
    data.frame(dataset = w$dataset,
               r = w[[paste0(metric, ".BT-AdEnet")]] / w[[paste0(metric, ".", other)]])
  }
  lim <- c(0.1, 10)
  panel <- function(z, main, ylab) {
    plot(NA, xlim = c(0.5, length(ds) + 0.5), ylim = lim, log = "y", xaxt = "n",
         xlab = "", ylab = ylab, main = main)
    axis(1, at = seq_along(ds), labels = short)
    grid(NA, NULL, col = GRID, lty = 1, lwd = 0.6)
    abline(h = 1, col = AXIS, lty = 2)
    for (k in seq_along(ds)) {
      v <- z$r[z$dataset == ds[k]]
      v <- v[is.finite(v)]
      if (!length(v)) next
      x <- k + 0.5 * (rank(v, ties.method = "first") / (length(v) + 1) - 0.5)
      clip <- v > lim[2] | v < lim[1]
      points(x, pmin(pmax(v, lim[1]), lim[2]), pch = ifelse(clip, 24, 21), col = cols[k],
             bg = "white", cex = 0.8)
      text(k, lim[2], sprintf("%.0f%% < 1", 100 * mean(v < 1)), cex = 0.72, pos = 1)
    }
  }
  pdf(file.path(out, "fig-realdata.pdf"), width = 9.4, height = 3.2, pointsize = 10)
  par(mfrow = c(1, 3))
  pan()
  panel(ratio_for("rtmspe", "Gaussian"), "(a) trimmed, against squared error", "error ratio, BT-AdEnet / comparator")
  panel(ratio_for("rmspe", "Gaussian"), "(b) untrimmed, against squared error", "error ratio")
  panel(ratio_for("rtmspe", "Tukey+SSVS"), "(c) trimmed, against spike-and-slab", "error ratio")
  dev.off()
}

## =============================== losses and priors at p = 100
if (have("comparators.rds")) {
  r <- readRDS(file.path(out, "comparators.rds"))
  a <- aggregate(cbind(me, leak) ~ arm + scenario + p, r, stats::median)
  a$loss <- sub(" .*", "", a$arm)
  a$prior <- sub(".* \\+ ", "", a$arm)
  pdf(file.path(out, "fig-comparators.pdf"), width = 9.4, height = 3.1, pointsize = 10)
  par(mfrow = c(1, 3))
  pan()
  lc <- c(Tukey = BLUE, Huber = ORANGE, Gauss = "#666666")
  pp <- c(AdEnet = 21, horseshoe = 24, SSVS = 22)
  ttl <- c(clean = "clean", vertical = "vertical outliers", both = "bad leverage")
  for (sc in names(ttl)) {
    z <- a[a$scenario == sc & a$p == 100, ]
    plot(z$leak, z$me, type = "n", log = "xy", xlab = "prior leakage", ylab = "model error",
         main = paste0("(", match(sc, names(ttl)), ") ", ttl[sc]))
    grid(col = GRID, lty = 1, lwd = 0.6)
    points(z$leak, z$me, col = lc[z$loss], pch = pp[z$prior], bg = "white", cex = 1.2, lwd = 1.6)
    if (sc == "clean") {
      legend("topleft", names(lc), col = lc, pch = 19, bty = "n", cex = 0.78)
      legend("bottomright", names(pp), pch = pp, col = AXIS, bty = "n", cex = 0.78, pt.bg = "white")
    }
  }
  dev.off()
}
## ================================ ingredients of the generalized posterior
## Analytic: loss, influence, prior, and exact grid posteriors of a location
## parameter, so this figure reads no .rds.
source(file.path(here, "bt_adenet.R"))
local({
  d <- 4.685
  u <- seq(-8, 8, length.out = 801)
  pdf(file.path(out, "fig-gibbs.pdf"), width = 7.4, height = 5.0, pointsize = 10)
  par(mfrow = c(2, 2))
  pan()
  plot(u, u^2 / 2, type = "n", ylim = c(0, 6), xlab = "standardized residual u",
       ylab = "loss", main = "(a) loss")
  grid_h()
  abline(v = c(-d, d), col = AXIS, lty = 3)
  lines(u, u^2 / 2, col = ORANGE, lwd = 2, lty = 2)
  lines(u, tukey_rho(u, d), col = BLUE, lwd = 2)
  legend("bottomright", c("Tukey biweight", "squared error"), col = c(BLUE, ORANGE),
         lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
  plot(u, u, type = "n", ylim = c(-3, 3), xlab = "standardized residual u",
       ylab = "derivative of the loss", main = "(b) influence")
  grid_h()
  abline(h = 0, v = c(-d, d), col = AXIS, lty = 3)
  lines(u, u, col = ORANGE, lwd = 2, lty = 2)
  lines(u, tukey_psi(u, d), col = BLUE, lwd = 2)
  legend("topleft", c(expression(psi[d](u)), "u"), col = c(BLUE, ORANGE),
         lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
  b <- seq(-4, 4, length.out = 801)
  dens <- function(a, c) {
    f <- exp(-a * abs(b) - c * b^2 / 2)
    f / sum(f * diff(b[1:2]))
  }
  plot(b, dens(5, 0.5), type = "n", log = "y", ylim = c(1e-4, 5), yaxt = "n",
       xlab = expression(beta[j]), ylab = "prior density", main = "(c) adaptive elastic-net prior")
  logax(2, c(1e-4, 1e-3, 1e-2, 0.1, 1))
  grid_h()
  lines(b, dens(0.5, 0.5), col = BLUE, lwd = 2)
  lines(b, dens(5, 0.5), col = ORANGE, lwd = 2, lty = 2)
  legend("topright", c(expression(paste("small weight, ", a[j] == 0.5)),
                     expression(paste("large weight, ", a[j] == 5))),
         col = c(BLUE, ORANGE), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
  set.seed(3)
  y <- c(stats::rnorm(36), rep(10, 4))
  m <- seq(-1.5, 3, length.out = 901)
  post <- function(w, loss) {
    lp <- -w * sapply(m, function(mu) sum(loss(y - mu)))
    f <- exp(lp - max(lp))
    f / sum(f * diff(m[1:2]))
  }
  pt_lo <- post(0.5, function(r) tukey_rho(r, d))
  pt_hi <- post(2, function(r) tukey_rho(r, d))
  pg <- post(1, function(r) r^2 / 2)
  plot(m, pt_hi, type = "n", ylim = c(0, 1.25 * max(pt_lo, pt_hi, pg)),
       xlab = expression(paste("location ", mu)), ylab = "posterior density",
       main = "(d) generalized posterior, 10% outliers")
  grid_h()
  abline(v = 0, col = AXIS, lty = 3)
  lines(m, pg, col = ORANGE, lwd = 2, lty = 2)
  lines(m, pt_lo, col = BLUE, lwd = 2, lty = 4)
  lines(m, pt_hi, col = BLUE, lwd = 2)
  legend("topright", c("biweight, w = 2", "biweight, w = 0.5", "squared error, w = 1"),
         col = c(BLUE, BLUE, ORANGE), lty = c(1, 4, 2), lwd = 2, bty = "n", cex = 0.85)
  dev.off()
})

## ============================== four chains on Boston and riboflavin
if (have("posterior-fits.rds")) {
  pf <- readRDS(file.path(out, "posterior-fits.rds"))
  chc <- c(BLUE, ORANGE, GREEN, PURPLE)
  chp <- c(21, 24, 22, 23)
  rhat1 <- function(m) diag_chains(lapply(seq_len(ncol(m)), function(k) m[, k, drop = FALSE]))[1, "rhat"]
  trace_panel <- function(lam, name, tag) {
    lm <- log10(lam)
    it <- seq(1, nrow(lm), by = 4)
    rg <- range(lm)
    plot(NA, xlim = c(1, nrow(lm)), ylim = rg + c(0, 0.3 * diff(rg)),
         xlab = "iteration after burn-in", ylab = expression(log[10] ~ lambda[1]),
         main = bquote(.(tag) ~ .(name) * ":  rank" ~ hat(R) == .(sprintf("%.2f", rhat1(lm)))))
    grid_h()
    for (k in seq_len(ncol(lm))) lines(it, lm[it, k], col = chc[k], lwd = 0.7)
    legend("top", paste("chain", seq_len(ncol(lm))), col = chc, lwd = 2, horiz = TRUE,
           bty = "n", cex = 0.72)
  }
  pdf(file.path(out, "fig-mcmc.pdf"), width = 9.4, height = 3.0, pointsize = 10)
  par(mfrow = c(1, 3))
  pan()
  trace_panel(pf$boston$trace$lam1, "Boston", "(a)")
  trace_panel(pf$ribo$trace$lam1, "riboflavin", "(b)")
  tr <- pf$ribo$trace
  th <- seq(1, nrow(tr$lam1), by = 10)
  plot(NA, xlim = range(log10(tr$lam1)), ylim = range(tr$l1), log = "y",
       xlab = expression(log[10] ~ lambda[1]), ylab = expression(paste("||", beta, "||")[1] ~ "(standardized)"),
       main = "(c) riboflavin: two shrinkage modes")
  grid_h()
  for (k in seq_len(ncol(tr$lam1)))
    points(log10(tr$lam1[th, k]), tr$l1[th, k], col = chc[k], pch = chp[k], bg = "white", cex = 0.6)
  legend("topright", paste("chain", seq_len(ncol(tr$lam1))), col = chc, pch = chp,
         pt.bg = "white", bty = "n", cex = 0.78)
  dev.off()

  ## ---------------------------------------- the Boston posterior as an object
  b <- pf$boston
  z <- sweep(b$beta, 2, b$sd_x, "*")
  q <- t(apply(z, 2, stats::quantile, c(0.025, 0.25, 0.5, 0.75, 0.975)))
  o <- order(q[, 3])
  mapz <- (b$map * b$sd_x)[o]
  q <- q[o, ]
  yy <- seq_along(o)
  cl <- ifelse(b$selected[o], BLUE, "#8C8C8C")
  pdf(file.path(out, "fig-posterior.pdf"), width = 9.4, height = 3.3, pointsize = 10)
  par(mfrow = c(1, 3))
  pan()
  par(mar = c(3.6, 4.4, 2.0, 0.8))
  plot(NA, xlim = range(q, mapz), ylim = c(0.5, length(o) + 0.5), yaxt = "n",
       xlab = "coefficient x standard deviation of x", ylab = "", main = "(a) coefficients")
  axis(2, at = yy, labels = b$names[o], cex.axis = 0.8)
  grid(NULL, NA, col = GRID, lty = 1, lwd = 0.6)
  abline(v = 0, col = AXIS, lty = 3)
  segments(q[, 1], yy, q[, 5], yy, col = cl, lwd = 1)
  segments(q[, 2], yy, q[, 4], yy, col = cl, lwd = 3.5)
  points(q[, 3], yy, pch = 21, col = cl, bg = "white", cex = 0.9)
  points(mapz, yy, pch = 4, col = ORANGE, cex = 0.9, lwd = 1.4)
  legend("bottomright", c("median, 50% and 95% intervals", "interval contains 0", "posterior mode"),
         col = c(BLUE, "#8C8C8C", ORANGE), pch = c(21, 21, 4), lty = c(1, 1, NA), lwd = c(3, 3, 1.4),
         pt.bg = "white", bty = "n", cex = 0.72)
  pan()
  l1 <- as.vector(b$trace$lam1)
  l2 <- as.vector(b$trace$lam2)
  th <- seq(1, length(l1), by = 4)
  plot(l1[th], l2[th], log = "xy", pch = 16, cex = 0.4, col = grDevices::adjustcolor(BLUE, 0.3),
       xlab = expression(lambda[1]), ylab = expression(lambda[2]), main = "(b) penalty parameters")
  grid(col = GRID, lty = 1, lwd = 0.6)
  abline(v = stats::median(l1), h = stats::median(l2), col = ORANGE, lty = 2)
  legend("bottomleft", c("joint posterior draws", "posterior medians"), col = c(BLUE, ORANGE),
         pch = c(16, NA), lty = c(NA, 2), bty = "n", cex = 0.78)
  hi <- b$p_reject > 0.5
  plot(b$fitted, b$p_reject, type = "n", ylim = c(0, 1.08),
       xlab = "posterior median fit, log(medv)", ylab = "posterior probability of rejection",
       main = "(c) observations beyond the rejection radius")
  grid_h()
  abline(h = 0.5, col = AXIS, lty = 3)
  points(b$fitted[!hi], b$p_reject[!hi], pch = 21, col = BLUE, bg = "white", cex = 0.6)
  points(b$fitted[hi], b$p_reject[hi], pch = 24, col = ORANGE, bg = "white", cex = 0.9)
  legend("right", c(sprintf("%d tracts above 0.5", sum(hi)), sprintf("%d at or below", sum(!hi))),
         col = c(ORANGE, BLUE), pch = c(24, 21), pt.bg = "white", bty = "n", cex = 0.78)
  dev.off()
}
cat("figures written\n")
