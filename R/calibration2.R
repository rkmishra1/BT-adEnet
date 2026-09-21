## Calibration study on a single design (REVIEW-R2 roadmap items 10 and 11).
##   Rscript R/calibration2.R [nrep] [nrep_freq] [B] [cores]
##
## Design: n = 200, p = 20, s = 6, rho = 0.6, 10% vertical outliers, MM pilot.
## One design serves E1 (coverage against w) and E4 (excess Tukey risk against w on
## a clean test set of 2000), so the two optimal learning rates are comparable.
## The calibrated posterior runs four chains from diverse starts with
## rank-normalized diagnostics.  Its coverage error is split into shrinkage bias
## and dispersion through z_j = (median_j - beta_j) / sd_j on the active set.  The
## first nrep_freq replications also get two frequentist intervals for the same
## T-AdEnet estimate: a percentile bootstrap with B resamples, and the sandwich
## J^-1 I J^-1 on the intercept and the selected support.
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
.libPaths(c(file.path(here, "..", "Rlib"), .libPaths()))
source(file.path(here, "bt_adenet.R"))
out <- file.path(here, "..")
args <- as.integer(commandArgs(TRUE))
arg <- function(k, def) if (length(args) >= k) args[k] else def
R <- arg(1, 200)
RF <- arg(2, 100)
B <- arg(3, 100)
cores <- arg(4, 16)

n <- 200
p <- 20
s <- 6
ws <- c(0.25, 0.5, 1, 1.5, 2, 3, 5)
gen <- function(nn, eps, seed) {
  set.seed(seed)
  S <- 0.6^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- matrix(rnorm(nn * p), nn, p) %*% chol(S)
  beta <- numeric(p)
  beta[seq_len(s)] <- (2 * rbinom(s, 1, .5) - 1) * runif(s, 1, 3)
  e <- rnorm(nn, 0, 2)
  if (eps > 0) {
    o <- sample(nn, round(eps * nn))
    e[o] <- rnorm(length(o), 15, 2)
  }
  list(X = X, y = as.vector(X %*% beta) + e, beta = beta)
}
covers <- function(lo, hi, b) mean(lo <= b & b <= hi)

one_rep <- function(r) {
  d <- gen(n, 0.10, 100 + r)
  te <- gen(2000, 0, 4000 + r)
  te$y <- as.vector(te$X %*% d$beta) + stats::rnorm(2000, 0, 2)
  A <- d$beta != 0

  sw <- do.call(rbind, lapply(ws, function(w) {
    g <- bt_adenet(d$X, d$y, learn = w, pilot = "lmrob", nsim = 4000, burn = 2000,
                   chains = 1, seed = r, verbose = FALSE)
    sg <- summary(g)
    idx <- round(seq(1, nrow(g$beta), length.out = 200))
    risk <- mean(vapply(idx, function(k)
      mean(tukey_rho((te$y - g$alpha[k] - te$X %*% g$beta[k, ]) / g$sigma)), 0))
    data.frame(rep = r, w = w, cov = covers(sg$lower[A], sg$upper[A], d$beta[A]),
               width = mean(sg$upper[A] - sg$lower[A]),
               excess = risk - mean(tukey_rho((te$y - te$X %*% d$beta) / g$sigma)))
  }))

  t0 <- proc.time()[3]
  f <- bt_adenet(d$X, d$y, pilot = "lmrob", nsim = 4000, burn = 2000, chains = 4,
                 starts = "diverse", diagnostics = TRUE, seed = r, verbose = FALSE)
  tb <- proc.time()[3] - t0
  sm <- summary(f)
  z <- (sm$median - d$beta) / apply(f$beta, 2, stats::sd)
  meth <- data.frame(rep = r, method = "BT-AdEnet", w = f$w,
    cov = covers(sm$lower[A], sm$upper[A], d$beta[A]),
    width = mean(sm$upper[A] - sm$lower[A]), tpr = mean(sm$selected[A]), secs = tb,
    rhat = max(f$diag[, "rhat"]), ess_bulk = min(f$diag[, "ess_bulk"]),
    ess_tail = min(f$diag[, "ess_tail"]), engine = f$pilot_engine)
  zz <- data.frame(rep = r, j = which(A), z = z[A], toward0 = sign(d$beta[A]) * z[A])

  if (r <= RF) {
    point <- function(i) bt_adenet(d$X[i, ], d$y[i], learn = 1, pilot = "lmrob",
                                   nsim = 1, burn = 0, chains = 1, verbose = FALSE)
    t0 <- proc.time()[3]
    m0 <- point(seq_len(n))
    tp <- proc.time()[3] - t0
    S <- which(m0$map$beta != 0)
    u <- as.vector(d$y - m0$map$alpha - d$X %*% m0$map$beta) / m0$sigma
    Z <- cbind(1, d$X[, S, drop = FALSE])
    J <- crossprod(Z, Z * tukey_psi1(u)) / n
    I <- crossprod(Z * tukey_psi(u)) / n
    Ji <- tryCatch(solve(J), error = function(e) NULL)
    se <- rep(0, p)
    ok <- FALSE
    if (!is.null(Ji)) {
      v <- diag(Ji %*% I %*% Ji) * m0$sigma^2 / n
      ok <- all(v > 0)
      if (ok) se[S] <- sqrt(v[-1])
    }
    t0 <- proc.time()[3]
    bs <- t(vapply(seq_len(B), function(b) {
      set.seed(1e6 * r + b)
      tryCatch(point(sample(n, replace = TRUE))$map$beta, error = function(e) rep(NA_real_, p))
    }, numeric(p)))
    tboot <- tp + proc.time()[3] - t0
    lo <- apply(bs, 2, stats::quantile, 0.025, na.rm = TRUE)
    hi <- apply(bs, 2, stats::quantile, 0.975, na.rm = TRUE)
    meth <- rbind(meth,
      data.frame(rep = r, method = "sandwich", w = NA,
        cov = if (ok) covers((m0$map$beta - 1.96 * se)[A], (m0$map$beta + 1.96 * se)[A], d$beta[A]) else NA,
        width = if (ok) mean(2 * 1.96 * se[A]) else NA, tpr = mean(m0$map$beta[A] != 0),
        secs = tp, rhat = NA, ess_bulk = NA, ess_tail = NA, engine = m0$pilot_engine),
      data.frame(rep = r, method = "bootstrap", w = NA,
        cov = covers(lo[A], hi[A], d$beta[A]), width = mean(hi[A] - lo[A]),
        tpr = mean(m0$map$beta[A] != 0), secs = tboot, rhat = NA, ess_bulk = NA,
        ess_tail = NA, engine = m0$pilot_engine))
  }
  list(sweep = sw, meth = meth, z = zz)
}

res <- par_lapply(R, one_rep, cores)
sw <- do.call(rbind, lapply(res, `[[`, "sweep"))
me <- do.call(rbind, lapply(res, `[[`, "meth"))
zz <- do.call(rbind, lapply(res, `[[`, "z"))
mcse <- function(v) stats::sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v)))
a_sw <- do.call(rbind, lapply(split(sw, sw$w), function(g) data.frame(w = g$w[1], R = nrow(g),
  cov = mean(g$cov), cov_se = mcse(g$cov), width = mean(g$width),
  excess = mean(g$excess), excess_se = mcse(g$excess))))
a_me <- do.call(rbind, lapply(split(me, me$method), function(g) data.frame(method = g$method[1],
  R = sum(!is.na(g$cov)), cov = mean(g$cov, na.rm = TRUE), cov_se = mcse(g$cov),
  width = mean(g$width, na.rm = TRUE), tpr = mean(g$tpr), secs = stats::median(g$secs),
  w = stats::median(g$w), rhat = stats::median(g$rhat), ess_bulk = stats::median(g$ess_bulk),
  ess_tail = stats::median(g$ess_tail))))
mu <- mean(zz$toward0)
sz <- stats::sd(zz$toward0)
a_z <- data.frame(nz = nrow(zz), toward0 = mu, toward0_se = mcse(zz$toward0), sd = sz,
  inside_196 = mean(abs(zz$z) <= 1.96),
  normal_fit = stats::pnorm((1.96 - mu) / sz) - stats::pnorm((-1.96 - mu) / sz),
  no_bias = 2 * stats::pnorm(1.96 / sz) - 1)
k <- which(a_sw$cov[-1] < 0.95 & a_sw$cov[-nrow(a_sw)] >= 0.95)[1]
cross <- if (is.na(k)) NA else a_sw$w[k] + (a_sw$cov[k] - 0.95) /
  (a_sw$cov[k] - a_sw$cov[k + 1]) * (a_sw$w[k + 1] - a_sw$w[k])
print(a_sw, row.names = FALSE)
print(a_me, row.names = FALSE)
print(a_z, row.names = FALSE)
cat(sprintf("coverage crosses 0.95 at w = %.3f, risk minimized at w = %.2f, engines: %s\n",
            cross, a_sw$w[which.min(a_sw$excess)], paste(unique(me$engine), collapse = ",")))
saveRDS(list(sweep = sw, meth = me, z = zz, agg_sweep = a_sw, agg_meth = a_me, agg_z = a_z,
             cross = cross, R = R, RF = RF, B = B), file.path(out, "calib2.rds"))
cat("\nwrote calib2.rds\n")
