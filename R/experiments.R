## Bayesian experiments for bayes-tadenet.tex.  These test claims that only exist
## in the Bayesian formulation; the selection-accuracy table (validate.R) is the
## baseline "does it work" check, not the point.
##
##   E1  calibration: does credible-interval coverage track the learning rate?
##   E2  Prop 2.3: is posterior sensitivity to one arbitrary outlier bounded?
##   E3  Cor 6.4: do the adaptive weights control the prior-leakage term?
##   E4  Thm 6.2: is excess Tukey risk U-shaped in w?
##
##   Rscript R/experiments.R
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
.libPaths(c(file.path(here, "..", "Rlib"), .libPaths()))
source(file.path(here, "bt_adenet.R"))
out <- file.path(here, "..")
## Rscript R/experiments.R [which, e.g. "E2"] [nrep] [cores].  E1, E3 and E4 are
## superseded by calibration2.R and ratio.R and kept for reproducibility.
args <- commandArgs(TRUE)
run <- if (length(args)) args[1] else "E1,E2,E3,E4"
nrep <- if (length(args) > 1) as.integer(args[2]) else NA
cores <- if (length(args) > 2) as.integer(args[3]) else 16

## A Gaussian-likelihood Bayesian adaptive elastic net comes free: T_d(u) -> u^2/2
## as d -> Inf, so the same sampler with a huge tuning constant IS the non-robust
## comparator.  One code path, so any difference is the loss and nothing else.
D_GAUSS <- Inf   # squared-error loss: the non-robust arm, same sampler

gen <- function(n, p, s, rho = 0.6, eps = 0.10, sdv = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  S <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- matrix(rnorm(n * p), n, p) %*% chol(S)
  beta <- numeric(p)
  beta[seq_len(s)] <- (2 * rbinom(s, 1, .5) - 1) * runif(s, 1, 3)
  e <- rnorm(n, 0, sdv)
  if (eps > 0) { o <- sample(n, round(eps * n)); e[o] <- rnorm(length(o), 15, sdv) }
  list(X = X, y = as.vector(X %*% beta) + e, beta = beta, Sigma = S)
}

## ============================================================== E1 calibration
## Prop 6.1 says credible sets are confidence sets only when w J = I.  Sweep w and
## watch marginal coverage cross the nominal level; check the calibrated values land
## near the crossing.
if (grepl("E1", run)) {
cat("\n=== E1: coverage as a function of the learning rate ===\n")
ws <- c(0.25, 0.5, 1, 1.5, 2, 3, 5)
R1 <- 20
e1 <- NULL
for (r in seq_len(R1)) {
  d <- gen(200, 20, 6, seed = 100 + r)
  A <- d$beta != 0
  ref <- bt_adenet(d$X, d$y, nsim = 4000, burn = 2000, chains = 1,
                   seed = r, verbose = FALSE)
  wcov <- calibrate_w(scale(d$X, apply(d$X, 2, median), apply(d$X, 2, mad)), d$y,
                      ref$map$alpha, ref$map$beta * 1, ref$sigma, method = "cov")
  for (w in ws) {
    f <- bt_adenet(d$X, d$y, learn = w, nsim = 4000, burn = 2000, chains = 1,
                   seed = r, verbose = FALSE)
    sm <- summary(f)
    e1 <- rbind(e1, data.frame(rep = r, w = w,
      cov = mean(sm$lower[A] <= d$beta[A] & d$beta[A] <= sm$upper[A]),
      width = mean(sm$upper[A] - sm$lower[A])))
  }
  e1 <- rbind(e1, data.frame(rep = r, w = NA, cov = NA, width = NA))  # spacer
  cat(sprintf("  rep %d/%d (w_hat = %.2f)\n", r, R1, ref$w)); flush.console()
}
e1 <- e1[!is.na(e1$w), ]
a1 <- aggregate(cbind(cov, width) ~ w, e1, mean)
print(a1, row.names = FALSE)
saveRDS(list(raw = e1, agg = a1), file.path(out, "e1-calibration.rds"))
}

## ============================================ E2 bounded posterior sensitivity
## Prop 2.3: shifting m = 1 observation by an arbitrary amount t moves the log
## posterior by at most w*d^2/6, uniformly.  Under a Gaussian likelihood the same
## quantity is unbounded.  Measure the posterior-mean displacement as t -> Inf.
if (grepl("E2", run)) {
cat("\n=== E2: posterior sensitivity to one arbitrary outlier ===\n")
ts <- c(0, 1, 5, 10, 50, 100, 1000, 1e5)
R2 <- if (is.na(nrep)) 10 else nrep
e2 <- do.call(rbind, par_lapply(R2, function(r) {
  d <- gen(150, 15, 6, eps = 0, seed = 200 + r)
  rows <- NULL
  for (dd in c(4.685, D_GAUSS)) {
    base <- NULL
    for (t in ts) {
      y2 <- d$y
      y2[1] <- y2[1] + t
      f <- try(bt_adenet(d$X, y2, d = dd, pilot = "lmrob", nsim = 3000, burn = 1500,
                         chains = 1, seed = r, verbose = FALSE), silent = TRUE)
      if (inherits(f, "try-error")) next
      pm <- colMeans(f$beta)
      if (is.null(base)) base <- pm
      ## the posterior mean is only part of the damage: a single outlier also
      ## destroys the scale estimate and, through it, the learning rate
      rows <- rbind(rows, data.frame(rep = r,
        loss = ifelse(is.finite(dd), "Tukey", "Gaussian"),
        t = t, shift = sqrt(sum((pm - base)^2)),
        err = sqrt(sum((pm - d$beta)^2)),
        sigma = f$sigma, w = f$w, size = sum(summary(f)$selected),
        engine = f$pilot_engine))
    }
  }
  rows
}, cores))
a2 <- aggregate(cbind(shift, err, sigma, w, size) ~ loss + t, e2, median)
print(a2[order(a2$loss, a2$t), ], row.names = FALSE)
saveRDS(list(raw = e2, agg = a2), file.path(out, "e2-sensitivity.rds"))
}

## ================================== E3 adaptive weights and the leakage term
## Cor 6.4: the prior-leakage term sum_{j notin A} 1/a_j is what decides whether the
## elastic-net prior behaves like a sparse prior.  iota = 0 switches the adaptive
## weighting off (w_j == 1), which is the Castillo et al. regime.
if (grepl("E3", run)) {
cat("\n=== E3: adaptive vs non-adaptive weights as p grows ===\n")
ps <- c(20, 50, 100, 200, 400)
R3 <- 10
e3 <- NULL
for (r in seq_len(R3)) for (p in ps) {
  d <- gen(200, p, 6, seed = 300 + r)
  A <- d$beta != 0
  for (io in c(1, 0)) {
    f <- bt_adenet(d$X, d$y, iota = io, nsim = 4000, burn = 2000, chains = 1,
                   seed = r, verbose = FALSE)
    sm <- summary(f)
    ## empirical leakage: posterior mass the inactive coordinates carry
    leak <- sum(abs(colMeans(f$beta))[!A])
    e3 <- rbind(e3, data.frame(rep = r, p = p,
      weights = ifelse(io > 0, "adaptive", "uniform"),
      leak = leak, l2 = sqrt(sum((colMeans(f$beta) - d$beta)^2)),
      tpr = sum(sm$selected & A) / sum(A),
      fpr = sum(sm$selected & !A) / sum(!A)))
  }
  cat(sprintf("  rep %d p=%d\n", r, p)); flush.console()
}
a3 <- aggregate(cbind(leak, l2, tpr, fpr) ~ weights + p, e3, median)
print(a3[order(a3$weights, a3$p), ], row.names = FALSE)
saveRDS(list(raw = e3, agg = a3), file.path(out, "e3-adaptive.rds"))
}

## ================================================ E4 excess risk against w
## Thm 6.2 predicts KL/(wn) + w B^2/4: decreasing then increasing in w.
if (grepl("E4", run)) {
cat("\n=== E4: excess Tukey risk against the learning rate ===\n")
R4 <- 10
e4 <- NULL
for (r in seq_len(R4)) {
  d <- gen(200, 30, 6, seed = 400 + r)
  te <- gen(2000, 30, 6, eps = 0, seed = 4000 + r)
  te$X <- te$X; te$y <- as.vector(te$X %*% d$beta) + rnorm(2000, 0, 2)
  for (w in ws) {
    f <- bt_adenet(d$X, d$y, learn = w, nsim = 4000, burn = 2000, chains = 1,
                   seed = r, verbose = FALSE)
    ## posterior-expected Tukey risk on a clean test sample
    idx <- sample(nrow(f$beta), 200)
    risk <- mean(apply(f$beta[idx, , drop = FALSE], 1, function(b)
      mean(tukey_rho((te$y - mean(f$alpha) - te$X %*% b) / f$sigma))))
    orc <- mean(tukey_rho((te$y - te$X %*% d$beta) / f$sigma))
    e4 <- rbind(e4, data.frame(rep = r, w = w, excess = risk - orc))
  }
  cat(sprintf("  rep %d/%d\n", r, R4)); flush.console()
}
a4 <- aggregate(excess ~ w, e4, median)
print(a4, row.names = FALSE)
saveRDS(list(raw = e4, agg = a4), file.path(out, "e4-risk.rds"))
}

cat("\nall experiments written to *.rds\n")
