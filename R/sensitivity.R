## Sensitivity to the constants fixed by convention (REVIEW-R2 roadmap item 10).
##   Rscript R/sensitivity.R [nrep] [cores]
## n = 200, p = 50, s = 6, rho = 0.6, 10% vertical outliers, MM pilot.  One knob is
## varied at a time with the others at their defaults.
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
.libPaths(c(file.path(here, "..", "Rlib"), .libPaths()))
source(file.path(here, "bt_adenet.R"))
out <- file.path(here, "..")
args <- as.integer(commandArgs(TRUE))
R <- if (length(args) >= 1) args[1] else 50
cores <- if (length(args) >= 2) args[2] else 16

gen <- function(seed, p = 50, s = 6, n = 200) {
  set.seed(seed)
  S <- 0.6^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- matrix(rnorm(n * p), n, p) %*% chol(S)
  b <- numeric(p)
  b[seq_len(s)] <- (2 * rbinom(s, 1, .5) - 1) * runif(s, 1, 3)
  e <- rnorm(n, 0, 2)
  o <- sample(n, 20)
  e[o] <- rnorm(20, 15, 2)
  list(X = X, y = as.vector(X %*% b) + e, beta = b, Sigma = S)
}

grid <- rbind(
  data.frame(knob = "d",     value = c(3.0, 4.0, 4.685, 5.5, 7.0), stringsAsFactors = FALSE),
  data.frame(knob = "iota",  value = c(0.5, 1, 1.5, 2)),
  data.frame(knob = "b_lam", value = c(1e-3, 1e-2, 1e-1, 1)),
  data.frame(knob = "floor", value = c(0, 1, 2)))    # 0 default, 1 p-scaled, 2 fixed 1e-3
tasks <- merge(grid, data.frame(rep = seq_len(R)), by = NULL)

run <- function(i) {
  g <- tasks[i, ]
  d <- gen(12000 + g$rep)
  a <- list(d = 4.685, iota = 1, b_lam = 1e-2, w_floor = NULL)
  if (g$knob == "d") a$d <- g$value
  if (g$knob == "iota") a$iota <- g$value
  if (g$knob == "b_lam") a$b_lam <- g$value
  if (g$knob == "floor") a$w_floor <- if (g$value == 1) "p" else if (g$value == 2) 1e-3 else NULL
  f <- bt_adenet(d$X, d$y, d = a$d, iota = a$iota, b_lam = a$b_lam,
                 w_floor = a$w_floor, pilot = "lmrob", nsim = 3000, burn = 1500,
                 chains = 2, diagnostics = TRUE, seed = g$rep, verbose = FALSE)
  sm <- summary(f)
  b <- coef(f)
  A <- d$beta != 0
  data.frame(knob = g$knob, value = g$value, rep = g$rep,
    me = as.numeric(t(b - d$beta) %*% d$Sigma %*% (b - d$beta)),
    tpr = sum(sm$selected & A) / sum(A), fpr = sum(sm$selected & !A) / sum(!A),
    cov = mean(sm$lower[A] <= d$beta[A] & d$beta[A] <= sm$upper[A]),
    sigma = f$sigma, w = f$w, size = sum(sm$selected),
    rhat_rank = max(f$diag[, "rhat"], na.rm = TRUE),
    ess_bulk = min(f$diag[, "ess_bulk"], na.rm = TRUE), engine = f$pilot_engine)
}

res <- do.call(rbind, par_lapply(nrow(tasks), run, cores))
saveRDS(res, file.path(out, "sensitivity.rds"))
agg <- aggregate(cbind(me, tpr, fpr, cov, w, rhat_rank) ~ knob + value, res, mean)
print(agg[order(agg$knob, agg$value), ], row.names = FALSE)
print(table(res$engine))
cat("\nwrote sensitivity.rds\n")
