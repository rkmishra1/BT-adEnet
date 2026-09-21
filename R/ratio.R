## E3 and E3b redone with the scale-free quantity of the revised Corollary 6.4,
## the separation ratio min_{j notin A} w_j / max_{j in A} w_j of the weights the
## sampler actually used (REVIEW-R2 roadmap item 3).
##   Rscript R/ratio.R [nrep] [cores]
## Design as E3: n = 200, s = 6, rho = 0.6, 10% vertical outliers.  The pilot is
## IRLS at every p, so that the pilot does not change along the p axis.
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
.libPaths(c(file.path(here, "..", "Rlib"), .libPaths()))
source(file.path(here, "bt_adenet.R"))
out <- file.path(here, "..")
args <- as.integer(commandArgs(TRUE))
R <- if (length(args) >= 1) args[1] else 20
cores <- if (length(args) >= 2) args[2] else 16

gen <- function(n, p, s, seed) {
  set.seed(seed)
  S <- 0.6^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- matrix(rnorm(n * p), n, p) %*% chol(S)
  beta <- numeric(p)
  beta[seq_len(s)] <- (2 * rbinom(s, 1, .5) - 1) * runif(s, 1, 3)
  e <- rnorm(n, 0, 2)
  o <- sample(n, round(0.10 * n))
  e[o] <- rnorm(length(o), 15, 2)
  list(X = X, y = as.vector(X %*% beta) + e, beta = beta)
}

ps <- c(20, 50, 100, 200, 400)
cfg <- rbind(
  expand.grid(p = ps, iota = c(0, 0.5, 1, 1.5, 2), floor = "default", rep = seq_len(R),
              stringsAsFactors = FALSE),
  expand.grid(p = ps, iota = 1, floor = "p", rep = seq_len(R), stringsAsFactors = FALSE))

run <- function(i) {
  k <- cfg[i, ]
  d <- gen(200, k$p, 6, seed = 300 + k$rep)
  A <- d$beta != 0
  f <- bt_adenet(d$X, d$y, iota = k$iota, w_floor = if (k$floor == "p") "p" else NULL,
                 pilot = "irls", nsim = 4000, burn = 2000, chains = 1, seed = k$rep,
                 verbose = FALSE)
  pm <- colMeans(f$beta)
  sm <- summary(f)
  lam1 <- stats::median(f$lambda[, 1])
  data.frame(k, leak = sum(abs(pm)[!A]), l2 = sqrt(sum((pm - d$beta)^2)),
    tpr = sum(sm$selected & A) / sum(A), fpr = sum(sm$selected & !A) / sum(!A),
    ratio = min(f$weights[!A]) / max(f$weights[A]),
    W = min(f$weights[!A]), wmaxA = max(f$weights[A]),
    prior_leak = sum(1 / (f$w * lam1 * f$weights[!A])), w = f$w, lam1 = lam1)
}

res <- do.call(rbind, par_lapply(nrow(cfg), run, cores))
saveRDS(res, file.path(out, "ratio.rds"))
res$grp <- paste0("iota=", res$iota, " floor=", res$floor)
sl <- do.call(rbind, lapply(split(res, res$grp), function(g) {
  m <- stats::lm(log(leak) ~ log(p), data = g)
  data.frame(group = g$grp[1], slope = unname(coef(m)[2]),
             se = unname(sqrt(diag(stats::vcov(m)))[2]))
}))
print(sl, row.names = FALSE)
agg <- aggregate(cbind(leak, prior_leak, ratio, W, l2, tpr, fpr) ~ iota + floor + p, res,
                 stats::median)
print(agg[order(agg$floor, agg$iota, agg$p), ], row.names = FALSE)
m <- stats::lm(log(leak) ~ log(ratio) + factor(p), data = res)
cat(sprintf("log leak on log ratio within p: slope %.3f (se %.3f)\n",
            coef(m)[2], sqrt(diag(stats::vcov(m)))[2]))
cat("\nwrote ratio.rds\n")
