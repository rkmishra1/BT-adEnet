## Comparison across losses and priors (REVIEW-R2 roadmap item 10).
##   Rscript R/comparators.R [nrep] [cores]
##
## Two factors, one sampler, so every difference is attributable:
##   loss  : Tukey biweight (redescending) / pseudo-Huber (convex, the Bayesian
##           Huberized lasso family) / squared error (d = Inf)
##   prior : adaptive elastic net / horseshoe (Makalic and Schmidt augmentation) /
##           continuous spike-and-slab (George and McCulloch)
## The MM pilot runs where p < n/2 (p = 35) and IRLS otherwise (p = 100).  Seeds
## match the first run, so its 20 replications are the first 20 here.
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
.libPaths(c(file.path(here, "..", "Rlib"), .libPaths()))
source(file.path(here, "bt_adenet.R"))
out <- file.path(here, "..")
args <- as.integer(commandArgs(TRUE))
R <- if (length(args) >= 1) args[1] else 50
cores <- if (length(args) >= 2) args[2] else 16

gen <- function(n, p, s, rho, scenario, sdv = 2, seed) {
  set.seed(seed)
  S <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- matrix(rnorm(n * p), n, p) %*% chol(S)
  beta <- numeric(p)
  beta[seq_len(s)] <- (2 * rbinom(s, 1, .5) - 1) * runif(s, 1, 3)
  e <- rnorm(n, 0, sdv)
  Xo <- X
  if (scenario != "clean") {
    o <- sample(n, round(0.10 * n))
    e[o] <- rnorm(length(o), 15, sdv)
    if (scenario == "both") Xo[o, ] <- rnorm(length(o) * p, 5, 1)
  }
  list(X = Xo, y = as.vector(X %*% beta) + e, beta = beta, Sigma = S)
}

arms <- list(
  list(id = "Tukey + AdEnet",    fam = "tukey", d = 4.685, prior = "adenet"),
  list(id = "Tukey + horseshoe", fam = "tukey", d = 4.685, prior = "horseshoe"),
  list(id = "Huber + AdEnet",    fam = "huber", d = 1.345, prior = "adenet"),
  list(id = "Huber + horseshoe", fam = "huber", d = 1.345, prior = "horseshoe"),
  list(id = "Gauss + AdEnet",    fam = "tukey", d = Inf,   prior = "adenet"),
  list(id = "Gauss + horseshoe", fam = "tukey", d = Inf,   prior = "horseshoe"),
  list(id = "Tukey + SSVS",      fam = "tukey", d = 4.685, prior = "ssvs"),
  list(id = "Huber + SSVS",      fam = "huber", d = 1.345, prior = "ssvs"),
  list(id = "Gauss + SSVS",      fam = "tukey", d = Inf,   prior = "ssvs"))

cfg <- expand.grid(scenario = c("clean", "vertical", "both"), p = c(35L, 100L),
                   rep = seq_len(R), stringsAsFactors = FALSE)
cfg$s <- 3 * floor(cfg$p^(1/3))
cfg$cell <- match(paste(cfg$scenario, cfg$p), unique(paste(cfg$scenario, cfg$p)))

run <- function(i) {
  k <- cfg[i, ]
  d <- gen(200, k$p, k$s, 0.6, k$scenario, seed = 9000 + 100 * k$cell + k$rep)
  A <- d$beta != 0
  do.call(rbind, lapply(arms, function(a) tryCatch({
    t0 <- proc.time()[3]
    f <- bt_adenet(d$X, d$y, fam = a$fam, d = a$d, prior = a$prior,
                   pilot = if (k$p < 100) "lmrob" else "irls",
                   nsim = 4000, burn = 2000, chains = 2, diagnostics = TRUE,
                   seed = k$rep, verbose = FALSE)
    sm <- summary(f)
    b <- coef(f)
    data.frame(arm = a$id, scenario = k$scenario, p = k$p, rep = k$rep,
      me = as.numeric(t(b - d$beta) %*% d$Sigma %*% (b - d$beta)),
      l2 = sqrt(sum((b - d$beta)^2)),
      tpr = sum(sm$selected & A) / sum(A),
      fpr = sum(sm$selected & !A) / sum(!A),
      size = sum(sm$selected),
      cov = mean(sm$lower[A] <= d$beta[A] & d$beta[A] <= sm$upper[A]),
      leak = sum(abs(b)[!A]),
      rhat = stats::median(f$rhat, na.rm = TRUE),
      rhat_rank = max(f$diag[, "rhat"], na.rm = TRUE),
      ess_bulk = min(f$diag[, "ess_bulk"], na.rm = TRUE),
      ess_tail = min(f$diag[, "ess_tail"], na.rm = TRUE),
      engine = f$pilot_engine, secs = proc.time()[3] - t0)
  }, error = function(e) NULL)))
}

res <- do.call(rbind, par_lapply(nrow(cfg), run, cores))
saveRDS(res, file.path(out, "comparators.rds"))
agg <- aggregate(cbind(me, tpr, fpr, cov, leak, rhat_rank, ess_bulk) ~ arm + scenario + p,
                 res, mean)
print(agg[order(agg$p, agg$scenario, agg$me), ], row.names = FALSE)
print(table(res$engine, res$p))
cat("\nwrote comparators.rds\n")
