## Why the prior ranking reverses between simulation and real data (REVIEW-R2
## roadmap item 9).  Candidate explanation: the spike-and-slab wins when the signal
## is sparse and strong and loses when it is spread thinly over many coefficients,
## as in the spectral and gene-expression designs.  Tukey loss throughout, three
## priors, three signal regimes, n = 200, p = 100, rho = 0.6, 10% vertical
## outliers in training, clean test set of 1000, IRLS pilot (p = n/2).
## Selection stability is the Jaccard index between fits on two random 3/4
## subsamples, with a pair of empty sets counted as missing rather than as 1.
##   Rscript R/diffuse.R [nrep] [cores]
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
.libPaths(c(file.path(here, "..", "Rlib"), .libPaths()))
source(file.path(here, "bt_adenet.R"))
out <- file.path(here, "..")
args <- as.integer(commandArgs(TRUE))
R <- if (length(args) >= 1) args[1] else 30
cores <- if (length(args) >= 2) args[2] else 16

n <- 200
p <- 100
regimes <- data.frame(regime = c("sparse", "moderate", "diffuse"), s = c(6, 20, 60),
                      lo = c(1, 0.5, 0.1), hi = c(3, 1.5, 0.5), stringsAsFactors = FALSE)
gen <- function(nn, reg, seed, beta = NULL, eps = 0.10) {
  set.seed(seed)
  S <- 0.6^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- matrix(rnorm(nn * p), nn, p) %*% chol(S)
  if (is.null(beta)) {
    beta <- numeric(p)
    beta[seq_len(reg$s)] <- (2 * rbinom(reg$s, 1, .5) - 1) * runif(reg$s, reg$lo, reg$hi)
  }
  e <- rnorm(nn, 0, 2)
  if (eps > 0) {
    o <- sample(nn, round(eps * nn))
    e[o] <- rnorm(length(o), 15, 2)
  }
  list(X = X, y = as.vector(X %*% beta) + e, beta = beta, Sigma = S)
}
jac <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0) NA_real_ else length(intersect(a, b)) / u
}

tasks <- merge(merge(regimes["regime"], data.frame(prior = c("adenet", "horseshoe", "ssvs"),
               stringsAsFactors = FALSE), by = NULL), data.frame(rep = seq_len(R)), by = NULL)

run <- function(i) {
  k <- tasks[i, ]
  cell <- match(k$regime, regimes$regime)
  reg <- regimes[cell, ]
  d <- gen(n, reg, seed = 20000 + 100 * cell + k$rep)
  te <- gen(1000, reg, seed = 30000 + 100 * cell + k$rep, beta = d$beta, eps = 0)
  A <- d$beta != 0
  fit <- function(idx) bt_adenet(d$X[idx, ], d$y[idx], prior = k$prior, pilot = "irls",
                                 nsim = 4000, burn = 2000, chains = 2, seed = k$rep,
                                 verbose = FALSE)
  f <- fit(seq_len(n))
  sm <- summary(f)
  b <- colMeans(f$beta)
  pred <- mean(f$alpha) + as.vector(te$X %*% b)
  set.seed(40000 + k$rep)
  i1 <- sample(n, 0.75 * n)
  i2 <- sample(n, 0.75 * n)
  s1 <- which(summary(fit(i1))$selected)
  s2 <- which(summary(fit(i2))$selected)
  data.frame(regime = k$regime, prior = k$prior, rep = k$rep,
    rmspe = sqrt(mean((te$y - pred)^2)),
    me = as.numeric(t(b - d$beta) %*% d$Sigma %*% (b - d$beta)),
    tpr = sum(sm$selected & A) / sum(A), fpr = sum(sm$selected & !A) / sum(!A),
    size = sum(sm$selected), cov = mean(sm$lower[A] <= d$beta[A] & d$beta[A] <= sm$upper[A]),
    stab = jac(s1, s2), null_pair = length(s1) + length(s2) == 0)
}

res <- do.call(rbind, par_lapply(nrow(tasks), run, cores))
saveRDS(res, file.path(out, "diffuse.rds"))
agg <- aggregate(cbind(rmspe, me, tpr, fpr, size, cov) ~ regime + prior, res, mean)
agg$stab <- aggregate(stab ~ regime + prior, res, mean, na.action = na.pass,
                      na.rm = TRUE)$stab
print(agg[order(agg$regime, agg$rmspe), ], row.names = FALSE)
w <- reshape(res[, c("regime", "prior", "rep", "rmspe")], idvar = c("regime", "rep"),
             timevar = "prior", direction = "wide")
print(aggregate(cbind(adenet_beats_ssvs = rmspe.adenet < rmspe.ssvs) ~ regime, w, mean))
cat("\nwrote diffuse.rds\n")
