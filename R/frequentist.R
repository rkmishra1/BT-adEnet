## A frequentist robust selection comparator on the design of comparators.R
## (REVIEW-R2 roadmap item 11): enetLTS (Kurnaz, Hoffmann and Filzmoser 2018)
## with its default tuning, reweighted coefficients.
##   Rscript R/frequentist.R [nrep] [cores]
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
.libPaths(c(file.path(here, "..", "Rlib"), .libPaths()))
source(file.path(here, "bt_adenet.R"))
suppressMessages(library(enetLTS))
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
cfg <- expand.grid(scenario = c("clean", "vertical", "both"), p = c(35L, 100L),
                   rep = seq_len(R), stringsAsFactors = FALSE)
cfg$s <- 3 * floor(cfg$p^(1/3))
cfg$cell <- match(paste(cfg$scenario, cfg$p), unique(paste(cfg$scenario, cfg$p)))

run <- function(i) {
  k <- cfg[i, ]
  d <- gen(200, k$p, k$s, 0.6, k$scenario, seed = 9000 + 100 * k$cell + k$rep)
  A <- d$beta != 0
  t0 <- proc.time()[3]
  f <- enetLTS(d$X, d$y, seed = k$rep)
  cf <- as.numeric(f$coefficients)
  b <- if (length(cf) == k$p + 1) cf[-1] else cf
  sel <- b != 0
  data.frame(arm = "enetLTS", scenario = k$scenario, p = k$p, rep = k$rep,
    me = as.numeric(t(b - d$beta) %*% d$Sigma %*% (b - d$beta)),
    tpr = sum(sel & A) / sum(A), fpr = sum(sel & !A) / sum(!A), size = sum(sel),
    alpha = f$alpha, ncoef = length(cf), secs = proc.time()[3] - t0)
}

res <- do.call(rbind, par_lapply(nrow(cfg), run, cores))
saveRDS(res, file.path(out, "frequentist.rds"))
print(aggregate(cbind(me, tpr, fpr, size, secs) ~ scenario + p, res, mean), row.names = FALSE)
cat("\nwrote frequentist.rds\n")
