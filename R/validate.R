## NOTE: the 2026-09-11 run requested the MM pilot, but robustbase was not on the
## library path and init_robust() fell back to IRLS without saying so.  Rerun
## 2026-09-14 with Rlib on the path, four chains from diverse starts, and
## rank-normalized diagnostics.  init_robust() now stops instead of falling back.
## Replicate study behind Table 1 of bayes-tadenet.tex.
##   Rscript R/validate.R [nrep]
## Writes ../sim-table.tex and ../sim-results.rds
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
.libPaths(c(file.path(here, "..", "Rlib"), .libPaths()))
source(file.path(here, "bt_adenet.R"))
##   Rscript R/validate.R [nrep] [cores] [largest p to run] [refit 1/0]
args <- as.numeric(commandArgs(TRUE))
R <- if (length(args) >= 1) args[1] else 50
cores <- if (length(args) >= 2) args[2] else 16
pmax_run <- if (length(args) >= 3) args[3] else Inf
refit <- if (length(args) >= 4) args[4] == 1 else TRUE
## a fifth argument pmin reruns only configurations with p >= pmin and keeps the
## other rows of the existing sim-results.rds
pmin_run <- if (length(args) >= 5) args[5] else 0

gen <- function(n, p, s, rho, scenario, sdv) {
  S <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- matrix(rnorm(n * p), n, p) %*% chol(S)
  beta <- numeric(p)
  beta[seq_len(s)] <- (2 * rbinom(s, 1, .5) - 1) * runif(s, 1, 3)
  e <- rnorm(n, 0, sdv)
  Xobs <- X
  if (scenario != "clean") {
    o <- sample(n, round(0.10 * n))
    e[o] <- rnorm(length(o), 15, sdv)
    if (scenario == "both") Xobs[o, ] <- rnorm(length(o) * p, 5, 1)  # bad leverage
  }
  list(X = Xobs, y = as.vector(X %*% beta) + e, beta = beta, Sigma = S)
}

cfg <- expand.grid(
  scenario = c("clean", "vertical", "both"),
  rho = c(0.6, 0.8), p = c(15L, 35L),
  stringsAsFactors = FALSE)
cfg <- rbind(cfg, data.frame(scenario = c("clean", "vertical", "both"),
                             rho = 0.6, p = 752L))
cfg$s <- 3 * floor(cfg$p^(1 / 3))
cfg$R <- ifelse(cfg$p > 100, max(10, R %/% 2), R)

me <- function(b, truth, S) as.numeric(t(b - truth) %*% S %*% (b - truth))
res <- list()

if (refit) {
tasks <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i)
  if (cfg$p[i] <= pmax_run && cfg$p[i] >= pmin_run) data.frame(i = i, r = seq_len(cfg$R[i])) else NULL))
fit_one <- function(t) {
  i <- tasks$i[t]
  r <- tasks$r[t]
  k <- cfg[i, ]
  set.seed(1000 * i + r)
  d <- gen(200, k$p, k$s, k$rho, k$scenario, 2)
  t0 <- proc.time()[3]
  ## Sec. 8 uses the MM pilot where it applies, and Sec. 7 must match, since the
  ## pilot fixes the adaptive weights (Cor 6.4).
  f <- bt_adenet(d$X, d$y, nsim = if (k$p > 100) 12000 else 5000,
                 burn = if (k$p > 100) 6000 else 2500,
                 pilot = if (k$p < 100) "lmrob" else "irls",
                 chains = 4, starts = "diverse", diagnostics = TRUE,
                 seed = r, verbose = FALSE)
  sm <- summary(f)
  A <- d$beta != 0
  cat(sprintf("[%d] p=%d rho=%.1f %-8s rep %d done\n", i, k$p, k$rho, k$scenario, r))
  cbind(k, rep = r, data.frame(
    me_map = me(f$map$beta, d$beta, d$Sigma),
    me_med = me(coef(f), d$beta, d$Sigma),
    me_sel = me(selected_coef(f), d$beta, d$Sigma),
    tpr = sum(sm$selected & A) / sum(A),
    fpr = sum(sm$selected & !A) / sum(!A),
    size = sum(sm$selected),
    exact = as.numeric(all(sm$selected == A)),
    cov = mean(sm$lower[A] <= d$beta[A] & d$beta[A] <= sm$upper[A]),
    w = f$w, sigma = f$sigma,
    rhat = max(f$rhat, na.rm = TRUE),
    rhat_rank = max(f$diag[, "rhat"], na.rm = TRUE),
    rank_bad = mean(f$diag[, "rhat"] > 1.01, na.rm = TRUE),
    ess_bulk = min(f$diag[, "ess_bulk"], na.rm = TRUE),
    ess_tail = min(f$diag[, "ess_tail"], na.rm = TRUE),
    engine = f$pilot_engine,
    secs = proc.time()[3] - t0))
}
res <- do.call(rbind, par_lapply(nrow(tasks), fit_one, cores))
if (pmin_run > 0 && file.exists(file.path(here, "..", "sim-results.rds"))) {
  old <- readRDS(file.path(here, "..", "sim-results.rds"))
  res <- rbind(old[old$p < pmin_run, names(res)], res)
}
saveRDS(res, file.path(here, "..", "sim-results.rds"))
}

res <- readRDS(file.path(here, "..", "sim-results.rds"))

agg <- do.call(rbind, lapply(split(res, list(res$p, res$rho, res$scenario), drop = TRUE),
  function(z) data.frame(p = z$p[1], rho = z$rho[1], scenario = z$scenario[1],
    s = z$s[1], R = nrow(z),
    me_map = median(z$me_map), me_med = median(z$me_med), me_sel = median(z$me_sel),
    tpr = mean(z$tpr), fpr = mean(z$fpr), size = mean(z$size),
    exact = mean(z$exact), cov = mean(z$cov), cov_se = sd(z$cov) / sqrt(nrow(z)),
    w = median(z$w), sigma = median(z$sigma), rhat = median(z$rhat_rank),
    ess_tail = median(z$ess_tail), engine = z$engine[1], secs = median(z$secs))))
agg <- agg[order(agg$p, agg$rho, match(agg$scenario, c("clean", "vertical", "both"))), ]

lab <- c(clean = "clean", vertical = "vertical", both = "leverage")
con <- file(file.path(here, "..", "sim-table.tex"), "w")
cap <- paste0(
  "\\caption{Estimation and selection, ", max(agg$R[agg$p < 100]), " replications per ",
  "configuration and ", min(agg$R), " at $p = 752$. The ME columns are the median model error of ",
  "the posterior mode, the posterior median and the summary \\eqref{eq:select}, which also defines ",
  "TPR, FPR, $|\\hat{\\mathcal{A}}|$ and exact recovery. Cov is the mean coverage of $95\\%$ ",
  "intervals on the active set, $\\hat R$ and ESS$_{\\text{tail}}$ are medians for the worst ",
  "coordinate, and $^{\\dagger}$ marks $\\hat R \\ge 1.05$.}")
writeLines(c(
"\\begin{table}[!htb]", "\\centering", cap,
"\\label{tab:sim}", "\\resizebox{\\textwidth}{!}{%",
"\\begin{tabular}{llrrrrrrrrlrrr}", "\\toprule",
"$p$ & scenario & $\\rho$ & ME$_{\\text{mode}}$ & ME$_{\\text{med}}$ & ME$_{\\text{sel}}$ & TPR & FPR & $|\\hat{\\mathcal{A}}|$ & exact & cov & $w$ & $\\hat R$ & ESS$_{\\text{tail}}$ \\\\",
"\\midrule"), con)
prev <- ""
for (i in seq_len(nrow(agg))) {
  z <- agg[i, ]
  pcol <- if (identical(as.character(z$p), prev)) "" else
    sprintf("%d ($s=%d$)", z$p, z$s)
  prev <- as.character(z$p)
  ## flag any cell whose chains did not meet the usual Rhat < 1.1 criterion:
  ## the numbers in that row describe the sampler, not the posterior
  flag <- if (z$rhat >= 1.05) "$^{\\dagger}$" else ""
  writeLines(sprintf(
    "%s & %s & %.1f & %.3f & %.3f & %.3f & %.3f & %.4f & %.1f & %.2f & %.2f (%.2f) & %.2f & %.2f%s & %.0f \\\\",
    pcol, lab[z$scenario], z$rho, z$me_map, z$me_med, z$me_sel,
    z$tpr, z$fpr, z$size, z$exact, z$cov, z$cov_se, z$w, z$rhat, flag, z$ess_tail), con)
}
writeLines(c("\\bottomrule", "\\end{tabular}}", "\\end{table}"), con)
close(con)
cat("\nwrote sim-table.tex and sim-results.rds\n")
print(agg, row.names = FALSE)
