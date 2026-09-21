## Full-data posterior fits behind the posterior and MCMC figures.
##   Rscript R/posterior-fits.R
## Boston and Riboflavin, four chains from the dispersed starts of Section 4.4.
## Writes posterior-fits.rds with the per-chain hyperparameter traces, the Boston
## draws and the per-coordinate rank diagnostics.
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
out <- file.path(here, "..")
## the dataset loaders live in realdata.R, above its main block
src <- readLines(file.path(here, "realdata.R"))
eval(parse(text = src[grep("^[.]libPaths", src)[1]:(grep("^## -+ main", src) - 1)]))

fit4 <- function(d) {
  t0 <- proc.time()[3]
  f <- bt_adenet(d$X, d$y, pilot = d$pilot, nsim = 4000, burn = 2000, chains = 4,
                 starts = "diverse", seed = 1, diagnostics = TRUE, verbose = FALSE)
  f$secs <- proc.time()[3] - t0
  f
}
## hyperparameter and L1-norm traces, one matrix column per chain
traces <- function(f) {
  list(lam1 = sapply(f$chains, function(o) o$lambda[, 1]),
       lam2 = sapply(f$chains, function(o) o$lambda[, 2]),
       l1 = sapply(f$chains, function(o) rowSums(abs(o$beta))))
}

b <- load_boston()
fb <- fit4(b)
mu <- sweep(b$X %*% t(fb$beta), 2, fb$alpha, "+")
res_b <- list(
  names = colnames(b$X), sd_x = apply(b$X, 2, stats::sd),
  beta = fb$beta, map = fb$map$beta, selected = summary(fb)$selected,
  ## posterior probability that observation i sits beyond the rejection radius
  p_reject = rowMeans(abs(b$y - mu) / fb$sigma > fb$d),
  fitted = apply(mu, 1, stats::median), y = b$y,
  trace = traces(fb), diag = fb$diag, w = fb$w, sigma = fb$sigma, secs = fb$secs)
cat(sprintf("Boston: w %.2f, max rank Rhat %.3f, %d selected, %.0fs\n", fb$w,
            max(fb$diag[, "rhat"], na.rm = TRUE), sum(res_b$selected), fb$secs))

r <- load_riboflavin()
fr <- fit4(r)
res_r <- list(trace = traces(fr), diag = fr$diag, w = fr$w, sigma = fr$sigma,
              map_size = sum(fr$map$beta != 0), selected = sum(summary(fr)$selected),
              secs = fr$secs)
cat(sprintf("Riboflavin: w %.2f, max rank Rhat %.3f, MAP %d, selected %d, %.0fs\n", fr$w,
            max(fr$diag[, "rhat"], na.rm = TRUE), res_r$map_size, res_r$selected, fr$secs))
cat("Riboflavin lambda1 median by chain:", signif(apply(res_r$trace$lam1, 2, median), 3), "\n")
cat("Riboflavin L1 median by chain:", signif(apply(res_r$trace$l1, 2, median), 3), "\n")

saveRDS(list(boston = res_b, ribo = res_r), file.path(out, "posterior-fits.rds"))
cat("wrote posterior-fits.rds\n")
