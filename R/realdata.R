## Real-data analysis for bayes-tadenet.tex.
##   Rscript R/realdata.R [nsplit_1 .. nsplit_4] [cores] [merge]
##   RD_PANEL=new Rscript R/realdata.R ...   second panel, written to realdata-new*.rds
##
## Protocol follows Section 7 of the frequentist paper (random 2:1 train/test
## splits, trimmed RMSPE, model size, Jaccard selection stability) so the two
## papers are comparable, and adds the measures only a posterior supplies:
## coverage and width of held-out predictive intervals.
##
## NOTE on the predictive interval.  A Gibbs posterior has no likelihood, so there
## is no posterior predictive distribution.  We report a PLUG-IN predictive: the
## posterior of the mean function x'beta convolved with the empirical distribution
## of the training residuals.  This is not a posterior predictive and is labelled
## as such in the manuscript.
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
.libPaths(c(file.path(here, "..", "Rlib"), .libPaths()))
source(file.path(here, "bt_adenet.R"))
out <- file.path(here, "..")
args <- commandArgs(TRUE)
panel <- Sys.getenv("RD_PANEL", "paper")
tag <- if (panel == "paper") "realdata" else "realdata-new"
## split counts scaled by measured cost per fit
nsplit <- if (panel == "paper")
  c("Boston housing" = 50, "Ames housing" = 25, "Glass" = 25, "Riboflavin" = 10) else
  c("Diabetes" = 50, "TopGear" = 50, "Octane" = 50, "NCI-60" = 25)
for (i in seq_len(min(4, length(args)))) nsplit[i] <- as.integer(args[i])

## ------------------------------------------------------------------ datasets
## Drop columns until the design (with intercept) has full column rank.  Real
## dummy-coded designs are routinely rank deficient -- Ames is deficient by 5 even
## under reference coding, because some factor levels are nested -- and an exactly
## unidentified direction is a flat ridge in the posterior that two chains explore
## differently, showing up as Rhat in the single digits.
full_rank <- function(X) {
  ## Only meaningful when p < n.  With p >= n the design is necessarily rank
  ## deficient and that is not a defect: Proposition 2.2 says the proper prior makes
  ## the posterior well defined regardless, and reducing to rank n-1 would throw away
  ## the p > n structure the analysis is about.
  if (ncol(X) >= nrow(X) - 1) return(X)
  q <- qr(cbind(1, X))
  keep <- sort(q$pivot[seq_len(q$rank)])
  keep <- keep[keep > 1] - 1                       # drop the intercept column
  X[, keep, drop = FALSE]
}

## The `black` column of MASS::Boston is 1000(Bk - 0.63)^2 with Bk the proportion
## of Black residents by town -- a covariate whose functional form encodes a racist
## premise, and the reason scikit-learn withdrew this dataset.  We drop it.
load_boston <- function() {
  b <- MASS::Boston
  X <- as.matrix(b[, setdiff(names(b), c("medv", "black"))])
  list(X = full_rank(X), y = log(b$medv),
       name = "Boston housing", pilot = "lmrob")
}

load_ames <- function() {
  a <- AmesHousing::make_ames()
  y <- log(a$Sale_Price); a$Sale_Price <- NULL
  X <- stats::model.matrix(~ ., data = a)[, -1, drop = FALSE]   # reference coding
  X <- X[, apply(X, 2, function(v) length(unique(v)) > 1 & mean(v != 0) > 0.02)]
  list(X = full_rank(X), y = y, name = "Ames housing", pilot = "irls")
}
## Archaeological glass (Lemberge et al. 2000).  The spectra and the chemical
## compositions ship in different packages, on the same 180 samples: cellWise has the
## 750 EPXMA frequencies, chemometrics the 13 oxide concentrations including PbO.
## Joining them reconstructs the design used by Alfons et al. (2013) and Kurnaz et al.
## (2018).  Preprocessing follows their rule: keep frequencies whose maximum intensity
## exceeds 0.1 per cent of the overall maximum.
load_glass <- function() {
  utils::data("data_glass", package = "cellWise", envir = environment())
  utils::data("glass", package = "chemometrics", envir = environment())
  S <- as.matrix(get("data_glass", envir = environment()))
  G <- as.matrix(get("glass", envir = environment()))
  stopifnot(nrow(S) == nrow(G))
  keep <- apply(S, 2, max) / max(S) > 0.001
  list(X = full_rank(S[, keep, drop = FALSE]), y = as.numeric(G[, "PbO"]),
       name = "Glass", pilot = "irls")
}

load_riboflavin <- function() {
  e <- new.env(); load(file.path(here, "..", "data", "riboflavin.RData"), e)
  r <- get("riboflavin", e)
  list(X = as.matrix(r$x), y = as.numeric(r$y),
       name = "Riboflavin", pilot = "irls")
}

## ------------------------------------------------ second panel (RD_PANEL=new)
## Each dataset asks one question the first panel leaves open.
## Diabetes (Efron, Hastie, Johnstone and Tibshirani 2004): the 10 baseline covariates
## with squares and interactions, p = 64.  Clean and n >> p: what does robustness cost?
load_diabetes <- function() {
  e <- new.env(); utils::data("diabetes", package = "lars", envir = e)
  list(X = full_rank(unclass(e$diabetes$x2)), y = e$diabetes$y,
       name = "Diabetes", pilot = "lmrob")
}
## TopGear (robustHD): log MPG on price, engine, body and equipment, p = 40 after dummy
## coding and dropping incomplete rows.  Three plug-in hybrids report 235 to 470 MPG
## under a test cycle that credits electric range, so they are outliers by construction.
load_topgear <- function() {
  e <- new.env(); utils::data("TopGear", package = "robustHD", envir = e)
  tg <- e$TopGear
  lab <- paste(tg$Maker, tg$Model)
  tg <- tg[, !names(tg) %in% c("Maker", "Model", "Type")]
  ok <- stats::complete.cases(tg)
  tg <- tg[ok, ]; lab <- lab[ok]
  y <- log(tg$MPG); tg$MPG <- NULL
  tg$Price <- log(tg$Price)
  X <- stats::model.matrix(~ ., data = tg)[, -1, drop = FALSE]
  list(X = full_rank(X), y = y, name = "TopGear", pilot = "lmrob", labels = lab,
       outliers = which(lab %in% c("BMW i3", "Chevrolet Volt", "Vauxhall Ampera")))
}
## Octane (rrcov): NIR spectra of 39 gasoline samples at 226 wavelengths.  Samples 25, 26
## and 36-39 contain added alcohol (Hubert, Rousseeuw and Vanden Branden 2005).  They are
## known outliers in the spectra, which is x, and need not be outlying in y.
load_octane <- function() {
  e <- new.env(); utils::data("octane", package = "rrcov", envir = e)
  list(X = as.matrix(e$octane[, -1]), y = e$octane$y, name = "Octane", pilot = "irls",
       outliers = c(25, 26, 36:39))
}
## NCI-60 (robustHD, the version of Alfons, Croux and Gelper 2013): protein KRT18, the
## protein with the largest MAD, on gene expression.  The 22283 genes are cut to the 2000
## with the largest variance, a filter that never sees y and so cannot leak test data.
load_nci60 <- function() {
  e <- new.env(); utils::data("nci60", package = "robustHD", envir = e)
  keep <- order(-apply(e$gene, 2, stats::var))[1:2000]
  X <- e$gene[, keep]
  colnames(X) <- make.unique(as.character(e$geneInfo$Symbol[keep]))
  list(X = X, y = as.numeric(e$protein[, 92]), name = "NCI-60", pilot = "irls")
}

## ------------------------------------------------------------------ measures
rtmspe <- function(r, trim = 0.9) {
  h <- floor(trim * length(r)); sqrt(mean(sort(r^2)[seq_len(h)]))
}
jaccard <- function(sets) {
  k <- length(sets); if (k < 2) return(NA_real_)
  v <- c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    u <- length(union(sets[[i]], sets[[j]]))
    ## a pair of empty selections says nothing about stability, so it is excluded
    ## (it used to count as 1, which made a method that selects nothing look stable)
    v <- c(v, if (u == 0) NA_real_ else length(intersect(sets[[i]], sets[[j]])) / u)
  }
  if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
}

## plug-in predictive interval (see header note)
pred_interval <- function(fit, Xte, res_train, level = 0.90, ndraw = 300) {
  idx <- sample(nrow(fit$beta), min(ndraw, nrow(fit$beta)))
  mu <- sweep(Xte %*% t(fit$beta[idx, , drop = FALSE]), 2, fit$alpha[idx], "+")
  e <- sample(res_train, length(idx), replace = TRUE)
  s <- sweep(mu, 2, e, "+")
  q <- c((1 - level) / 2, 1 - (1 - level) / 2)
  t(apply(s, 1, stats::quantile, probs = q))
}

## ------------------------------------------------------------------ one split
run_split <- function(d, tr, te, arm, seed, nsim, burn) {
  t0 <- proc.time()[3]
  f <- try(bt_adenet(d$X[tr, , drop = FALSE], d$y[tr], d = arm$d, prior = arm$prior,
                     pilot = d$pilot, nsim = nsim, burn = burn, chains = 2, seed = seed,
                     learn = if (is.null(arm$learn)) "calibrate" else arm$learn,
                     diagnostics = TRUE, verbose = FALSE), silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)
  ## Predict from the draws JOINTLY, then summarise.  Pairing the marginal median
  ## of alpha with the marginal medians of beta is not a coherent fit: alpha absorbs
  ## sum_j beta_j * xc_j / xs_j from the internal standardisation, so the two are
  ## strongly dependent once p is large.  On Ames (p = 176) that pairing put the
  ## intercept at -77 against a response mean of 12, leaving residuals with mean
  ## -1.74 and standard deviation 0.14 -- a fit correct up to a constant offset,
  ## and predictive intervals that covered nothing.  p = 12 on Boston hid it.
  fitted_med <- function(Xm) {
    mu <- sweep(Xm %*% t(f$beta), 2, f$alpha, "+")
    apply(mu, 1, stats::median)
  }
  res_tr <- d$y[tr] - fitted_med(d$X[tr, , drop = FALSE])
  r <- d$y[te] - fitted_med(d$X[te, , drop = FALSE])
  pi <- pred_interval(f, d$X[te, , drop = FALSE], res_tr)
  covd <- d$y[te] >= pi[, 1] & d$y[te] <= pi[, 2]
  bulk <- order(abs(r))[seq_len(floor(0.9 * length(r)))]   # the 90% best-fit test points
  sm <- summary(f)
  ## Second-panel measures.  They come after every random draw above, so the first
  ## panel's numbers do not move.  MAP is the frequentist Tukey adaptive elastic net,
  ## the lasso and the training median are non-robust and null references.
  xte <- d$X[te, , drop = FALSE]
  r_map <- d$y[te] - f$map$alpha - as.vector(xte %*% f$map$beta)
  rej <- abs(res_tr) / f$sigma > f$d
  tro <- tr %in% d$outliers
  lasso <- if (arm$id == "BT-AdEnet") {
    g <- glmnet::cv.glmnet(d$X[tr, , drop = FALSE], d$y[tr], nfolds = 5)
    rtmspe(d$y[te] - as.vector(stats::predict(g, xte, s = "lambda.min")))
  } else NA_real_
  list(row = data.frame(
      rtmspe = rtmspe(r), rmspe = sqrt(mean(r^2)),
      cover_all = mean(covd), cover_bulk = mean(covd[bulk]),
      width = stats::median(pi[, 2] - pi[, 1]),
      size = sum(sm$selected), rhat = stats::median(f$rhat, na.rm = TRUE),
      ## for p in the hundreds, max(Rhat) is set by a handful of coordinates in
      ## correlated dummy blocks; the fraction failing the criterion is the
      ## informative summary
      rhat_bad = mean(f$rhat > 1.1, na.rm = TRUE),
      rhat_rank = max(f$diag[, "rhat"], na.rm = TRUE),
      ess_bulk = min(f$diag[, "ess_bulk"], na.rm = TRUE),
      ess_tail = min(f$diag[, "ess_tail"], na.rm = TRUE),
      ## calibrate_w() returns w = 1 when the mode selects nothing
      map_size = sum(f$map$beta != 0), engine = f$pilot_engine,
      sigma = f$sigma, w = f$w,
      map_rtmspe = rtmspe(r_map), lasso_rtmspe = lasso,
      null_rtmspe = rtmspe(d$y[te] - stats::median(d$y[tr])),
      reject = mean(rej),
      out_flag = if (any(tro)) mean(rej[tro]) else NA_real_,
      clean_flag = if (length(d$outliers)) mean(rej[!tro]) else NA_real_,
      secs = proc.time()[3] - t0),
    sel = which(sm$selected))
}

## ---------------------------------------------------------------------- main
## Splits run in parallel.  Seeds are per split, so results do not depend on the
## worker, but wall-clock seconds are measured under contention and are not the
## cost figures reported in the paper.
cores <- if (length(args) > 4) as.integer(args[5]) else 16
## "merge" as a sixth argument keeps the datasets of the existing realdata.rds that
## this run skips, so one dataset can be rerun without refitting the others.
merge_old <- length(args) > 5 && args[6] == "merge" && file.exists(file.path(out, paste0(tag, ".rds")))
arms <- list(list(id = "BT-AdEnet", d = 4.685, prior = "adenet"),
             list(id = "Tukey+SSVS", d = 4.685, prior = "ssvs"),
             list(id = "Gaussian",   d = Inf,   prior = "adenet"))
## the second panel adds w = 1, to separate what calibration does from what the loss does
if (panel != "paper")
  arms <- c(arms, list(list(id = "BT w=1", d = 4.685, prior = "adenet", learn = 1)))
loaders <- if (panel == "paper") list(load_boston, load_ames, load_glass, load_riboflavin) else
  list(load_diabetes, load_topgear, load_octane, load_nci60)
res <- NULL
sets <- list()
full <- list()
for (ld in loaders) {
  d <- ld()
  n <- nrow(d$X)
  R <- nsplit[[d$name]]
  if (is.null(R) || is.na(R) || R < 1) next
  cat(sprintf("\n== %s: n=%d p=%d, %d splits, pilot=%s ==\n",
              d$name, n, ncol(d$X), R, d$pilot))
  tasks <- expand.grid(a = seq_along(arms), r = seq_len(R))
  out_d <- par_lapply(nrow(tasks), function(t) {
    aa <- arms[[tasks$a[t]]]
    r <- tasks$r[t]
    set.seed(7000 + r)
    tr <- sample(n, floor(2 * n / 3))
    te <- setdiff(seq_len(n), tr)
    o <- run_split(d, tr, te, aa, r, nsim = 4000, burn = 2000)
    if (is.null(o)) return(NULL)
    cat(sprintf("  %-10s split %2d/%d  rtmspe=%.3f size=%d (%.0fs)\n",
                aa$id, r, R, o$row$rtmspe, o$row$size, o$row$secs))
    list(row = cbind(data.frame(dataset = d$name, arm = aa$id, rep = r), o$row), sel = o$sel)
  }, cores)
  res <- rbind(res, do.call(rbind, lapply(out_d, `[[`, "row")))
  for (o in out_d) sets[[paste(o$row$dataset, o$row$arm, o$row$rep)]] <- o$sel
  if (panel != "paper") {
    ## one fit on all the data: which variables it keeps and which points it rejects
    f <- bt_adenet(d$X, d$y, pilot = d$pilot, nsim = 4000, burn = 2000, chains = 2,
                   seed = 1, diagnostics = TRUE, verbose = FALSE)
    mu <- apply(sweep(d$X %*% t(f$beta), 2, f$alpha, "+"), 1, stats::median)
    z <- (d$y - mu) / f$sigma
    nm <- colnames(d$X)
    if (is.null(nm)) nm <- paste0("x", seq_len(ncol(d$X)))
    full[[d$name]] <- list(selected = nm[summary(f)$selected], map = nm[f$map$beta != 0],
      z = z, rejected = which(abs(z) > f$d), outliers = d$outliers, labels = d$labels,
      w = f$w, sigma = f$sigma, rhat = max(f$diag[, "rhat"], na.rm = TRUE))
    cat(sprintf("  full fit: %d selected, %d at MAP, rejected rows %s\n",
                length(full[[d$name]]$selected), length(full[[d$name]]$map),
                paste(full[[d$name]]$rejected, collapse = " ")))
  }
}
if (merge_old) {
  old <- readRDS(file.path(out, paste0(tag, ".rds")))
  redone <- unique(res$dataset)
  res <- rbind(old$res[!old$res$dataset %in% redone, ], res)
  keep <- !vapply(names(old$sets), function(k) any(startsWith(k, paste0(redone, " "))), TRUE)
  sets <- c(old$sets[keep], sets)
  full <- c(old$full[setdiff(names(old$full), redone)], full)
}
saveRDS(list(res = res, sets = sets, full = full), file.path(out, paste0(tag, ".rds")))
agg <- do.call(rbind, lapply(split(res, list(res$dataset, res$arm), drop = TRUE),
  function(z) data.frame(dataset = z$dataset[1], arm = z$arm[1], R = nrow(z),
    rtmspe = median(z$rtmspe), rmspe = median(z$rmspe),
    cover_all = mean(z$cover_all), cover_bulk = mean(z$cover_bulk),
    width = median(z$width), size = median(z$size), null = sum(z$size == 0),
    sigma = median(z$sigma), w = median(z$w), w_fallback = sum(z$map_size == 0),
    rhat = median(z$rhat_rank), ess_bulk = median(z$ess_bulk),
    ess_tail = median(z$ess_tail), rhat_bad = mean(z$rhat_bad), secs = median(z$secs),
    stability = jaccard(sets[paste(z$dataset, z$arm, z$rep)]),
    map_rtmspe = if (is.null(z$map_rtmspe)) NA else median(z$map_rtmspe),
    lasso_rtmspe = if (is.null(z$lasso_rtmspe)) NA else median(z$lasso_rtmspe, na.rm = TRUE),
    null_rtmspe = if (is.null(z$null_rtmspe)) NA else median(z$null_rtmspe),
    reject = if (is.null(z$reject)) NA else median(z$reject),
    out_flag = if (is.null(z$out_flag)) NA else mean(z$out_flag, na.rm = TRUE),
    clean_flag = if (is.null(z$clean_flag)) NA else mean(z$clean_flag, na.rm = TRUE))))
print(agg, row.names = FALSE)
saveRDS(agg, file.path(out, paste0(tag, "-agg.rds")))
cat(sprintf("\nwrote %s.rds\n", tag))
