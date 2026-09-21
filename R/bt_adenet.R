## Bayesian Tukey Adaptive Elasticnet (BT-AdEnet)
##
## Gibbs posterior built on the Tukey biweight loss of the T-AdEnet paper,
## with the adaptive elastic-net penalty represented as a proper scale-mixture
## prior.  MAP of this posterior == the frequentist T-AdEnet estimator.
##
## Base R only (MASS optional).  See ../bayes-tadenet.tex for the derivation.

## ---------------------------------------------------------------- Tukey loss

## d = Inf selects the squared-error loss.  This is not a curiosity: it is the
## non-robust comparator used in the sensitivity experiment, obtained through the
## same sampler so that any difference is attributable to the loss alone.  A large
## finite d will NOT do -- once |r|/sigma exceeds it the "Gaussian" arm redescends
## too, which is precisely the regime the experiment probes.
tukey_rho <- function(u, d = 4.685) {
  if (is.infinite(d)) return(u^2 / 2)
  s <- (u / d)^2
  ifelse(s <= 1, (d^2 / 6) * (1 - (1 - s)^3), d^2 / 6)
}

## psi = rho'
tukey_psi <- function(u, d = 4.685) {
  if (is.infinite(d)) return(u)
  s <- (u / d)^2
  ifelse(s <= 1, u * (1 - s)^2, 0)
}

## IRLS weight psi(u)/u -- always >= 0, used for the Gauss-Newton preconditioner
tukey_w <- function(u, d = 4.685) {
  if (is.infinite(d)) return(rep(1, length(u)))
  s <- (u / d)^2
  ifelse(s <= 1, (1 - s)^2, 0)
}

## psi' -- used only for the sandwich calibration (can be negative: rho is non-convex)
tukey_psi1 <- function(u, d = 4.685) {
  if (is.infinite(d)) return(rep(1, length(u)))
  s <- (u / d)^2
  ifelse(s <= 1, (1 - s) * (1 - 5 * s), 0)
}

## ------------------------------------------------------- robust initial fit

## MM-style pilot.  Uses robustbase::lmrob when available (that is what the
## paper uses); otherwise falls back to a self-contained ridge -> Huber -> Tukey
## IRLS.  The ridge amount is the calibration knob for p >= n.
## pilot = "auto" uses robustbase::lmrob when it is installed and p < n/2, and the
## built-in IRLS otherwise.  Which one ran materially affects the adaptive weights,
## so the choice is recorded in the fit object rather than left implicit: installing
## robustbase silently changes results otherwise.
init_robust <- function(x, y, d = 4.685, ridge = NULL, maxit = 200,
                        pilot = c("auto", "lmrob", "irls"), fam = "tukey") {
  pilot <- match.arg(pilot)
  n <- nrow(x); p <- ncol(x)
  ## In p >= n a weak ridge interpolates the data and MAD(residual) collapses to
  ## zero, which then poisons sigma-hat and the adaptive weights.  Pick the ridge
  ## by targeting an effective degrees of freedom, and correct the scale for it.
  edf <- NA_real_
  if (is.null(ridge)) {
    if (p >= n / 2) {
      sv2 <- svd(x, nu = 0, nv = 0)$d^2
      target <- n / 4
      lo <- log(1e-6 * sum(sv2)); hi <- log(1e6 * sum(sv2))
      for (b in 1:60) {
        mid <- (lo + hi) / 2
        if (sum(sv2 / (sv2 + exp(mid))) > target) lo <- mid else hi <- mid
      }
      ridge <- exp((lo + hi) / 2)
      edf <- sum(sv2 / (sv2 + ridge))
    } else {
      ridge <- 1e-6 * n
      edf <- p
    }
  }
  df_corr <- sqrt(max(1 - edf / n, 0.25))

  ## An explicit request for the MM pilot must not degrade silently.  Before this
  ## check, every simulation script that asked for "lmrob" without robustbase on
  ## the library path ran the IRLS pilot instead.
  if (pilot == "lmrob" && p < n / 2 && !requireNamespace("robustbase", quietly = TRUE))
    stop("pilot = \"lmrob\" requested but robustbase is not on .libPaths()")
  if (pilot != "irls" && p < n / 2 && requireNamespace("robustbase", quietly = TRUE)) {
    fit <- try(robustbase::lmrob(y ~ x, method = "MM"), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      cf <- stats::coef(fit)
      return(list(alpha = unname(cf[1]), beta = unname(cf[-1]),
                  sigma = fit$scale, edf = p, ridge = 0, engine = "lmrob"))
    }
  }

  ## ridge start, then IRLS: Huber (convex, gets us close) then Tukey.
  ## When p > n the ridge normal equations are solved through the n x n
  ## push-through identity (X'OX + rI)^-1 X'O z = Xt'(Xt Xt' + rI)^-1 O^(1/2) z with
  ## Xt = O^(1/2) X.  At p = 4088 the p x p solve is O(p^3) per IRLS step and was
  ## most of the cost of a fit.  The p <= n branch keeps the original expressions.
  push <- function(om, z) {
    xt <- x * sqrt(om)
    as.vector(crossprod(xt, solve(tcrossprod(xt) + ridge * diag(n), sqrt(om) * z)))
  }
  alpha <- stats::median(y)
  if (p <= n) {
    xtx <- crossprod(x)
    xty <- crossprod(x, y)
    beta <- as.vector(solve(xtx + ridge * diag(p), xty - alpha * colSums(x)))
  } else {
    beta <- push(rep(1, n), y - alpha)
  }
  r <- y - alpha - x %*% beta
  sigma <- max(stats::mad(r), 1e-8) / df_corr

  irls <- function(beta, alpha, wfun, iters) {
    for (it in seq_len(iters)) {
      r <- as.vector(y - alpha - x %*% beta)
      sigma <<- max(stats::mad(r), 1e-8) / df_corr
      om <- wfun(r / sigma)
      om <- pmax(om, 1e-6)                     # keep the weighted normal eqs PD
      alpha <- sum(om * (y - x %*% beta)) / sum(om)
      xw <- x * om
      beta <- if (p <= n) as.vector(solve(crossprod(x, xw) + ridge * diag(p),
                                          crossprod(xw, y - alpha)))
              else push(om, y - alpha)
    }
    list(alpha = alpha, beta = beta)
  }
  hw <- function(u, k = 1.345) pmin(1, k / pmax(abs(u), 1e-12))
  fit <- irls(beta, alpha, hw, min(50, maxit))
  fit <- irls(fit$beta, fit$alpha, function(u) w_(u, d, fam), min(150, maxit))
  r <- as.vector(y - fit$alpha - x %*% fit$beta)
  list(alpha = fit$alpha, beta = fit$beta,
       sigma = max(stats::mad(r), 1e-8) / df_corr, edf = edf, ridge = ridge,
       engine = "irls")
}

## ---------------------------------------- frequentist T-AdEnet (Algorithm 1)

## Proximal AdaGrad for
##   sum_i rho_d((y_i - alpha - x_i'b)/sigma) + lam1 sum_j w_j|b_j| + lam2/2 ||b||^2
## AdaGrad's cumulative travel along a coordinate with a consistent gradient sign
## is O(eta*sqrt(maxit)), so a fixed eta = 0.01 needs ~1e5 iterations to reach a
## coefficient of size 3 and silently under-converges at the 1e4 cap.  We set
## eta so that the reachable range covers the pilot coefficients.
tadenet_fit <- function(x, y, lam1, lam2, w, sigma, d = 4.685,
                        beta0 = NULL, alpha0 = NULL, fam = "tukey",
                        eta = NULL, eps = 1e-8, maxit = 1e4, tol = 1e-6) {
  n <- nrow(x); p <- ncol(x)
  beta <- if (is.null(beta0)) numeric(p) else beta0
  alpha <- if (is.null(alpha0)) stats::median(y) else alpha0
  if (is.null(eta))
    eta <- max(abs(beta), stats::mad(y) / 2, 1) / sqrt(maxit)
  G <- numeric(p); Ga <- 0
  obj_prev <- Inf
  for (l in seq_len(maxit)) {
    r <- as.vector(y - alpha - x %*% beta)
    ps <- psi_(r / sigma, d, fam)
    g <- -as.vector(crossprod(x, ps)) / sigma
    ga <- -sum(ps) / sigma
    G <- G + g^2; Ga <- Ga + ga^2
    st <- eta / (sqrt(G) + eps)
    u <- beta - st * g
    beta <- sign(u) * pmax(abs(u) - st * lam1 * w, 0) / (1 + st * lam2)
    alpha <- alpha - (eta / (sqrt(Ga) + eps)) * ga
    obj <- sum(rho_(as.vector(y - alpha - x %*% beta) / sigma, d, fam))
    if (abs(obj - obj_prev) / n < tol) break
    obj_prev <- obj
  }
  r <- as.vector(y - alpha - x %*% beta)
  list(alpha = alpha, beta = beta, iters = l,
       loss = sum(rho_(r / sigma, d, fam)),
       objective = sum(rho_(r / sigma, d, fam)) +
         lam1 * sum(w * abs(beta)) + lam2 / 2 * sum(beta^2))
}

## ---------------------------------------------------- pseudo-Huber family
## The convex comparator.  Kawakami and Hashimoto (2022) build a Bayesian Huberized
## lasso on this loss; it is smooth, convex, and has bounded psi, so the posterior is
## log-concave.  Note rho is UNBOUNDED (it grows like d|u|), so neither Prop 2.3 nor
## Thm 6.2 applies to it -- which is exactly the comparison of interest.
ph_rho  <- function(u, d) d^2 * (sqrt(1 + (u / d)^2) - 1)
ph_psi  <- function(u, d) u / sqrt(1 + (u / d)^2)
ph_w    <- function(u, d) 1 / sqrt(1 + (u / d)^2)
ph_psi1 <- function(u, d) (1 + (u / d)^2)^(-1.5)

## Family dispatchers.  Default "tukey" so the S-scale of s_scale() below, which must
## stay Tukey in every arm, is unaffected.
rho_  <- function(u, d, fam = "tukey") if (fam == "huber") ph_rho(u, d)  else tukey_rho(u, d)
psi_  <- function(u, d, fam = "tukey") if (fam == "huber") ph_psi(u, d)  else tukey_psi(u, d)
w_    <- function(u, d, fam = "tukey") if (fam == "huber") ph_w(u, d)    else tukey_w(u, d)
psi1_ <- function(u, d, fam = "tukey") if (fam == "huber") ph_psi1(u, d) else tukey_psi1(u, d)

## Tuning constant for the SCALE, not for the regression loss.  As in
## MM-estimation (Yohai 1987) the two are different: d = 4.685 gives the loss 95%
## Gaussian efficiency, while d0 = 1.5476 is the value at which
## E_Phi[rho_d0] = (1/2) max rho_d0, i.e. the S-scale has 50% breakdown and is
## still Fisher consistent at the normal.  Using 4.685 for both puts the scale
## breakdown at 12%, and a 10% contaminated sample then inflates sigma-hat
## by a factor of two.
D_SCALE <- 1.5476

## Fisher-consistency constant for the S-scale: kappa = E_Phi[rho_d(Z)]
tukey_kappa <- function(d = D_SCALE)
  stats::integrate(function(u) tukey_rho(u, d) * stats::dnorm(u), -Inf, Inf)$value


## S-scale of a residual vector: the unique root of mean(rho_d0(r/s)) = kappa.
## mean(rho) decreases strictly in s from d0^2/6 to 0, so the root is unique.
s_scale <- function(r, d = D_SCALE, kappa = tukey_kappa(d), df = 0, n = length(r),
                    floor = 0) {
  g <- function(s) mean(tukey_rho(r / s, d)) - kappa
  lo <- 1e-8; hi <- max(stats::mad(r), stats::sd(r), 1e-8)
  while (g(hi) > 0 && hi < 1e8) hi <- hi * 2
  s <- if (g(hi) > 0) hi else stats::uniroot(g, c(lo, hi), tol = 1e-10)$root
  ## p >> n lets a dense fit interpolate, driving the S-scale to zero; without a
  ## floor the next lambda_max overflows and the path dies.
  max(s / sqrt(max(1 - df / n, 0.25)), floor)
}

## RBIC path over lam1 with a concomitant scale.  Holding sigma fixed at a
## pilot value fails when p >> n: the pilot underfits, sigma-hat is inflated,
## the penalty is then too strong relative to the loss, and the path is trapped
## at an over-sparse fixed point.  Re-estimating sigma at every lambda (as in
## scaled sparse regression, Sun & Zhang 2012) breaks that loop, and the model
## selection criterion becomes the classical BIC on the S-scale,
##      RBIC(lambda) = 2 n log sigma-hat(lambda) + log(n) |A(lambda)|.
tadenet_path <- function(x, y, w, sigma, lam2 = 1e-2, d = 4.685,
                         nlam = if (ncol(x) > 200) 20L else 30L,
                         lam_min_ratio = 1e-3,
                         beta0 = NULL, alpha0 = NULL, concomitant = TRUE,
                         inner = 2, dfmax = NULL, path_maxit = 2000,
                         fam = "tukey", ...) {
  n <- nrow(x); kap <- tukey_kappa(D_SCALE)
  ## A concomitant scale rewards interpolation: with p >> n and no cap on the
  ## model size, lambda -> 0 drives sigma-hat -> 0 and 2n log sigma-hat -> -Inf
  ## faster than log(n)|A| grows, so RBIC is unbounded below.  Capping the model
  ## size at n/2 (as glmnet's dfmax does) restores identifiability.
  if (is.null(dfmax)) dfmax <- floor(n / 2)
  sfloor <- 1e-3 * max(stats::mad(y), 1e-8)
  sigma <- max(sigma, sfloor)
  ## lambda_max is the smallest penalty giving the null fit, so it must be computed
  ## from the residuals of the INTERCEPT-ONLY model.  Using y itself silently fails
  ## whenever the response has a non-zero location: for log(medv) on the Boston data
  ## y/sigma ~ 28 > d, every observation is rejected, psi == 0 and lambda_max == 0.
  a0 <- if (is.null(alpha0)) stats::median(y) else alpha0
  r0 <- y - a0 - if (is.null(beta0)) 0 else as.vector(x %*% beta0)
  lmax <- max(abs(crossprod(x, psi_(r0 / sigma, d, fam))) / (sigma * w))
  if (!is.finite(lmax) || lmax <= 0)
    lmax <- max(abs(crossprod(x, y - a0))) / sigma
  grid <- exp(seq(log(lmax), log(lam_min_ratio * lmax), length.out = nlam))
  best <- NULL; b <- beta0; a <- alpha0; sg <- sigma
  for (lam1 in grid) {
    for (k in seq_len(if (concomitant) inner else 1L)) {
      f <- tadenet_fit(x, y, lam1, lam2, w, sg, d, beta0 = b, alpha0 = a,
                       fam = fam, maxit = path_maxit)
      b <- f$beta; a <- f$alpha
      if (!concomitant) break
      r <- as.vector(y - a - x %*% b)
      sg_new <- s_scale(r, D_SCALE, kap, df = sum(b != 0), n = n, floor = sfloor)
      if (abs(sg_new - sg) / sg < 1e-3) { sg <- sg_new; break }
      sg <- sg_new
    }
    df <- sum(b != 0)
    rbic <- if (concomitant) 2 * n * log(sg) + log(n) * df
            else 2 * f$loss + log(n) * df
    if (df <= dfmax && (is.null(best) || rbic < best$rbic))
      best <- list(rbic = rbic, lam1 = lam1, lam2 = lam2, sigma = sg,
                   beta = b, alpha = a, loss = f$loss)
    if (df > dfmax) break
  }
  if (is.null(best))
    stop("no fit on the path satisfied the dfmax cap; raise dfmax or lam_min_ratio")
  ## full-precision refit at the selected lambda.  The refit can activate more
  ## coordinates than the path did, and it must not take the model past dfmax: at
  ## p > n an oversized support makes J singular, the calibrated learning rate
  ## collapses toward zero, and the posterior reduces to its nearly flat prior.
  f <- tadenet_fit(x, y, best$lam1, lam2, w, best$sigma, d,
                   beta0 = best$beta, alpha0 = best$alpha, fam = fam, ...)
  best$capped <- sum(f$beta != 0) > dfmax
  if (!best$capped) {
    best$beta <- f$beta
    best$alpha <- f$alpha
    best$loss <- f$loss
  }
  best$sigma <- s_scale(as.vector(y - f$alpha - x %*% f$beta), D_SCALE, kap,
                        df = sum(f$beta != 0), n = n, floor = sfloor)
  best
}

## ------------------------------------------------- learning-rate calibration

## Trace-matched Gibbs-posterior learning rate.
##   method "cov" (default): match traces of the COVARIANCES,
##       w = tr(J^-1) / tr(J^-1 I J^-1)
##     -- targets average credible-interval width, which is what we report.
##   method "lhw": Lyddon, Holmes & Walker (2019, Lemma 1) match traces of the
##     PRECISIONS (the Fisher information number),
##       w = tr(J I^-1 J') / tr(J)
## The two agree whenever I is proportional to J (e.g. at the Gaussian model) and
## differ under contamination.
calibrate_w <- function(x, y, alpha, beta, sigma, d = 4.685, support = NULL,
                        method = c("cov", "lhw"), fam = "tukey") {
  method <- match.arg(method)
  if (is.null(support)) support <- which(beta != 0)
  if (!length(support)) return(1)
  xs <- x[, support, drop = FALSE]
  u <- as.vector(y - alpha - x %*% beta) / sigma
  n <- nrow(x)
  I <- crossprod(xs * psi_(u, d, fam)) / (n * sigma^2)
  J <- crossprod(xs, xs * psi1_(u, d, fam)) / (n * sigma^2)
  ev <- eigen((J + t(J)) / 2, symmetric = TRUE, only.values = TRUE)$values
  ## Non-convexity can make the empirical J indefinite.  Shifting it by its most
  ## negative eigenvalue sends that direction's eigenvalue to zero and the
  ## calibrated rate with it (w near 1e-8 on a riboflavin split), so the posterior
  ## reverts to its nearly flat prior.  Calibrate instead on the eigen-directions
  ## with positive eigenvalue, where the loss is locally convex.  A positive
  ## definite J takes the unchanged branch below.
  if (min(ev) <= 1e-8) {
    es <- eigen((J + t(J)) / 2, symmetric = TRUE)
    keep <- es$values > 1e-8 * max(es$values)
    if (!any(keep)) return(1)
    V <- es$vectors[, keep, drop = FALSE]
    lam <- es$values[keep]
    Ip <- crossprod(V, I %*% V)
    wt <- if (method == "cov") sum(1 / lam) / sum(diag(Ip) / lam^2) else {
      Ipi <- tryCatch(solve(Ip), error = function(e) NULL)
      if (is.null(Ipi)) return(1)
      sum(lam^2 * diag(Ipi)) / sum(lam)
    }
    return(if (!is.finite(wt) || wt <= 0) 1 else wt)
  }
  Ji <- tryCatch(solve(J), error = function(e) NULL)
  if (is.null(Ji)) return(1)
  wt <- if (method == "cov") {
    sum(diag(Ji)) / sum(diag(Ji %*% I %*% Ji))
  } else {
    Ii <- tryCatch(solve(I), error = function(e) NULL)
    if (is.null(Ii)) return(1)
    sum(diag(J %*% Ii %*% t(J))) / sum(diag(J))
  }
  if (!is.finite(wt) || wt <= 0) 1 else wt
}

## ----------------------------------------------------------------- samplers

## Inverse Gaussian (Michael, Schucany & Haas 1976), vectorised
rinvgauss <- function(n, mu, lambda) {
  nu <- stats::rnorm(n)^2
  xx <- mu + mu^2 * nu / (2 * lambda) -
    (mu / (2 * lambda)) * sqrt(4 * mu * lambda * nu + mu^2 * nu^2)
  z <- stats::runif(n)
  ifelse(z <= mu / (mu + xx), xx, mu^2 / xx)
}

## log normaliser of  exp(-a|b| - c b^2/2)  over the real line
log_enet_norm <- function(a, c) {
  s <- sqrt(c)
  log(2) + 0.5 * log(2 * pi / c) + a^2 / (2 * c) +
    stats::pnorm(-a / s, log.p = TRUE)
}

## ------------------------------------------------------------------- BT-AdEnet

#' Bayesian Tukey Adaptive Elasticnet
#'
#' @param x,y      design matrix (no intercept column) and response
#' @param d        Tukey tuning constant
#' @param iota     adaptive-weight exponent, w_j = |beta_j^MM|^(-iota)
#' @param learn    Gibbs-posterior learning rate: "calibrate", or a number
#' @param lambda   NULL to sample (lam1, lam2); or c(lam1, lam2) to hold fixed
#' @param nsim,burn,thin,chains  MCMC settings
bt_adenet <- function(x, y, d = 4.685, iota = 1,
                      learn = "calibrate", lambda = NULL,
                      nsim = 4000, burn = 2000, thin = 1, chains = 2,
                      a_lam = 1, b_lam = 1e-2, seed = NULL, w_floor = NULL,
                      pilot = c("auto", "lmrob", "irls"), fam = c("tukey", "huber"),
                      prior = c("adenet", "horseshoe", "ssvs"),
                      standardize = TRUE, starts = c("local", "diverse"),
                      diagnostics = FALSE, verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  x <- as.matrix(x); y <- as.vector(y)
  n <- nrow(x); p <- ncol(x)
  stopifnot(length(y) == n)

  ## --- standardise (robustly: median / mad, so contamination cannot set the scale)
  if (standardize) {
    xc <- apply(x, 2, stats::median)
    xs <- apply(x, 2, function(v) {
      s <- stats::mad(v)                       # binary / heavily-tied columns give 0
      if (s > 0) return(s)
      s <- stats::sd(v); if (s > 0) s else 1
    })
  } else {
    xc <- numeric(p); xs <- rep(1, p)
  }
  X <- sweep(sweep(x, 2, xc, "-"), 2, xs, "/")

  ## --- pilot.  The dense pilot underfits badly when p >= n, which inflates
  ## MAD(residual); one sparse refinement pass fixes the scale and the weights
  ## before they are handed to the sampler.
  fam <- match.arg(fam); prior <- match.arg(prior)
  starts <- match.arg(starts)
  pilot <- init_robust(X, y, d, pilot = match.arg(pilot), fam = fam)
  sigma <- pilot$sigma
  ## The floor sets W = min_{j notin A} w_j, which Corollary 6.4 requires to grow
  ## with p.  The default sigma/sqrt(n) does not depend on p at all, so the leakage
  ## term is only suppressed by a constant factor; w_floor = "p" makes the floor
  ## p-dependent, as the corollary asks.
  fl <- if (is.null(w_floor)) sigma / sqrt(n)
        else if (identical(w_floor, "p")) sigma / (sqrt(n) * p^(1 / iota))
        else as.numeric(w_floor)
  mk_w <- function(b, sig) {
    ww <- 1 / pmax(abs(b), fl)^iota
    ww <- ww / mean(ww)
    pmax(ww, 1e-6)          # bounded dynamic range keeps lambda_max finite
  }
  wj <- mk_w(pilot$beta, sigma)
  map <- tadenet_path(X, y, wj, sigma, d = d, fam = fam,
                      beta0 = pilot$beta, alpha0 = pilot$alpha)
  capped_any <- isTRUE(map$capped)
  ## one adaptive re-weighting pass: the pilot weights come from a dense fit,
  ## the second set from the sparse fit the path just produced
  for (pass in 1:2) {
    wj_new <- mk_w(map$beta, map$sigma)
    cand <- tadenet_path(X, y, wj_new, map$sigma, d = d, fam = fam,
                         beta0 = map$beta, alpha0 = map$alpha)
    capped_any <- capped_any || isTRUE(cand$capped)
    if (cand$rbic >= map$rbic) break          # keep the better of the two passes
    map <- cand; wj <- wj_new
  }
  sigma <- map$sigma
  if (verbose) message(sprintf(
    "pilot: engine=%s sigma=%.3f | MAP: lam1=%.4g support=%d",
    pilot$engine, sigma, map$lam1, sum(map$beta != 0)))

  ## --- learning rate
  wlr <- if (identical(learn, "calibrate"))
    calibrate_w(X, y, map$alpha, map$beta, sigma, d, method = "cov", fam = fam)
  else if (identical(learn, "lhw"))
    calibrate_w(X, y, map$alpha, map$beta, sigma, d, method = "lhw", fam = fam)
  else as.numeric(learn)
  if (verbose) message(sprintf("learning rate w = %.4f", wlr))

  ## Gauss-Newton preconditioner pieces, frozen at the MAP (Sec. 4.3 of the tex)
  u0 <- as.vector(y - map$alpha - X %*% map$beta) / sigma
  om0 <- pmax(w_(u0, d, fam), 1e-4)
  Sj <- as.vector(crossprod(X^2, om0)) * wlr / sigma^2
  Sa <- sum(om0) * wlr / sigma^2

  ## Proper prior on the intercept, uniform on [-K, K].  Under a flat prior the
  ## posterior is improper: once every residual is past d the bounded loss stops
  ## penalising alpha, and the integrand tends to a positive constant.  K is far
  ## outside anything the chain can reach, so seeded draws equal those of the
  ## flat-prior sampler unless a proposal leaves the interval.
  a_half <- 1e6 * (1 + max(abs(y)) + abs(map$alpha))

  ## A ridge least-squares start for starts = "diverse".  It sits in the basin a
  ## contaminated sample creates if one exists, which local jitter cannot reach.
  ls_start <- function() {
    a0 <- stats::median(y)
    lam <- 1e-3 * sum(X^2) / max(n, p)
    b0 <- if (p <= n) solve(crossprod(X) + lam * diag(p), crossprod(X, y - a0))
          else crossprod(X, solve(tcrossprod(X) + lam * diag(n), y - a0))
    list(alpha = a0, beta = as.vector(b0))
  }

  fixed_lam <- !is.null(lambda)
  keep <- length(seq(burn + 1, nsim, by = thin))
  out <- vector("list", chains)

  for (ch in seq_len(chains)) {
    beta <- map$beta; alpha <- map$alpha
    ## Dispersed start for the Rhat diagnostic.  The perturbation must be measured
    ## in posterior standard deviations, not in units of sigma: a fixed per-
    ## coordinate jitter has norm growing like sqrt(p), so at p = 752 it throws
    ## chain 2 clear of the mode and, the target being non-convex, it never
    ## returns -- which shows up as Rhat in the single digits rather than as
    ## honest multimodality.
    if (ch > 1) {
      post_sd <- 1 / sqrt(Sj + 1 + wlr * max(map$lam2, 1e-4))
      beta <- beta + 2 * post_sd * stats::rnorm(p)
    }
    if (starts == "diverse" && ch %% 4 == 3) {
      beta <- pilot$beta
      alpha <- pilot$alpha
    }
    if (starts == "diverse" && ch %% 4 == 0) {
      s0 <- ls_start()
      beta <- s0$beta
      alpha <- s0$alpha
    }
    lam1 <- if (fixed_lam) lambda[1] else max(map$lam1, 1e-4)
    lam2 <- if (fixed_lam) lambda[2] else max(map$lam2, 1e-4)
    tau2 <- rep(1, p)
    ## horseshoe state (Makalic and Schmidt 2016 inverse-gamma augmentation)
    hs_lam2 <- rep(1, p); hs_nu <- rep(1, p); hs_tau2 <- 1; hs_xi <- 1
    ## SSVS state (George and McCulloch continuous spike and slab).  gamma_j is a
    ## genuine inclusion indicator, so this is the one prior here that yields
    ## posterior inclusion probabilities rather than requiring a summary rule.
    ss_g <- rep(1, p); ss_pi <- 0.1
    ss_spike <- (sigma / (10 * sqrt(n)))^2; ss_slab <- (2 * stats::mad(y))^2
    G <- matrix(0, 0, p)
    step <- 1e-2; acc <- 0; acc_h <- 0; sd_h <- 0.2

    B <- matrix(NA_real_, keep, p); A <- numeric(keep); L <- matrix(NA_real_, keep, 2)
    k <- 0L

    logpost_beta <- function(alpha, beta, prec) {
      if (abs(alpha) > a_half) return(-Inf)
      r <- as.vector(y - alpha - X %*% beta)
      -wlr * sum(rho_(r / sigma, d, fam)) - 0.5 * sum(prec * beta^2)
    }
    grad_beta <- function(alpha, beta, prec) {
      r <- as.vector(y - alpha - X %*% beta)
      ps <- psi_(r / sigma, d, fam)
      list(b = wlr * as.vector(crossprod(X, ps)) / sigma - prec * beta,
           a = wlr * sum(ps) / sigma)
    }

    for (it in seq_len(nsim)) {
      ## ---- 1. beta, alpha | tau2, lam2 : preconditioned MALA
      cc <- if (prior == "adenet") wlr * lam2 else 0
      prec <- switch(prior,
        adenet    = 1 / tau2 + cc,
        horseshoe = 1 / (hs_tau2 * hs_lam2),
        ssvs      = 1 / ifelse(ss_g == 1, ss_slab, ss_spike))
      Pb <- Sj + prec; Pa <- Sa                       # diagonal mass matrix
      g <- grad_beta(alpha, beta, prec)
      mb <- beta + 0.5 * step^2 * g$b / Pb
      ma <- alpha + 0.5 * step^2 * g$a / Pa
      bp <- mb + step * stats::rnorm(p) / sqrt(Pb)
      ap <- ma + step * stats::rnorm(1) / sqrt(Pa)
      gp <- grad_beta(ap, bp, prec)
      mb2 <- bp + 0.5 * step^2 * gp$b / Pb
      ma2 <- ap + 0.5 * step^2 * gp$a / Pa
      lq_fwd <- -sum(Pb * (bp - mb)^2) / (2 * step^2) - Pa * (ap - ma)^2 / (2 * step^2)
      lq_bwd <- -sum(Pb * (beta - mb2)^2) / (2 * step^2) - Pa * (alpha - ma2)^2 / (2 * step^2)
      lacc <- logpost_beta(ap, bp, prec) - logpost_beta(alpha, beta, prec) +
        lq_bwd - lq_fwd
      if (is.finite(lacc) && log(stats::runif(1)) < lacc) {
        beta <- bp; alpha <- ap; acc <- acc + 1
      }
      if (it <= burn) step <- step * exp((min(exp(lacc), 1) - 0.574) / sqrt(it))

      if (prior == "adenet") {
        ## ---- 2. tau2 | beta, lam1 : 1/tau_j^2 ~ InvGauss(a_j/|beta_j|, a_j^2)
        aj <- wlr * lam1 * wj
        mu <- pmin(aj / pmax(abs(beta), 1e-8), 1e8)
        tau2 <- 1 / pmax(rinvgauss(p, mu, aj^2), 1e-12)
      } else if (prior == "horseshoe") {
        ## ---- 2'. horseshoe: beta_j | lam_j, tau ~ N(0, tau^2 lam_j^2) with
        ## half-Cauchy hyperpriors, sampled through the auxiliary inverse-gamma
        ## representation, which is conditionally conjugate throughout.
        rig <- function(n, shape, rate) 1 / stats::rgamma(n, shape, rate = rate)
        hs_lam2 <- rig(p, 1, 1 / hs_nu + beta^2 / (2 * hs_tau2))
        hs_nu   <- rig(p, 1, 1 + 1 / hs_lam2)
        hs_tau2 <- rig(1, (p + 1) / 2, 1 / hs_xi + sum(beta^2 / hs_lam2) / 2)
        hs_xi   <- rig(1, 1, 1 + 1 / hs_tau2)
        hs_lam2 <- pmax(hs_lam2, 1e-12); hs_tau2 <- max(hs_tau2, 1e-12)
      } else {
        ## ---- 2''. SSVS: gamma_j | beta_j is Bernoulli with odds given by the ratio
        ## of slab to spike density at beta_j, and pi carries a Beta(1,1) prior.
        l1_ <- stats::dnorm(beta, 0, sqrt(ss_slab),  log = TRUE) + log(ss_pi)
        l0_ <- stats::dnorm(beta, 0, sqrt(ss_spike), log = TRUE) + log1p(-ss_pi)
        ss_g <- stats::rbinom(p, 1, 1 / (1 + exp(pmin(l0_ - l1_, 700))))
        ss_pi <- stats::rbeta(1, 1 + sum(ss_g), 1 + p - sum(ss_g))
      }

      ## ---- 3. (lam1, lam2) | beta, tau2 : random-walk MH on the log scale
      if (!fixed_lam && prior == "adenet") {
        lp <- function(l1, l2) {
          a <- wlr * l1 * wj; c2 <- wlr * l2
          sum(log(a) - a^2 * tau2 / 2 - log_enet_norm(a, c2)) -
            c2 * sum(beta^2) / 2 +
            (a_lam - 1) * log(l1) - b_lam * l1 +
            (a_lam - 1) * log(l2) - b_lam * l2 +
            log(l1) + log(l2)                          # Jacobian
        }
        l1p <- lam1 * exp(stats::rnorm(1, 0, sd_h))
        l2p <- lam2 * exp(stats::rnorm(1, 0, sd_h))
        la <- lp(l1p, l2p) - lp(lam1, lam2)
        if (is.finite(la) && log(stats::runif(1)) < la) {
          lam1 <- l1p; lam2 <- l2p; acc_h <- acc_h + 1
        }
        if (it <= burn) sd_h <- sd_h * exp((min(exp(la), 1) - 0.3) / sqrt(it))
      }

      if (it > burn && ((it - burn) %% thin == 0)) {
        k <- k + 1L; B[k, ] <- beta; A[k] <- alpha; L[k, ] <- c(lam1, lam2)
        if (prior == "ssvs") G <- rbind(G, ss_g)
      }
    }
    out[[ch]] <- list(gamma = G, beta = B[seq_len(k), , drop = FALSE], alpha = A[seq_len(k)],
                      lambda = L[seq_len(k), , drop = FALSE],
                      acc_beta = acc / nsim, acc_lambda = acc_h / nsim, step = step)
    if (verbose) message(sprintf("chain %d: acc(beta)=%.2f acc(lambda)=%.2f",
                                 ch, out[[ch]]$acc_beta, out[[ch]]$acc_lambda))
  }

  ## --- back-transform to the original scale and summarise
  Braw <- do.call(rbind, lapply(out, `[[`, "beta"))
  Araw <- unlist(lapply(out, `[[`, "alpha"))
  B <- sweep(Braw, 2, xs, "/")
  A <- Araw - as.vector(Braw %*% (xc / xs))

  res <- list(
    beta = B, alpha = A,
    lambda = do.call(rbind, lapply(out, `[[`, "lambda")),
    chains = out, w = wlr, sigma = sigma, weights = wj, d = d, fam = fam,
    prior = prior,
    ## posterior inclusion probabilities, available only under the SSVS prior
    pip = if (prior == "ssvs")
      colMeans(do.call(rbind, lapply(out, `[[`, "gamma"))) else NULL,
    pilot_engine = pilot$engine,
    map = list(alpha = map$alpha - sum(map$beta * xc / xs),
               beta = map$beta / xs, lam1 = map$lam1, lam2 = map$lam2,
               capped = capped_any),
    rhat = rhat_chains(lapply(out, function(o) sweep(o$beta, 2, xs, "/"))),
    diag = if (diagnostics && chains > 1)
      diag_chains(lapply(out, function(o) sweep(o$beta, 2, xs, "/"))) else NULL,
    alpha_bound = a_half, starts = starts,
    call = match.call())
  class(res) <- "bt_adenet"
  res
}

## split-free Gelman-Rubin across chains (NA with a single chain)
rhat_chains <- function(mats) {
  if (length(mats) < 2) return(rep(NA_real_, ncol(mats[[1]])))
  m <- length(mats); nn <- min(sapply(mats, nrow))
  mats <- lapply(mats, function(z) z[seq_len(nn), , drop = FALSE])
  means <- t(sapply(mats, colMeans))
  vars <- t(sapply(mats, function(z) apply(z, 2, stats::var)))
  W <- colMeans(vars)
  Bn <- apply(means, 2, stats::var)
  vhat <- (nn - 1) / nn * W + Bn
  ifelse(W > 0, sqrt(vhat / W), NA_real_)
}

## Rank-normalized split R-hat with bulk and tail ESS (Vehtari, Gelman, Simpson,
## Carpenter and Burkner 2021).  rhat_chains() above is kept so that earlier
## tables stay reproducible.  mats is a list of draws x parameters matrices.
acov_fft <- function(v) {
  n <- length(v)
  f <- stats::fft(c(v - mean(v), numeric(n)))
  Re(stats::fft(Mod(f)^2, inverse = TRUE))[seq_len(n)] / (2 * n) / n
}
diag_chains <- function(mats) {
  h <- min(sapply(mats, nrow)) %/% 2
  if (h < 4) return(NULL)
  halves <- unlist(lapply(mats, function(z)
    list(z[seq_len(h), , drop = FALSE], z[h + seq_len(h), , drop = FALSE])),
    recursive = FALSE)
  m <- length(halves)
  zrank <- function(a) matrix(stats::qnorm((rank(a) - 3 / 8) / (length(a) + 1 / 4)), h)
  rhat1 <- function(a) {
    W <- mean(apply(a, 2, stats::var))
    if (!(W > 0)) return(NA_real_)
    sqrt(((h - 1) / h * W + stats::var(colMeans(a))) / W)
  }
  ess1 <- function(a) {
    W <- mean(apply(a, 2, stats::var))
    if (!(W > 0)) return(NA_real_)
    vplus <- (h - 1) / h * W + stats::var(colMeans(a))
    rho <- 1 - (W - rowMeans(apply(a, 2, acov_fft))) / vplus
    tau <- -1
    t <- 1
    while (t < h) {
      pr <- rho[t] + rho[t + 1]
      if (pr < 0) break
      tau <- tau + 2 * pr
      t <- t + 2
    }
    m * h / max(tau, 1 / log10(m * h))
  }
  one <- function(v) {
    q <- stats::quantile(v, c(0.05, 0.95))
    c(rhat = max(rhat1(zrank(v)), rhat1(zrank(abs(v - stats::median(v))))),
      ess_bulk = ess1(zrank(v)),
      ess_tail = suppressWarnings(min(ess1(1 * (v <= q[1])), ess1(1 * (v <= q[2])),
                                      na.rm = TRUE)))
  }
  t(sapply(seq_len(ncol(mats[[1]])), function(j) one(sapply(halves, function(z) z[, j]))))
}

## ------------------------------------------------------------------ summary

#' Posterior summary and credible-interval variable selection
summary.bt_adenet <- function(object, level = 0.95, ...) {
  q <- c((1 - level) / 2, 1 - (1 - level) / 2)
  ci <- t(apply(object$beta, 2, stats::quantile, probs = q))
  ## Under SSVS selection is by posterior inclusion probability at the median
  ## model, which is a genuine posterior statement. Under the continuous priors it
  ## has to be a summary rule, since they put no mass on sparse vectors.
  sel <- if (!is.null(object$pip)) object$pip > 0.5 else (ci[, 1] > 0 | ci[, 2] < 0)
  data.frame(
    mean = colMeans(object$beta),
    median = apply(object$beta, 2, stats::median),
    lower = ci[, 1], upper = ci[, 2],
    pip = if (is.null(object$pip)) NA_real_ else object$pip,
    selected = sel,
    rhat = object$rhat)
}

coef.bt_adenet <- function(object, type = c("median", "mean", "map"), ...) {
  type <- match.arg(type)
  switch(type,
         median = apply(object$beta, 2, stats::median),
         mean = colMeans(object$beta),
         map = object$map$beta)
}

#' Sparsified point estimate: posterior median zeroed outside the credible set
selected_coef <- function(object, level = 0.95) {
  s <- summary(object, level)
  ifelse(s$selected, s$median, 0)
}

## Experiment helper: parallel lapply over task indices that logs and drops failed
## tasks.  Every task sets its own seeds, so results do not depend on the worker.
par_lapply <- function(n, f, cores) {
  r <- parallel::mclapply(seq_len(n), function(i) tryCatch(f(i), error = function(e) {
    message("task ", i, " failed: ", conditionMessage(e))
    NULL
  }), mc.cores = cores, mc.preschedule = FALSE)
  Filter(Negate(is.null), r)
}
