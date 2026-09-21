## LaTeX tables for bayes-tadenet.tex from the experiment .rds files.
##   Rscript R/tables.R
## Captions describe what is reported and never what it shows, so a rerun cannot
## leave a caption asserting something the new numbers contradict.
here <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
out <- file.path(here, "..")
tb <- function(name, lines) writeLines(lines, file.path(out, name))
fmt <- function(x, k = 3) ifelse(is.finite(x), formatC(x, format = "f", digits = k), "--")
mcse <- function(v) stats::sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v)))
pm <- function(m, s, k = 2) paste0(fmt(m, k), " (", fmt(s, k), ")")
have <- function(f) file.exists(file.path(out, f))
lab <- c(clean = "clean", vertical = "vertical", both = "leverage")

## ---- E1 and E4 on one design, and the frequentist intervals --------------
if (have("calib2.rds")) {
  cb <- readRDS(file.path(out, "calib2.rds"))
  a <- cb$agg_sweep[order(cb$agg_sweep$w), ]
  m <- cb$agg_meth
  wcal <- m$w[m$method == "BT-AdEnet"]
  tb("e1-table.tex", c(
  "\\begin{table}[!htb]\\centering",
  sprintf(paste0("\\caption{Coverage, width and excess risk against the learning rate, %d ",
    "replications. Coverage and width are means over the active set, excess risk is the Tukey risk ",
    "on a clean test sample minus the risk at $\\beta^{*}$, and Monte Carlo standard errors are in ",
    "parentheses.}"), a$R[1]),
  "\\label{tab:e1}", "\\resizebox{\\textwidth}{!}{%",
  paste0("\\begin{tabular}{l", strrep("r", nrow(a)), "}"), "\\toprule",
  paste0("$w$ & ", paste(fmt(a$w, 2), collapse = " & "), " \\\\"), "\\midrule",
  paste0("coverage & ", paste(pm(a$cov, a$cov_se, 3), collapse = " & "), " \\\\"),
  paste0("width & ", paste(fmt(a$width, 3), collapse = " & "), " \\\\"),
  paste0("excess risk & ", paste(pm(a$excess, a$excess_se, 4), collapse = " & "), " \\\\"),
  "\\bottomrule", "\\end{tabular}}", "\\end{table}"))

  mm <- m[match(c("BT-AdEnet", "bootstrap", "sandwich"), m$method), ]
  tb("methods-table.tex", c(
  "\\begin{table}[!htb]\\centering",
  sprintf(paste0("\\caption{Posterior, bootstrap and sandwich intervals around the same estimate, ",
    "design of Table~\\ref{tab:e1}. Coverage and width are over the active set, TPR is the share of ",
    "active coefficients selected, and seconds, $\\hat R$ and effective sample sizes are medians, the ",
    "last two for the worst coordinate.}")),
  "\\label{tab:methods}",
  "\\begin{tabular}{lrrrrrrr}", "\\toprule",
  "method & coverage & width & TPR & seconds & $\\hat R$ & ESS$_{\\text{bulk}}$ & ESS$_{\\text{tail}}$ \\\\",
  "\\midrule",
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\", mm$method, pm(mm$cov, mm$cov_se, 3),
          fmt(mm$width, 3), fmt(mm$tpr, 2), fmt(mm$secs, 1), fmt(mm$rhat, 3),
          fmt(mm$ess_bulk, 0), fmt(mm$ess_tail, 0)),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"))
  cat("calibration: bias split\n")
  print(cb$agg_z, row.names = FALSE)
}

## ---- E2 bounded posterior sensitivity ------------------------------------
if (have("e2-sensitivity.rds")) {
  e2 <- readRDS(file.path(out, "e2-sensitivity.rds"))
  a <- e2$agg[order(e2$agg$loss, e2$agg$t), ]
  ts <- sort(unique(a$t))
  row <- function(ls, col, k = 3) {
    z <- a[a$loss == ls, ]
    paste(sapply(ts, function(tt) {
      v <- z[[col]][z$t == tt]
      if (!length(v)) "--" else fmt(v, k)
    }), collapse = " & ")
  }
  tb("e2-table.tex", c(
  "\\begin{table}[!htb]\\centering",
  sprintf(paste0("\\caption{Posterior mean shift, scale and learning rate when $y_1$ is replaced ",
    "by $y_1+t$, medians over %d replications. The Gaussian arm is the same sampler with ",
    "$d=\\infty$.}"), length(unique(e2$raw$rep))),
  "\\label{tab:e2}",
  paste0("\\begin{tabular}{ll", strrep("r", length(ts)), "}"), "\\toprule",
  paste0("loss & $t$ & ", paste(formatC(ts, format = "g"), collapse = " & "), " \\\\"),
  "\\midrule",
  paste0("\\multirow{3}{*}{Tukey} & $\\|\\Delta\\hat\\beta\\|$ & ", row("Tukey", "shift"), " \\\\"),
  paste0(" & $\\hat\\sigma$ & ", row("Tukey", "sigma", 2), " \\\\"),
  paste0(" & $w$ & ", row("Tukey", "w", 2), " \\\\"),
  "\\midrule",
  paste0("\\multirow{3}{*}{Gaussian} & $\\|\\Delta\\hat\\beta\\|$ & ", row("Gaussian", "shift"), " \\\\"),
  paste0(" & $\\hat\\sigma$ & ", row("Gaussian", "sigma", 2), " \\\\"),
  paste0(" & $w$ & ", row("Gaussian", "w", 3), " \\\\"),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"))
}

## ---- leakage and the separation ratio (E3 and E3b redone) ---------------
if (have("ratio.rds")) {
  r <- readRDS(file.path(out, "ratio.rds"))
  r$Q <- (r$p - 6) / (r$prior_leak * r$w * r$lam1) / r$wmaxA
  r$cfg <- ifelse(r$floor == "p", "pfloor", paste0("iota", r$iota))
  cfgs <- c("iota0", "iota0.5", "iota1", "iota1.5", "iota2", "pfloor")
  nm <- c(iota0 = "$\\iota=0$ (uniform)", iota0.5 = "$\\iota=0.5$", iota1 = "$\\iota=1$",
          iota1.5 = "$\\iota=1.5$", iota2 = "$\\iota=2$", pfloor = "$\\iota=1$, $p$-scaled floor")
  ps <- sort(unique(r$p))
  slope <- function(g, col) {
    if (stats::sd(log(g[[col]])) == 0) return(c(0, 0))
    m <- stats::lm(log(g[[col]]) ~ log(g$p))
    c(unname(coef(m)[2]), unname(sqrt(diag(stats::vcov(m)))[2]))
  }
  block <- function(col, k) unlist(lapply(cfgs, function(cf) {
    g <- r[r$cfg == cf, ]
    sl <- slope(g, col)
    sprintf("%s & %s & %s \\\\", nm[cf],
      paste(fmt(sapply(ps, function(pp) stats::median(g[[col]][g$p == pp])), k), collapse = " & "),
      pm(sl[1], sl[2], 2))
  }))
  tb("e3-table.tex", c(
  "\\begin{table}[!htb]\\centering",
  sprintf(paste0("\\caption{Posterior leakage and separation ratio $Q$ as $p$ grows, medians over ",
    "%d replications. The last column is the slope of the log quantity on $\\log p$ with its ",
    "standard error.}"),
    length(unique(r$rep))),
  "\\label{tab:e3}", "\\resizebox{\\textwidth}{!}{%",
  paste0("\\begin{tabular}{l", strrep("r", length(ps)), "r}"), "\\toprule",
  paste0("weights & ", paste0("$p=", ps, "$", collapse = " & "), " & slope \\\\"), "\\midrule",
  paste0("\\multicolumn{", length(ps) + 2, "}{l}{\\emph{posterior leakage}} \\\\"),
  block("leak", 2), "\\midrule",
  paste0("\\multicolumn{", length(ps) + 2, "}{l}{\\emph{separation ratio $Q$}} \\\\"),
  block("Q", 1),
  "\\bottomrule", "\\end{tabular}}", "\\end{table}"))
  m <- stats::lm(log(leak) ~ log(Q) + factor(p), data = r)
  cat(sprintf("ratio: log leak on log Q within p, slope %.3f (se %.3f)\n",
              coef(m)[2], sqrt(diag(stats::vcov(m)))[2]))
}

## ---- sparse against diffuse signal ---------------------------------------
if (have("diffuse.rds")) {
  r <- readRDS(file.path(out, "diffuse.rds"))
  regs <- c("sparse", "moderate", "diffuse")
  pri <- c(adenet = "AdEnet", horseshoe = "horseshoe", ssvs = "SSVS")
  beats <- sapply(regs, function(g) {
    w <- merge(r[r$regime == g & r$prior == "adenet", c("rep", "rmspe")],
               r[r$regime == g & r$prior == "ssvs", c("rep", "rmspe")], by = "rep")
    100 * mean(w$rmspe.x < w$rmspe.y)
  })
  lines <- c("\\begin{table}[!htb]\\centering",
  sprintf(paste0("\\caption{Priors under sparse, moderate and diffuse signals, %d replications. ",
    "RMSPE is on a clean test set, ME is the median model error, cov the mean coverage on the ",
    "active set, and stability the Jaccard index between selections on two subsamples, with the ",
    "number of replications where both were empty in parentheses.}"),
    length(unique(r$rep))),
  "\\label{tab:diffuse}",
  "\\begin{tabular}{llrrrrrr}", "\\toprule",
  "signal & prior & RMSPE & ME & TPR & FPR & cov & stability \\\\", "\\midrule")
  for (g in regs) {
    if (g != "sparse") lines <- c(lines, "\\midrule")
    for (pr in names(pri)) {
      z <- r[r$regime == g & r$prior == pr, ]
      lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s & %s (%d) \\\\",
        if (pr == "adenet") g else "", pri[pr], pm(mean(z$rmspe), mcse(z$rmspe), 3),
        fmt(stats::median(z$me), 2), fmt(mean(z$tpr), 2), fmt(mean(z$fpr), 3),
        fmt(mean(z$cov), 2), fmt(mean(z$stab, na.rm = TRUE), 2), sum(z$null_pair)))
    }
  }
  tb("diffuse-table.tex", c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}"))
}

## ---- comparators, one table per p, with the enetLTS row ------------------
if (have("comparators.rds")) {
  r <- readRDS(file.path(out, "comparators.rds"))
  fq <- if (have("frequentist.rds")) readRDS(file.path(out, "frequentist.rds")) else NULL
  R <- length(unique(r$rep))
  one_table <- function(pp, cap, label) {
    lines <- c("\\begin{table}[!htb]\\centering", cap, label,
      "\\begin{tabular}{llrrrrrr}", "\\toprule",
      "setting & arm & ME & TPR & FPR & cov & leak & ESS$_{\\text{bulk}}$ \\\\", "\\midrule")
    for (sc in c("clean", "vertical", "both")) {
      z <- r[r$p == pp & r$scenario == sc, ]
      s <- do.call(rbind, lapply(unique(z$arm), function(am) {
        g <- z[z$arm == am, ]
        data.frame(arm = am, me = stats::median(g$me), tpr = mean(g$tpr), fpr = mean(g$fpr),
                   cov = mean(g$cov), cse = mcse(g$cov), leak = stats::median(g$leak),
                   ess = stats::median(g$ess_bulk))
      }))
      s <- s[order(s$me), ]
      if (sc != "clean") lines <- c(lines, "\\midrule")
      lines <- c(lines, sprintf("%s & \\texttt{%s} & %s & %s & %s & %s & %s & %s \\\\",
        ifelse(seq_len(nrow(s)) == 1, lab[sc], ""), s$arm, fmt(s$me, 3), fmt(s$tpr, 2),
        fmt(s$fpr, 3), pm(s$cov, s$cse, 2), fmt(s$leak, 2), fmt(s$ess, 0)))
      g <- if (is.null(fq)) NULL else fq[fq$p == pp & fq$scenario == sc, ]
      if (length(g) && nrow(g))
        lines <- c(lines, sprintf(" & \\texttt{enetLTS} & %s & %s & %s & -- & -- & -- \\\\",
          fmt(stats::median(g$me), 3), fmt(mean(g$tpr), 2), fmt(mean(g$fpr), 3)))
    }
    c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  }
  cap35 <- sprintf(paste0("\\caption{Losses and priors at $p=35$ ($s=9$), %d replications, ordered ",
    "by model error within each block. ME is the median model error, TPR and FPR are means, cov is ",
    "the mean coverage of $95\\%%$ intervals on the active set, leak the median prior leakage, and ",
    "ESS$_{\\text{bulk}}$ the median smallest bulk effective sample size. The last row of each block ",
    "is enetLTS \\parencite{Kurnaz2018}.}"), R)
  cap100 <- paste0("\\caption{Losses and priors at $p=100$ ($s=12$), otherwise as ",
    "Table~\\ref{tab:comparators}.}")
  tb("comparators-table.tex", c(one_table(35, cap35, "\\label{tab:comparators}"), "",
                                one_table(100, cap100, "\\label{tab:comparators100}")))

  a <- aggregate(secs ~ arm + p, r, stats::median)
  ps <- sort(unique(a$p))
  lines <- c("\\begin{table}[!htb]\\centering",
  paste0("\\caption{Median seconds per fit in the runs of Tables~\\ref{tab:comparators} ",
    "and~\\ref{tab:comparators100}. Fits shared a machine, so compare rows rather than absolute ",
    "times.}"),
  "\\label{tab:runtime}",
  paste0("\\begin{tabular}{l", strrep("r", length(ps)), "}"), "\\toprule",
  paste0("arm & ", paste0("$p=", ps, "$", collapse = " & "), " \\\\"), "\\midrule")
  for (am in sort(unique(a$arm)))
    lines <- c(lines, sprintf("\\texttt{%s} & %s \\\\", am,
      paste(sapply(ps, function(pp) fmt(a$secs[a$arm == am & a$p == pp], 1)), collapse = " & ")))
  if (!is.null(fq)) {
    b <- aggregate(secs ~ p, fq, stats::median)
    lines <- c(lines, sprintf("\\texttt{enetLTS} & %s \\\\",
      paste(sapply(ps, function(pp) fmt(b$secs[b$p == pp], 1)), collapse = " & ")))
  }
  tb("runtime-table.tex", c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}"))
}

## ---- sensitivity -----------------------------------------------------------
if (have("sensitivity.rds")) {
  r <- readRDS(file.path(out, "sensitivity.rds"))
  nm <- c(d = "$d$", iota = "$\\iota$", b_lam = "hyperprior rate", floor = "weight floor")
  fl <- c("0" = "default", "1" = "$p$-scaled", "2" = "$10^{-3}$")
  lines <- c("\\begin{table}[!htb]\\centering",
  sprintf(paste0("\\caption{Sensitivity to the fixed constants, %d replications at $p=50$. One ",
    "constant varies at a time, with defaults in bold, and the columns are as in ",
    "Table~\\ref{tab:comparators}.}"),
    length(unique(r$rep))),
  "\\label{tab:sensitivity}",
  "\\begin{tabular}{llrrrr}", "\\toprule",
  "constant & value & ME & TPR & FPR & cov \\\\", "\\midrule")
  for (k in c("d", "iota", "b_lam", "floor")) {
    if (k != "d") lines <- c(lines, "\\midrule")
    vals <- sort(unique(r$value[r$knob == k]))
    dflt <- c(d = 4.685, iota = 1, b_lam = 1e-2, floor = 0)[k]
    for (i in seq_along(vals)) {
      g <- r[r$knob == k & r$value == vals[i], ]
      v <- if (k == "floor") fl[as.character(vals[i])] else formatC(vals[i], format = "g")
      if (isTRUE(all.equal(vals[i], unname(dflt)))) v <- paste0("\\textbf{", v, "}")
      lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s \\\\",
        if (i == 1) nm[k] else "", v, fmt(stats::median(g$me), 3), fmt(mean(g$tpr), 2),
        fmt(mean(g$fpr), 4), pm(mean(g$cov), mcse(g$cov), 2)))
    }
  }
  tb("sensitivity-table.tex", c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}"))
}

## ---- real data -------------------------------------------------------------
if (have("realdata-agg.rds")) {
  a <- readRDS(file.path(out, "realdata-agg.rds"))
  ord <- c("BT-AdEnet" = 1, "Tukey+SSVS" = 2, "Gaussian" = 3)
  a <- a[order(match(a$dataset, c("Boston housing", "Ames housing", "Glass", "Riboflavin")),
               ord[a$arm]), ]
  lines <- c("\\begin{table}[!htb]\\centering",
  paste0("\\caption{Real data, medians over $2{:}1$ training and test splits. RTMSPE is the ",
    "prediction error trimmed at $90\\%$ and RMSPE the untrimmed one. Cover is the mean coverage of ",
    "$90\\%$ plug-in predictive intervals, null the number of splits selecting nothing, and the count ",
    "after $w$ the splits where $w=1$ was used. $\\hat R$ and ESS$_{\\text{tail}}$ are for the worst ",
    "coordinate.}"),
  "\\label{tab:realdata}", "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{llrrrrrrrrrrr}", "\\toprule",
  paste("dataset & arm & $R$ & RTMSPE & RMSPE & cover & width & $|\\hat{\\mathcal{A}}|$ & null",
        "& stability & $w$ & $\\hat R$ & ESS$_{\\text{tail}}$ \\\\"),
  "\\midrule")
  prev <- ""
  for (i in seq_len(nrow(a))) {
    z <- a[i, ]
    dcol <- if (identical(as.character(z$dataset), prev)) "" else as.character(z$dataset)
    if (dcol != "" && i > 1) lines <- c(lines, "\\midrule")
    prev <- as.character(z$dataset)
    lines <- c(lines, sprintf(
      "%s & \\texttt{%s} & %d & %s & %s & %s & %s & %.0f & %d & %s & %s (%d) & %s & %s \\\\",
      dcol, z$arm, z$R, fmt(z$rtmspe, 4), fmt(z$rmspe, 4), fmt(z$cover_all, 3), fmt(z$width, 3),
      z$size, z$null, fmt(z$stability, 3), fmt(z$w, 2), z$w_fallback, fmt(z$rhat, 2),
      fmt(z$ess_tail, 0)))
  }
  tb("realdata-table.tex", c(lines, "\\bottomrule", "\\end{tabular}}", "\\end{table}"))
}
cat("tables written\n")
