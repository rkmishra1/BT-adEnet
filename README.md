# BT-adEnet

A Gibbs Posterior for Robust Sparse Regression under a Bounded Loss: Exact
Poisson Augmentation and a Conjugate Sampler.

The Tukey biweight criterion is bounded, so `exp(-Tukey)` is not integrable and
there is no likelihood. The model is therefore a **generalized (Gibbs)
posterior**: the biweight enters as a loss with a calibrated learning rate `w`,
and the adaptive elastic-net penalty is reinterpreted as a proper prior via a
scale mixture of normals. The T-AdEnet estimator is the exact posterior mode, so
every posterior summary describes uncertainty around the frequentist solution.

## Code

| file | produces | needed for |
|---|---|---|
| `R/bt_adenet.R` | — | Algorithm 1: loss, pilot, T-AdEnet solver, sampler. Every other script `source()`s it |
| `R/validate.R` | `sim-results.rds` | Table 1, `sim*` macros |
| `R/calibration2.R` | `calib2.rds` | Tables 2–3, Figure 4, `cal*` macros |
| `R/experiments.R` | `e2-sensitivity.rds` | Figure 5, `es*` macros (run arm E2; E1/E3/E4 are in the file but superseded) |
| `R/ratio.R` | `ratio.rds` | Table 4, Figure 6 (leakage + weight floor), `ra*` macros |
| `R/diffuse.R` | `diffuse.rds` | diffuse-prior table, `df*` macros |
| `R/comparators.R` | `comparators.rds` | comparator + runtime tables, Figure 7, `cp*` macros |
| `R/frequentist.R` | `frequentist.rds` | the enetLTS column of the comparator table |
| `R/sensitivity.R` | `sensitivity.rds` | sensitivity table, `se*` macros |
| `R/realdata.R` | `realdata.rds`, `realdata-agg.rds` | real-data table, Figure 8, `rd*` macros |
| `R/posterior-fits.R` | `posterior-fits.rds` | Figures 9–10, `pf*` macros |
| `R/tables.R` | the 10 `\input` files | every table |
| `R/figures.R` | the 7 generated PDFs | every figure except the TikZ factor graph, which is inline in the `.tex` |

The experiment scripts write their `.rds` outputs to the repository root;
`R/tables.R` and `R/figures.R` read those files and produce the LaTeX tables
and figure PDFs used in the manuscript.

## Running

All scripts are run from the repository root, e.g.

```bash
Rscript R/validate.R [nrep] [cores] [largest p to run] [refit 1/0]
Rscript R/calibration2.R [nrep] [nrep_freq] [B] [cores]
Rscript R/experiments.R [which, e.g. "E2"] [nrep] [cores]
Rscript R/ratio.R [nrep] [cores]
Rscript R/diffuse.R [nrep] [cores]
Rscript R/comparators.R [nrep] [cores]
Rscript R/frequentist.R [nrep] [cores]
Rscript R/sensitivity.R [nrep] [cores]
Rscript R/realdata.R [nsplit_1 .. nsplit_4] [cores] [merge]
Rscript R/posterior-fits.R
Rscript R/tables.R
Rscript R/figures.R
```

The sampler in `R/bt_adenet.R` is base R only; `robustbase` is used for the MM
pilot when installed and `p < n/2`, otherwise a self-contained
ridge → Huber → Tukey IRLS pilot is used. The comparator arm additionally uses
`enetLTS`.
