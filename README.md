# Bayesian Reconciliation of Hierarchical Forecasts

Companion code for **["Probabilistic Reconciliation of Hierarchical Forecast via Bayes' Rule"](https://link.springer.com/chapter/10.1007/978-3-030-67664-3_13)**, Corani, Azzimonti, Augusto & Zaffalon, ECML-PKDD 2020 (published in LNAI vol. 12459, Springer 2021).

This repository implements **pMinT** and **LG**, two Bayesian methods that reconcile probabilistic forecasts produced independently across the levels of a hierarchical time series. Both methods recover the well-known MinT reconciliation as their point-forecast limit, and additionally return a coherent joint posterior — not just the means — for use in downstream decision-making.

If you just want to reproduce the paper's experiments, jump to [Reproducing the paper](#reproducing-the-paper).

---

## What problem does this solve?

When you forecast a quantity at multiple aggregation levels — country → region → city, or company → division → product — the level-by-level forecasts almost never sum up consistently. The country forecast disagrees with the sum of its regions; the regions disagree with the sum of their cities. *Reconciliation* is the process of adjusting these incoherent forecasts into a single coherent set.

The standard frequentist solution is **MinT** (Wickramasuriya et al., 2019), which finds the minimum-variance unbiased reconciliation by solving a generalized least-squares problem. MinT gives you reconciled point forecasts and a reconciled covariance, but it is fundamentally a regression-style adjustment.

This repository takes a **Bayesian view of the same problem**. The bottom-level forecasts are treated as the prior; the upper-level forecasts are treated as noisy observations that the bottom levels must remain consistent with through the summation matrix. Applying Bayes' rule yields a posterior over the bottom-level series, from which a coherent joint sample for the entire hierarchy can be drawn. Two variants are implemented:

- **LG** (Linear Gaussian) — assumes upper and bottom forecast errors are independent. Recovers MinT exactly when forecast errors are jointly Gaussian.
- **pMinT** (probabilistic MinT) — drops the independence assumption and models the full cross-covariance between upper and bottom errors. Generalizes LG and is the more accurate method when residuals are correlated across levels (which is the empirically common case).

Both methods produce, in one pass, a coherent **mean**, **median**, **full joint covariance**, and a **probabilistic sample** of arbitrary size — usable for energy-score evaluation, quantile reporting, or any downstream decision under uncertainty.

## Methods at a glance

| Method | Upper/bottom correlation | Output | Recovers MinT? |
|---|---|---|---|
| BottomUp | n/a | Point + cov | No (different estimator) |
| MinT (frequentist baseline) | implicit | Point + cov | — |
| **LG** (this repo) | independent | Point + cov + joint sample | Yes, in the Gaussian limit |
| **pMinT** (this repo) | full cross-covariance | Point + cov + joint sample | Yes, in the Gaussian limit |

The core update is implemented in `bayesRecon` in [`reconciliationMethods.R`](reconciliationMethods.R). Matrix-heavy operations (`cppMatMult`, `cppInverse`) are dispatched to C++ via Rcpp in [`cppAlgebra.cpp`](cppAlgebra.cpp) to keep large-hierarchy experiments tractable.

## Repository layout

```
.
├── reconciliationMethods.R      # bayesRecon — the core update; covariance estimators; energy score
├── hierRecBayesianExperiment.R  # one experiment: fit base, reconcile, score
├── batchHier.R                  # runs the full experiment grid (parallelizable)
├── loadTourism.R                # tourism dataset loader → hts object
├── drawLargeHierarchy.R         # synthetic 4-2-1 hierarchy generator
├── cppAlgebra.cpp               # Rcpp matrix multiplication / inversion
├── TourismData_v3.csv           # Australian tourism data (Hyndman et al.)
└── result_analysis.ipynb        # post-hoc analysis of result CSVs
```

## Requirements

R ≥ 3.6 with the following packages (auto-installed on first run by `reconciliationMethods.R`):

- `hts` — hierarchical time series, MinT baselines
- `forecast` — base forecasts (`auto.arima`, `ets`)
- `huge`, `SHIP`, `corpcor` — shrinkage and graphical-lasso covariance estimation
- `MASS`, `tmvtnorm`, `truncnorm` — (truncated) multivariate normal sampling
- `matlib`, `matrixStats`, `wordspace` — matrix utilities
- `Rcpp` — C++ glue for the inner loops
- `optparse` — CLI parsing for `batchHier.R`

The Python analysis notebook (`result_analysis.ipynb`) needs `pandas`, `numpy`, and `matplotlib`.

## Quick start

Run a single experiment from R:

```r
source("hierRecBayesianExperiment.R")

# Reconcile 1-step-ahead ETS forecasts on the infantgts hierarchy
hierRecBayesianExperiment(dset = "infantgts", fmethod = "ets",  h = 1, iTest = 1)

# Same on tourism with auto.arima
hierRecBayesianExperiment(dset = "tourism",   fmethod = "arima", h = 1, iTest = 2)

# Synthetic 4-bottom / 2-mid / 1-top hierarchy, 3-step ahead
hierRecBayesianExperiment(dset = "syntheticLarge", h = 3, synth_n = 300)
```

Arguments:

| Arg | Meaning | Default |
|---|---|---|
| `dset` | `infantgts`, `tourism`, or `syntheticLarge` | required |
| `fmethod` | base forecasting method: `ets` or `arima` | `ets` |
| `h` | forecast horizon (1–4 in the paper) | `1` |
| `iTest` | rolling train/test split index (1–45); enables parallelization | `1` |
| `seed` | RNG seed for the synthetic case | `0` |
| `synth_n` | length of the synthetic series | `100` |
| `sampleSize` | posterior sample size for energy-score / quantiles | `100000` |
| `enforceKhOne` | `TRUE` = kh=1, `FALSE` = kh=h (see paper § 3) | `FALSE` |
| `savePredictions` | dump posterior samples to `prediction/` | `TRUE` |

Each run writes a CSV under `results/` named `{h}_{dset}_{fmethod}_{from_split}to{to_split}.csv`, containing **RMSE, MAE, and energy score** for every (group × reconciliation scheme) combination — `All`, `Upper`, `Bottom`, and each hierarchy depth, scored against `BottomUp`, `pMinT`, `LG`, and the raw base forecasts.

## Using `bayesRecon` directly

If you already have base forecasts and a residual covariance estimate from your own pipeline, you can call the reconciliation step in isolation:

```r
source("reconciliationMethods.R")

# preds:        named vector of base forecasts, length = n_total nodes
# mSumMatrix:   hierarchy summation matrix S (from hts::smatrix)
# mCovar:       residual covariance estimate (see estimateCovariance below)
# reconType:    "pmint" or "lg"
# sampleSize:   posterior sample size

out <- bayesRecon(preds       = preds,
                  mSumMatrix  = S,
                  mCovar      = Sigma_hat,
                  reconType   = "pmint",
                  sampleSize  = 100000,
                  kh          = 1)

# Coherent point forecast for the full hierarchy
out$coherentPreds

# Bottom-level posterior (mean, covariance)
out$posteriorMean
out$posteriorVariance

# Full coherent joint sample (sampleSize × n_total) — use this for any
# downstream probabilistic metric (energy score, quantiles, intervals).
out$sample
```

### Covariance estimators

`estimateCovariance(residuals, method, ...)` exposes the covariance estimators used throughout the paper:

- `"diagonal"` — independent errors
- `"sam"` — sample covariance
- `"sammint"` — uncentered sample covariance (MinT-compatible)
- `"glasso"` — sparse precision via graphical lasso (`huge`)
- `"shr"` — Schäfer-Strimmer shrinkage to a diagonal target
- `"shrmint"` — the shrinkage estimator used by MinT

Mixing different covariance estimators with `bayesRecon` is how the paper's ablation table is generated.

## Reproducing the paper

The full experimental grid is in [`batchHier.R`](batchHier.R). Example:

```bash
Rscript batchHier.R -d infantgts -m arima -p "results/" -k true
Rscript batchHier.R --help    # for all options
```

The grid covers two real datasets (`infantgts`, `tourism`) plus the synthetic `syntheticLarge` hierarchy, both base methods (`ets`, `arima`), horizons 1–4, and all reconciliation schemes. Splits are independent and trivially parallelizable across `iTest`. Aggregated tables and figures are produced by `result_analysis.ipynb`.

## Implementation notes

A few details worth highlighting:

- **Vectorized energy score.** `energyScore()` implements the sample energy score from Gneiting & Raftery using a single random permutation of the sample matrix instead of splitting it in half — produces results numerically identical to `scoringRules::es_sample` at roughly half the cost.
- **Truncated sampling.** `sampleMVN(..., positivity = TRUE)` supports non-negativity constraints via `tmvtnorm::rtmvnorm` (Gibbs) for the joint case, or fast independent `truncnorm::rtruncnorm` draws when the marginals are sufficient.
- **Numerical safety.** The posterior covariance is forced symmetric (`makeSymmetric`) and, if needed, nudged to the nearest positive-definite matrix at decreasing tolerance levels (`make.positive.definite`) before any sampling step. This makes the pipeline robust to the near-singular covariances that show up with shrinkage estimators on deep hierarchies.

## Citation

If you use this code, please cite:

> Corani, G., Azzimonti, D., Augusto, J. P. S. C., & Zaffalon, M. (2021). **Probabilistic Reconciliation of Hierarchical Forecast via Bayes' Rule.** In Hutter, F., Kersting, K., Lijffijt, J., Valera, I. (eds), *Machine Learning and Knowledge Discovery in Databases — ECML PKDD 2020*. Lecture Notes in Computer Science (LNAI), vol. 12459, pp. 211–226. Springer, Cham. [https://doi.org/10.1007/978-3-030-67664-3_13](https://doi.org/10.1007/978-3-030-67664-3_13)

```bibtex
@inproceedings{corani2021probabilistic,
  title     = {Probabilistic Reconciliation of Hierarchical Forecast via Bayes' Rule},
  author    = {Corani, Giorgio and Azzimonti, Dario and Augusto, Jo{\~a}o P. S. C. and Zaffalon, Marco},
  booktitle = {Machine Learning and Knowledge Discovery in Databases (ECML PKDD 2020)},
  editor    = {Hutter, Frank and Kersting, Kristian and Lijffijt, Jefrey and Valera, Isabel},
  series    = {Lecture Notes in Computer Science},
  volume    = {12459},
  pages     = {211--226},
  publisher = {Springer, Cham},
  year      = {2021},
  doi       = {10.1007/978-3-030-67664-3_13}
}
```

## License

Released for academic reproducibility. If you need a different license for downstream use, open an issue.
