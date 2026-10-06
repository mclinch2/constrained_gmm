# Gaussian Finite-Mixture Simulation Code

This folder contains the code for two simulation studies accompanying the manuscript **“Bayesian Analysis Using a Constrained Mixture of Normal-Inverse-Gamma Models.”** One study uses an intercept-only Gaussian mixture (`y`); the other uses a Gaussian mixture of regressions (`xy`). The scripts compare the Direct Sampler (DS), constrained direct samplers, DS-ML, MC-MCMC, and a conventional MCMC baseline.

## Files

| File | Purpose |
| --- | --- |
| `sim_study_y.R` | Main driver for the **intercept-only** simulation. It generates replicated Gaussian-mixture data, runs all competing samplers, calculates clustering, mixing, density-accuracy, and timing metrics, and writes time-stamped tab-delimited results plus an `.RData` workspace. |
| `sim_study_y_functions.R` | Functions used by the intercept-only driver. It implements the full MCMC and direct/constrained samplers, builds candidate partitions from subsets or machine-learning clustering methods, draws model parameters conditional on allocations, and provides candidate diagnostics. |
| `sim_study_xy.R` | Main driver for the **mixture-of-regressions** simulation. It generates data from mixtures of simple linear regressions, runs the competing samplers, evaluates ARI, ESS, conditional KS/L2 accuracy, and run time, and writes time-stamped results plus an `.RData` workspace. |
| `sim_xy_functions.R` | Functions used by the regression driver. It implements collapsed posterior calculations, candidate-set construction, DS/DS-Const/DS-ML, full MCMC, and MC-MCMC for the regression model. |
| `check_results.R` | Validates completed result files before analysis. It checks for duplicate cell/replicate rows, incorrect replicate counts, incomplete method output, and invalid or out-of-range metrics. It exits with an error status if a check fails. |
| `summarize_results.R` | Reads either study’s result files and prints manuscript-style summaries by design cell and method, including replicate count, ARI, two ESS measures, KS, L2, and CPU time. It summarizes results but does not perform the integrity checks above. |
| `plot_sim_y_paper.R` | Combines intercept-only result files, averages metrics over replicates, converts internal method names to manuscript labels, and saves one faceted PNG figure per simulation scenario. |
| `plot_sim_xy_paper.R` | Performs the corresponding aggregation and plotting for the regression study, using a layout and color scheme comparable to the intercept-only figures. |

## Typical workflow

1. Run the appropriate simulation driver. Each driver sources its matching functions file.
2. Run `check_results.R` on the resulting text files.
3. Use `summarize_results.R` to review manuscript-style numerical summaries.
4. Run the matching plotting script to create the figures.


## Inputs and outputs

Both simulation drivers accept the following optional positional arguments:

```text
N  k  clust_sep  sig_const  prior_M  reps  cores
```

Comma-separated values can be supplied for the design parameters. `clust_sep` codes the degree of component overlap, and `prior_M = 2` uses the manuscript setting of 25 available mixture components. Exact DS enumeration is run only when `N <= 13`.

The drivers use these environment variables:

- `FMM_SRC`: path to the matching functions file.
- `FMM_OUT`: directory in which simulation output is saved.
- `FMM_TAG`: optional suffix that prevents filename collisions across concurrent jobs.
- `FMM_SUB_NITER`, `FMM_SUB_NBURN`, and `FMM_SUB_NTHIN`: optional controls for the regression study’s subset fit.

The plotting scripts read all matching result files from `RESULTS_DIR`. The intercept-only figures are saved under `plots_metric_byN_M25_int_only/`; regression figures are saved under `plots_metric_byN_xy_sims_with10000/`.

## R packages

The simulation scripts use `salso`, `coda`, `parallel`, `miscPack`, `mclust`, `kernlab`, `cluster`, `dbscan`, and `e1071`; the regression study additionally uses `mvtnorm` and `matrixStats`. Plotting uses packages from the tidyverse ecosystem, principally `dplyr`, `tidyr`, `purrr`, `stringr`, and `ggplot2` (the intercept-only plotting script also loads `xtable`). The checking and summary scripts use base R.

