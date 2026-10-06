# Bayesian Analysis Using a Constrained Mixture of Normal-Inverse-Gamma Models

This repository contains the R code used for two Gaussian finite-mixture simulation studies:

1. an intercept-only (Y-only) Gaussian mixture study; and
2. a Gaussian mixture-of-regressions (XY) study.

The code compares standard MCMC, an exhaustive direct sampler for very small samples, constrained direct samplers, a machine-learning candidate-set direct sampler, and a collapsed-label MCMC/composition method.

## Repository structure

The scripts are expected to be organized as follows:

```text
constrained_gmm/
├── code/
│   ├── sim_study_y.R
│   ├── sim_study_y_functions.R
│   ├── sim_study_xy.R
│   ├── sim_xy_functions.R
│   ├── check_results.R
│   ├── summarize_results.R
│   ├── plot_sim_y_paper.R
│   ├── plot_sim_xy_paper.R
│   ├── run_sim_y.sh
│   ├── run_sim_xy.sh
│   ├── run_analysis_y.sh
│   ├── run_analysis_xy.sh
│   └── run_plots.sh
├── logs/
├── test_results_y/
└── test_results_xy/
```

The result and plot directories are created as needed. The `logs/` directory must exist **before** submitting a SLURM job because SLURM opens its output and error files before executing the shell script.

```bash
mkdir -p logs
```

## Software requirements

The simulations were developed and tested with R 4.4.0 on a Linux SLURM cluster.

The main R package dependencies are:

```r
cran_packages <- c(
  "salso",
  "coda",
  "mvtnorm",
  "mclust",
  "kernlab",
  "cluster",
  "dbscan",
  "e1071",
  "matrixStats",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "stringr"
)

missing_packages <- cran_packages[
  !vapply(cran_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  install.packages(missing_packages)
}
```

### Installing `miscPack` from GitHub

The simulations require Garritt Page's [`miscPack`](https://github.com/gpage2990/miscPack) package. This is an R package named `miscPack`; there is no function named `misc_pack()`. After installation, load it with `library(miscPack)` or call its mixture function explicitly as `miscPack::gaussian_mixture()`.

`miscPack` is not currently installed by the generic CRAN block above. Its package metadata lists `MASS`, `spam`, and `scaleGMRF` as dependencies. Because `scaleGMRF` is also installed from GitHub and imports `INLA`, use the following installation block for a new R library:

```r
options(
  repos = c(
    CRAN = "https://cloud.r-project.org",
    INLA = "https://inla.r-inla-download.org/R/stable"
  )
)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

if (!requireNamespace("scaleGMRF", quietly = TRUE)) {
  remotes::install_github(
    "LFerrariIt/scaleGMRF",
    dependencies = NA,
    upgrade = "never"
  )
}

remotes::install_github(
  "gpage2990/miscPack",
  dependencies = NA,
  upgrade = "never"
)
```

On a cluster, first load the intended version of R and point `R_LIBS_USER` to the library in which the packages should be installed. For example:

```bash
module purge
module load R/4.4.0

export R_LIBS_USER="/path/to/your/R/library"
unset R_LIBS
mkdir -p "${R_LIBS_USER}"

Rscript --vanilla
```

Then paste the preceding R installation block at the R prompt. If the cluster does not permit internet access from compute nodes, perform the installation on an internet-enabled login or development node.

Verify both the installed location and the required function with:

```r
if (!requireNamespace("miscPack", quietly = TRUE)) {
  stop("miscPack is not installed")
}

cat("miscPack version:", as.character(packageVersion("miscPack")), "\n")
cat("miscPack location:", find.package("miscPack"), "\n")

if (!exists(
  "gaussian_mixture",
  where = asNamespace("miscPack"),
  inherits = FALSE
)) {
  stop("The installed miscPack does not contain gaussian_mixture()")
}

cat("miscPack::gaussian_mixture() is available.\n")
```

For a strictly reproducible release, pin the GitHub installation to the exact tested commit by using `"gpage2990/miscPack@<commit-sha>"` and record that commit in this README.

The current Y-only plotting script may also load `xtable`. It is not used by the plotting code and can be removed from `plot_sim_y_paper.R`; otherwise install it with:

```r
install.packages("xtable")
```

## Clone the repository

```bash
git clone https://github.com/mclinch2/constrained_gmm.git
cd constrained_gmm
mkdir -p logs
```

All SLURM commands below should be submitted from the repository root.

## Simulation design

Both studies use the following default design:

- sample sizes: `N = 10, 500, 1000, 10000`;
- true number of components: `k = 3`;
- cluster-separation settings: `1, 2, 3`;
- component standard deviation: `0.25`;
- default fitted-model option: `PRIOR_M = 2`;
- five replicates per cell for the initial test;
- 50 replicates per cell for the full study.

For `k = 3`, the separation codes are:

| Code | Description | Component means/intercepts |
|---:|---|---|
| 1 | Moderate overlap | `(-1, 0, 1)` |
| 2 | Large overlap | `(-2/3, 0, 2/3)` |
| 3 | No overlap | `(-2, 0, 2)` |

`sig_const = 0.25` is a standard deviation, not a variance.

### Meaning of `PRIOR_M`

`PRIOR_M` is an option code rather than the literal number of fitted components:

| `PRIOR_M` | Fitted model |
|---:|---|
| 1 | `M = k + 2` |
| 2 | `M = 25` |
| Any other value | `M = k` |

Thus, the default `PRIOR_M=2` means that the fitted model contains $M=25$ components.

## Running a small test with SLURM

The recommended first step is the five-replicate test. Each simulation shell script submits an array of 12 jobs, one for every combination of four sample sizes and three separation settings.

### Y-only study

```bash
sbatch code/run_sim_y.sh
```

### XY regression study

```bash
sbatch code/run_sim_xy.sh
```

Each study runs:

```text
12 design cells × 5 replicates = 60 simulated datasets
```

The default test results are written to:

```text
test_results_y/
test_results_xy/
```

If a user-specific R library is needed, provide it without editing the script:

```bash
sbatch \
  --export=ALL,FMM_R_LIB=/path/to/R/library \
  code/run_sim_y.sh
```

Use the analogous command for `run_sim_xy.sh`.

## Running the full 50-replicate study with SLURM

Use separate output directories for the final study so the five-replicate test files are not mixed with the final results.

### Y-only study

```bash
sbatch \
  --export=ALL,FMM_REPS=50,FMM_OUT="${PWD}/sim_results_y" \
  code/run_sim_y.sh
```

### XY regression study

```bash
sbatch \
  --export=ALL,FMM_REPS=50,FMM_OUT="${PWD}/sim_results_xy" \
  code/run_sim_xy.sh
```

To use a specific R library as well:

```bash
sbatch \
  --export=ALL,FMM_REPS=50,FMM_OUT="${PWD}/sim_results_y",FMM_R_LIB=/path/to/R/library \
  code/run_sim_y.sh
```

Each complete study runs:

```text
12 design cells × 50 replicates = 600 simulated datasets
```

These are computationally intensive jobs. Do not run the full simulations on a login node.

## Running the R simulation drivers directly

The R drivers accept seven positional arguments:

```text
1. N
2. k
3. clust_sep
4. sig_const
5. prior_M
6. replicates
7. cores
```

Comma-separated values can be supplied for design parameters.

### One five-replicate Y-only cell

```bash
FMM_OUT="${PWD}/test_results_y" \
Rscript --vanilla code/sim_study_y.R \
  10 3 1 0.25 2 5 1
```

### One five-replicate XY cell

```bash
FMM_OUT="${PWD}/test_results_xy" \
Rscript --vanilla code/sim_study_xy.R \
  10 3 1 0.25 2 5 1
```

### All Y-only design cells

```bash
FMM_OUT="${PWD}/sim_results_y" \
Rscript --vanilla code/sim_study_y.R \
  10,500,1000,10000 3 1,2,3 0.25 2 50 1
```

### All XY design cells

```bash
FMM_OUT="${PWD}/sim_results_xy" \
Rscript --vanilla code/sim_study_xy.R \
  10,500,1000,10000 3 1,2,3 0.25 2 50 1
```

Direct execution is useful for a small local test. SLURM is recommended for the full study.

## Simulation output

The Y-only driver writes files matching:

```text
study1_saveK_results*.txt
MixSimSaveK_y_*.RData
```

The XY driver writes files matching:

```text
study_xy_saveK_results*.txt
MixSimSaveK_xy_*.RData
```

The tab-delimited text files contain one row per design-cell replicate. Columns contain the design identifiers and method-specific metrics.

The principal reported metrics are:

- mean and standard deviation of adjusted Rand index (ARI);
- effective sample size for the first mixture moment;
- effective sample size for the second mixture moment;
- Kolmogorov-Smirnov discrepancy;
- L2 density discrepancy; and
- elapsed time in seconds.

The `.RData` files retain additional R objects, including posterior draws of the number of occupied components.

## Checking and analyzing the five-replicate results

The analysis scripts perform three operations:

1. check the result files for duplicate or incomplete design cells;
2. print a table-style numerical summary; and
3. create the study figures.

### Y-only analysis

```bash
sbatch code/run_analysis_y.sh
```

### XY analysis

```bash
sbatch code/run_analysis_xy.sh
```

The table-style summary is printed to the corresponding SLURM log file under `logs/`.

## Analyzing the final 50-replicate results

### Y-only analysis

```bash
sbatch \
  --export=ALL,FMM_REPS=50,FMM_OUT="${PWD}/sim_results_y" \
  code/run_analysis_y.sh
```

### XY analysis

```bash
sbatch \
  --export=ALL,FMM_REPS=50,FMM_OUT="${PWD}/sim_results_xy" \
  code/run_analysis_xy.sh
```

`FMM_REPS` must match the number of replicates requested when the simulations were run.

## Running the result check directly

For Y-only results:

```bash
Rscript --vanilla code/check_results.R \
  5 test_results_y/study1_saveK_results*.txt
```

For XY results:

```bash
Rscript --vanilla code/check_results.R \
  5 test_results_xy/study_xy_saveK_results*.txt
```

The revised check exits with a nonzero status if it finds duplicate rows, incorrect replicate counts, incomplete method results, or invalid metric values. The analysis shell scripts use `set -e`, so a failed check prevents summaries and figures from being created.

## Printing summaries directly

For Y-only test results:

```bash
Rscript --vanilla code/summarize_results.R \
  test_results_y --reps 5
```

For XY test results:

```bash
Rscript --vanilla code/summarize_results.R \
  test_results_xy --reps 5
```

For a full study, replace the directory and use `--reps 50`.

## Creating figures

The combined plotting shell script checks the results and creates figures for both studies:

```bash
sbatch code/run_plots.sh
```

Plot only one study with:

```bash
sbatch --export=ALL,FMM_STUDY=y  code/run_plots.sh
sbatch --export=ALL,FMM_STUDY=xy code/run_plots.sh
```

For the final results:

```bash
sbatch \
  --export=ALL,FMM_REPS=50,FMM_OUT_Y="${PWD}/sim_results_y",FMM_OUT_XY="${PWD}/sim_results_xy" \
  code/run_plots.sh
```

Figures are written to:

```text
plots_metric_byN_M25_int_only/
plots_metric_byN_xy_sims_with10000/
```

Under the default design, each study creates three PNG files:

```text
k3_sep1_sig0.25.png
k3_sep2_sig0.25.png
k3_sep3_sig0.25.png
```

Each figure contains:

- one column for each sample size;
- one row for each reported metric;
- one bar for each applicable method; and
- the mean across simulation replicates printed above each bar.

The exhaustive direct sampler is only run for the small-sample setting, so it is absent from the larger-(N) panels.

## Methods represented in the output

The principal method labels are:

| Output label | Description |
|---|---|
| `MCMC` | Standard data-augmented Gibbs sampler |
| `DS` | Exhaustive direct sampler for very small (N) |
| `DS-Const 0.1` | Constrained direct sampler using a 10% subset |
| `DS-Const 0.25` | Constrained direct sampler using a 25% subset |
| `DS-Const 0.5` | Constrained direct sampler using a 50% subset |
| `DS-Const-MAP` | MAP-refined constrained sampler; Y-only study |
| `DS-ML` | Candidate set generated using clustering algorithms |
| `MC-MCMC` | Collapsed label chain followed by composition draws |

An additional blocked Gibbs reference implementation is retained in the Y-only simulation code but is excluded from the manuscript plots.

## Reproducibility notes

- The simulation drivers call `sessionInfo()` so the R version and loaded package versions are recorded in the log.
- Result filenames include a timestamp and, under the SLURM scripts, a job-specific tag to avoid array tasks overwriting one another.
- For a fixed `N`, `k`, and replicate number, the current seed construction reuses the same underlying labels and standardized random errors across separation settings. This is a common-random-numbers design rather than independent simulation across separation settings.
- The large-overlap setting uses `(-2/3, 0, 2/3)` for $k=3$.
- The Y-only DS-ML implementation switches to scalable clustering algorithms for the largest sample size. The XY DS-ML implementation uses scalable mclust and spectral-clustering settings at large (N).
- Simulation results should be checked before files from multiple jobs are combined. Placing reruns of the same design cell in one result directory creates duplicate replicates, which `check_results.R` will flag.
- Each simulation task writes its result file after all replicates assigned to that task finish. For long full-study jobs, interruption before completion can therefore lose the unfinished task's in-memory progress.

## Recommended workflow

1. Install and verify all packages, including the required `miscPack::gaussian_mixture()` implementation.
2. Run the five-replicate Y and XY tests.
3. Run `run_analysis_y.sh` and `run_analysis_xy.sh` and resolve every failed integrity check.
4. Inspect the numerical summaries and figures.
5. Run the full 50-replicate studies into new output directories.
6. Run the analysis scripts with `FMM_REPS=50` and the corresponding final result directories.

## Citation

If you use this code, please cite the accompanying manuscript:

> Clinch, Madelyn, Jonathan R. Bradley, Andrés F. Barrientos, and Garritt L. Page (2026). “Bayesian Analysis Using a Constrained Mixture of Normal-Inverse-Gamma Models.” *arXiv preprint* [arXiv:2606.23435](https://arxiv.org/abs/2606.23435).

BibTeX:

```bibtex
@misc{clinch2026bayesian,
  title         = {Bayesian Analysis Using a Constrained Mixture of
                   Normal-Inverse-Gamma Models},
  author        = {Clinch, Madelyn and Bradley, Jonathan R. and
                   Barrientos, Andr\'{e}s F. and Page, Garritt L.},
  year          = {2026},
  eprint        = {2606.23435},
  archivePrefix = {arXiv},
  url           = {https://arxiv.org/abs/2606.23435}
}
```
