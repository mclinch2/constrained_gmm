#!/usr/bin/env Rscript

## ============================================================================
## Simulation study: Gaussian mixture of regressions with a NIG prior
## (Section 3.4 and Appendix E.2 of the manuscript).
##
## Design.  Sample sizes N in {10, 500, 1000, 10000} are crossed with three
## component-separation settings, with 50 replicate data sets per cell.  Each
## replicate draws from an equally weighted three-component mixture of simple
## linear regressions -- the components share the slope and differ in the
## intercept -- and fits every competing sampler to the same data.  The
## covariate matrix is held fixed within a cell, so inference is conditional
## on X and only the response is replicated.
##
## Samplers, with the name each one carries in the manuscript:
##   mcmc              MCMC            full_mcmc_sampling_regression_dependent()
##   ds                Direct Sampler  direct_sampling_xy(Z_sampling = "ds"),
##                                     N <= 13 only
##   ds-gibbs          MC-MCMC         direct_sampling_gibbs_regression()
##   ds{10,25,50}-obs  DS-Const        direct_sampling_xy(Z_sampling = "ds_obs")
##   ds-ml             DS-ML           direct_sampling_xy(Z_sampling = "ml")
##
## Metrics.  Mean and standard deviation of the adjusted Rand index against the
## true labels, effective sample size for the first two posterior moments of
## the conditional mixture, the conditional Kolmogorov-Smirnov and L2 distances
## of Appendix E.2, and wall-clock time in seconds.  Column names match the
## intercept-only study: KSmean_* and L2mean_* are the distances of the
## posterior MEAN curve, sdKS_* and sdL2_* the spread across draws.
##
## Usage.  All arguments are optional and positional; each list is comma
## separated.
##   Rscript sim_study_xy.R [N] [k] [clust_sep] [sig_const] [prior_M] [reps] [cores]
## The built-in default for N is 500 alone, so reproducing the published study
## means passing the full set explicitly:
##   Rscript sim_study_xy.R 10,500,1000,10000
##
## Output.  One tab-delimited row per (cell, replicate), time-stamped, written
## to the directory named in the CONFIG block below.
## ============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(salso)        # ARI(), enumerate.partitions()
  library(coda)         # effectiveSize()
  library(mvtnorm)      # rmvnorm() for the beta draws
  library(parallel)     # mclapply() over parameter cells
  library(miscPack)     # gaussian_mixture()
  library(mclust)       # DS-ML: model-based clustering
  library(kernlab)      # DS-ML: spectral clustering
  library(cluster)      # DS-ML: pam()
  library(dbscan)       # DS-ML: density-based clustering (optional generator)
  library(e1071)
  library(matrixStats)
})

## ---------------------------------------------------------------- CONFIG ----
## By default the companion functions file is taken from this script's own
## directory and results are written to ./sim_results_xy.  Either path may be
## relocated without editing the script, which is how the study was run on a
## cluster:
##   export FMM_SRC=/path/to/sim_xy_functions.R
##   export FMM_OUT=/path/to/sim_results_xy
.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
src_file <- Sys.getenv("FMM_SRC", unset = file.path(.script_dir, "sim_xy_functions.R"))
out_dir  <- Sys.getenv("FMM_OUT", unset = file.path(getwd(), "sim_results_xy"))
if (!file.exists(src_file))
  stop("functions file not found: ", src_file, "  (set FMM_SRC)")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
## Override with FMM_SUB_NITER / FMM_SUB_NBURN / FMM_SUB_NTHIN.
sub_niter <- as.integer(Sys.getenv("FMM_SUB_NITER", unset = "110000"))
sub_nburn <- as.integer(Sys.getenv("FMM_SUB_NBURN", unset = "10000"))
sub_nthin <- as.integer(Sys.getenv("FMM_SUB_NTHIN", unset = "100"))
cat("sourcing: ", src_file, "\nwriting to: ", out_dir,
    "\nsubset fit: niter=", sub_niter, " nburn=", sub_nburn, " nthin=", sub_nthin,
    "\n", sep = "")
source(src_file)
sessionInfo()   # recorded in the log so the run is reproducible

## ---------------------------------------------------------- design settings --
## Positional command-line arguments override any of the defaults below; see
## the usage line in the file header.  clust_sep and prior_M are codes whose
## meanings are given in sim_data_xy() and ss_study_xy() respectively.
args <- commandArgs(TRUE)
N_set <- if (length(args) >= 1) as.numeric(strsplit(args[1], ",")[[1]]) else c(500)
k_set <- if (length(args) >= 2) as.numeric(strsplit(args[2], ",")[[1]]) else c(3)
clust_sep_set <- if (length(args) >= 3) as.numeric(strsplit(args[3], ",")[[1]]) else c(1, 2, 3)
sig_const_set <- if (length(args) >= 4) as.numeric(strsplit(args[4], ",")[[1]]) else c(0.25)
prior_M_set <- if (length(args) >= 5) as.numeric(strsplit(args[5], ",")[[1]]) else c(2)
n_reps <- if (length(args) >= 6) as.numeric(args[6]) else 50
mc_cores <- if (length(args) >= 7) as.numeric(args[7]) else 8

## One row per design cell; cells are distributed across cores by mclapply()
## at the bottom of the script.
param_grid <- expand.grid(
  k = k_set,
  N = N_set,
  clust_sep = clust_sep_set,
  sig_const = sig_const_set,
  prior_M = prior_M_set,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
param_list <- split(param_grid, seq(nrow(param_grid)))
param_list <- lapply(param_list, as.list)

cat("Ns:      ", paste(N_set, collapse = ", "), "\n")
cat("reps:    ", n_reps, "\n")
cat("cores:   ", mc_cores, "\n")

## Draw one replicate from an equally weighted k-component mixture of simple
## linear regressions.  The components share the slope (fixed at 1) and are
## separated only through their intercepts, so clust_sep controls how far apart
## the parallel regression lines lie.
##
## @param k          number of true components.
## @param N          sample size.
## @param clust_sep  separation code, reported in the manuscript as
##                   1 = moderate overlap, 2 = large overlap, 3 = no overlap;
##                   the intercepts are (-1, 0, 1), (-2/3, 0, 2/3) and
##                   (-2, 0, 2) respectively at k = 3. 
## @param sig_const  component standard deviation sigma_k 
## @param xi         N x 2 covariate matrix (intercept and one N(0,1) covariate),
##                   supplied by the caller so it stays fixed across replicates.
## @return data.frame(id, xi.1, xi.2, yi, si), where si holds the true labels.
sim_data_xy <- function(k=3, N=10, clust_sep = 1, sig_const = 0.25, xi=NULL){
  if (is.null(xi)) xi <- cbind(1, rnorm(N))
  if(clust_sep == 1){
    beta0k <- seq(-floor(k/2), floor(k/2))            # moderate: (-1, 0, 1)
  }else if(clust_sep == 2){
    beta0k <- seq(-floor(k/2), floor(k/2)) / (k/2)    # large:  (-2/3, 0, 2/3)
  }else if(clust_sep == 3){
    beta0k <- seq(-2 * floor(k/2), 2 * floor(k/2), by = 2)  # none: (-2, 0, 2)
  }else{
    beta0k <- seq(-floor(k/2), floor(k/2))            # fall back to moderate
  }
  si <- sample(1:k, N, replace=TRUE, prob=rep(1/k,k)) # true component labels
  sig <- rep(sig_const, k)                           # component sd's
  beta1k <- rep(1, k)                                # common slope
  betak <- cbind(beta0k, beta1k)
  yi <- rowSums(sapply(1:k, function(j) {
    rnorm(N, xi %*% betak[j,], sig[j]) * (si == j)
  }))
  
  df <- data.frame(id = 1:N, xi = xi, yi = yi, si = si)
  return(df)
}

## ---------------------------------------------------------------------------
## Truth and goodness-of-fit utilities
## ---------------------------------------------------------------------------

true_params_xy <- function(k, clust_sep, sig_const) {
  if (clust_sep == 1) {
    beta0k <- seq(-floor(k/2), floor(k/2))
  } else if (clust_sep == 2) {
    beta0k <- seq(-floor(k/2), floor(k/2)) / (k/2)
  } else if (clust_sep == 3) {
    beta0k <- seq(-2 * floor(k/2), 2 * floor(k/2), by = 2)
  } else {
    beta0k <- seq(-floor(k/2), floor(k/2))
  }
  beta1k <- rep(1, k)
  beta  <- cbind(beta0k, beta1k)                 # k x p (p=2 here)
  sigma2 <- rep(sig_const^2, k)
  pi <- rep(1/k, k)                              # weights
  list(pi = pi, beta = beta, sigma2 = sigma2)
}

## Probability integral transform u_i = F(Y_i | x_i) under a normal mixture
## regression.  Retained as a diagnostic; the reported metrics use
## calc_obs_KS_L2() below.
mix_pit <- function(y, X, pi, beta, sigma2) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  
  pi <- as.numeric(pi)
  pi <- pi / sum(pi)
  
  beta <- as.matrix(beta)                        # K x p
  K <- nrow(beta)
  
  sigma2 <- as.numeric(sigma2)
  sigma2 <- pmax(sigma2, 1e-12)
  sdvec <- sqrt(sigma2)
  
  Mu <- X %*% t(beta)                            # n x K
  
  n <- length(y)
  ymat  <- matrix(y, nrow = n, ncol = K)
  sdmat <- matrix(sdvec, nrow = n, ncol = K, byrow = TRUE)
  
  P <- pnorm((ymat - Mu) / sdmat)                # n x K
  u <- as.numeric(P %*% pi)                      # length n
  pmin(pmax(u, 0), 1)
}

## One-sample Kolmogorov-Smirnov statistic of u against Uniform(0, 1),
## the companion to mix_pit().
ks_unif <- function(u) {
  u <- sort(pmin(pmax(as.numeric(u), 0), 1))
  n <- length(u)
  i <- seq_len(n)
  max( max(i/n - u), max(u - (i - 1)/n) )
}

## Extract (pi, beta, sigma2) for a single posterior draw.  beta is a 3-D
## array whose draw index may be the first or the last margin depending on the
## sampler, so the layout is inferred from which margin has length nsamp.
##
## @return list(pi, beta [K x p], sigma2).
extract_draw_params <- function(obj, t, nsamp) {
  if (is.null(obj$pi) || is.null(obj$beta) || is.null(obj$sigma2)) {
    stop("Object must contain $pi, $beta, $sigma2.")
  }
  pi <- as.numeric(obj$pi[t, ])
  s2 <- as.numeric(obj$sigma2[t, ])
  
  bA <- obj$beta
  d <- dim(bA)
  if (length(d) != 3) stop("beta must be a 3D array.")
  
  if (d[1] == nsamp) {          # [draw, K, p]
    beta <- bA[t, , , drop = FALSE]
    beta <- matrix(beta, nrow = d[2], ncol = d[3])   # K x p
  } else if (d[3] == nsamp) {   # [K, p, draw]
    beta <- bA[, , t, drop = FALSE]
    beta <- matrix(beta, nrow = d[1], ncol = d[2])   # K x p
  } else {
    stop("Cannot infer beta layout: expected dim(beta)[1]==nsamp or dim(beta)[3]==nsamp.")
  }
  
  list(pi = pi, beta = beta, sigma2 = s2)
}

## ---------------------------------------------------------------------------
## Appendix E.2: the conditional KS and L2 distances, evaluated at the OBSERVED
## design points (x_i, Y_i) rather than on a grid.
##
##   F^[j](Y_i|x_i) = sum_k pi_k^[j] Phi( (Y_i - x_i' beta_k^[j]) / sigma_k^[j] )
##   Fhat(Y_i|x_i)  = (1/nsamp) sum_j F^[j](Y_i|x_i)
##   KS = max_{i=1..N} | F_true(Y_i|x_i) - Fhat(Y_i|x_i) |
##   L2 = { (1/N) sum_i [ f_true(Y_i|x_i) - fhat(Y_i|x_i) ]^2 }^(1/2)
##
## Evaluating at the observations rather than on a grid is what makes this
## affordable at the largest sample size: at N = 10,000 with M = 25 it uses all
## N observations and all 1,000 draws in roughly 28 seconds.
##
## Note that KS is a maximum over i, an extreme order statistic, so it drifts
## upward with N and is not directly comparable across the N facets.
## ---------------------------------------------------------------------------

## Conditional CDF and density of a normal mixture regression, evaluated at
## each observed (x_i, Y_i) for one parameter value.
mix_F_f_obs <- function(y, X, pi, beta, sigma2) {
  X    <- as.matrix(X); N <- length(y)
  pi   <- as.numeric(pi); pi <- pi / sum(pi)
  beta <- as.matrix(beta)
  sdv  <- sqrt(pmax(as.numeric(sigma2), 1e-12))
  Mu   <- X %*% t(beta)                       # N x K
  Z    <- (y - Mu) / rep(sdv, each = N)
  list(F = as.vector(pnorm(Z) %*% pi),
       f = as.vector((dnorm(Z) / rep(sdv, each = N)) %*% pi))
}

## KS and L2 of the posterior mean curve against the truth, plus the standard
## deviation of the corresponding per-draw distances.
##
## @param obj        a sampler's output; needs $pi, $beta and $sigma2.
## @param truth_par  output of true_params_xy().
## @return list(KS, L2, KS_draw_sd, L2_draw_sd).
calc_obs_KS_L2 <- function(obj, truth_par, y, X, nsamp) {
  N  <- length(y)
  tr <- mix_F_f_obs(y, X, truth_par$pi, truth_par$beta, truth_par$sigma2)

  Fsum <- fsum <- numeric(N)
  KS_draw <- L2_draw <- numeric(nsamp)

  for (t in seq_len(nsamp)) {
    par_t <- extract_draw_params(obj, t, nsamp = nsamp)
    cur   <- mix_F_f_obs(y, X, par_t$pi, par_t$beta, par_t$sigma2)
    Fsum  <- Fsum + cur$F
    fsum  <- fsum + cur$f
    KS_draw[t] <- max(abs(tr$F - cur$F))
    L2_draw[t] <- sqrt(mean((tr$f - cur$f)^2))
  }

  list(KS         = max(abs(tr$F - Fsum / nsamp)),
       L2         = sqrt(mean((tr$f - fsum / nsamp)^2)),
       KS_draw_sd = sd(KS_draw),
       L2_draw_sd = sd(L2_draw))
}

## Run every sampler on n_reps replicate data sets from one design cell.
##
## @param k,N,clust_sep,sig_const  passed through to sim_data_xy().
## @param prior_M  code for the number of mixture components M carried by the
##                 fitted model: 1 gives M = k + 2, 2 gives M = 25 (the setting
##                 used in the manuscript), anything else gives M = k.
## @return list(summary, K_draws), where summary is an n_reps x length(coln)
##         matrix of metrics and K_draws holds the per-draw number of occupied
##         components for each sampler.
ss_study_xy <- function(k, N, clust_sep, sig_const, prior_M){
  ndata <- n_reps
  ## Posterior sample size retained per sampler.  
  nsamp <- 1000
  methods <- c("mcmc","ds","ds-gibbs","ds10-obs","ds25-obs","ds50-obs","ds-ml")
  ## KSmean_/L2mean_ are the distances of the posterior mean curve, matching
  ## the column names used by the intercept-only study.
  metrics <- c("mnARI","sdARI","speed",
               "essMean","essM2",
               "KSmean","sdKS",
               "L2mean","sdL2")
  
  coln <- c(outer(metrics, methods, paste, sep = "_"))
  coln <- c(coln, "N","k","clust_sep","sig_const","M","rep")
  
  out <- matrix(NA_real_, nrow = ndata, ncol = length(coln),
                dimnames = list(NULL, coln))
  ## Number of occupied components K per draw, kept for the posterior-of-K
  ## summaries.  The exhaustive direct sampler is only feasible for N <= 13.
  K_store <- list(
    mcmc        = matrix(NA_integer_, ndata, nsamp),
    ds          = if (N <= 13) matrix(NA_integer_, ndata, nsamp) else NULL,
    ds_gibbs    = matrix(NA_integer_, ndata, nsamp),
    ds_obs10    = matrix(NA_integer_, ndata, nsamp),
    ds_obs25    = matrix(NA_integer_, ndata, nsamp),
    ds_obs50    = matrix(NA_integer_, ndata, nsamp),
    ds_ml       = matrix(NA_integer_, ndata, nsamp)
  )
  ## The covariate matrix is held fixed across replicates -- only the response
  ## is replicated -- so inference is conditional on X.  Seeding on (N, k) alone
  ## means every clust_sep / sig_const / prior_M cell shares the same X, and X
  ## is reproducible from the recorded design.
  set.seed(90210 + 10000*N + 100*k)
  xi <- cbind(rep(1, N), rnorm(N))
  ml_large_cutoff <- 1000   
  
  for(d in 1:ndata){
    ## A failure in any one sampler costs only this replicate: the row is
    ## blanked, its design columns are restored, and the loop continues.
    tryCatch({
      ## Seed depends on the cell and the replicate, so every method sees the
      ## same data and reruns of a single cell are reproducible.
      set.seed(100000 + 10000*N + 100*k + d)
      cat("N =", N, " d =", d, "clust_sep = ", clust_sep, "sig_const = ", sig_const, "prior_M = ", prior_M, "k = ", k, "\n")
      df <- sim_data_xy(k=k, N=N, clust_sep = clust_sep, sig_const = sig_const, xi)
      si <- df$si
      yi <- df$yi
      
      ## Prior hyperparameters of Equation (1): symmetric Dirichlet weight
      ## alpha, prior scale sig2_beta on the regression coefficients, and the
      ## inverse-gamma shape and scale (omega, kappa) on the variances.
      omega <- 2; kappa <-1
      alpha <- 1; sig2_beta <- 10
      if(prior_M == 1){
        M <- k + 2 # slightly more components than the truth
      }else if(prior_M == 2){
        M <- 25    # far more components than the truth (manuscript setting)
      }else{
        M <- k     # components fixed at the truth
      }
      prior_for_M <- rep(1,M)/M   # uniform prior on K for the direct sampler

      ## ------------------------- fit every sampler to this replicate -------
      ## Timings bracket each call individually and are reported in seconds.

      ## DS-ML -- k-means, EM and spectral clustering at every N, over
      ## K = 2,...,M, matching Section 2.4.2 and Section 3.4 (3 algorithms,
      ## (M-1) K values).  
      ml_big <- N > ml_large_cutoff
      time1 <- Sys.time()
      ds_ml <- direct_sampling_xy(y = yi, X = xi, M = M, Z_sampling = "ml",
                                  n_sub = 1,
                                  ml_algorithms = c("gmm", "spectral", "kmeans"),
                                  ml_features = "yx", ml_scale = TRUE,
                                  ml_fast = ml_big,
                                  spectral_nystrom_sample = if (ml_big) 400 else NULL,
                                  seed = d*11111, nsamp = nsamp,
                                  sig2_beta = sig2_beta, omega = omega, kappa = kappa,
                                  ml_K_seq = 2:M)
      time2 <- Sys.time()
      out[d,"speed_ds-ml"] <- as.numeric(time2 - time1, units = "secs")

      ## DS-Const with a subset of 10% of the observations.
      time1 <- Sys.time()
      ds_postobs10 <- direct_sampling_xy(y=yi, X=xi, M=M, alpha=alpha,
                                         sig2_beta=sig2_beta, omega=omega, kappa=kappa,
                                         Z_sampling = "ds_obs", nsamp=nsamp, n_sub=0.1,
                                         sub_niter=sub_niter, sub_nburn=sub_nburn, sub_nthin=sub_nthin)
      time2 <- Sys.time()
      out[d,"speed_ds10-obs"] <- as.numeric(time2 - time1, units = "secs")

      ## DS-Const, 25% subset.
      time1 <- Sys.time()
      ds_postobs25 <- direct_sampling_xy(y=yi, X=xi, M=M, alpha=alpha,
                                         sig2_beta=sig2_beta, omega=omega, kappa=kappa,
                                         Z_sampling = "ds_obs", nsamp=nsamp, n_sub=0.25,
                                         sub_niter=sub_niter, sub_nburn=sub_nburn, sub_nthin=sub_nthin)
      time2 <- Sys.time()
      out[d,"speed_ds25-obs"] <- as.numeric(time2 - time1, units = "secs")

      ## DS-Const, 50% subset.
      time1 <- Sys.time()
      ds_postobs50 <- direct_sampling_xy(y=yi, X=xi, M=M, alpha=alpha,
                                         sig2_beta=sig2_beta, omega=omega, kappa=kappa,
                                         Z_sampling = "ds_obs", nsamp=nsamp, n_sub=0.5,
                                         sub_niter=sub_niter, sub_nburn=sub_nburn, sub_nthin=sub_nthin)
      time2 <- Sys.time()
      out[d,"speed_ds50-obs"] <- as.numeric(time2 - time1, units = "secs")

      ## MCMC: nthin = (niter - nburn) / nsamp = (40000 - 10000) / 30 = 1000 draws
      time1 <- Sys.time()
      fit <- full_mcmc_sampling_regression_dependent(y=yi, X=xi, M=M, alpha=alpha,
                                                     omega=omega, kappa=kappa, sig2_beta=sig2_beta,
                                                     niter=40000, nburn=10000, nthin=30)
      time2 <- Sys.time()
      out[d,"speed_mcmc"] <- as.numeric(time2 - time1, units = "secs")

      ## Direct Sampler: enumerates every set partition of the N observations,
      ## so it is only tractable for small N.  Section 3.1 places a discrete
      ## uniform prior on K in {1,...,M}, with M capped at N.
      if(N <= 13){
        Mcap <- min(M, N)
        time1 <- Sys.time()
        ds <- direct_sampling_xy(y=yi, X=xi, M=Mcap, alpha=alpha,
                                 sig2_beta=sig2_beta, omega=omega, kappa=kappa,
                                 Z_sampling = "ds", nsamp=nsamp,
                                 prior_for_M = rep(1, Mcap)/Mcap)
        time2 <- Sys.time()
        out[d,"speed_ds"] <- as.numeric(time2 - time1, units = "secs")
      }

      ## MC-MCMC: collapsed label Gibbs sampler of Equation (8), followed by
      ## composition draws of the parameters.
      time1 <- Sys.time()
      ds_MCMC <- direct_sampling_gibbs_regression(yi, xi, M = M, alpha = alpha,
                                                  sig2_beta = sig2_beta, omega = omega, kappa = kappa, nsamp = nsamp)
      time2 <- Sys.time()
      out[d,"speed_ds-gibbs"] <- as.numeric(time2 - time1, units = "secs")

      ## ------------------------------- summarize this replicate ----------

      ## Adjusted Rand index of each posterior label draw against the truth.
      ARI.fit <- sapply(1:nsamp, function(x) ARI(fit$z[x,], si))
      if(N <= 13) ARI.ds <- sapply(1:nsamp, function(x) ARI(ds$Z[x,], si))
      ARI.ds_MCMC <- sapply(1:nsamp, function(x) ARI(ds_MCMC$Z[x,], si))
      ARI.ds_postobs10 <- sapply(1:nsamp, function(x) ARI(ds_postobs10$Z[x,], si))
      ARI.ds_postobs25 <- sapply(1:nsamp, function(x) ARI(ds_postobs25$Z[x,], si))
      ARI.ds_postobs50 <- sapply(1:nsamp, function(x) ARI(ds_postobs50$Z[x,], si))
      ARI.ds_ml <- sapply(1:nsamp, function(x) ARI(ds_ml$Z[x,], si))

      K_store$mcmc[d, ]     <- as.integer(fit$kpost)
      if (N <= 13) K_store$ds[d, ] <- as.integer(ds$K)
      K_store$ds_gibbs[d, ]  <- as.integer(ds_MCMC$K)
      K_store$ds_obs10[d, ] <- as.integer(ds_postobs10$K)
      K_store$ds_obs25[d, ] <- as.integer(ds_postobs25$K)
      K_store$ds_obs50[d, ] <- as.integer(ds_postobs50$K)
      K_store$ds_ml[d, ] <- as.integer(ds_ml$K)

      out[d, "mnARI_mcmc"] <- mean(ARI.fit); out[d, "sdARI_mcmc"] <- sd(ARI.fit)
      if (N <= 13) { out[d, "mnARI_ds"] <- mean(ARI.ds); out[d, "sdARI_ds"] <- sd(ARI.ds) }
      out[d, "mnARI_ds-gibbs"]  <- mean(ARI.ds_MCMC);    out[d, "sdARI_ds-gibbs"]  <- sd(ARI.ds_MCMC)
      out[d, "mnARI_ds10-obs"]   <- mean(ARI.ds_postobs10); out[d, "sdARI_ds10-obs"]   <- sd(ARI.ds_postobs10)
      out[d, "mnARI_ds25-obs"]   <- mean(ARI.ds_postobs25); out[d, "sdARI_ds25-obs"]   <- sd(ARI.ds_postobs25)
      out[d, "mnARI_ds50-obs"]   <- mean(ARI.ds_postobs50); out[d, "sdARI_ds50-obs"]   <- sd(ARI.ds_postobs50)
      out[d, "mnARI_ds-ml"]   <- mean(ARI.ds_ml); out[d, "sdARI_ds-ml"]   <- sd(ARI.ds_ml)

      ## Effective sample size of the conditional mixture mean, averaged over
      ## the design points: mean_i sum_k pi_k x_i' beta_k.  Uses the spectral
      ## estimator in coda::effectiveSize().
      out[d, "essMean_mcmc"] <- effectiveSize(vapply(seq_len(nsamp), function(t) {
        B <- fit$beta[t, , ]
        w <- fit$pi[t, ]
        w <- w / sum(w)
        Mu <- xi %*% t(B)
        m  <- as.vector(Mu %*% w)
        mean(m)
      }, numeric(1)))
      
      
      if (N <= 13) out[d, "essMean_ds"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B  <- ds$beta[t,,]            # M x p
        Mu <- xi %*% t(B)             # N x M
        w  <- as.numeric(ds$pi[t, ])
        w  <- w / sum(w)
        m  <- as.vector(Mu %*% w)     # N
        mean(m)
      }))

      out[d, "essMean_ds-gibbs"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B  <- ds_MCMC$beta[t,,]
        Mu <- xi %*% t(B)
        w  <- as.numeric(ds_MCMC$pi[t, ])
        w  <- w / sum(w)
        m  <- as.vector(Mu %*% w)
        mean(m)
      }))

      out[d, "essMean_ds10-obs"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B  <- ds_postobs10$beta[t,,]
        Mu <- xi %*% t(B)
        w  <- as.numeric(ds_postobs10$pi[t,])
        w  <- w / sum(w)
        m  <- as.vector(Mu %*% w)
        mean(m)
      }))

      out[d, "essMean_ds25-obs"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B  <- ds_postobs25$beta[t,,]
        Mu <- xi %*% t(B)
        w  <- as.numeric(ds_postobs25$pi[t,]); w <- w / sum(w)
        mean(as.vector(Mu %*% w))
      }))
      
      out[d, "essMean_ds50-obs"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B  <- ds_postobs50$beta[t,,]
        Mu <- xi %*% t(B)
        w  <- as.numeric(ds_postobs50$pi[t,]); w <- w / sum(w)
        mean(as.vector(Mu %*% w))
      }))

      out[d, "essMean_ds-ml"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B  <- ds_ml$beta[t,,]
        Mu <- xi %*% t(B)
        w  <- as.numeric(ds_ml$pi[t, ])
        w  <- w / sum(w)
        m  <- as.vector(Mu %*% w)
        mean(m)
      }))

      ## Effective sample size of the conditional second moment,
      ## mean_i sum_k pi_k ((x_i' beta_k)^2 + sigma2_k).
      out[d, "essM2_mcmc"] <- effectiveSize(vapply(seq_len(nsamp), function(t) {
        B  <- fit$beta[t, , ]
        w  <- fit$pi[t, ];  w <- w / sum(w)
        s2 <- fit$sigma2[t, ]
        MU   <- xi %*% t(B)
        term <- MU^2 + matrix(s2, nrow = N, ncol = length(s2), byrow = TRUE)
        m2_i <- as.vector(term %*% w)
        mean(m2_i)
      }, numeric(1)))
      
      if (N <= 13) out[d, "essM2_ds"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B   <- ds$beta[t,,]
        Mu  <- xi %*% t(B)
        w   <- as.numeric(ds$pi[t, ])
        w   <- w / sum(w)
        m2_i <- (Mu^2) %*% w + sum(w * as.numeric(ds$sigma2[t, ]))
        mean(as.vector(m2_i))
      }))
      
      out[d, "essM2_ds-gibbs"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B   <- ds_MCMC$beta[t,,]
        Mu  <- xi %*% t(B)
        w   <- as.numeric(ds_MCMC$pi[t, ])
        w   <- w / sum(w)
        m2_i <- (Mu^2) %*% w + sum(w * as.numeric(ds_MCMC$sigma2[t, ]))
        mean(m2_i)
      }))
      
      out[d, "essM2_ds10-obs"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B   <- ds_postobs10$beta[t,,]
        Mu  <- xi %*% t(B)
        w   <- as.numeric(ds_postobs10$pi[t, ])
        w   <- w / sum(w)
        m2_i <- (Mu^2) %*% w + sum(w * as.numeric(ds_postobs10$sigma2[t, ]))
        mean(as.vector(m2_i))
      }))
      
      out[d, "essM2_ds25-obs"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B   <- ds_postobs25$beta[t,,]
        Mu  <- xi %*% t(B)
        w   <- as.numeric(ds_postobs25$pi[t, ])
        w   <- w / sum(w)
        m2_i <- (Mu^2) %*% w + sum(w * as.numeric(ds_postobs25$sigma2[t, ]))
        mean(as.vector(m2_i))
      }))
      
      out[d, "essM2_ds50-obs"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B   <- ds_postobs50$beta[t,,]
        Mu  <- xi %*% t(B)
        w   <- as.numeric(ds_postobs50$pi[t, ])
        w   <- w / sum(w)
        m2_i <- (Mu^2) %*% w + sum(w * as.numeric(ds_postobs50$sigma2[t, ]))
        mean(as.vector(m2_i))
      }))

      out[d, "essM2_ds-ml"] <- effectiveSize(sapply(seq_len(nsamp), function(t) {
        B   <- ds_ml$beta[t,,]
        Mu  <- xi %*% t(B)
        w   <- as.numeric(ds_ml$pi[t, ])
        w   <- w / sum(w)
        m2_i <- (Mu^2) %*% w + sum(w * as.numeric(ds_ml$sigma2[t, ]))
        mean(as.vector(m2_i))
      }))
      
      
      ## Conditional KS and L2 against the truth at the observed points
      ## (Appendix E.2), using all N observations and all nsamp draws.
      truth_par <- true_params_xy(k = k, clust_sep = clust_sep, sig_const = sig_const)

      store_cond_metrics <- function(obj, method_name) {
        rr <- calc_obs_KS_L2(obj = obj, truth_par = truth_par,
                             y = yi, X = xi, nsamp = nsamp)
        out[d, paste0("KSmean_", method_name)] <<- rr$KS
        out[d, paste0("L2mean_", method_name)] <<- rr$L2
        out[d, paste0("sdKS_",   method_name)] <<- rr$KS_draw_sd
        out[d, paste0("sdL2_",   method_name)] <<- rr$L2_draw_sd
      }

      store_cond_metrics(list(pi = fit$pi, beta = fit$beta, sigma2 = fit$sigma2), "mcmc")
      if (N <= 13) {
        store_cond_metrics(list(pi = ds$pi, beta = ds$beta, sigma2 = ds$sigma2), "ds")
      }
      store_cond_metrics(list(pi = ds_MCMC$pi, beta = ds_MCMC$beta, sigma2 = ds_MCMC$sigma2), "ds-gibbs")
      store_cond_metrics(list(pi = ds_postobs10$pi, beta = ds_postobs10$beta, sigma2 = ds_postobs10$sigma2), "ds10-obs")
      store_cond_metrics(list(pi = ds_postobs25$pi, beta = ds_postobs25$beta, sigma2 = ds_postobs25$sigma2), "ds25-obs")
      store_cond_metrics(list(pi = ds_postobs50$pi, beta = ds_postobs50$beta, sigma2 = ds_postobs50$sigma2), "ds50-obs")
      store_cond_metrics(list(pi = ds_ml$pi, beta = ds_ml$beta, sigma2 = ds_ml$sigma2), "ds-ml")

      ## Design columns identifying this row.
      out[d, c("N","k","clust_sep","sig_const","M","rep")] <- c(N, k, clust_sep, sig_const, M, d)
    },
    error = function(e){
      message("ERROR on d=", d, ": ", conditionMessage(e))
      cl <- conditionCall(e)
      if (!is.null(cl))
        message("  in call: ", paste(deparse(cl), collapse = " "))
      out[d, ] <<- NA_real_
      out[d, c("N","k","clust_sep","sig_const","M","rep")] <<- c(N, k, clust_sep, sig_const, M, d)
    })
  }
  list(summary = out, K_draws = K_store)
}

## ------------------------------------------------------------------- run ----
## Cells are independent, so they are spread across cores; mc.set.seed gives
## each worker its own stream, and the seeds set inside ss_study_xy() keep the
## generated data reproducible regardless of the scheduling.
set.seed(1998)
ss_out <- mclapply(param_list, function(p) {
  ss_study_xy(k = p$k, N = p$N, clust_sep = p$clust_sep, sig_const = p$sig_const, prior_M = p$prior_M)
}, mc.cores = mc_cores, mc.preschedule = FALSE, mc.set.seed = TRUE)

## A worker that dies outright (segfault, out of memory, or an error raised
## outside the per-replicate tryCatch) comes back as a try-error.  Report those
## cells loudly and keep the ones that finished, rather than losing the whole
## run to a failed `[[`.
worker_ok <- !vapply(ss_out, inherits, logical(1), "try-error")
if (any(!worker_ok)) {
  message("WARNING: ", sum(!worker_ok), " of ", length(ss_out),
          " parameter cells failed entirely and are missing from the output:")
  for (i in which(!worker_ok)) {
    cnd <- attr(ss_out[[i]], "condition")
    message("  cell ", i, " (",
            paste(names(param_list[[i]]), unlist(param_list[[i]]), sep = "=", collapse = ", "),
            "): ", if (is.null(cnd)) "unknown" else conditionMessage(cnd))
  }
}
if (!any(worker_ok)) stop("every parameter cell failed; nothing to write.")
summaries <- do.call(rbind, lapply(ss_out[worker_ok], `[[`, "summary"))
res <- summaries

## ---------------------------------------------------------------- output ----
## One tab-delimited row per (cell, replicate).
##
## The filename is time-stamped to the second, which collides if several jobs
## finish together -- as array tasks do.  Set FMM_TAG (the SLURM scripts use the
## job and array-task id) to append a unique suffix.  The plotting script's glob
## is study_xy_saveK_results*.txt, so tagged files are still picked up.
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag   <- Sys.getenv("FMM_TAG", unset = "")
if (nzchar(tag)) stamp <- paste0(stamp, "_", gsub("[^A-Za-z0-9._-]", "-", tag))
outfile <- file.path(out_dir, sprintf("study_xy_saveK_results%s.txt", stamp))
write.table(res, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Wrote results to:", outfile, "\n")
save.image(file = file.path(out_dir, sprintf("MixSimSaveK_xy_%s.RData", stamp)))