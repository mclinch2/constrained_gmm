#!/usr/bin/env Rscript

## ============================================================================
## Simulation study: intercept-only Gaussian finite mixture with a NIG prior
## (Section 3.3 and Appendix E.1 of the manuscript).
##
## Design.  Sample sizes N in {10, 500, 1000, 10000} are crossed with three
## component-separation settings, with 50 replicate data sets per cell.  Each
## replicate draws from an equally weighted three-component normal mixture,
## fits every competing sampler to the same data, and records clustering,
## mixing, and density-estimation summaries.
##
## Samplers, with the name each one carries in the manuscript:
##   mcmc                 MCMC            full_mcmc_sampling_dependent()
##   ds                   Direct Sampler  direct_sampling()       (N <= 13 only)
##   ds-gibbs             MC-MCMC         direct_sampling_gibbs() 
##   ds{10,25,50}-obs     DS-Const        direct_sampling_postpred()
##   ds{10,25,50}-obsmap  DS-Const-MAP    direct_sampling_obs_refine()
##   ds-ml                DS-ML           direct_sampling_postpred_ML()
##   ds-mcmc              (not reported)  direct_sampling_MCMC(), a blocked
##                        Gibbs sampler retained for reference.  It is NOT the
##                        MC-MCMC of Equation (8) and is dropped downstream by
##                        plot_sim_y_paper.R.
##
## Metrics.  Mean and standard deviation of the adjusted Rand index against the
## true labels, effective sample size for the first two posterior moments of
## the mixture, the Kolmogorov-Smirnov distance to the empirical CDF, the L2
## distance to the true density, and wall-clock time in seconds.
##
## Two variants of each density metric are stored under similar names.  The
## manuscript reports the first pair; see compute_ks_ise() below.
##   KSmean_* / L2mean_*  KS and L2 of the posterior MEAN curve
##   mnKS_*   / mnL2_*    mean over draws of the per-draw KS and L2
##
## Usage.  All arguments are optional and positional; each list is comma
## separated.  Defaults reproduce the published study.
##   Rscript sim_study_y.R [N] [k] [clust_sep] [sig_const] [prior_M] [reps] [cores]
##
## Output.  One tab-delimited row per (cell, replicate) plus a workspace image,
## both time-stamped, written to the directory named in the CONFIG block below.
## ============================================================================

rm(list=ls())

suppressPackageStartupMessages({
  library(salso)      # ARI(), enumerate.partitions()
  library(coda)       # effectiveSize()
  library(parallel)   # mclapply() over parameter cells
  library(miscPack)   # gaussian_mixture(), the candidate generator
  library(mclust)     # DS-ML: model-based clustering
  library(kernlab)    # DS-ML: spectral clustering
  library(cluster)    # DS-ML: pam()/clara()
  library(dbscan)     # DS-ML: density-based clustering (optional generator)
  library(e1071)      # DS-ML: fuzzy c-means (optional generator)
})

## ---------------------------------------------------------------- CONFIG ----
## By default the companion functions file is taken from this script's own
## directory and results are written to ./sim_results_y.  Either path may be
## relocated without editing the script, which is how the study was run on a
## cluster:
##   export FMM_SRC=/path/to/sim_study_y_functions.R
##   export FMM_OUT=/path/to/sim_results_y
.script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
src_file <- Sys.getenv("FMM_SRC", unset = file.path(.script_dir, "sim_study_y_functions.R"))
out_dir  <- Sys.getenv("FMM_OUT", unset = file.path(getwd(), "sim_results_y"))
if (!file.exists(src_file))
  stop("functions file not found: ", src_file, "  (set FMM_SRC)")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("sourcing: ", src_file, "\nwriting to: ", out_dir, "\n", sep = "")
source(src_file)
sessionInfo()   # recorded in the log so the run is reproducible

## ---------------------------------------------------------- design settings --
## Positional command-line arguments override any of the defaults below; see
## the usage line in the file header.  clust_sep and prior_M are codes whose
## meanings are given in sim_data() and ss_study() respectively.
args <- commandArgs(TRUE)
N_set <- if (length(args) >= 1) as.numeric(strsplit(args[1], ",")[[1]]) else c(10, 500, 1000, 10000)
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

## Draw one replicate data set from an equally weighted k-component normal
## mixture with common component standard deviation.
##
## @param k          number of true components.
## @param N          sample size.
## @param clust_sep  separation code, reported in the manuscript as
##                   1 = moderate overlap, 2 = large overlap, 3 = no overlap.
##                   The large-overlap means are (-2/3, 0, 2/3)
## @param sig_const  component standard deviation sigma_k 
##                   The study uses 0.25, so sigma^2_k = 0.0625 
## @return list(data.frame(id, yi, si), true_pi, true_mu, true_sig), where si
##         holds the true component labels used for the ARI.
sim_data <- function(k=3, N=10, clust_sep = 1, sig_const = 0.25){
  if(clust_sep == 1){
    muk <- seq(-floor(k/2), floor(k/2))                 # moderate: (-1, 0, 1)
  }else if(clust_sep == 2){
    muk <- seq(-floor(k/2), floor(k/2)) / (k/2)         # large:  (-2/3, 0, 2/3)
  }else if(clust_sep == 3){
    muk <- seq(-2 * floor(k/2), 2 * floor(k/2), by = 2) # none:   (-2, 0, 2)
  }else{
    muk <- seq(-floor(k/2), floor(k/2))                 # fall back to moderate
  }
  si <- sample(1:k, N, replace=TRUE, prob=rep(1/k,k)) # true component labels
  sig <- rep(sig_const, k)                            # component sd's
  yi <- rowSums(sapply(1:k, function(j) {
    rnorm(N, muk[j], sig[j]) * (si == j)
  }))
  true_pi <- rep(1/k,k)
  true_mu <- muk
  true_sig <- sig
  df <- data.frame(id = 1:N, yi = yi, si = si)
  return(list(df, true_pi, true_mu, true_sig))
}

## Trapezoidal integral of y over the grid x; used for the integrated squared
## error between an estimated and the true density.
trapz <- function(x, y) { sum((y[-1] + y[-length(y)]) * diff(x)) * 0.5 }

## Mixture CDF at the points x for a single posterior draw (w, mu, sd),
## evaluated as a length(x) by length(w) matrix product.
mix_cdf_draw <- function(x, w, mu, sd) {
  rowSums(
    matrix(w,  nrow = length(x), ncol = length(w), byrow = TRUE) *
      pnorm(
        matrix(x,  nrow = length(x), ncol = length(w)),
        matrix(mu, nrow = length(x), ncol = length(w), byrow = TRUE),
        matrix(sd, nrow = length(x), ncol = length(w), byrow = TRUE)
      )
  )
}

## Density-estimation summaries for one fitted sampler.
##
## Each posterior draw (w, mu, sd2) defines a mixture CDF and density on the
## grid t_grid.  Two families of summary are returned, and they answer
## different questions:
##
##   * "_mean" / "_median": the KS and L2 distances of the POSTERIOR MEAN (or
##     median) curve.  These are the quantities defined in Appendix E.1 and
##     reported in the manuscript; they measure the accuracy of the Bayes
##     estimate of the density.
##   * "_avg" / "_sd": the mean and standard deviation of the per-draw
##     distances.  These describe the spread of the posterior over curves and
##     are retained for diagnostics only.
##
## @param w,mu,sd2  nsamp x M matrices of posterior weights, means and
##                  variances (vectors are promoted to one-row matrices).
##                  Weights are renormalised before use.
## @param t_grid    evaluation grid, spanning the observed range of y.
## @param F_emp     empirical CDF of the observed data on t_grid (KS target).
## @param f_true    true mixture density on t_grid (L2 target).
## @param do_plots,main_prefix  optional CDF and density overlays for
##                  inspecting a single fit; off in production runs.
## @return named list of KS, ISE and L2 summaries.
compute_ks_ise <- function(w, mu, sd2, t_grid, F_emp, f_true,
                           do_plots = FALSE, main_prefix = "") {
  sd <- sqrt(sd2)
  if (is.null(dim(mu))) mu <- matrix(mu, nrow = 1)
  if (is.null(dim(w)))  w  <- matrix(w,  nrow = nrow(mu))
  if (is.null(dim(sd))) sd <- matrix(sd, nrow = nrow(mu))
  ns <- nrow(mu); K <- ncol(mu); G <- length(t_grid)
  w <- w / pmax(rowSums(w), 1e-12)
  
  F_draws <- matrix(NA_real_, nrow = ns, ncol = G)
  f_draws <- matrix(NA_real_, nrow = ns, ncol = G)
  
  for (i in seq_len(ns)) {
    F_draws[i, ] <- mix_cdf_draw(t_grid, w[i, ], mu[i, ], sd[i, ])
    f_draws[i, ] <- rowSums(
      matrix(w[i, ], nrow = G, ncol = K, byrow = TRUE) *
        dnorm(
          matrix(t_grid, nrow = G, ncol = K),
          matrix(mu[i, ], nrow = G, ncol = K, byrow = TRUE),
          matrix(sd[i, ], nrow = G, ncol = K, byrow = TRUE)
        )
    )
  }
  
  ## Posterior mean and median curves, taken pointwise over the grid.
  F_mean <- colMeans(F_draws)
  f_mean <- colMeans(f_draws)
  F_median <- apply(F_draws, 2, median, na.rm = TRUE)
  f_median <- apply(f_draws, 2, median, na.rm = TRUE)
  
  ## KS distance: of the mean curve, of the median curve, then per draw.
  KS_mean   <- max(abs(F_mean   - F_emp))
  KS_median <- max(abs(F_median - F_emp))
  KS_by_draw <- apply(F_draws, 1, function(Fi) max(abs(Fi - F_emp)))
  KS_avg <- mean(KS_by_draw)
  KS_sd  <- stats::sd(KS_by_draw)
  
  ## Integrated squared error and its square root (L2), same three variants.
  ISE_mean   <- trapz(t_grid, (f_mean   - f_true)^2)
  ISE_median <- trapz(t_grid, (f_median - f_true)^2)
  L2_mean    <- sqrt(ISE_mean)
  L2_median  <- sqrt(ISE_median)
  
  ISE_by_draw <- apply(f_draws, 1, function(fi) trapz(t_grid, (fi - f_true)^2))
  L2_by_draw  <- sqrt(ISE_by_draw)
  ISE_avg <- mean(ISE_by_draw)
  ISE_sd  <- stats::sd(ISE_by_draw)
  L2_avg  <- mean(L2_by_draw)
  L2_sd   <- stats::sd(L2_by_draw)
  
  if (isTRUE(do_plots)) {
    ## CDF: empirical against the posterior mean and median
    plot(t_grid, F_emp, type = "s", lwd = 2, ylim = c(0, 1),
         xlab = "t", ylab = "CDF",
         main = paste0(main_prefix,
                       if (nzchar(main_prefix)) " " else "",
                       "Empirical vs Posterior Mean/Median CDF"))
    lines(t_grid, F_mean,   lwd = 2, lty = 2)
    lines(t_grid, F_median, lwd = 2, lty = 3)
    legend("bottomright",
           c("F_emp", "F_mean", "F_median"),
           lwd = 2, lty = c(1, 2, 3), bty = "n")
    
    ## Density: truth against the posterior mean and median
    ymax <- max(f_true, f_mean, f_median, na.rm = TRUE)
    plot(t_grid, f_true, type = "l", lwd = 2, ylim = c(0, ymax * 1.05),
         xlab = "t", ylab = "Density",
         main = paste0(main_prefix,
                       if (nzchar(main_prefix)) " " else "",
                       "True vs Posterior Mean/Median PDF"))
    lines(t_grid, f_mean,   lwd = 2, lty = 2)
    lines(t_grid, f_median, lwd = 2, lty = 3)
    legend("topright",
           c("f_true", "f_mean", "f_median"),
           lwd = 2, lty = c(1, 2, 3), bty = "n")
  }
  
  list(
    KS_mean   = KS_mean,
    KS_median = KS_median,
    KS_avg    = KS_avg,
    KS_sd     = KS_sd,
    
    ISE_mean   = ISE_mean,
    ISE_median = ISE_median,
    L2_mean    = L2_mean,
    L2_median  = L2_median,
    
    ISE_avg = ISE_avg,
    ISE_sd  = ISE_sd,
    L2_avg  = L2_avg,
    L2_sd   = L2_sd
  )
}


## Run every sampler on n_reps replicate data sets from one design cell.
##
## @param k,N,clust_sep,sig_const  passed straight to sim_data().
## @param prior_M  code for the number of mixture components M carried by the
##                 fitted model: 1 gives M = k + 2, 2 gives M = 25 (the setting
##                 used in the manuscript), anything else gives M = k.
## @return list(summary, K_draws), where summary is an n_reps x length(coln)
##         matrix of metrics and K_draws holds the per-draw number of occupied
##         components for each sampler.
ss_study <- function(k, N, clust_sep, sig_const, prior_M){
  ndata <- n_reps
  methods <- c("mcmc","ds","ds-mcmc","ds10-obs","ds25-obs","ds50-obs",
               "ds10-obsmap","ds25-obsmap","ds50-obsmap","ds-gibbs","ds-ml")
  metrics <- c("mnARI","sdARI","speed",
               "essMean","essM2",
               "postMean","postMed",
               "mnKS","sdKS","KSmed","KSmean",   
               "mnISE","sdISE","ISEmed","ISEmean",
               "mnL2","sdL2","L2med","L2mean")
  coln <- c(outer(metrics, methods, paste, sep = "_"))
  coln <- c(coln, "N","k","clust_sep","sig_const","M","rep")
  out <- matrix(NA_real_, nrow = ndata, ncol = length(coln), dimnames = list(NULL, coln))

  ## Posterior sample size retained per sampler.  Every method below is tuned
  ## to return exactly this many draws so that ESS and the density metrics are
  ## computed from equally sized samples.
  nsamp <- 1000

  ## Number of occupied components K per draw, kept for the posterior-of-K
  ## summaries.  The exhaustive direct sampler is only feasible for N <= 13.
  K_store <- list(
    mcmc        = matrix(NA_integer_, ndata, nsamp),
    ds          = if (N <= 13) matrix(NA_integer_, ndata, nsamp) else NULL,
    ds_mcmc     = matrix(NA_integer_, ndata, nsamp),
    ds_obs10    = matrix(NA_integer_, ndata, nsamp),
    ds_obs25    = matrix(NA_integer_, ndata, nsamp),
    ds_obs50    = matrix(NA_integer_, ndata, nsamp), 
    ds_obsmap10  = matrix(NA_integer_, ndata, nsamp), 
    ds_obsmap25  = matrix(NA_integer_, ndata, nsamp), 
    ds_obsmap50  = matrix(NA_integer_, ndata, nsamp),
    ds_gibbs    = matrix(NA_integer_, ndata, nsamp), 
    ds_ml       = matrix(NA_integer_, ndata, nsamp)
  )
  for(d in 1:ndata){
    ## A failure in any one sampler costs only this replicate: the row is
    ## blanked, its design columns are restored, and the loop continues.
    tryCatch({
      ## Seed depends on the cell and the replicate, so every method sees the
      ## same data and reruns of a single cell are reproducible.
      set.seed(100000 + 10000*N + 100*k + d)
      cat("N =", N, " d =", d, "clust_sep = ", clust_sep, "sig_const = ", sig_const, "prior_M = ", prior_M, "k = ", k, "\n")
      df.out <- sim_data(k=k, N=N, clust_sep = clust_sep, sig_const = sig_const)
      df <- df.out[[1]]
      true_pi <- df.out[[2]]
      true_mu <- df.out[[3]]
      true_sig <- df.out[[4]]
      si <- df$si
      yi <- df$yi
      
      ## Common evaluation grid and targets for the density metrics: the
      ## empirical CDF of this replicate and the true mixture density.
      t_grid <- seq(min(yi), max(yi), length.out = 1000)
      F_emp  <- ecdf(yi)(t_grid)
      
      f_true <- rowSums(matrix(true_pi, nrow = length(t_grid), ncol = length(true_pi), byrow = TRUE) *
                          dnorm(matrix(t_grid, nrow = length(t_grid), ncol = length(true_pi)),
                                matrix(true_mu, nrow = length(t_grid), ncol = length(true_pi), byrow = TRUE),
                                matrix(true_sig, nrow = length(t_grid), ncol = length(true_pi), byrow = TRUE)))
      
      ## Prior hyperparameters of Equation (1): symmetric Dirichlet weight
      ## alpha, prior variance sig2_mu on the component means, and the
      ## inverse-gamma shape and scale (omega, kappa) on the variances.
      alpha <- 1
      sig2_mu <- 100
      omega   <- 2
      kappa   <- 1
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

      ## MCMC: nthin = (niter - nburn) / nsamp = (5000 - 1000) / 4 = 1000 draws.
      time1 <- Sys.time()
      fit <- full_mcmc_sampling_dependent(y=yi, M=M, sig2_mu = sig2_mu, omega=omega, 
                                          kappa=kappa, alpha=alpha, niter=5000, nburn=1000, nthin=4)
      time2 <- Sys.time()
      out[d,"speed_mcmc"] <- as.numeric(time2 - time1, units = "secs") 
      
      ## Direct Sampler: enumerates every set partition of the N observations,
      ## so it is only tractable for small N.  M is capped at N because a
      ## partition cannot occupy more than N components.
      if(N <= 13){
        Mcap <- min(M, N)
        prior_for_Mcap <- rep(1,Mcap)/Mcap
        time1 <- Sys.time()
        ds <- direct_sampling(y=yi, sig2_mu=sig2_mu, omega=omega, kappa=kappa, M=Mcap,
                              nsamp=nsamp, prior_for_M = prior_for_Mcap)
        time2 <- Sys.time()
        out[d,"speed_ds"] <- as.numeric(time2 - time1, units = "secs") 
      }
      
      ## Blocked Gibbs reference sampler (not reported in the manuscript).
      time1 <- Sys.time()
      ds_MCMC <- direct_sampling_MCMC(y=yi, sig2_mu=sig2_mu, omega=omega, kappa=kappa, M=M,
                                      niter=1500, nburn=500, nthin=1)
      time2 <- Sys.time()
      out[d,"speed_ds-mcmc"] <- as.numeric(time2 - time1, units = "secs")
      
      
      ## DS-Const with a subset of 10% of the observations.
      time1 <- Sys.time()
      ds_postobs10 <- direct_sampling_postpred(y=yi, sig2_mu=sig2_mu, omega=omega, kappa=kappa, M=M,
                                               nsamp=nsamp, n_sub=0.1, d=d, z_true = si, 
                                               plot_candidates = FALSE, max_plot_candidates = 5, 
                                               plot_pZ_K = FALSE, plot_diagn = FALSE, add_truth = FALSE)
      time2 <- Sys.time()
      out[d,"speed_ds10-obs"] <- as.numeric(time2 - time1, units = "secs")
      
      ## DS-Const, 25% subset.
      time1 <- Sys.time()
      ds_postobs25 <- direct_sampling_postpred(y=yi, sig2_mu=sig2_mu, omega=omega, kappa=kappa, M=M,
                                               nsamp=nsamp, n_sub=0.25, d=d, z_true = si, 
                                               plot_candidates = FALSE, max_plot_candidates = 5, 
                                               plot_pZ_K = FALSE, plot_diagn = FALSE, add_truth = FALSE)
      time2 <- Sys.time()
      out[d,"speed_ds25-obs"] <- as.numeric(time2 - time1, units = "secs")
      
      ## DS-Const, 50% subset.
      time1 <- Sys.time()
      ds_postobs50 <- direct_sampling_postpred(y=yi, sig2_mu=sig2_mu, omega=omega, kappa=kappa, M=M,
                                               nsamp=nsamp, n_sub=0.5, d=d, z_true = si, 
                                               plot_candidates = FALSE, max_plot_candidates = 5, 
                                               plot_pZ_K = FALSE, plot_diagn = FALSE, add_truth = FALSE)
      time2 <- Sys.time()
      out[d,"speed_ds50-obs"] <- as.numeric(time2 - time1, units = "secs")
      
      ## DS-ML: the machine-learning clustering algorithms are fit to the FULL
      ## data set.  
      time1 <- Sys.time()
      ds_ml <- direct_sampling_postpred_ML(y=yi, sig2_mu=sig2_mu, omega=omega, kappa=kappa, M=M,
                                           nsamp=nsamp, d=d, z_true = si, plot_candidates = FALSE,
                                           max_plot_candidates = 5, plot_pZ_K = FALSE, plot_diagn = FALSE,
                                           add_truth = FALSE)
      time2 <- Sys.time()
      out[d,"speed_ds-ml"] <- as.numeric(time2 - time1, units = "secs")
      
      
      ## DS-Const-MAP, 10% subset: as DS-Const but the candidate set is grown
      ## from a MAP initialisation of the subset labels.
      time1 <- Sys.time()
      ds_postobsmap10 <- direct_sampling_obs_refine(y=yi, sig2_mu=sig2_mu, omega=omega, kappa=kappa, M=M,
                                            nsamp=nsamp, n_sub=0.1, d=d, plot_candidates = FALSE, max_plot_candidates = 5, 
                                            plot_pZ_K = FALSE, plot_diagn = FALSE, z_true = si)
      time2 <- Sys.time()
      out[d,"speed_ds10-obsmap"] <- as.numeric(time2 - time1, units = "secs")
      
      ## DS-Const-MAP, 25% subset.
      time1 <- Sys.time()
      ds_postobsmap25 <- direct_sampling_obs_refine(y=yi, sig2_mu=sig2_mu, omega=omega, kappa=kappa, M=M,
                                            nsamp=nsamp, n_sub=0.25, d=d, plot_candidates = FALSE, max_plot_candidates = 5, 
                                            plot_pZ_K = FALSE, plot_diagn = FALSE, z_true = si)
      time2 <- Sys.time()
      out[d,"speed_ds25-obsmap"] <- as.numeric(time2 - time1, units = "secs")
      
      ## DS-Const-MAP, 50% subset.
      time1 <- Sys.time()
      ds_postobsmap50 <- direct_sampling_obs_refine(y=yi, sig2_mu=sig2_mu, omega=omega, kappa=kappa, M=M,
                                            nsamp=nsamp, n_sub=0.5, d=d, plot_candidates = FALSE, max_plot_candidates = 5, 
                                            plot_pZ_K = FALSE, plot_diagn = FALSE, z_true = si)
      time2 <- Sys.time()
      out[d,"speed_ds50-obsmap"] <- as.numeric(time2 - time1, units = "secs")
      
      ## MC-MCMC: collapsed label Gibbs sampler of Equation (8), followed by
      ## composition draws of the parameters.
      time1 <- Sys.time()
      ds_gibbs <- direct_sampling_gibbs(y=yi, sig2_mu=sig2_mu, omega=omega, kappa=kappa, M=M,
                                        nsamp=nsamp)
      time2 <- Sys.time()
      out[d,"speed_ds-gibbs"] <- as.numeric(time2 - time1, units = "secs")
      
      
      ## ------------------------------- summarize this replicate ----------
      K_store$mcmc[d, ]     <- as.integer(fit$kpost)         
      if (N <= 13) K_store$ds[d, ] <- as.integer(ds$K)
      K_store$ds_mcmc[d, ]  <- as.integer(ds_MCMC$K)
      K_store$ds_gibbs[d, ]  <- as.integer(ds_gibbs$K)
      K_store$ds_obs10[d, ] <- as.integer(ds_postobs10$K)
      K_store$ds_obs25[d, ] <- as.integer(ds_postobs25$K)
      K_store$ds_obs50[d, ] <- as.integer(ds_postobs50$K)
      K_store$ds_ml[d, ] <- as.integer(ds_ml$K)
      K_store$ds_obsmap10[d, ] <- as.integer(ds_postobsmap10$K)
      K_store$ds_obsmap25[d, ] <- as.integer(ds_postobsmap25$K)
      K_store$ds_obsmap50[d, ] <- as.integer(ds_postobsmap50$K)
      
      ## Adjusted Rand index of each posterior label draw against the truth.
      ARI.fit <- sapply(1:nsamp, function(x) ARI(fit$z[x,], si))
      if(N <= 13) ARI.ds <- sapply(1:nsamp, function(x) ARI(ds$Z[x,], si))
      ARI.ds_MCMC <- sapply(1:nsamp, function(x) ARI(ds_MCMC$Z[x,], si))
      ARI.ds_gibbs <- sapply(1:nsamp, function(x) ARI(ds_gibbs$Z[x,], si))
      ARI.ds_postobs10 <- sapply(1:nsamp, function(x) ARI(ds_postobs10$Z[x,], si))
      ARI.ds_postobs25 <- sapply(1:nsamp, function(x) ARI(ds_postobs25$Z[x,], si))
      ARI.ds_postobs50 <- sapply(1:nsamp, function(x) ARI(ds_postobs50$Z[x,], si))
      ARI.ds_ml <- sapply(1:nsamp, function(x) ARI(ds_ml$Z[x,], si))
      ARI.ds_postobsmap10 <- sapply(1:nsamp, function(x) ARI(ds_postobsmap10$Z[x,], si))
      ARI.ds_postobsmap25 <- sapply(1:nsamp, function(x) ARI(ds_postobsmap25$Z[x,], si))
      ARI.ds_postobsmap50 <- sapply(1:nsamp, function(x) ARI(ds_postobsmap50$Z[x,], si))
      
      ## compute mean ARI and sd ARI
      out[d, "mnARI_mcmc"] <- mean(ARI.fit); out[d, "sdARI_mcmc"] <- sd(ARI.fit)
      if (N <= 13) { out[d, "mnARI_ds"] <- mean(ARI.ds); out[d, "sdARI_ds"] <- sd(ARI.ds) }
      out[d, "mnARI_ds-mcmc"]  <- mean(ARI.ds_MCMC);    out[d, "sdARI_ds-mcmc"]  <- sd(ARI.ds_MCMC)
      out[d, "mnARI_ds-gibbs"]  <- mean(ARI.ds_gibbs);    out[d, "sdARI_ds-gibbs"]  <- sd(ARI.ds_gibbs)
      out[d, "mnARI_ds10-obs"]   <- mean(ARI.ds_postobs10); out[d, "sdARI_ds10-obs"]   <- sd(ARI.ds_postobs10)
      out[d, "mnARI_ds25-obs"]   <- mean(ARI.ds_postobs25); out[d, "sdARI_ds25-obs"]   <- sd(ARI.ds_postobs25)
      out[d, "mnARI_ds50-obs"]   <- mean(ARI.ds_postobs50); out[d, "sdARI_ds50-obs"]   <- sd(ARI.ds_postobs50)
      out[d, "mnARI_ds-ml"]   <- mean(ARI.ds_ml); out[d, "sdARI_ds-ml"]   <- sd(ARI.ds_ml)
      out[d, "mnARI_ds10-obsmap"]   <- mean(ARI.ds_postobsmap10); out[d, "sdARI_ds10-obsmap"]   <- sd(ARI.ds_postobsmap10)
      out[d, "mnARI_ds25-obsmap"]   <- mean(ARI.ds_postobsmap25); out[d, "sdARI_ds25-obsmap"]   <- sd(ARI.ds_postobsmap25)
      out[d, "mnARI_ds50-obsmap"]   <- mean(ARI.ds_postobsmap50); out[d, "sdARI_ds50-obsmap"]   <- sd(ARI.ds_postobsmap50)
      
      ## Effective sample size of the mixture mean sum_k pi_k mu_k, using the
      ## spectral estimator in coda::effectiveSize().
      out[d, "essMean_mcmc"] <- effectiveSize(sapply(1:nrow(fit$mu), function(x) sum(fit$mu[x,]*fit$pi[x,])))
      if (N <= 13) out[d, "essMean_ds"] <- effectiveSize(sapply(1:nsamp, function(x) sum(ds$mu[x,]*ds$pi[x,])))
      out[d, "essMean_ds-mcmc"] <- effectiveSize(sapply(1:nrow(ds_MCMC$mu), function(x) sum(ds_MCMC$mu[x,]*ds_MCMC$pi[x,])))
      out[d, "essMean_ds-gibbs"] <- effectiveSize(sapply(1:nsamp, function(x) sum(ds_gibbs$mu[x,]*ds_gibbs$pi[x,])))
      out[d, "essMean_ds10-obs"] <- effectiveSize(sapply(1:nsamp, function(x) sum(ds_postobs10$mu[x,]*ds_postobs10$pi[x,])))
      out[d, "essMean_ds25-obs"] <- effectiveSize(sapply(1:nsamp, function(x) sum(ds_postobs25$mu[x,]*ds_postobs25$pi[x,])))
      out[d, "essMean_ds50-obs"] <- effectiveSize(sapply(1:nsamp, function(x) sum(ds_postobs50$mu[x,]*ds_postobs50$pi[x,])))
      out[d, "essMean_ds-ml"] <- effectiveSize(sapply(1:nsamp, function(x) sum(ds_ml$mu[x,]*ds_ml$pi[x,])))
      out[d, "essMean_ds10-obsmap"] <- effectiveSize(sapply(1:nsamp, function(x) sum(ds_postobsmap10$mu[x,]*ds_postobsmap10$pi[x,])))
      out[d, "essMean_ds25-obsmap"] <- effectiveSize(sapply(1:nsamp, function(x) sum(ds_postobsmap25$mu[x,]*ds_postobsmap25$pi[x,])))
      out[d, "essMean_ds50-obsmap"] <- effectiveSize(sapply(1:nsamp, function(x) sum(ds_postobsmap50$mu[x,]*ds_postobsmap50$pi[x,])))
      
      ## Posterior mean of the mixture mean.
      out[d, "postMean_mcmc"] <- mean(sapply(1:nrow(fit$mu), function(x) sum(fit$mu[x,]*fit$pi[x,])))
      if (N <= 13) out[d, "postMean_ds"] <- mean(sapply(1:nsamp, function(x) sum(ds$mu[x,]*ds$pi[x,])))
      out[d, "postMean_ds-mcmc"] <- mean(sapply(1:nrow(ds_MCMC$mu), function(x) sum(ds_MCMC$mu[x,]*ds_MCMC$pi[x,])))
      out[d, "postMean_ds-gibbs"] <- mean(sapply(1:nsamp, function(x) sum(ds_gibbs$mu[x,]*ds_gibbs$pi[x,])))
      out[d, "postMean_ds10-obs"] <- mean(sapply(1:nsamp, function(x) sum(ds_postobs10$mu[x,]*ds_postobs10$pi[x,])))
      out[d, "postMean_ds25-obs"] <- mean(sapply(1:nsamp, function(x) sum(ds_postobs25$mu[x,]*ds_postobs25$pi[x,])))
      out[d, "postMean_ds50-obs"] <- mean(sapply(1:nsamp, function(x) sum(ds_postobs50$mu[x,]*ds_postobs50$pi[x,])))
      out[d, "postMean_ds-ml"] <- mean(sapply(1:nsamp, function(x) sum(ds_ml$mu[x,]*ds_ml$pi[x,])))
      out[d, "postMean_ds10-obsmap"] <- mean(sapply(1:nsamp, function(x) sum(ds_postobsmap10$mu[x,]*ds_postobsmap10$pi[x,])))
      out[d, "postMean_ds25-obsmap"] <- mean(sapply(1:nsamp, function(x) sum(ds_postobsmap25$mu[x,]*ds_postobsmap25$pi[x,])))
      out[d, "postMean_ds50-obsmap"] <- mean(sapply(1:nsamp, function(x) sum(ds_postobsmap50$mu[x,]*ds_postobsmap50$pi[x,])))

      ## Posterior median of the mixture mean.
      out[d, "postMed_mcmc"] <- median(sapply(1:nrow(fit$mu), function(x) sum(fit$mu[x,]*fit$pi[x,])))
      if (N <= 13) out[d, "postMed_ds"] <- median(sapply(1:nsamp, function(x) sum(ds$mu[x,]*ds$pi[x,])))
      out[d, "postMed_ds-mcmc"] <- median(sapply(1:nrow(ds_MCMC$mu), function(x) sum(ds_MCMC$mu[x,]*ds_MCMC$pi[x,])))
      out[d, "postMed_ds-gibbs"] <- median(sapply(1:nsamp, function(x) sum(ds_gibbs$mu[x,]*ds_gibbs$pi[x,])))
      out[d, "postMed_ds10-obs"] <- median(sapply(1:nsamp, function(x) sum(ds_postobs10$mu[x,]*ds_postobs10$pi[x,])))
      out[d, "postMed_ds25-obs"] <- median(sapply(1:nsamp, function(x) sum(ds_postobs25$mu[x,]*ds_postobs25$pi[x,])))
      out[d, "postMed_ds50-obs"] <- median(sapply(1:nsamp, function(x) sum(ds_postobs50$mu[x,]*ds_postobs50$pi[x,])))
      out[d, "postMed_ds-ml"] <- median(sapply(1:nsamp, function(x) sum(ds_ml$mu[x,]*ds_ml$pi[x,])))
      out[d, "postMed_ds10-obsmap"] <- median(sapply(1:nsamp, function(x) sum(ds_postobsmap10$mu[x,]*ds_postobsmap10$pi[x,])))
      out[d, "postMed_ds25-obsmap"] <- median(sapply(1:nsamp, function(x) sum(ds_postobsmap25$mu[x,]*ds_postobsmap25$pi[x,])))
      out[d, "postMed_ds50-obsmap"] <- median(sapply(1:nsamp, function(x) sum(ds_postobsmap50$mu[x,]*ds_postobsmap50$pi[x,])))

      ## Effective sample size of the mixture second moment,
      ## sum_k pi_k (sigma^2_k + mu_k^2).
      out[d, "essM2_mcmc"] <- effectiveSize(sapply(1:nrow(fit$mu), function(x) {
        sum(fit$pi[x,] * (fit$sigma2[x,] + fit$mu[x,]^2))
      }))
      if (N <= 13) out[d, "essM2_ds"] <- effectiveSize(sapply(1:nsamp, function(x) {
        sum(ds$pi[x,] * (ds$sigma2[x,] + ds$mu[x,]^2))
      }))
      out[d, "essM2_ds-mcmc"] <- effectiveSize(sapply(1:nrow(ds_MCMC$mu), function(x) {
        sum(ds_MCMC$pi[x,] * (ds_MCMC$sigma2[x,] + ds_MCMC$mu[x,]^2))
      }))
      out[d, "essM2_ds-gibbs"] <- effectiveSize(sapply(1:nrow(ds_gibbs$mu), function(x) {
        sum(ds_gibbs$pi[x,] * (ds_gibbs$sigma2[x,] + ds_gibbs$mu[x,]^2))
      }))
      out[d, "essM2_ds10-obs"] <- effectiveSize(sapply(1:nsamp, function(x) {
        sum(ds_postobs10$pi[x,] * (ds_postobs10$sigma2[x,] + ds_postobs10$mu[x,]^2))
      }))
      out[d, "essM2_ds25-obs"] <- effectiveSize(sapply(1:nsamp, function(x) {
        sum(ds_postobs25$pi[x,] * (ds_postobs25$sigma2[x,] + ds_postobs25$mu[x,]^2))
      }))
      out[d, "essM2_ds50-obs"] <- effectiveSize(sapply(1:nsamp, function(x) {
        sum(ds_postobs50$pi[x,] * (ds_postobs50$sigma2[x,] + ds_postobs50$mu[x,]^2))
      }))
      out[d, "essM2_ds-ml"] <- effectiveSize(sapply(1:nsamp, function(x) {
        sum(ds_ml$pi[x,] * (ds_ml$sigma2[x,] + ds_ml$mu[x,]^2))
      }))
      out[d, "essM2_ds10-obsmap"] <- effectiveSize(sapply(1:nsamp, function(x) {
        sum(ds_postobsmap10$pi[x,] * (ds_postobsmap10$sigma2[x,] + ds_postobsmap10$mu[x,]^2))
      }))
      out[d, "essM2_ds25-obsmap"] <- effectiveSize(sapply(1:nsamp, function(x) {
        sum(ds_postobsmap25$pi[x,] * (ds_postobsmap25$sigma2[x,] + ds_postobsmap25$mu[x,]^2))
      }))
      out[d, "essM2_ds50-obsmap"] <- effectiveSize(sapply(1:nsamp, function(x) {
        sum(ds_postobsmap50$pi[x,] * (ds_postobsmap50$sigma2[x,] + ds_postobsmap50$mu[x,]^2))
      }))
      
      ## Density-estimation metrics, one call per sampler.
      density.mcmc <- compute_ks_ise(fit$pi, fit$mu, fit$sigma2, t_grid, F_emp, f_true)
      density.dsmcmc <- compute_ks_ise(ds_MCMC$pi, ds_MCMC$mu, ds_MCMC$sigma2, t_grid, F_emp, f_true)
      density.dsgibbs <- compute_ks_ise(ds_gibbs$pi, ds_gibbs$mu, ds_gibbs$sigma2, t_grid, F_emp, f_true)
      density.dsobsmap10 <- compute_ks_ise(ds_postobsmap10$pi, ds_postobsmap10$mu, ds_postobsmap10$sigma2, t_grid, F_emp, f_true)
      density.dsobsmap25 <- compute_ks_ise(ds_postobsmap25$pi, ds_postobsmap25$mu, ds_postobsmap25$sigma2, t_grid, F_emp, f_true)
      density.dsobsmap50 <- compute_ks_ise(ds_postobsmap50$pi, ds_postobsmap50$mu, ds_postobsmap50$sigma2, t_grid, F_emp, f_true)
      density.dsml<- compute_ks_ise(ds_ml$pi, ds_ml$mu, ds_ml$sigma2, t_grid, F_emp, f_true)
      density.dsobs10 <- compute_ks_ise(ds_postobs10$pi, ds_postobs10$mu, ds_postobs10$sigma2, t_grid, F_emp, f_true)
      density.dsobs25 <- compute_ks_ise(ds_postobs25$pi, ds_postobs25$mu, ds_postobs25$sigma2, t_grid, F_emp, f_true)
      density.dsobs50 <- compute_ks_ise(ds_postobs50$pi, ds_postobs50$mu, ds_postobs50$sigma2, t_grid, F_emp, f_true)
      
      
      if(N <= 13){density.ds <- compute_ks_ise(ds$pi, ds$mu, ds$sigma2, t_grid, F_emp, f_true)}
      
      ## Unpack the KS/ISE/L2 summaries into the result row.
      out[d, "mnKS_mcmc"] <- density.mcmc$KS_avg; out[d, "sdKS_mcmc"] <- density.mcmc$KS_sd
      out[d, "KSmed_mcmc"] <- density.mcmc$KS_median; out[d, "KSmean_mcmc"] <- density.mcmc$KS_mean
      out[d, "mnISE_mcmc"] <- density.mcmc$ISE_avg; out[d, "sdISE_mcmc"] <- density.mcmc$ISE_sd
      out[d, "ISEmed_mcmc"] <- density.mcmc$ISE_median; out[d, "ISEmean_mcmc"] <- density.mcmc$ISE_mean
      out[d, "mnL2_mcmc"] <- density.mcmc$L2_avg; out[d, "sdL2_mcmc"] <- density.mcmc$L2_sd
      out[d, "L2med_mcmc"] <- density.mcmc$L2_median; out[d, "L2mean_mcmc"] <- density.mcmc$L2_mean
      
      if(N <= 13){
        out[d, "mnKS_ds"] <- density.ds$KS_avg; out[d, "sdKS_ds"] <- density.ds$KS_sd
        out[d, "KSmed_ds"] <- density.ds$KS_median; out[d, "KSmean_ds"] <- density.ds$KS_mean
        out[d, "mnISE_ds"] <- density.ds$ISE_avg; out[d, "sdISE_ds"] <- density.ds$ISE_sd
        out[d, "ISEmed_ds"] <- density.ds$ISE_median; out[d, "ISEmean_ds"] <- density.ds$ISE_mean
        out[d, "mnL2_ds"] <- density.ds$L2_avg; out[d, "sdL2_ds"] <- density.ds$L2_sd
        out[d, "L2med_ds"] <- density.ds$L2_median; out[d, "L2mean_ds"] <- density.ds$L2_mean
      }
      
      out[d, "mnKS_ds-mcmc"] <- density.dsmcmc$KS_avg; out[d, "sdKS_ds-mcmc"] <- density.dsmcmc$KS_sd
      out[d, "KSmed_ds-mcmc"] <- density.dsmcmc$KS_median; out[d, "KSmean_ds-mcmc"] <- density.dsmcmc$KS_mean
      out[d, "mnISE_ds-mcmc"] <- density.dsmcmc$ISE_avg; out[d, "sdISE_ds-mcmc"] <- density.dsmcmc$ISE_sd
      out[d, "ISEmed_ds-mcmc"] <- density.dsmcmc$ISE_median; out[d, "ISEmean_ds-mcmc"] <- density.dsmcmc$ISE_mean
      out[d, "mnL2_ds-mcmc"] <- density.dsmcmc$L2_avg; out[d, "sdL2_ds-mcmc"] <- density.dsmcmc$L2_sd
      out[d, "L2med_ds-mcmc"] <- density.dsmcmc$L2_median; out[d, "L2mean_ds-mcmc"] <- density.dsmcmc$L2_mean
      
      out[d, "mnKS_ds-gibbs"] <- density.dsgibbs$KS_avg; out[d, "sdKS_ds-gibbs"] <- density.dsgibbs$KS_sd
      out[d, "KSmed_ds-gibbs"] <- density.dsgibbs$KS_median; out[d, "KSmean_ds-gibbs"] <- density.dsgibbs$KS_mean
      out[d, "mnISE_ds-gibbs"] <- density.dsgibbs$ISE_avg; out[d, "sdISE_ds-gibbs"] <- density.dsgibbs$ISE_sd
      out[d, "ISEmed_ds-gibbs"] <- density.dsgibbs$ISE_median; out[d, "ISEmean_ds-gibbs"] <- density.dsgibbs$ISE_mean
      out[d, "mnL2_ds-gibbs"] <- density.dsgibbs$L2_avg; out[d, "sdL2_ds-gibbs"] <- density.dsgibbs$L2_sd
      out[d, "L2med_ds-gibbs"] <- density.dsgibbs$L2_median; out[d, "L2mean_ds-gibbs"] <- density.dsgibbs$L2_mean
      
      out[d, "mnKS_ds10-obs"] <- density.dsobs10$KS_avg; out[d, "sdKS_ds10-obs"] <- density.dsobs10$KS_sd
      out[d, "KSmed_ds10-obs"] <- density.dsobs10$KS_median; out[d, "KSmean_ds10-obs"] <- density.dsobs10$KS_mean
      out[d, "mnISE_ds10-obs"] <- density.dsobs10$ISE_avg; out[d, "sdISE_ds10-obs"] <- density.dsobs10$ISE_sd
      out[d, "ISEmed_ds10-obs"] <- density.dsobs10$ISE_median; out[d, "ISEmean_ds10-obs"] <- density.dsobs10$ISE_mean
      out[d, "mnL2_ds10-obs"] <- density.dsobs10$L2_avg; out[d, "sdL2_ds10-obs"] <- density.dsobs10$L2_sd
      out[d, "L2med_ds10-obs"] <- density.dsobs10$L2_median; out[d, "L2mean_ds10-obs"] <- density.dsobs10$L2_mean
      
      out[d, "mnKS_ds25-obs"] <- density.dsobs25$KS_avg; out[d, "sdKS_ds25-obs"] <- density.dsobs25$KS_sd
      out[d, "KSmed_ds25-obs"] <- density.dsobs25$KS_median; out[d, "KSmean_ds25-obs"] <- density.dsobs25$KS_mean
      out[d, "mnISE_ds25-obs"] <- density.dsobs25$ISE_avg; out[d, "sdISE_ds25-obs"] <- density.dsobs25$ISE_sd
      out[d, "ISEmed_ds25-obs"] <- density.dsobs25$ISE_median; out[d, "ISEmean_ds25-obs"] <- density.dsobs25$ISE_mean
      out[d, "mnL2_ds25-obs"] <- density.dsobs25$L2_avg; out[d, "sdL2_ds25-obs"] <- density.dsobs25$L2_sd
      out[d, "L2med_ds25-obs"] <- density.dsobs25$L2_median; out[d, "L2mean_ds25-obs"] <- density.dsobs25$L2_mean
      
      out[d, "mnKS_ds50-obs"] <- density.dsobs50$KS_avg; out[d, "sdKS_ds50-obs"] <- density.dsobs50$KS_sd
      out[d, "KSmed_ds50-obs"] <- density.dsobs50$KS_median; out[d, "KSmean_ds50-obs"] <- density.dsobs50$KS_mean
      out[d, "mnISE_ds50-obs"] <- density.dsobs50$ISE_avg; out[d, "sdISE_ds50-obs"] <- density.dsobs50$ISE_sd
      out[d, "ISEmed_ds50-obs"] <- density.dsobs50$ISE_median; out[d, "ISEmean_ds50-obs"] <- density.dsobs50$ISE_mean
      out[d, "mnL2_ds50-obs"] <- density.dsobs50$L2_avg; out[d, "sdL2_ds50-obs"] <- density.dsobs50$L2_sd
      out[d, "L2med_ds50-obs"] <- density.dsobs50$L2_median; out[d, "L2mean_ds50-obs"] <- density.dsobs50$L2_mean
      
      out[d, "mnKS_ds-ml"] <- density.dsml$KS_avg; out[d, "sdKS_ds-ml"] <- density.dsml$KS_sd
      out[d, "KSmed_ds-ml"] <- density.dsml$KS_median; out[d, "KSmean_ds-ml"] <- density.dsml$KS_mean
      out[d, "mnISE_ds-ml"] <- density.dsml$ISE_avg; out[d, "sdISE_ds-ml"] <- density.dsml$ISE_sd
      out[d, "ISEmed_ds-ml"] <- density.dsml$ISE_median; out[d, "ISEmean_ds-ml"] <- density.dsml$ISE_mean
      out[d, "mnL2_ds-ml"] <- density.dsml$L2_avg; out[d, "sdL2_ds-ml"] <- density.dsml$L2_sd
      out[d, "L2med_ds-ml"] <- density.dsml$L2_median; out[d, "L2mean_ds-ml"] <- density.dsml$L2_mean
      
      out[d, "mnKS_ds10-obsmap"] <- density.dsobsmap10$KS_avg; out[d, "sdKS_ds10-obsmap"] <- density.dsobsmap10$KS_sd
      out[d, "KSmed_ds10-obsmap"] <- density.dsobsmap10$KS_median; out[d, "KSmean_ds10-obsmap"] <- density.dsobsmap10$KS_mean
      out[d, "mnISE_ds10-obsmap"] <- density.dsobsmap10$ISE_avg; out[d, "sdISE_ds10-obsmap"] <- density.dsobsmap10$ISE_sd
      out[d, "ISEmed_ds10-obsmap"] <- density.dsobsmap10$ISE_median; out[d, "ISEmean_ds10-obsmap"] <- density.dsobsmap10$ISE_mean
      out[d, "mnL2_ds10-obsmap"] <- density.dsobsmap10$L2_avg; out[d, "sdL2_ds10-obsmap"] <- density.dsobsmap10$L2_sd
      out[d, "L2med_ds10-obsmap"] <- density.dsobsmap10$L2_median; out[d, "L2mean_ds10-obsmap"] <- density.dsobsmap10$L2_mean
      
      out[d, "mnKS_ds25-obsmap"] <- density.dsobsmap25$KS_avg; out[d, "sdKS_ds25-obsmap"] <- density.dsobsmap25$KS_sd
      out[d, "KSmed_ds25-obsmap"] <- density.dsobsmap25$KS_median; out[d, "KSmean_ds25-obsmap"] <- density.dsobsmap25$KS_mean
      out[d, "mnISE_ds25-obsmap"] <- density.dsobsmap25$ISE_avg; out[d, "sdISE_ds25-obsmap"] <- density.dsobsmap25$ISE_sd
      out[d, "ISEmed_ds25-obsmap"] <- density.dsobsmap25$ISE_median; out[d, "ISEmean_ds25-obsmap"] <- density.dsobsmap25$ISE_mean
      out[d, "mnL2_ds25-obsmap"] <- density.dsobsmap25$L2_avg; out[d, "sdL2_ds25-obsmap"] <- density.dsobsmap25$L2_sd
      out[d, "L2med_ds25-obsmap"] <- density.dsobsmap25$L2_median; out[d, "L2mean_ds25-obsmap"] <- density.dsobsmap25$L2_mean
      
      out[d, "mnKS_ds50-obsmap"] <- density.dsobsmap50$KS_avg; out[d, "sdKS_ds50-obsmap"] <- density.dsobsmap50$KS_sd
      out[d, "KSmed_ds50-obsmap"] <- density.dsobsmap50$KS_median; out[d, "KSmean_ds50-obsmap"] <- density.dsobsmap50$KS_mean
      out[d, "mnISE_ds50-obsmap"] <- density.dsobsmap50$ISE_avg; out[d, "sdISE_ds50-obsmap"] <- density.dsobsmap50$ISE_sd
      out[d, "ISEmed_ds50-obsmap"] <- density.dsobsmap50$ISE_median; out[d, "ISEmean_ds50-obsmap"] <- density.dsobsmap50$ISE_mean
      out[d, "mnL2_ds50-obsmap"] <- density.dsobsmap50$L2_avg; out[d, "sdL2_ds50-obsmap"] <- density.dsobsmap50$L2_sd
      out[d, "L2med_ds50-obsmap"] <- density.dsobsmap50$L2_median; out[d, "L2mean_ds50-obsmap"] <- density.dsobsmap50$L2_mean
      
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
## each worker its own stream, and the per-replicate seed set inside ss_study()
## keeps the generated data reproducible regardless of the scheduling.
set.seed(1998)
ss_out <- mclapply(param_list, function(p) {
  ss_study(k = p$k, N = p$N, clust_sep = p$clust_sep, sig_const = p$sig_const, prior_M = p$prior_M)
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
## One tab-delimited row per (cell, replicate).  The workspace image preserves
## K_draws and the design objects for any follow-up summaries.
##
## The filename is time-stamped to the second, which collides if several jobs
## finish together -- as array tasks do.  Set FMM_TAG (the SLURM scripts use the
## job and array-task id) to append a unique suffix.  The plotting script's glob
## is study1_saveK_results*.txt, so tagged files are still picked up.
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
tag   <- Sys.getenv("FMM_TAG", unset = "")
if (nzchar(tag)) stamp <- paste0(stamp, "_", gsub("[^A-Za-z0-9._-]", "-", tag))
outfile <- file.path(out_dir, sprintf("study1_saveK_results%s.txt", stamp))
write.table(res, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Wrote results to:", outfile, "\n")
save.image(file = file.path(out_dir, sprintf("MixSimSaveK_y_%s.RData", stamp)))