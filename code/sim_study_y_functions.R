## ===========================================================================
##  Posterior samplers for the intercept-only Gaussian finite mixture model
##  with a Normal-Inverse-Gamma prior.
##
##  Companion to: "Bayesian analysis using a constrained mixture of
##  normal-inverse-gamma models", Section 3.3 and Appendix E.1.
##  Driver: sim_study_y.R
##
##  ---------------------------------------------------------------------------
##  Model (manuscript Equation 1, with p = 1 and x_i = 1)
##
##      Y_i | Z_ik = 1        ~  N(mu_k, sigma2_k)
##      Z_i | pi              ~  Categorical(pi_1, ..., pi_K)
##      pi                    ~  Dirichlet(alpha_1, ..., alpha_K)
##      mu_k | sigma2_k       ~  N(0, sigma2_k * sigma2_beta)      (scale dependent)
##      sigma2_k              ~  InverseGamma(omega, kappa)
##
##  ---------------------------------------------------------------------------
##  Argument naming
##
##      sig2_mu   sigma2_beta, the prior scale on the component mean
##      omega     inverse-gamma shape
##      kappa     inverse-gamma scale
##      alpha     symmetric Dirichlet concentration
##      M         K, the number of available mixture components
##      nsamp     number of posterior draws returned
##
##  Simulation-study values (Section 3.1): alpha = 1, sigma2_beta = 100,
##  omega = 2, kappa = 1, M = 25 (M = 10 when N = 10).
##
##  ---------------------------------------------------------------------------
##  Contents
##
##  Shared building blocks
##      relabel_labels            canonical 1..k relabelling of an allocation
##      partition_log_terms       partition-only pieces of Equation (6)
##      log_partition_weight      Equation (6) with K fixed at M
##      candidate_probs           Equation (12) normalization
##      draw_params_given_Z       Equations (3)-(5)
##      compose_from_candidates   method of composition over a candidate set
##
##  Samplers (manuscript name in quotes)
##      full_mcmc_sampling_dependent   "MCMC", the unconstrained Gibbs baseline
##      direct_sampling                "Direct Sampler", full enumeration, N <= 13
##      gibbs_labels_only_intercept    label chain of Equation (8)
##      direct_sampling_gibbs          "MC-MCMC"
##      direct_sampling_postpred       "DS-Const"
##      direct_sampling_postpred_ML    "DS-ML"
##      direct_sampling_obs_refine     "DS-Const-MAP"
##      direct_sampling_MCMC           blocked Gibbs; not used in the manuscript
##
##  Candidate-set construction and diagnostics
##      get_ml_candidates         partitions from ML clustering algorithms
##      map_kmeans_candidates     k-means initialized at a MAP partition
##      summarize_candidates, plot_top_candidates, plot_pZ_vs_K,
##      plot_candidate_diagnostics
##
##  Requires: salso (enumerate.partitions, ARI), miscPack (gaussian_mixture);
##  optionally mclust, cluster, kernlab, dbscan, e1071, ggplot2 for the ML
##  candidate generators and diagnostics.  The Dirichlet and inverse-gamma
##  draws are defined below rather than taken from MCMCpack.
## ===========================================================================

## Draw n vectors from Dirichlet(alpha), returned as an n x length(alpha) matrix.
rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  x / as.vector(x %*% rep(1, l))
}

## Draw n values from InverseGamma(shape, scale); shape may be a vector.
rinvgamma <- function(n, shape, scale = 1) {
  1 / rgamma(n = n, shape = shape, rate = scale)
}


## Relabel an allocation canonically as 1, 2, ... in order of first appearance,
## so that allocations inducing the same partition compare equal.
relabel_labels <- function(z) {
  match(z, unique(z))
}

## ===========================================================================
##  Shared building blocks
## ===========================================================================

## Partition-dependent factors of Equation (6), for the intercept-only case
## (p = 1, X_k = 1_{n_k}).  Two identities make this O(1) per component:
##
##   |I_{n_k} + sigma2_beta 1 1'|          = 1 + sigma2_beta n_k
##   y_k'(sigma2_beta 1 1' + I)^{-1} y_k   = sum(y^2)
##                                           - sigma2_beta (sum y)^2
##                                             / (1 + sigma2_beta n_k)
##
## the second by Sherman-Morrison.  Unoccupied components contribute exactly
## zero, since at n_k = 0 both the numerator and denominator reduce to
## omega*log(kappa); this reproduces the product over {k : n_j > 0} in (6).
##
## None of the returned quantities depend on K, so callers that vary K may
## compute them once and reuse them (see direct_sampling).
##
## @param parts  candidates x N matrix of allocations with entries in 1..M
## @return list: n_j (block counts), k (occupied blocks), log_dir_occ
##   (sum of lgamma(n_j + alpha) over occupied blocks), log_nig (the
##   normal-inverse-gamma product)
partition_log_terms <- function(parts, y, M, alpha, sig2_mu, omega, kappa) {
  n_j   <- t(apply(parts, 1L, function(z) tabulate(as.integer(z), nbins = M)))
  y2sum <- t(apply(parts, 1L, function(z)
    vapply(1:M, function(k) sum(y[z == k]^2), numeric(1))))
  sum_y <- t(apply(parts, 1L, function(z)
    vapply(1:M, function(k) sum(y[z == k]), numeric(1))))
  sumy2 <- sum_y^2

  occ  <- n_j > 0
  lnum <- lden <- matrix(0, nrow = nrow(n_j), ncol = ncol(n_j))

  lnum[occ] <- n_j[occ] * log(1 / sqrt(2 * pi)) +
    (-0.5 * log1p(sig2_mu * n_j[occ])) +
    (omega * log(kappa) - lgamma(omega)) +
    lgamma(0.5 * n_j[occ] + omega)

  quad_term <- 0.5 * (y2sum[occ] -
                        sig2_mu / (sig2_mu * n_j[occ] + 1) * sumy2[occ]) + kappa
  lden[occ] <- (0.5 * n_j[occ] + omega) * log(quad_term)

  list(n_j         = n_j,
       k           = rowSums(occ),
       log_dir_occ = rowSums(lgamma(n_j + alpha) * occ),
       log_nig     = rowSums(lnum - lden))
}

## Log marginal label weight f({Z_ik}, K, y) of Equation (6), with K fixed at M.
## This is the case for every constrained sampler, where C(y) is an explicit
## list of allocations.  Factors that do not depend on the partition -- f(K),
## Gamma(sum alpha_k) and prod Gamma(alpha_k) -- are omitted; they cancel on
## normalization.
##
## Not applicable when K varies across candidates: those factors then become
## K-dependent.  See direct_sampling() for that case.
log_partition_weight <- function(parts, y, M, alpha, sig2_mu, omega, kappa) {
  tm  <- partition_log_terms(parts, y, M, alpha, sig2_mu, omega, kappa)
  n_j <- tm$n_j
  rowSums(lgamma(n_j + alpha)) - lgamma(rowSums(n_j + alpha)) + tm$log_nig
}

## Equation (12): sampling probabilities over a candidate set, normalized on
## the log scale for numerical stability.
candidate_probs <- function(log_w) {
  p <- exp(log_w - max(log_w))
  p / sum(p)
}

## Draw (pi, sigma2, mu) | Z, y from the conditionals of Equations (3)-(5),
## intercept-only with prior mean m = 0.  Components with n_k = 0 revert to
## their priors -- Dirichlet(alpha), IG(omega, kappa) and
## N(0, sigma2_k sigma2_beta).
##
## @return list with elements pi, sigma2, mu, each of length M
draw_params_given_Z <- function(Z, y, M, alpha, sig2_mu, omega, kappa, m = 0) {
  n_k    <- tabulate(as.integer(Z), nbins = M)
  y2_sum <- vapply(1:M, function(x) sum(y[Z == x]^2), numeric(1))
  sum_y  <- vapply(1:M, function(x) sum(y[Z == x]),   numeric(1))

  pii    <- as.numeric(rdirichlet(1, n_k + alpha))

  astar  <- 0.5 * n_k + omega
  bstar  <- 0.5 * (y2_sum - sig2_mu / (n_k * sig2_mu + 1) * sum_y^2) + kappa
  sigma2 <- rinvgamma(M, shape = astar, scale = bstar)

  mns <- (sum_y + m / sig2_mu) / (n_k + 1 / sig2_mu)
  vs  <- sigma2 / (n_k + 1 / sig2_mu)
  mu  <- rnorm(M, mns, sqrt(vs))

  list(pi = pii, sigma2 = sigma2, mu = mu)
}

## Method of composition: draw nsamp independent (Z, pi, sigma2, mu) tuples by
## sampling Z from the candidate set with probabilities pZ, then drawing the
## continuous parameters from their conditionals given Z.
##
## @return list: mu, sigma2, pi (nsamp x M), Z (nsamp x N), K (occupied
##   clusters per draw), fitted.values
compose_from_candidates <- function(unique_parts, pZ, y, M, alpha,
                                    sig2_mu, omega, kappa, nsamp, m = 0) {
  N <- length(y)
  mus <- sigma2s <- pis <- matrix(NA_real_, nrow = nsamp, ncol = M)
  Zs  <- matrix(NA_integer_, nrow = nsamp, ncol = N)
  Ys  <- matrix(NA_real_,    nrow = nsamp, ncol = N)
  K   <- numeric(nsamp)

  n_cand <- nrow(unique_parts)
  for (i in seq_len(nsamp)) {
    Z <- unique_parts[sample.int(n_cand, 1L, prob = pZ), ]
    K[i] <- length(unique(Z))

    draws <- draw_params_given_Z(Z, y, M, alpha, sig2_mu, omega, kappa, m)

    mus[i, ]     <- draws$mu
    sigma2s[i, ] <- draws$sigma2
    pis[i, ]     <- draws$pi
    Zs[i, ]      <- Z
    Ys[i, ]      <- draws$mu[Z]
  }

  list(mu = mus, sigma2 = sigma2s, pi = pis, Z = Zs, K = K, fitted.values = Ys)
}

## ===========================================================================
##  "MCMC": the unconstrained Gibbs sampler used as the baseline
## ===========================================================================
## Data-augmented Gibbs sampler over (pi, sigma2, mu, z). 
## @return list: mu, sigma2, pi (ntotal x M), z (ntotal x N), kpost (occupied
##   clusters per retained draw), where ntotal = floor((niter - nburn)/nthin)
full_mcmc_sampling_dependent <- function(y, M=3, alpha =1, omega=2, kappa=1, sig2_mu = 100,
                                         niter=1000, nburn=0, nthin=1){
  N <- length(y)
  ntotal <- floor((niter - nburn) / nthin)
  # initialize parameters
  z <- sample(1:M, N, replace=TRUE)
  mu <- rep(0, M)
  sig2 <- rep(1, M)
  wi <- rep(1/M, M)

  # stored values
  z_keep <- matrix(NA, ntotal, N)
  mu_keep <- sig2_keep <- pi_keep <- matrix(NA, ntotal, M)
  k_keep <- matrix(NA, ntotal, 1)

  ii <- 1
  for(i in 1:niter){
    # update w
    nk <- sapply(1:M, function(jj) sum(z==jj))
    alpha_star <- alpha + nk
    w <- as.numeric(rdirichlet(1, alpha_star))

    sum_y   <- sapply(1:M, function(k) sum(y[z == k]))
    ssq <- sapply(1:M, function(x) sum((y[z==x]-mu[x])^2))

    # sample sigma2 -- IG(omega + (n_k + p)/2, kappa + SS/2 + mu^2/(2 sigma^2_beta))
    astar <-  omega + 0.5*(nk + 1)
    bstar <- kappa + 0.5*ssq + 0.5*mu^2/sig2_mu
    sig2 <- rinvgamma(M, shape=astar, scale=bstar)

    # sample mu
    mns <- (sum_y)/(nk + 1/sig2_mu)
    vs  <- (sig2)/(nk + 1/sig2_mu)
    mu <- rnorm(M, mns, sqrt(vs))

    # update zi
    un_norm <- sapply(1:N, function(x) dnorm(y[x], mu, sqrt(sig2), log=TRUE) + log(w))
    pZ <- exp(t(un_norm) - apply(un_norm, 2, max))/apply(exp(t(un_norm) - apply(un_norm, 2, max)),1,sum)
    z <- sapply(1:N, function(x) sample(1:M, 1, prob=pZ[x,]))

    # Count number of occupied components
    krep <- sum(sapply(1:M, function(jj) sum(z==jj)) >0)

    # save iterates
    if((i > nburn) & (((i - nburn) %% nthin) == 0)){
      z_keep[ii,] <- z
      mu_keep[ii,] <- mu
      sig2_keep[ii,] <- sig2
      pi_keep[ii,] <- w
      k_keep[ii] <- krep
      ii <- ii + 1
      if(ii > ntotal) break
    }

  }
  list(mu=mu_keep, sigma2 = sig2_keep, pi = pi_keep, z = z_keep, kpost=k_keep)
}

## ===========================================================================
##  "Direct Sampler": enumerate the label space (feasible for N <= 13 only)
## ===========================================================================
## Samples exactly from f({Z_ik}, K | y) by enumeration, with K random under the
## discrete uniform prior of Section 3.1.
##
## With K random, four factors of Equation (6) that are constant at fixed K
## become K-dependent -- f(K) = rho_K, Gamma(K*alpha), 1/Gamma(N + K*alpha) and
## 1/Gamma(alpha)^K -- and all are retained.  The log weight of partition P
## under K is
##
##   log w(P, K) = log rho_K
##               + lgamma(K + 1) - lgamma(K - k + 1)        # labellings
##               + lgamma(K*alpha) - lgamma(N + K*alpha)    # Dirichlet normalizer
##               + sum_{j: n_j>0} lgamma(n_j + alpha) - k*lgamma(alpha)
##               + [normal-inverse-gamma product, independent of K]
##
## Draws are padded to M columns: a draw with K components fills columns 1..K,
## and columns K+1..M carry pi = 0 with prior draws of sigma2 and mu, so they
## contribute nothing to the mixture density or to any moment summary.
##
## @param prior_for_M  length-M probability vector rho over K = 1..M; a vector
##   of ones instead fixes K at M
## @return list: mu, sigma2, pi (nsamp x M), Z (nsamp x N), K (occupied
##   clusters), Kcomp (sampled number of components), pKgY (analytic f(K | y))
direct_sampling <- function(y, M = 3, alpha = 1,
                            sig2_mu = 100, omega = 2, kappa = 1,
                            nsamp = 1000, prior_for_M) {
  # prior_for_M:
  #   - if a vector of ones, K is treated as fixed at M
  #   - otherwise, a length-M probability vector over K = 1..M (sums to 1)

  N <- length(y)

  ## ---- CAP M UP FRONT (and align prior) ----
  M_in <- M
  M <- min(M, N)

  if (length(prior_for_M) != M) {
    if (all(prior_for_M == 1)) {
      prior_for_M <- rep(1, M)
    } else if (length(prior_for_M) == M_in) {
      prior_for_M <- prior_for_M[seq_len(M)]
    } else {
      stop("prior_for_M must be either all ones (fixed-M) or a length-M vector that sums to 1.")
    }
  }
  fixed_K <- (sum(prior_for_M) == M)
  if (!fixed_K) {
    if (any(prior_for_M < 0) || abs(sum(prior_for_M) - 1) > 1e-8) {
      stop("In the variable-K case, prior_for_M must be a nonnegative vector of length M that sums to 1.")
    }
  }

  ## ---- Enumerate partitions ----
  all_parts <- enumerate.partitions(N)  # rows = partitions, columns = 1..N labels
  all_parts_M <- all_parts[apply(all_parts, 1, function(x) length(unique(x))) <= M, , drop = FALSE]

  if (min(all_parts_M) < 1L || max(all_parts_M) > M)
    stop("enumerate.partitions() returned labels outside 1..", M, " (range ",
         min(all_parts_M), "..", max(all_parts_M),
         "); direct_sampling() assumes 1-based contiguous cluster labels.")

  ## ---- Partition-only terms (independent of K, so computed once) ----
  tm    <- partition_log_terms(all_parts_M, y, M, alpha, sig2_mu, omega, kappa)
  kP    <- tm$k
  base  <- tm$log_dir_occ - kP * lgamma(alpha) + tm$log_nig

  ## ---- K-dependent terms ----
  K_grid <- if (fixed_K) M else seq_len(M)
  log_w  <- matrix(-Inf, nrow = nrow(all_parts_M), ncol = length(K_grid))

  for (kk in seq_along(K_grid)) {
    K  <- K_grid[kk]
    ok <- which(kP <= K)
    if (!length(ok)) next
    log_w[ok, kk] <-
      (if (fixed_K) 0 else log(prior_for_M[K])) +
      (lgamma(K + 1) - lgamma(K - kP[ok] + 1)) +
      (lgamma(K * alpha) - lgamma(N + K * alpha)) +
      base[ok]
  }

  w    <- exp(log_w - max(log_w[is.finite(log_w)]))
  w[!is.finite(log_w)] <- 0
  pKgY <- colSums(w) / sum(w)

  mus      <- matrix(NA_real_, nrow = nsamp, ncol = M)
  sigma2s  <- matrix(NA_real_, nrow = nsamp, ncol = M)
  pis      <- matrix(NA_real_, nrow = nsamp, ncol = M)
  Zs       <- matrix(NA_integer_, nrow = nsamp, ncol = N)
  Ys       <- matrix(NA_real_, nrow = nsamp, ncol = N)
  Ks       <- numeric(nsamp)   # number of OCCUPIED clusters, comparable across methods
  Kcomp    <- numeric(nsamp)   # number of COMPONENTS, the model parameter K

  for (i in seq_len(nsamp)) {
    kk <- sample.int(length(K_grid), 1L, prob = pKgY)
    K  <- K_grid[kk]

    # Sample partition Z | K
    wK   <- w[, kk]
    pick <- sample.int(length(wK), 1L, prob = wK)
    Z    <- as.integer(all_parts_M[pick, ])

    # Draw the K component parameters, then pad out to M columns with pi = 0
    draws <- draw_params_given_Z(Z, y, K, alpha, sig2_mu, omega, kappa)
    pad   <- M - K
    if (pad > 0L) {
      s2_pad <- rinvgamma(pad, shape = omega, scale = kappa)
      mu_all <- c(draws$mu,     rnorm(pad, 0, sqrt(s2_pad * sig2_mu)))
      s2_all <- c(draws$sigma2, s2_pad)
      pi_all <- c(draws$pi,     rep(0, pad))
    } else {
      mu_all <- draws$mu; s2_all <- draws$sigma2; pi_all <- draws$pi
    }

    mus[i, ]     <- mu_all
    sigma2s[i, ] <- s2_all
    pis[i, ]     <- pi_all
    Zs[i, ]      <- Z
    Ys[i, ]      <- mu_all[Z]
    Ks[i]        <- length(unique(Z))
    Kcomp[i]     <- K
  }

  list(mu = mus, sigma2 = sigma2s, pi = pis, Z = Zs, K = Ks,
       Kcomp = Kcomp, pKgY = setNames(pKgY, K_grid), fitted.values = Ys)
}

## ===========================================================================
##  Blocked Gibbs sampler (not used in the manuscript)
## ===========================================================================
## Updates labels from f(z_i | mu, sigma2, pi) rather than from the collapsed
## conditional of Equation (8), so this is NOT the manuscript's MC-MCMC.
## Retained for reference and reported in the driver as "ds-mcmc"; the plotting
## script excludes it.
direct_sampling_MCMC <- function(y, M = 3, alpha = 1, sig2_mu = 100, omega = 2, kappa = 1,
                                 niter = 1000, nburn = 0, nthin = 1){
  N <- length(y)
  ntotal <- floor((niter - nburn) / nthin)

  mus <- sigma2s <- pis <- matrix(NA, nrow=ntotal, ncol=M)
  Zs <- matrix(NA, nrow=ntotal, ncol=N)
  K_keep <- numeric(ntotal)

  mu <- rep(0, M); sigma2 <- rep(1, M); pii <- rep(1/M, M)

  ii <- 1
  for (iter in 1:niter) {
    un_norm <- sapply(1:N, function(x) dnorm(y[x], mu, sqrt(sigma2), log=TRUE) + log(pii))
    pZ <- exp(t(un_norm) - apply(un_norm, 2, max))/apply(exp(t(un_norm) - apply(un_norm, 2, max)),1,sum)
    Z <- sapply(1:N, function(x) sample(1:M, 1, prob=pZ[x,]))

    krep <- length(unique(Z))

    draws  <- draw_params_given_Z(Z, y, M, alpha, sig2_mu, omega, kappa)
    pii    <- draws$pi
    sigma2 <- draws$sigma2
    mu     <- draws$mu

    if ((iter > nburn) && (((iter - nburn) %% nthin) == 0)) {
      mus[ii, ]     <- mu
      sigma2s[ii, ] <- sigma2
      pis[ii, ]     <- pii
      Zs[ii, ]      <- Z
      K_keep[ii]    <- krep
      ii <- ii + 1
      if (ii > ntotal) break
    }
  }

  list(mu = mus, sigma2 = sigma2s, pi = pis, Z = Zs, K = K_keep)
}

## ===========================================================================
##  Candidate-set diagnostics
## ===========================================================================
## All of these run only when the corresponding plot_* flag is TRUE, so they
## stay out of the timed path of the samplers.

plot_candidate_diagnostics <- function(unique_parts, pZ, max_top = 30, n_sub) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; skipping plots.")
    return(invisible(NULL))
  }
  library(ggplot2)
  library(dplyr)

  stopifnot(nrow(unique_parts) == length(pZ))

  N <- ncol(unique_parts)
  title_info <- paste0("(N = ", N, ", n_sub = ", n_sub, ")")

  K_cand <- apply(unique_parts, 1, function(z) length(unique(z)))

  ord <- order(pZ, decreasing = TRUE)
  df <- data.frame(
    candidate = seq_along(pZ),
    rank      = seq_along(pZ),
    pZ        = pZ[ord],
    K         = K_cand[ord]
  )

  ## 1) Barplot of top candidates
  df_top <- head(df, max_top)

  p1 <- ggplot(df_top, aes(x = factor(rank), y = pZ, fill = factor(K))) +
    geom_col() +
    geom_text(aes(label = K), vjust = -0.3, size = 3) +
    labs(
      x = "Candidate rank (by pZ)",
      y = "pZ",
      fill = "K",
      title = paste("Top candidate partitions", title_info),
      subtitle = "Numbers above bars = K (no. of occupied clusters)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  print(p1)

  ## 2) Total probability mass by K
  df_K <- df |>
    group_by(K) |>
    summarize(total_pZ = sum(pZ), max_pZ = max(pZ), n_cand = dplyr::n(), .groups = "drop")

  p2 <- ggplot(df_K, aes(x = factor(K), y = total_pZ)) +
    geom_col() +
    geom_text(aes(label = paste0("n=", n_cand)), vjust = -0.3, size = 3) +
    labs(
      x = "K (number of occupied clusters)",
      y = "Total probability mass",
      title = paste("Total pZ by K", title_info),
      subtitle = "Labels show how many candidate partitions had that K"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  print(p2)

  ## 3) Cumulative probability vs rank
  df$cum_pZ <- cumsum(df$pZ)

  p3 <- ggplot(df, aes(x = rank, y = cum_pZ)) +
    geom_step() +
    geom_hline(yintercept = 0.9, linetype = "dashed") +
    labs(
      x = "Candidate rank (sorted by pZ)",
      y = "Cumulative probability",
      title = paste("Cumulative probability mass across candidates", title_info),
      subtitle = "Dashed line at 0.9; shows how many candidates carry 90% of mass"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  print(p3)

  invisible(list(df = df, df_K = df_K))
}

## Console summary of a candidate set: singleton counts, ARI against a known
## truth where supplied, and the leading candidates by posterior probability.
summarize_candidates <- function(unique_parts, pZ, unnorm_pZ, M, z_true = NULL,
                                 cand_names = NULL, top = 15) {
  print(summary(unnorm_pZ))
  print(summary(pZ))
  print(summary(log(pZ)))

  L <- nrow(unique_parts)

  n_singletons <- apply(unique_parts, 1, function(z_row) {
    sum(tabulate(as.integer(z_row), nbins = M) == 1)
  })

  cat("Unweighted summary of singleton counts over candidate partitions:\n")
  print(summary(n_singletons))
  cat("Posterior mean number of singleton clusters:", sum(pZ * n_singletons), "\n")
  cat("Posterior Pr( >= 1 singleton cluster ):", sum(pZ[n_singletons > 0]), "\n")

  ari_vec <- NULL
  if (!is.null(z_true)) {
    ari_vec <- apply(unique_parts, 1, function(z_row) ARI(z_row, z_true))
    cat("\nUnweighted summary of ARI over candidate partitions:\n")
    print(summary(ari_vec))
    cat("Posterior mean ARI:", sum(pZ * ari_vec), "\n")
  }

  k <- min(top, L)
  top_idx <- order(pZ, decreasing = TRUE)[1:k]

  cat("\nTop", k, "candidates by posterior probability:\n\n")
  for (j in seq_along(top_idx)) {
    idx <- top_idx[j]
    sizes_j <- tabulate(as.integer(unique_parts[idx, ]), nbins = M)

    cat("Candidate", idx)
    if (!is.null(cand_names) && nzchar(cand_names[idx])) cat(" (", cand_names[idx], ")", sep = "")
    cat("\n")
    cat("  pZ               :", pZ[idx], "\n")
    cat("  singleton count  :", sum(sizes_j == 1), "\n")
    cat("  cluster sizes    :", paste(sizes_j, collapse = " "), "\n")
    if (!is.null(ari_vec)) cat("  ARI w.r.t. truth :", ari_vec[idx], "\n")
    cat("\n")
  }

  invisible(ari_vec)
}

plot_top_candidates <- function(unique_parts, pZ, y, z_true = NULL,
                                max_plot_candidates = 5, n_sub = NA,
                                cand_names = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; skipping candidate plot.")
    return(invisible(NULL))
  }
  N <- length(y)
  title_info <- sprintf("(N = %d, n_sub = %s)", N, n_sub)

  ord      <- order(pZ, decreasing = TRUE)
  keep_ids <- ord[seq_len(min(length(pZ), max_plot_candidates))]

  df_list <- list()
  panel_levels <- character(0)
  idx <- 1L

  if (!is.null(z_true)) {
    if (length(z_true) != N) stop("z_true must have length N = ", N)
    panel_levels <- c(panel_levels, "Truth")
    df_list[[idx]] <- data.frame(panel = "Truth", obs_idx = seq_len(N),
                                 y = y, cluster = as.integer(z_true))
    idx <- idx + 1L
  }

  for (j in seq_along(keep_ids)) {
    g   <- keep_ids[j]
    nm  <- if (!is.null(cand_names) && nzchar(cand_names[g])) paste0(" [", cand_names[g], "]") else ""
    lab <- paste0("Cand ", g, nm, "\n p=", sprintf("%.3f", pZ[g]))
    panel_levels <- c(panel_levels, lab)
    df_list[[idx]] <- data.frame(panel = lab, obs_idx = seq_len(N),
                                 y = y, cluster = as.integer(unique_parts[g, ]))
    idx <- idx + 1L
  }

  df_long <- do.call(rbind, df_list)
  df_long$panel <- factor(df_long$panel, levels = panel_levels)

  print(
    ggplot2::ggplot(df_long, ggplot2::aes(x = obs_idx, y = y, colour = factor(cluster))) +
      ggplot2::geom_point(size = 1.8, alpha = 0.8) +
      ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
      ggplot2::labs(x = "Observation index", y = "y", colour = "Cluster",
                    title = paste("True labels and top candidate partitions", title_info)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
                     plot.title = ggplot2::element_text(hjust = 0.5))
  )
  invisible(NULL)
}

plot_pZ_vs_K <- function(unique_parts, pZ, unnorm_pZ, N, n_sub = NA) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; skipping pZ vs K plot.")
    return(invisible(NULL))
  }
  title_info <- sprintf("(N = %d, n_sub = %s)", N, n_sub)

  df_pZ <- data.frame(
    candidate = seq_along(pZ),
    K         = apply(unique_parts, 1, function(z) length(unique(z))),
    pZ        = pZ,
    logw      = unnorm_pZ,
    log_pZ    = log(pZ)
  )
  best_id <- which.max(df_pZ$pZ)

  one_panel <- function(yvar, ylab, title) {
    ggplot2::ggplot(df_pZ, ggplot2::aes(x = K, y = .data[[yvar]])) +
      ggplot2::geom_jitter(data = subset(df_pZ, candidate != best_id),
                           width = 0.15, height = 0, alpha = 0.5) +
      ggplot2::geom_point(data = df_pZ[best_id, , drop = FALSE],
                          ggplot2::aes(x = K, y = .data[[yvar]]),
                          inherit.aes = FALSE, color = "red", size = 4, shape = 19) +
      ggplot2::labs(x = "Number of occupied clusters (K)", y = ylab,
                    title = paste(title, title_info), subtitle = "Red point = max p(Z)") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))
  }

  print(one_panel("pZ",     "p(Z | y) for candidate partition",     "Candidate partition weights vs K"))
  print(one_panel("log_pZ", "log p(Z | y)",                          "Log posterior probabilities vs K"))
  print(one_panel("logw",   "unnormalized log weight  log w(Z | y)", "Unnormalized log weights vs K"))
  invisible(NULL)
}

## ===========================================================================
##  "DS-Const": candidate set from a Gibbs fit to a subset of the data
## ===========================================================================
## Appendix B.  Fit the unconstrained model by Gibbs sampling to n_sub * N
## observations; for each retained draw, predict allocations for all N
## observations; the distinct full-data allocations form C(y).  Labels are then
## sampled from C(y) via Equation (12) and the continuous parameters drawn from
## their conditionals.
##
## The subset chain must supply at least nsamp draws, which is checked.
##
## @param n_sub  subset fraction; the study uses 0.1, 0.25 and 0.5
## @param d      replicate index, used to seed the subset draw so that
##   DS-Const-MAP reuses the same subset
## @return the list from compose_from_candidates(), plus unique_parts and pZ so
##   that DS-Const-MAP can reuse the candidate set
direct_sampling_postpred <- function(y, M=3, alpha=1,
                                     sig2_mu=100, omega=2, kappa=1,
                                     nsamp=1000, n_sub=0.1, d,
                                     plot_candidates = FALSE, max_plot_candidates = 5,
                                     plot_pZ_K = FALSE, plot_diagn = FALSE, z_true = NULL,
                                     add_truth = FALSE,
                                     sub_niter = 11000, sub_nburn = 1000, sub_nthin = 10){

  N <- length(y)
  ys <- y
  m <- 0
  set.seed(d*11111)
  y_sub <- y[sample(1:N, floor(n_sub*N), replace=FALSE)]
  # The prediction step below indexes nsamp draws of the subset fit.
  if (floor((sub_niter - sub_nburn) / sub_nthin) < nsamp)
    stop("subset Gibbs yields ", floor((sub_niter - sub_nburn) / sub_nthin),
         " draws but nsamp = ", nsamp, ".  Increase sub_niter or reduce nsamp/sub_nthin.")

  fit <- gaussian_mixture(y=y_sub, N=M, m=0, v=sig2_mu, a=omega, b=kappa, alpha=alpha,
                          niter=sub_niter, nburn=sub_nburn, nthin=sub_nthin)

  partitions <- matrix(NA, nrow=nsamp, ncol=N)
  for(jj in 1:nsamp){
    un_norm <- sapply(1:N, function(x) dnorm(ys[x], fit$mu[jj,], sqrt(fit$sigma2[jj,]), log=TRUE) + log(fit$w[jj,]))
    pZ <- exp(t(un_norm) - apply(un_norm, 2, max))/apply(exp(t(un_norm) - apply(un_norm, 2, max)),1,sum)
    partitions[jj,] <- sapply(1:N, function(x) sample(1:M, 1, prob=pZ[x,]))
    partitions[jj,] <- relabel_labels(partitions[jj,])
  }

  unique_parts <- partitions[!duplicated(partitions), , drop = FALSE]
  if(add_truth == TRUE){
    unique_parts <- rbind(z_true, unique_parts)
  }

  unnorm_pZ <- log_partition_weight(unique_parts, y, M, alpha, sig2_mu, omega, kappa)
  pZ        <- candidate_probs(unnorm_pZ)

  ## Diagnostics only -- gated so they stay out of the timed path.
  if (plot_pZ_K || plot_diagn) {
    summarize_candidates(unique_parts, pZ, unnorm_pZ, M, z_true = z_true)
  }
  if (plot_candidates) {
    plot_top_candidates(unique_parts, pZ, y, z_true, max_plot_candidates, n_sub)
  }
  if (plot_pZ_K) {
    plot_pZ_vs_K(unique_parts, pZ, unnorm_pZ, N, n_sub)
  }
  if (plot_diagn) {
    plot_candidate_diagnostics(unique_parts, pZ, max_top = 30, n_sub)
  }

  out <- compose_from_candidates(unique_parts, pZ, y, M, alpha,
                                 sig2_mu, omega, kappa, nsamp, m)
  out$unique_parts <- unique_parts
  out$pZ           <- pZ
  out
}


## ---------------------------------------------------------------------------
## Candidate partitions from machine-learning clustering algorithms
## ---------------------------------------------------------------------------
## Each requested algorithm is fit over K_seq and the resulting partitions are
## deduplicated.  Algorithms are fit to a training subset (all observations by
## default) and every observation is then labelled by its nearest center or
## medoid, so the returned partitions always span all N observations.
##
## @param y            numeric vector, or N x p matrix of features
## @param M            maximum number of components
## @param algorithms   any of "kmeans", "gmm", "hclust_ward", "pam", "clara",
##                     "spectral", "dbscan", "fuzzyc".  "clara" is a large-N
##                     k-medoids and is the scalable stand-in for "pam"
## @param K_seq        values of K to fit; defaults to 2:min(M, N)
## @param use_mclust_autoK  if TRUE, "gmm" chooses G by BIC and contributes a
##                     single candidate; FALSE (the manuscript's specification)
##                     fits one model per K
## @param subset_frac  fraction of observations used to FIT the algorithms;
##                     labels are always produced for all N
## @param relabel      apply relabel_labels() row-wise
## @return matrix of distinct candidate partitions, one row per candidate, with
##   row names recording every algorithm that produced that partition
get_ml_candidates <- function(y, M, algorithms = c("kmeans","gmm","hclust_ward","pam"),
                              K_seq = NULL, use_mclust_autoK = FALSE, kmeans_nstart = 5,
                              seed = NULL, relabel = FALSE, subset_frac = 1, subset_index = NULL,
                              verbose = TRUE) {
  if (is.vector(y)) {
    y_mat <- matrix(y, ncol = 1L)
  } else {
    y_mat <- as.matrix(y)
  }

  N <- nrow(y_mat)
  all_idx <- seq_len(N)

  if (!is.null(seed)) set.seed(seed)

  ## ---- choose training subset ----
  if (!is.null(subset_index)) {
    subset_index <- sort(unique(subset_index))
    subset_index <- subset_index[subset_index %in% all_idx]
    if (length(subset_index) < 2L) {
      stop("subset_index must contain at least 2 valid indices.")
    }
    idx_train <- subset_index
  } else if (subset_frac < 1) {
    if (subset_frac <= 0 || subset_frac > 1) {
      stop("subset_frac must be in (0, 1].")
    }
    n_train <- max(2L, floor(subset_frac * N))
    idx_train <- sort(sample(all_idx, n_train, replace = FALSE))
  } else {
    idx_train <- all_idx
  }

  y_train  <- y_mat[idx_train, , drop = FALSE]
  N_train  <- nrow(y_train)

  M_eff <- min(M, N_train)

  if (is.null(K_seq)) {
    K_seq <- if (M_eff > 1L) 2:M_eff else 1L
  }

  if (N <= 13) {
    K_seq <- 2:(N - 1)
  }

  K_seq <- sort(unique(K_seq))
  K_seq <- K_seq[K_seq >= 1 & K_seq <= M_eff]
  if (length(K_seq) == 0L) {
    stop("K_seq has no valid values in [1, M_eff].")
  }

  cand_list <- list()

  assign_by_centers <- function(centers, X_full) {
    C <- as.matrix(centers)
    X <- as.matrix(X_full)
    X_sq <- rowSums(X^2)
    C_sq <- rowSums(C^2)
    d_sq <- outer(X_sq, C_sq, "+") - 2 * (X %*% t(C))
    max.col(-d_sq)  # index of smallest distance
  }

  centers_from_labels <- function(z_train, K) {
    centers <- matrix(NA_real_, nrow = K, ncol = ncol(y_mat))
    for (k in 1:K) {
      centers[k, ] <- colMeans(y_train[z_train == k, , drop = FALSE])
    }
    centers
  }

  for (alg in algorithms) {

    if (alg == "kmeans") {
      for (K in K_seq) {
        km_fit <- try(stats::kmeans(y_train, centers = K, nstart = kmeans_nstart), silent = TRUE)
        if (inherits(km_fit, "try-error")) {
          warning("kmeans failed for K = ", K, "; skipping this candidate.")
        } else {
          cand_list[[paste0("kmeans_K", K)]] <- as.integer(assign_by_centers(km_fit$centers, y_mat))
        }
      }

    } else if (alg == "gmm") {
      if (!requireNamespace("mclust", quietly = TRUE)) {
        warning("Skipping 'gmm': package 'mclust' not installed.")
      } else if (use_mclust_autoK) {
        # Let Mclust choose G in 1..M_eff via BIC -- ONE candidate only.
        fit <- try(mclust::Mclust(y_train, G = 1:M_eff, verbose = FALSE), silent = TRUE)
        if (inherits(fit, "try-error")) {
          warning("Mclust auto-K failed; skipping 'gmm' candidate.")
        } else {
          pred <- try(mclust::predict.Mclust(fit, newdata = y_mat), silent = TRUE)
          if (inherits(pred, "try-error")) {
            warning("predict.Mclust failed; skipping 'gmm' candidate.")
          } else {
            cand_list[[paste0("gmm_autoK_G", fit$G)]] <- as.integer(pred$classification)
          }
        }
      } else {
        for (K in K_seq) {
          fit <- try(mclust::Mclust(y_train, G = K, verbose = FALSE), silent = TRUE)
          if (inherits(fit, "try-error")) {
            warning("Mclust failed for G = ", K, "; skipping this candidate.")
          } else {
            pred <- try(mclust::predict.Mclust(fit, newdata = y_mat), silent = TRUE)
            if (inherits(pred, "try-error")) {
              warning("predict.Mclust failed for G = ", K, "; skipping.")
            } else {
              cand_list[[paste0("gmm_G", K)]] <- as.integer(pred$classification)
            }
          }
        }
      }

    } else if (alg == "hclust_ward") {
      d_train <- stats::dist(y_train)
      if (length(d_train) == 0L || max(d_train) == 0) {
        warning("Data (subset) are (almost) constant; skipping 'hclust_ward'.")
      } else {
        hc <- stats::hclust(d_train, method = "ward.D2")
        for (K in K_seq) {
          z_train <- as.integer(stats::cutree(hc, k = K))
          cand_list[[paste0("hclust_ward_K", K)]] <-
            as.integer(assign_by_centers(centers_from_labels(z_train, K), y_mat))
        }
      }

    } else if (alg == "spectral") {
      if (!requireNamespace("kernlab", quietly = TRUE)) {
        warning("Skipping 'spectral': package 'kernlab' not installed.")
      } else {

        if (nrow(y_train) > 5000L) {
          warning("Skipping 'spectral': ", nrow(y_train), " points would need a ",
                  nrow(y_train), " x ", nrow(y_train), " affinity matrix. ",
                  "This branch does not use the Nystrom approximation; use the ",
                  "regression code path, or subsample, if spectral is needed here.")
          next
        }
        d_train <- stats::dist(y_train)
        if (length(d_train) == 0L || max(d_train) == 0) {
          warning("Data (subset) are (almost) constant; skipping 'spectral'.")
        } else {
          sigma_val <- 1 / (2 * (stats::median(d_train)^2))
          for (K in K_seq) {
            fit <- try(kernlab::specc(x = y_train, centers = K, kernel = "rbfdot",
                                      kpar = list(sigma = sigma_val)), silent = TRUE)
            if (inherits(fit, "try-error")) {
              warning("kernlab::specc failed for centers = ", K, ": ",
                      conditionMessage(attr(fit, "condition")), " Skipping.")
            } else {
              cand_list[[paste0("spectral_K", K)]] <-
                as.integer(assign_by_centers(centers_from_labels(as.integer(fit), K), y_mat))
            }
          }
        }
      }

    } else if (alg == "pam") {
      if (!requireNamespace("cluster", quietly = TRUE)) {
        warning("Skipping 'pam': package 'cluster' not installed.")
      } else {
        for (K in K_seq) {
          fit <- try(cluster::pam(y_train, k = K), silent = TRUE)
          if (inherits(fit, "try-error")) {
            warning("cluster::pam failed for k = ", K, "; skipping.")
          } else {
            cand_list[[paste0("pam_K", K)]] <-
              as.integer(assign_by_centers(as.matrix(fit$medoids), y_mat))
          }
        }
      }

    } else if (alg == "clara") {
      # Large-N k-medoids: CLARA applies PAM to repeated subsamples, so it is
      # O(N) rather than O(N^2) in memory and scales to N = 10,000+.
      if (!requireNamespace("cluster", quietly = TRUE)) {
        warning("Skipping 'clara': package 'cluster' not installed.")
      } else {
        for (K in K_seq) {
          fit <- try(cluster::clara(y_train, k = K, samples = 15,
                                    sampsize = min(N_train, max(100 + 10 * K,
                                                       ceiling(0.05 * N_train))),
                                    pamLike = TRUE), silent = TRUE)
          if (inherits(fit, "try-error")) {
            warning("cluster::clara failed for k = ", K, "; skipping.")
          } else {
            cand_list[[paste0("clara_K", K)]] <-
              as.integer(assign_by_centers(as.matrix(fit$medoids), y_mat))
          }
        }
      }

    } else if (alg == "dbscan") {
      if (!requireNamespace("dbscan", quietly = TRUE)) {
        warning("Skipping 'dbscan': package 'dbscan' not installed.")
      } else {
        d_train <- stats::dist(y_train)
        if (length(d_train) == 0L || max(d_train) == 0) {
          warning("Data (subset) are (almost) constant; skipping 'dbscan'.")
        } else {
          eps_seq <- as.numeric(stats::quantile(d_train, probs = c(0.01, 0.03, 0.05), na.rm = TRUE))
          eps_seq <- unique(eps_seq[eps_seq > 0])
          if (length(eps_seq) == 0L) {
            warning("Could not find positive eps values; skipping 'dbscan'.")
          } else {
            minPts <- max(3L, ncol(y_train) + 1L)
            for (eps_val in eps_seq) {
              fit <- try(dbscan::dbscan(y_train, eps = eps_val, minPts = minPts), silent = TRUE)
              if (inherits(fit, "try-error")) {
                warning("dbscan failed for eps = ", eps_val, "; skipping.")
              } else {
                z_raw <- fit$cluster  # 0 = noise
                if (all(z_raw == 0L)) next
                z_train <- as.integer(factor(z_raw))
                K_found <- max(z_train)
                if (K_found > M) next
                cand_list[[paste0("dbscan_eps", signif(eps_val, 3), "_K", K_found)]] <-
                  as.integer(assign_by_centers(centers_from_labels(z_train, K_found), y_mat))
              }
            }
          }
        }
      }

    } else if (alg == "fuzzyc") {
      if (!requireNamespace("e1071", quietly = TRUE)) {
        warning("Skipping 'fuzzyc': package 'e1071' not installed.")
      } else {
        for (K in K_seq) {
          fit <- try(e1071::cmeans(y_train, centers = K, m = 2), silent = TRUE)
          if (inherits(fit, "try-error")) {
            warning("e1071::cmeans failed for K = ", K, "; skipping.")
          } else {
            cand_list[[paste0("fuzzyc_K", K)]] <-
              as.integer(assign_by_centers(as.matrix(fit$centers), y_mat))
          }
        }
      }

    } else {
      warning("Unknown algorithm name: '", alg, "'. Skipping.")
    }
  } # end for(alg)

  if (length(cand_list) == 0L) {
    stop("No candidate partitions generated. Check installed packages / algorithms / K_seq.")
  }

  ## ---- bind into matrix (candidates x N) ----
  parts <- do.call(rbind, cand_list)
  rownames(parts) <- names(cand_list)

  if (relabel) {
    for (i in seq_len(nrow(parts))) {
      parts[i, ] <- relabel_labels(parts[i, ])
    }
  }

  ## ---- drop duplicate partitions, but keep all algorithm names ----
  keys         <- apply(parts, 1L, paste, collapse = "_")
  keep_idx     <- !duplicated(keys)
  unique_parts <- parts[keep_idx, , drop = FALSE]
  keys_kept    <- keys[keep_idx]

  all_names    <- rownames(parts)
  rownames(unique_parts) <- vapply(
    keys_kept, FUN.VALUE = character(1L),
    FUN = function(k) paste(all_names[keys == k], collapse = "; ")
  )

  if (verbose) {
    cat("Generated", nrow(parts), "raw candidates from ML clustering.\n")
    cat("After removing duplicates:", nrow(unique_parts), "unique candidates remain.\n")
    cat("Training subset size:", N_train, "of", N, "total observations.\n")
  }

  unique_parts
}


## ===========================================================================
##  "DS-ML": candidate set from ML clustering algorithms on the full data
## ===========================================================================
## Section 2.4.2.  Each algorithm is fit for K = 2, ..., M, giving (M - 1)
## candidates per algorithm, and the distinct partitions form C(y). 
##
## @return the list from compose_from_candidates(), plus unique_parts, pZ and
##   the algorithm set actually used
direct_sampling_postpred_ML <- function(y, M = 3, alpha = 1,
                                        sig2_mu = 100, omega = 2, kappa = 1,
                                        nsamp = 1000, n_sub = 1, d,
                                        plot_candidates = FALSE, max_plot_candidates = 5,
                                        plot_pZ_K = FALSE, plot_diagn = FALSE, z_true = NULL,
                                        add_truth = FALSE,
                                        algorithms = NULL, large_n_cutoff = 999,
                                        use_mclust_autoK = FALSE) {

  N <- length(y)
  m <- 0

  if (is.null(algorithms)) {
    algorithms <- if (N > large_n_cutoff) {
      c("kmeans", "gmm", "clara")                     
    } else {
      c("hclust_ward", "pam", "kmeans", "gmm")        
    }
  }

  unique_parts <- get_ml_candidates(
    y                = y,
    M                = M,
    algorithms       = algorithms,
    K_seq            = 2:M,
    use_mclust_autoK = use_mclust_autoK,
    kmeans_nstart    = 10,
    seed             = d * 11111,
    relabel          = TRUE,
    subset_frac      = n_sub,
    verbose          = FALSE
  )

  cand_names <- rownames(unique_parts)
  if (is.null(cand_names)) {
    cand_names <- paste0("cand_", seq_len(nrow(unique_parts)))
    rownames(unique_parts) <- cand_names
  }

  if (add_truth == TRUE) {
    unique_parts <- rbind(z_true, unique_parts)
    cand_names   <- c("Truth", cand_names)
    rownames(unique_parts) <- cand_names
  }

  unnorm_pZ <- log_partition_weight(unique_parts, y, M, alpha, sig2_mu, omega, kappa)
  pZ        <- candidate_probs(unnorm_pZ)

  ## Diagnostics only -- gated so they stay out of the timed path.
  if (plot_pZ_K || plot_diagn) {
    summarize_candidates(unique_parts, pZ, unnorm_pZ, M, z_true = z_true,
                         cand_names = cand_names)
  }
  if (plot_candidates) {
    plot_top_candidates(unique_parts, pZ, y, z_true, max_plot_candidates, n_sub, cand_names)
  }
  if (plot_pZ_K) {
    plot_pZ_vs_K(unique_parts, pZ, unnorm_pZ, N, n_sub)
  }
  if (plot_diagn) {
    plot_candidate_diagnostics(unique_parts, pZ, max_top = 30, n_sub)
  }

  out <- compose_from_candidates(unique_parts, pZ, y, M, alpha,
                                 sig2_mu, omega, kappa, nsamp, m)
  out$unique_parts <- unique_parts
  out$pZ           <- pZ
  out$algorithms   <- algorithms
  out
}


## ===========================================================================
##  "DS-Const-MAP": k-means candidate set initialized at the DS-Const mode
## ===========================================================================

## Build one k-means partition per K in K_seq, each initialized from the
## marginal posterior mode Z_map of the DS-Const candidate set (Section 2.4.2).
## Initial centers for a given K are
##   K <= k_map : the means of the K largest MAP clusters;
##   K >  k_map : all k_map MAP means, then farthest-point augmentation --
##                repeatedly adding the observation furthest from the current
##                center set -- until K distinct center exist.
## Values of K admitting fewer than K distinct center are skipped.
##
## @return matrix of candidate partitions, one row per K, or NULL if none
map_kmeans_candidates <- function(y, Z_map, K_seq, iter_max = 100) {
  N <- length(y)
  Z_map <- as.integer(Z_map)

  sizes  <- tabulate(Z_map, nbins = max(Z_map))
  occ    <- which(sizes > 0)
  mu_map <- vapply(occ, function(k) mean(y[Z_map == k]), numeric(1))
  # order the MAP clusters largest-first so that "the K largest" is well defined
  mu_map <- mu_map[order(sizes[occ], decreasing = TRUE)]
  k_map  <- length(mu_map)

  n_distinct <- length(unique(y))
  cand_list  <- list()

  for (K in K_seq) {
    if (K < 2L || K > n_distinct) next

    if (K <= k_map) {
      centers <- mu_map[seq_len(K)]
    } else {
      centers <- mu_map
      while (length(unique(centers)) < K) {
        d_near <- vapply(y, function(yy) min(abs(yy - centers)), numeric(1))
        if (max(d_near) <= 0) break          # no further distinct center exists
        centers <- c(centers, y[which.max(d_near)])
      }
    }

    centers <- unique(centers)
    if (length(centers) < 2L) next

    km <- try(stats::kmeans(y, centers = matrix(centers, ncol = 1L),
                            iter.max = iter_max, nstart = 1L), silent = TRUE)
    if (inherits(km, "try-error")) {
      warning("MAP-initialized kmeans failed for K = ", K, "; skipping.")
      next
    }
    cand_list[[paste0("mapkm_K", length(centers))]] <- relabel_labels(km$cluster)
  }

  if (length(cand_list) == 0L) return(NULL)

  parts <- do.call(rbind, cand_list)
  rownames(parts) <- names(cand_list)
  parts
}

## DS-Const-MAP.  Run the DS-Const candidate construction, take its marginal
## posterior mode, use that mode to initialize k-means for K = 2, ..., M, and
## sample labels from the resulting candidate set via Equation (12).
##
## The subset draw is seeded on d exactly as in direct_sampling_postpred(), so
## both methods see the same subset and the mode is the DS-Const mode.
##
## @param include_map  also admit the DS-Const mode itself as a candidate;
##   FALSE matches Section 2.4.2, where C(y) is the set of k-means partitions
## @return the list from compose_from_candidates(), plus unique_parts, pZ and
##   the MAP partition Z_map
direct_sampling_obs_refine <- function(y, M=3, alpha=1, sig2_mu=100, omega=2, kappa=1,
                                       nsamp=1000, n_sub=0.1, d,
                                       plot_candidates = FALSE, max_plot_candidates = 5,
                                       plot_pZ_K = FALSE, plot_diagn = FALSE, z_true = NULL,
                                       plot_truth_lab = FALSE, include_map = FALSE,
                                       sub_niter = 11000, sub_nburn = 1000, sub_nthin = 10){

  N <- length(y)
  m <- 0

  ## ---- Step 1: DS-Const candidate set (same seed => same subset as DS-Const)
  set.seed(d*11111)
  y_sub <- y[sample(1:N, floor(n_sub*N), replace=FALSE)]
  # The prediction step below indexes nsamp draws of the subset fit.
  if (floor((sub_niter - sub_nburn) / sub_nthin) < nsamp)
    stop("subset Gibbs yields ", floor((sub_niter - sub_nburn) / sub_nthin),
         " draws but nsamp = ", nsamp, ".  Increase sub_niter or reduce nsamp/sub_nthin.")

  fit <- gaussian_mixture(y=y_sub, N=M, m=0, v=sig2_mu, a=omega, b=kappa, alpha=alpha,
                          niter=sub_niter, nburn=sub_nburn, nthin=sub_nthin)
  ys <- y

  partitions <- matrix(NA, nrow=nsamp, ncol=N)
  for(jj in 1:nsamp){
    un_norm <- sapply(1:N, function(x) dnorm(ys[x], fit$mu[jj,], sqrt(fit$sigma2[jj,]), log=TRUE) + log(fit$w[jj,]))
    pZ <- exp(t(un_norm) - apply(un_norm, 2, max))/apply(exp(t(un_norm) - apply(un_norm, 2, max)),1,sum)
    partitions[jj,] <- sapply(1:N, function(x) sample(1:M, 1, prob=pZ[x,]))
    partitions[jj,] <- relabel_labels(partitions[jj,])
  }

  const_parts <- partitions[!duplicated(partitions), , drop = FALSE]
  const_logw  <- log_partition_weight(const_parts, y, M, alpha, sig2_mu, omega, kappa)

  ## ---- Step 2: marginal posterior mode of the DS-Const candidate set
  Z_map <- as.integer(const_parts[which.max(const_logw), ])

  ## ---- Step 3: MAP-initialized k-means candidate set over K = 2, ..., M
  K_seq        <- 2:min(M, N - 1L)
  unique_parts <- map_kmeans_candidates(y, Z_map, K_seq)

  if (is.null(unique_parts)) {
    warning("No MAP-initialized k-means candidates were produced; ",
            "falling back to the DS-Const MAP partition.")
    unique_parts <- matrix(relabel_labels(Z_map), nrow = 1L)
    rownames(unique_parts) <- "map"
  } else if (include_map) {
    unique_parts <- rbind(map = relabel_labels(Z_map), unique_parts)
  }

  unique_parts <- unique_parts[!duplicated(unique_parts), , drop = FALSE]
  cand_names   <- rownames(unique_parts)

  ## ---- Step 4: Equation (12) over the candidate set
  unnorm_pZ <- log_partition_weight(unique_parts, y, M, alpha, sig2_mu, omega, kappa)
  pZ        <- candidate_probs(unnorm_pZ)

  ## Diagnostics only -- gated so they stay out of the timed path.
  if (plot_pZ_K || plot_diagn) {
    summarize_candidates(unique_parts, pZ, unnorm_pZ, M, z_true = z_true,
                         cand_names = cand_names)
  }
  if (plot_candidates) {
    plot_top_candidates(unique_parts, pZ, y, z_true, max_plot_candidates, n_sub, cand_names)
  }
  if (plot_pZ_K) {
    plot_pZ_vs_K(unique_parts, pZ, unnorm_pZ, N, n_sub)
  }
  if (plot_diagn) {
    plot_candidate_diagnostics(unique_parts, pZ, max_top = 30, n_sub)
  }
  if (!is.null(z_true) && plot_truth_lab) {
    if (length(z_true) != N) {
      warning("Length of z_true does not match N; skipping truth vs k-means plot.")
    } else if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 not installed; skipping truth vs k-means plot.")
    } else {
      Z_best <- unique_parts[which.max(pZ), ]
      print(
        ggplot2::ggplot(data.frame(truth = factor(z_true), kmeans = factor(Z_best)),
                        ggplot2::aes(x = truth, y = kmeans)) +
          ggplot2::geom_jitter(width = 0.15, height = 0.15, alpha = 0.6) +
          ggplot2::labs(x = "True cluster label", y = "k-means cluster label",
                        title = "Truth vs modal candidate (MAP-initialized k-means)") +
          ggplot2::theme_minimal()
      )
    }
  }

  ## ---- Step 5: method of composition over the candidate set
  out <- compose_from_candidates(unique_parts, pZ, y, M, alpha,
                                 sig2_mu, omega, kappa, nsamp, m)
  out$unique_parts <- unique_parts
  out$pZ           <- pZ
  out$Z_map        <- Z_map
  out
}


## ===========================================================================
##  "MC-MCMC": collapsed Gibbs sampler for the marginal posterior of the labels
## ===========================================================================
## Implements the full conditional of Equation (8).  Observation i is removed
## from its component, and the log weight for reallocating it to component j is
## the log ratio of the Equation (6) joint before and after the move. 
## @return list: z (ntotal x N label draws), kpost (occupied clusters)
gibbs_labels_only_intercept <- function(
    y, M = 25, alpha = 1, omega = 2, kappa = 1, sig2_mu = 100,
    niter = 2500, nburn = 500, nthin = 2, seed = NULL
){
  if (!is.null(seed)) set.seed(seed)

  N <- length(y)
  ntotal <- floor((niter - nburn) / nthin)
  stopifnot(ntotal > 0L)

  if (length(alpha) == 1L) alpha <- rep(alpha, M)
  stopifnot(length(alpha) == M)

  y2 <- y * y
  n_grid <- 0:N
  lg_w_n   <- lgamma(omega + 0.5 * n_grid)     # lgamma(omega + n/2)
  log_1psn <- log1p(sig2_mu * n_grid)          # log(1 + sigma^2_beta n)

  lg_w <- function(n) lg_w_n[n + 1L]
  l1psn <- function(n) log_1psn[n + 1L]

  nj <- integer(M)
  S1 <- numeric(M)
  S2 <- numeric(M)

  z <- sample.int(M, N, replace = TRUE)
  for (j in 1:M) {
    idx <- (z == j)
    nj[j] <- sum(idx)
    if (nj[j] > 0L) {
      yj <- y[idx]
      S1[j] <- sum(yj)
      S2[j] <- sum(yj * yj)
    }
  }

  z_keep <- matrix(NA_integer_, nrow = ntotal, ncol = N)
  k_keep <- integer(ntotal)

  logw <- numeric(M)
  pj   <- numeric(M)

  ii <- 1L
  const_norm <- -0.5 * log(2 * pi)

  for (iter in 1:niter) {
    for (i in sample.int(N)) {
      j_old <- z[i]

      nj[j_old] <- nj[j_old] - 1L
      S1[j_old] <- S1[j_old] - y[i]
      S2[j_old] <- S2[j_old] - y2[i]

      log_prior <- log(nj + alpha)

      n  <- nj
      s1 <- S1
      s2 <- S2

      Q_before <- s2 - (sig2_mu * s1 * s1) / (1 + sig2_mu * n)

      n_after  <- n + 1L
      s1_after <- s1 + y[i]
      s2_after <- s2 + y2[i]
      Q_after  <- s2_after - (sig2_mu * s1_after * s1_after) / (1 + sig2_mu * n_after)

      dm <- const_norm +
        (lg_w(n_after) - lg_w(n)) +
        (omega + 0.5 * n)  * log(kappa + 0.5 * Q_before) -
        (omega + 0.5 * n_after) * log(kappa + 0.5 * Q_after) -
        0.5 * (l1psn(n_after) - l1psn(n))

      logw[] <- log_prior + dm

      mmax <- max(logw)
      pj[] <- exp(logw - mmax)
      s <- sum(pj)
      if (!is.finite(s) || s <= 0) {
        j_new <- which.max(logw)
      } else {
        u <- runif(1) * s
        cs <- 0.0
        j_new <- 1L
        for (j in 1:M) { cs <- cs + pj[j]; if (u <= cs) { j_new <- j; break } }
      }

      z[i] <- j_new
      nj[j_new] <- nj[j_new] + 1L
      S1[j_new] <- S1[j_new] + y[i]
      S2[j_new] <- S2[j_new] + y2[i]
    }

    krep <- sum(nj > 0L)
    if ((iter > nburn) && (((iter - nburn) %% nthin) == 0)) {
      z_keep[ii, ] <- z
      k_keep[ii]   <- krep
      ii <- ii + 1L
      if (ii > ntotal) break
    }
  }

  list(z = z_keep, kpost = k_keep)
}


## MC-MCMC.  Method of composition on top of the Equation (8) label chain: take
## z^[t] from the chain, then draw pi, sigma2 and mu from their conditionals
## given z^[t].  The chain must supply at least nsamp draws, which is checked.
direct_sampling_gibbs <- function(y, M=3, alpha=1, sig2_mu=100,
                                  omega=2, kappa=1, nsamp=1000,
                                  niter = 3000, nburn = 1000, nthin = 2){

  N <- length(y)
  m <- 0
  # The composition step below indexes nsamp draws of the label chain.
  ntotal <- floor((niter - nburn) / nthin)
  if (ntotal < nsamp)
    stop("direct_sampling_gibbs: label chain yields ", ntotal,
         " draws but nsamp = ", nsamp,
         ".  Increase niter or reduce nsamp/nthin.")

  gibbs_fit <- gibbs_labels_only_intercept(y, M, alpha, omega,
                                           kappa, sig2_mu, niter = niter,
                                           nburn = nburn, nthin = nthin, seed = NULL)

  mus <- sigma2s <- pis <- matrix(NA, nrow=nsamp, ncol=M)
  Zs <- matrix(NA, nrow=nsamp, ncol=N)
  Ys <- matrix(NA, nrow=nsamp, ncol=N)
  K <- numeric(nsamp)

  for(i in 1:nsamp){
    Z <- gibbs_fit$z[i,]
    K[i] <- length(unique(Z))

    draws <- draw_params_given_Z(Z, y, M, alpha, sig2_mu, omega, kappa, m)

    mus[i,]     <- draws$mu
    sigma2s[i,] <- draws$sigma2
    pis[i,]     <- draws$pi
    Zs[i,]      <- Z
    Ys[i,]      <- draws$mu[Z]
  }

  list(mu=mus, sigma2=sigma2s, pi=pis, Z=Zs, K=K, fitted.values = Ys)
}
