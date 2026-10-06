## ===========================================================================
##  Posterior samplers for the Gaussian mixture-of-regressions model with a
##  Normal-Inverse-Gamma prior.
##
##  Companion to: "Bayesian analysis using a constrained mixture of
##  normal-inverse-gamma models", Section 3.4 and Appendix E.2.
##  Driver: sim_study_xy.R
##
##  ---------------------------------------------------------------------------
##  Model (manuscript Equation 1)
##
##      Y_i | Z_ik = 1        ~  N(x_i' beta_k, sigma2_k)
##      Z_i | pi              ~  Categorical(pi_1, ..., pi_K)
##      pi                    ~  Dirichlet(alpha_1, ..., alpha_K)
##      beta_k | sigma2_k     ~  MVN(0_p, sigma2_k * sigma2_beta * I_p)
##      sigma2_k              ~  InverseGamma(omega, kappa)
##      K                     ~  Categorical(rho_1, ..., rho_M)
##
##  The intercept-only samplers in sim_study_y_functions.R are the special case
##  p = 1, x_i = 1; this file keeps the general p x 1 covariate vector, so the
##  scalar identities used there become p x p Cholesky solves here.
##
##  ---------------------------------------------------------------------------
##  Argument naming
##
##      X           N x p covariate matrix, intercept column included
##      sig2_beta   sigma2_beta, the prior scale on the regression coefficients
##      omega       inverse-gamma shape
##      kappa       inverse-gamma scale
##      alpha       symmetric Dirichlet concentration
##      M           K, the number of available mixture components
##      nsamp       number of posterior draws returned
##      Z_sampling  which candidate generator to use; see direct_sampling_xy()
##
##  Simulation-study values (Section 3.1): alpha = 1, sigma2_beta = 10,
##  omega = 2, kappa = 1, M = 25 (M = min(M, N) when N = 10).
##
##  ---------------------------------------------------------------------------
##  Contents
##
##  Shared building blocks
##      relabel_labels           canonical 1..k relabelling of an allocation
##      dedup_rows               drop repeated candidate partitions
##      log_sum_exp              stable log of a sum of exponentials
##      cluster_log_marginal     the k-th collapsed factor of Equation (6)
##      cluster_Q                quadratic form Q_k of Equation (4)
##      partition_log_terms_xy   partition-only pieces of Equation (6)
##      compute_pZ_candidates    Equation (12) normalisation, K fixed at M
##      compute_pZK_candidates   the same when K is random
##
##  Samplers (manuscript name in quotes)
##      direct_sampling_xy                       one entry point for the
##                                               "Direct Sampler", "DS-Const",
##                                               "DS-ML" and an MCMC baseline,
##                                               selected by Z_sampling
##      full_mcmc_sampling_regression_dependent  "MCMC", unconstrained Gibbs
##      gibbs_labels_only_regression             label chain of Equation (8)
##      direct_sampling_gibbs_regression         "MC-MCMC"
##
##  Candidate-set construction and diagnostics
##      make_candidates            dispatches to the requested generator
##      get_ml_candidates          partitions from ML clustering algorithms
##      print_candidate_summaries, plot_truth_and_candidates
##
##  Requires: mvtnorm (rmvnorm), salso (enumerate.partitions, ARI),
##  miscPack (gaussian_mixture); optionally mclust, cluster, kernlab, dbscan,
##  ggplot2 for the ML candidate generators and diagnostics.  The Dirichlet and
##  inverse-gamma draws are defined below rather than taken from MCMCpack.
## ===========================================================================

## ---------------------------------------------------------------------------
## The inverse-gamma parameterisation is the manuscript's: for
## sigma2 ~ IG(omega, kappa) the density is proportional to
## x^(-omega-1) exp(-kappa/x), with mean kappa/(omega-1) when omega > 1.
## ---------------------------------------------------------------------------

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


## ===========================================================================
##  Shared building blocks
## ===========================================================================

## Relabel positive integer labels canonically as 1, 2, ... in order of first
## appearance, so that allocations inducing the same partition compare equal.
## keep_zero = TRUE preserves 0 as a "noise" label, as dbscan returns.
relabel_labels <- function(z, keep_zero = FALSE) {
  z <- as.integer(z)
  if (keep_zero) {
    pos <- z > 0L
    if (!any(pos)) return(z)
    zpos <- z[pos]
    z[pos] <- as.integer(factor(zpos, levels = unique(zpos)))
    return(z)
  } else {
    return(as.integer(factor(z, levels = unique(z))))
  }
}

## Drop repeated rows of a candidate matrix.  .
dedup_rows <- function(mat) {
  mat <- as.matrix(mat)
  keys <- apply(mat, 1, paste, collapse = "_")
  mat[!duplicated(keys), , drop = FALSE]
}

## log(sum(exp(v))) computed by shifting out the maximum.
log_sum_exp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

## The k-th factor of the collapsed marginal in Equation (6): the density of
## the observations in one cluster with beta_k and sigma2_k integrated out,
##
##   -n_k/2 log(2 pi) - 1/2 log|I_p + sigma2_beta X_k'X_k|
##   + omega log(kappa) - lgamma(omega)
##   + lgamma(omega + n_k/2) - (omega + n_k/2) log(kappa + Q_k/2),
##
## where A = I_p + sigma2_beta X_k'X_k, b = X_k'y_k and
## Q_k = y_k'y_k - sigma2_beta b'A^{-1}b.
## An empty cluster contributes 0, matching the product over 
## occupied components in Equation (6).
##
## @return scalar log marginal.
cluster_log_marginal <- function(yk, Xk, sig2_beta, omega, kappa) {
  nk <- length(yk)
  if (nk == 0L) return(0)
  
  p <- ncol(Xk)
  A <- diag(p) + sig2_beta * crossprod(Xk)       # p x p
  R <- chol(A)
  logdetA <- 2 * sum(log(diag(R)))
  
  b <- crossprod(Xk, yk)                         # p x 1
  Ainv_b <- backsolve(R, forwardsolve(t(R), b))  # A^{-1} b
  Q <- sum(yk^2) - sig2_beta * sum(b * Ainv_b)
  Q <- max(Q, 0)  
  
  out <- (-0.5 * nk * log(2*pi)) -
    0.5 * logdetA +
    omega * log(kappa) - lgamma(omega) +
    lgamma(omega + 0.5 * nk) -
    (omega + 0.5 * nk) * log(kappa + 0.5 * Q)
  
  out
}

## The quadratic form Q_k above on its own, which is the only data-dependent
## part of the inverse-gamma scale in Equation (4), sigma2_k | Z ~
## IG(omega + n_k/2, kappa + Q_k/2).
cluster_Q <- function(yk, Xk, sig2_beta) {
  nk <- length(yk)
  if (nk == 0L) return(0)
  
  p <- ncol(Xk)
  A <- diag(p) + sig2_beta * crossprod(Xk)
  R <- chol(A)
  
  b <- crossprod(Xk, yk)
  Ainv_b <- backsolve(R, forwardsolve(t(R), b))
  Q <- sum(yk^2) - sig2_beta * sum(b * Ainv_b)
  max(Q, 0)
}

## Normalized p(Z | y, X) over a candidate set, with K fixed at M -- the
## constrained posterior of Equation (12) restricted to the candidate label
## space C(y).  Each unnormalized weight is the Dirichlet-multinomial term of
## Equation (6) times the collapsed cluster marginals above.
##
## @param parts  candidates x N matrix of allocations with labels in 1..M.
## @return vector of probabilities, one per candidate, summing to 1.
compute_pZ_candidates <- function(parts, y, X, M, alpha, sig2_beta, omega, kappa) {
  parts <- as.matrix(parts)
  N <- length(y)
  stopifnot(ncol(parts) == N)
  
  if (length(alpha) == 1L) alpha <- rep(alpha, M)
  stopifnot(length(alpha) == M)
  
  logw <- numeric(nrow(parts))
  
  for (r in seq_len(nrow(parts))) {
    z <- parts[r, ]
    nk <- tabulate(z, nbins = M)
    
    log_dm <- lgamma(sum(alpha)) - lgamma(N + sum(alpha)) +
      sum(lgamma(nk + alpha)) - sum(lgamma(alpha))
    
    log_like <- 0
    for (k in 1:M) {
      idx <- which(z == k)
      if (length(idx) == 0L) next
      log_like <- log_like + cluster_log_marginal(
        yk = y[idx],
        Xk = X[idx, , drop = FALSE],
        sig2_beta = sig2_beta,
        omega = omega,
        kappa = kappa
      )
    }
    
    logw[r] <- log_dm + log_like
  }
  
  lse <- log_sum_exp(logw)
  p <- exp(logw - lse)
  p / sum(p)
}

###############################################
## Partition-only pieces of Equation (6), reusable across K
###############################################
## Returns, for each candidate partition: the number of occupied blocks k, the
## occupied-block Dirichlet contribution sum_{j: n_j>0} lgamma(n_j + alpha), and
## the product of collapsed cluster marginals.  None depend on K.
partition_log_terms_xy <- function(parts, y, X, M, alpha, sig2_beta, omega, kappa) {
  parts <- as.matrix(parts)
  if (length(alpha) != 1L) stop("partition_log_terms_xy assumes a scalar (symmetric) alpha.")
  R <- nrow(parts)
  k_occ <- integer(R); log_dir_occ <- numeric(R); log_lik <- numeric(R)

  for (r in seq_len(R)) {
    z  <- parts[r, ]
    nk <- tabulate(z, nbins = M)
    occ <- which(nk > 0)
    k_occ[r]       <- length(occ)
    log_dir_occ[r] <- sum(lgamma(nk[occ] + alpha))
    ll <- 0
    for (kk in occ) {
      idx <- which(z == kk)
      ll <- ll + cluster_log_marginal(y[idx], X[idx, , drop = FALSE],
                                      sig2_beta, omega, kappa)
    }
    log_lik[r] <- ll
  }
  list(k = k_occ, log_dir_occ = log_dir_occ, log_lik = log_lik)
}

## Joint weights p(Z, K | y, X) over a candidate set when K is random, used by
## the enumerating direct sampler of Section 3.1.
##
##   log w(P, K) = log rho_K
##               + lgamma(K + 1) - lgamma(K - k + 1)        # labellings
##               + lgamma(K*alpha) - lgamma(N + K*alpha)    # Dirichlet normalizer
##               + sum_{j: n_j>0} lgamma(n_j + alpha) - k*lgamma(alpha)
##               + [collapsed cluster marginals, independent of K]
##
## Candidates with more occupied blocks than K are impossible and are left at
## weight zero.  This mirrors direct_sampling() in the intercept-only code.
##
## @param prior_for_M  rho_1, ..., rho_M, the prior on K.
## @return candidates x M matrix of unnormalized weights, one column per K.
compute_pZK_candidates <- function(parts, y, X, M, alpha, sig2_beta, omega, kappa,
                                   prior_for_M) {
  N  <- length(y)
  tm <- partition_log_terms_xy(parts, y, X, M, alpha, sig2_beta, omega, kappa)
  kP <- tm$k
  base <- tm$log_dir_occ - kP * lgamma(alpha) + tm$log_lik

  log_w <- matrix(-Inf, nrow = nrow(parts), ncol = M)
  for (K in seq_len(M)) {
    ok <- which(kP <= K)
    if (!length(ok)) next
    log_w[ok, K] <- log(prior_for_M[K]) +
      (lgamma(K + 1) - lgamma(K - kP[ok] + 1)) +
      (lgamma(K * alpha) - lgamma(N + K * alpha)) +
      base[ok]
  }
  w <- exp(log_w - max(log_w[is.finite(log_w)]))
  w[!is.finite(log_w)] <- 0
  w
}

## ===========================================================================
##  Candidate-set construction
## ===========================================================================

## Build the candidate label space C(y) of Section 2.4, dispatching on the
## requested strategy:
##
##   "ds"      every set partition of 1..N with at most M blocks.  Exhaustive,
##             so only feasible for N < 13; this is the "Direct Sampler".
##   "ds_obs"  allocations drawn from a mixture-of-regressions fit to a random
##             subset of n_sub * N observations, then applied to all N.  This is
##             "DS-Const", and requires miscPack::gaussian_mixture().
##   "ds_pred" as "ds_obs", but observations are allocated by their posterior
##             predictive mean rather than their observed response.
##   "ml"      partitions returned by classical clustering algorithms run at a
##             range of K.  This is "DS-ML".
##
##
## @param n_sub    subset fraction for "ds_obs"/"ds_pred".
## @param nsamp    number of allocations drawn before deduplication.
## @param sub_a,sub_b  inverse-gamma hyperparameters for the subset fit;
##                 default to omega and kappa so the generator matches the
##                 inferential model.
## @param ml_features  what the ML algorithms cluster on: the response and
##                 covariates jointly ("yx"), the covariates alone ("x"),
##                 residuals from a single OLS fit ("resid", "resid_only"),
##                 fitted values with residuals ("yhat_resid"), or the response
##                 alone ("y").
## @return unique candidates x N matrix of allocations with labels in 1..M.
make_candidates <- function(y, X, M,
                            Z_sampling = c("ds_obs","ds_pred","ml","ds"),
                            nsamp = 1000, n_sub = 0.1,
                            alpha = 1, sig2_beta = 100, omega = 2, kappa = 1,
                            sub_a = NULL, sub_b = NULL,
                            ml_algorithms = c("hclust_ward","pam","kmeans","spectral","dbscan","gmm"),
                            ml_K_seq = NULL,
                            ml_features = c("yx","x","resid","yhat_resid","y","resid_only"),
                            ml_scale = FALSE,
                            ml_fast = FALSE,
                            spectral_nystrom_sample = NULL,
                            sub_niter = 110000, sub_nburn = 10000, sub_nthin = 100,
                            seed = NULL,
                            verbose = TRUE) {
  
  Z_sampling <- match.arg(Z_sampling)
  ml_features <- match.arg(ml_features)
  
  if (!is.null(seed)) set.seed(seed)
  N <- length(y)
  X <- as.matrix(X)
  
  scale_safe <- function(A) {
    A <- as.matrix(A)
    mu <- colMeans(A)
    sdv <- apply(A, 2, sd)
    sdv[sdv == 0] <- 1
    sweep(sweep(A, 2, mu, "-"), 2, sdv, "/")
  }
  
  ## Full enumeration: Bell(N) set partitions, restricted to at most M blocks.
  if (Z_sampling == "ds") {
    if (N >= 13) stop("Z_sampling='ds' only implemented for N < 13.")
    all_parts <- enumerate.partitions(N)
    parts <- all_parts[apply(all_parts, 1, function(z) length(table(z))) <= M, , drop = FALSE]
    return(dedup_rows(parts))
  }
  
  if (Z_sampling %in% c("ds_obs","ds_pred")) {
    if (!exists("gaussian_mixture")) stop("gaussian_mixture() not found for ds_obs/ds_pred.")
    samp  <- sample.int(N, max(1, floor(n_sub * N)), replace = FALSE)
    y_sub <- y[samp]
    X_sub <- X[samp, , drop = FALSE]
    
    ## Inverse-gamma hyperparameters for the subset fit.  Defaulting to the
    ## omega/kappa of the inferential model keeps the generator consistent
    ## with Equation (1); pass sub_a/sub_b to override.
    if (is.null(sub_a)) sub_a <- omega
    if (is.null(sub_b)) sub_b <- kappa

    ## The subset fit must supply at least nsamp draws for the allocation loop
    ## below.  sub_niter is settable from the driver (FMM_SUB_NITER), so the
    ## coupling is checked rather than assumed.
    if (floor((sub_niter - sub_nburn) / sub_nthin) < nsamp)
      stop("subset Gibbs yields ", floor((sub_niter - sub_nburn) / sub_nthin),
           " draws but nsamp = ", nsamp,
           ".  Increase sub_niter or reduce nsamp/sub_nthin.")

    fit <- gaussian_mixture(
      y = y_sub, Xmat = X_sub, N = M,
      m = 0, v = sig2_beta, a = sub_a, b = sub_b, alpha = alpha,
      niter = sub_niter, nburn = sub_nburn, nthin = sub_nthin
    )
    
    if (Z_sampling == "ds_pred") {
      ## Posterior predictive mean for every observation, averaged over the
      ## subset fit's draws and its mixture components.
      p <- ncol(X)
      ys_all <- sapply(seq_len(N), function(i) {
        sapply(seq_len(nsamp), function(t) {
          Bt <- matrix(fit$beta[t, ], nrow = M, byrow = TRUE)
          sum((X[i, ] %*% t(Bt)) * fit$w[t, ])
        })
      })
      ys <- colMeans(ys_all)  
      Xs <- X
    } else {
      ys <- y
      Xs <- X
    }
    
    ## One allocation per retained draw: assign each observation to a component
    ## with probability proportional to its weighted regression likelihood.
    p <- ncol(X)
    parts <- matrix(NA_integer_, nrow = nsamp, ncol = N)
    for (jj in seq_len(nsamp)) {
      Bt <- matrix(fit$beta[jj, ], nrow = M, byrow = TRUE)
      un_norm <- sapply(seq_len(N), function(i) {
        dnorm(ys[i],
              mean = as.numeric(Xs[i, , drop = FALSE] %*% t(Bt)),
              sd = sqrt(fit$sigma2[jj, ]),  
              log = TRUE) + log(fit$w[jj, ])
      })
      un_norm <- matrix(un_norm, nrow = M, ncol = N)
      
      un_norm <- sweep(un_norm, 2, apply(un_norm, 2, max), "-")
      P <- exp(un_norm)
      P <- sweep(P, 2, colSums(P), "/")
      z <- apply(P, 2, function(pr) sample.int(M, 1, prob = pr))
      parts[jj, ] <- relabel_labels(z)
    }
    return(dedup_rows(parts))
  }
  
  if (Z_sampling == "ml") {
    if (!exists("get_ml_candidates")) stop("get_ml_candidates() not found for Z_sampling='ml'.")

    ## Feature matrix the clustering algorithms see (Section 2.4.2).  Constant
    ## columns are dropped from "yx" because the intercept carries no signal.
    if (ml_features == "yx") {
      const_cols <- apply(X, 2, sd) < 1e-8
      feat <- cbind(y = y, X = X[, !const_cols, drop = FALSE])
    } else if (ml_features == "x") {
      feat <- X
    } else if (ml_features == "resid") {
      bhat <- lm.fit(X, y)$coefficients
      r <- y - as.numeric(X %*% bhat)
      feat <- cbind(resid = r, X)
    } else if (ml_features == "yhat_resid") {
      bhat <- lm.fit(X, y)$coefficients
      yhat <- as.numeric(X %*% bhat)
      r <- y - yhat
      feat <- cbind(yhat = yhat, resid = r)
    } else if (ml_features == "y") {
      feat <- matrix(y, ncol = 1L)
    } else { # resid_only
      bhat <- lm.fit(X, y)$coefficients
      r <- y - as.numeric(X %*% bhat)
      feat <- matrix(r, ncol = 1L)
    }
    
    if (ml_scale) feat <- scale_safe(feat)
    
    parts <- get_ml_candidates(
      y = feat, M = M,
      algorithms = ml_algorithms,
      K_seq = ml_K_seq,
      seed = seed,
      relabel = TRUE,
      ml_fast = ml_fast,
      spectral_nystrom_sample = spectral_nystrom_sample,
      verbose = verbose
    )
    return(dedup_rows(parts))
  }
  
  stop("Unknown Z_sampling.")
}


## ===========================================================================
##  Samplers
## ===========================================================================

## Direct (MCMC-free) posterior sampler for the mixture of regressions, and the
## single entry point for four of the compared methods.  Z_sampling selects the
## candidate generator and hence the manuscript method:
##
##   "ds"       "Direct Sampler"  full enumeration, N < 13
##   "ds_obs"   "DS-Const"        subset-fit allocations
##   "ds_pred"  --                as "ds_obs", allocating on predicted y
##   "ml"       "DS-ML"           ML clustering partitions
##   "mcmc"     --                Gibbs label updates instead of a candidate
##                                set, for reference within this function
##
## Each of the nsamp draws is independent for the candidate-based settings: a
## label configuration is drawn from the constrained posterior of Equation (12)
## and the parameters follow by the method of composition, Equations (3)-(5).
##
## @param prior_for_M  a prior rho_1..rho_M on K.  Honored only when
##                     Z_sampling = "ds": the constrained methods enumerate a
##                     finite candidate set at K = M, so K is not random there.
## @return list(beta [nsamp x M x p], sigma2, pi, Z, K, fitted.values), plus
##         candidate_Z and either candidate_pZ or pKgY for the non-MCMC paths.
direct_sampling_xy <- function(
    y, X, M = 3,
    alpha = 1,
    sig2_beta = 100, omega = 2, kappa = 1,
    Z_sampling = c("ds_obs","ds_pred","ml","ds","mcmc"),
    nsamp = 1000, n_sub = 0.1,
    ml_algorithms = c("hclust_ward","pam","kmeans","spectral","dbscan","gmm"),
    ml_K_seq = NULL,
    ml_features = c("yx","x","resid","yhat_resid","y","resid_only"),
    ml_scale = FALSE,
    ml_fast = FALSE,
    spectral_nystrom_sample = NULL,
    sub_niter = 110000, sub_nburn = 10000, sub_nthin = 100,
    seed = NULL,
    verbose = TRUE,
    prior_for_M = NULL
) {
  Z_sampling <- match.arg(Z_sampling)
  X <- as.matrix(X)
  N <- length(y)
  p <- ncol(X)
  stopifnot(nrow(X) == N)
  if (!is.null(seed)) set.seed(seed)
  
  ## Candidate set and its posterior weights, once, before the draw loop.
  if (Z_sampling != "mcmc") {
    candZ <- make_candidates(
      y, X, M,
      Z_sampling = Z_sampling,
      nsamp = nsamp, n_sub = n_sub,
      alpha = alpha, sig2_beta = sig2_beta, omega = omega, kappa = kappa,
      ml_algorithms = ml_algorithms, ml_K_seq = ml_K_seq,
      ml_features = ml_features, ml_scale = ml_scale, ml_fast = ml_fast,
      spectral_nystrom_sample = spectral_nystrom_sample,
      sub_niter = sub_niter, sub_nburn = sub_nburn, sub_nthin = sub_nthin,
      seed = seed, verbose = verbose
    )
    ## Random K applies only to the enumerating direct sampler (Section 3.1);
    ## every other method fixes K = M.  A malformed prior_for_M falls back to
    ## the fixed-K path rather than failing.
    var_K <- (Z_sampling == "ds") && !is.null(prior_for_M) &&
             abs(sum(prior_for_M) - 1) < 1e-8 && length(prior_for_M) == M
    if (var_K) {
      wK   <- compute_pZK_candidates(candZ, y, X, M, alpha, sig2_beta, omega,
                                     kappa, prior_for_M)
      pKgY <- colSums(wK) / sum(wK)
    } else {
      pZ <- compute_pZ_candidates(candZ, y, X, M, alpha, sig2_beta, omega, kappa)
    }
  } else var_K <- FALSE
  
  betas   <- array(NA_real_, dim = c(nsamp, M, p))
  sigma2s <- matrix(NA_real_, nrow = nsamp, ncol = M)
  pis     <- matrix(NA_real_, nrow = nsamp, ncol = M)
  Zs      <- matrix(NA_integer_, nrow = nsamp, ncol = N)
  Ys      <- matrix(NA_real_, nrow = nsamp, ncol = N)
  K       <- integer(nsamp)
  
  z <- sample.int(M, N, replace = TRUE)
  beta <- matrix(0, M, p)
  sigma2 <- rep(1, M)
  pii <- rep(1/M, M)
  
  if (length(alpha) == 1L) alpha_vec <- rep(alpha, M) else alpha_vec <- alpha
  
  for (s in seq_len(nsamp)) {
    
    ## Labels: a Gibbs update for "mcmc", otherwise an independent draw from
    ## the constrained posterior over the candidate set.
    if (Z_sampling == "mcmc") {
      logpk <- matrix(NA_real_, M, N)
      for (k in 1:M) {
        mu <- as.numeric(X %*% beta[k, ])
        logpk[k, ] <- dnorm(y, mu, sqrt(sigma2[k]), log = TRUE) + log(pii[k])
      }
      logpk <- sweep(logpk, 2, apply(logpk, 2, max), "-")
      pk <- exp(logpk)
      pk <- sweep(pk, 2, colSums(pk), "/")
      z <- apply(pk, 2, function(pr) sample.int(M, 1, prob = pr))
      z <- relabel_labels(z)
    } else if (var_K) {
      K_s <- sample.int(M, 1, prob = pKgY)
      z   <- candZ[sample.int(nrow(candZ), 1, prob = wK[, K_s]), ]
    } else {
      ridx <- sample.int(nrow(candZ), 1, prob = pZ)
      z <- candZ[ridx, ]
    }
    K_active <- if (var_K) K_s else M
    
    Zs[s, ] <- z
    K[s] <- length(unique(z))
    
    nk <- tabulate(z, nbins = M)
    ## pi | Z ~ Dirichlet(alpha + n), Equation (5).  With random K only the
    ## first K_active components exist; the remainder are padded with zero
    ## weight 
    pii <- numeric(M)
    pii[seq_len(K_active)] <-
      as.numeric(rdirichlet(1, (alpha_vec + nk)[seq_len(K_active)]))
    pis[s, ] <- pii
    
    ## sigma2_k | Z ~ IG(omega + n_k/2, kappa + Q_k/2), Equation (4), with beta
    ## integrated out.  An empty component reverts to the IG(omega, kappa) prior.
    for (k in 1:M) {
      idx <- which(z == k)
      nk_k <- length(idx)
      Qk <- if (nk_k == 0L) 0 else cluster_Q(y[idx], X[idx, , drop = FALSE], sig2_beta)
      sigma2[k] <- rinvgamma(
        1,
        shape = omega + 0.5 * nk_k,
        scale = kappa + 0.5 * Qk
      )
      sigma2[k] <- max(sigma2[k], 1e-12)
    }
    sigma2s[s, ] <- sigma2
    
    ## beta_k | sigma2_k, Z ~ MVN(A^{-1} X_k'y_k, sigma2_k A^{-1}) with
    ## A = X_k'X_k + I_p/sigma2_beta, Equation (3).  An empty component reverts
    ## to the MVN(0_p, sigma2_k sigma2_beta I_p) prior, which is the same
    ## expression with X_k'X_k = 0.
    for (k in 1:M) {
      idx <- which(z == k)
      if (length(idx) == 0L) {
        Prec <- (1 / sig2_beta) * diag(p)
        R <- chol(Prec)
        mu <- rep(0, p)
      } else {
        Xk <- X[idx, , drop = FALSE]
        yk <- y[idx]
        Prec <- crossprod(Xk) + (1 / sig2_beta) * diag(p)
        R <- chol(Prec)
        mu <- backsolve(R, forwardsolve(t(R), crossprod(Xk, yk)))
      }
      V <- sigma2[k] * chol2inv(R)
      beta[k, ] <- as.numeric(mvtnorm::rmvnorm(1, mean = mu, sigma = V))
    }
    betas[s, , ] <- beta
    
    ## Fitted values under this draw's own allocation, used by the density and
    ## KS metrics in the driver.
    Ys[s, ] <- vapply(seq_len(N), function(i) as.numeric(X[i, ] %*% beta[z[i], ]), numeric(1))
  }
  
  out <- list(beta = betas, sigma2 = sigma2s, pi = pis, Z = Zs, K = K, fitted.values = Ys)
  if (Z_sampling != "mcmc") {
    out$candidate_Z <- candZ
    if (var_K) out$pKgY <- pKgY else out$candidate_pZ <- pZ
  }
  out
}

## Candidate partitions from classical clustering algorithms (Section 2.4.2).
## Every algorithm is run at each K in K_seq, and the resulting partitions are
## pooled and deduplicated.  Each is wrapped so that a failure or a missing
## optional package drops that algorithm rather than the whole run.
##
## @param y        N x q feature matrix built by make_candidates(); constant
##                 columns (an intercept) are dropped before clustering.
## @param M        upper bound on the number of components.
## @param K_seq    values of K to run each algorithm at; defaults to 2..min(M, N).
## @param use_mclust_autoK  if TRUE, mclust chooses K by BIC and contributes a
##                 single partition instead of one per K.  FALSE in the study,
##                 so the GMM is represented at every K like the others.
## @param relabel  canonicalize each partition to 1..K before deduplication.
## @return unique candidates x N integer matrix, rows named by algorithm and K.
get_ml_candidates <- function(
    y,
    M,
    algorithms = c("kmeans", "gmm", "hclust_ward", "pam",
                   "spectral", "dbscan"),
    K_seq = NULL,
    use_mclust_autoK = FALSE,
    kmeans_nstart = 3,
    seed = NULL,
    relabel = FALSE,
    ml_scale = FALSE,
    ml_fast = FALSE,
    spectral_nystrom_sample = NULL,
    verbose = TRUE
) {
  if (is.vector(y)) {
    y_mat <- matrix(as.numeric(y), ncol = 1L)
  } else {
    y_mat <- as.matrix(y)
  }
  storage.mode(y_mat) <- "double"
  
  N <- nrow(y_mat)
  if (N < 2L) stop("get_ml_candidates: need at least 2 observations.")
  
  if (any(!is.finite(y_mat))) {
    stop("get_ml_candidates: y contains NA/NaN/Inf; clean data before clustering.")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  ## Optionally scale
  if (ml_scale) {
    scale_safe <- function(A) {
      A <- as.matrix(A)
      mu <- colMeans(A)
      sdv <- apply(A, 2, sd)
      sdv[sdv == 0] <- 1
      sweep(sweep(A, 2, mu, "-"), 2, sdv, "/")
    }
    y_mat <- scale_safe(y_mat)
  }
  
  ## Drop constant columns: an intercept carries no clustering signal, and a
  ## zero-variance column breaks distance-based algorithms.
  if (ncol(y_mat) > 1L) {
    sds <- apply(y_mat, 2, sd)
    const_cols <- sds < .Machine$double.eps
    if (any(const_cols)) {
      if (verbose) {
        nm <- colnames(y_mat)
        if (is.null(nm)) nm <- which(const_cols)
        else nm <- nm[const_cols]
        message(
          "get_ml_candidates: dropping ", sum(const_cols),
          " constant column(s): ", paste(nm, collapse = ", ")
        )
      }
      y_train <- y_mat[, !const_cols, drop = FALSE]
    } else {
      y_train <- y_mat
    }
  } else {
    y_train <- y_mat
  }
  
  if (ncol(y_train) == 0L) {
    stop("get_ml_candidates: all columns are constant; cannot cluster.")
  }
  
  ## Default grid of K, capped at N since no partition can have more blocks
  ## than observations.
  if (is.null(K_seq)) {
    M_eff <- min(M, N)
    K_seq <- if (M_eff <= 1L) 1L else 2:M_eff
  }
  K_seq <- as.integer(K_seq)
  K_seq <- K_seq[K_seq >= 1L & K_seq <= N]
  
  cand_list <- list()
  
  relabel_partition <- function(z) {
    z <- as.integer(z)
    pos <- z > 0L
    if (!any(pos)) return(z)
    labs <- sort(unique(z[pos]))
    map  <- seq_along(labs)
    names(map) <- labs
    z[pos] <- map[as.character(z[pos])]
    z
  }
  
  ## ---------------------------------------------------------------------
  ## One block per algorithm.  Each contributes zero or more named candidates.
  ## ---------------------------------------------------------------------
  for (alg in algorithms) {
    if (verbose) message("Running algorithm: ", alg)
    
    if (alg == "kmeans") {
      for (K in K_seq) {
        if (K > N) next
        fit <- try(
          stats::kmeans(y_train, centers = K, nstart = kmeans_nstart),
          silent = TRUE
        )
        if (inherits(fit, "try-error")) {
          if (verbose) message("  kmeans_K", K, " failed; skipping.")
          next
        }
        cand_list[[paste0("kmeans_K", K)]] <- as.integer(fit$cluster)
      }
      
    } else if (alg == "gmm") {
      if (!requireNamespace("mclust", quietly = TRUE)) {
        warning("Skipping 'gmm': package 'mclust' not installed.")
        next
      }
      
      if (isTRUE(use_mclust_autoK)) {
        fit <- tryCatch(
          mclust::Mclust(y_train, verbose = FALSE),
          error = function(e) {
            if (verbose) message("  gmm_autoK error: ", conditionMessage(e))
            NULL
          }
        )
        if (!is.null(fit)) {
          cand_list[["gmm_autoK"]] <- as.integer(fit$classification)
        }
      } else {
        for (K in K_seq) {
          if (K > N) next
          fit <- tryCatch(
            mclust::Mclust(y_train, G = K, verbose = FALSE,
                           modelNames = if (!ml_fast) NULL
                                        else if (ncol(y_train) == 1L) "V" else "VVV"),
            error = function(e) {
              if (verbose) message("  gmm_K", K, " error: ", conditionMessage(e))
              NULL
            }
          )
          if (!is.null(fit)) {
            cand_list[[paste0("gmm_K", K)]] <- as.integer(fit$classification)
          }
        }
      }
      
    } else if (alg == "hclust_ward") {
      dmat <- try(stats::dist(y_train), silent = TRUE)
      if (inherits(dmat, "try-error")) {
        if (verbose) message("  hclust: dist() failed; skipping.")
        next
      }
      hc <- try(stats::hclust(dmat, method = "ward.D2"), silent = TRUE)
      if (inherits(hc, "try-error")) {
        if (verbose) message("  hclust: hclust() failed; skipping.")
        next
      }
      for (K in K_seq) {
        if (K > N) next
        z <- stats::cutree(hc, k = K)
        cand_list[[paste0("hclust_ward_K", K)]] <- as.integer(z)
      }
      
    } else if (alg == "pam") {
      if (!requireNamespace("cluster", quietly = TRUE)) {
        warning("Skipping 'pam': package 'cluster' not installed.")
        next
      }
      for (K in K_seq) {
        if (K > N) next
        fit <- try(cluster::pam(y_train, k = K), silent = TRUE)
        if (inherits(fit, "try-error")) {
          if (verbose) message("  pam_K", K, " failed; skipping.")
          next
        }
        cand_list[[paste0("pam_K", K)]] <- as.integer(fit$clustering)
      }
      
    }else if (alg == "dbscan") {
      if (!requireNamespace("dbscan", quietly = TRUE)) {
        warning("Skipping 'dbscan': package 'dbscan' not installed.")
        next
      }
      if (N < 2L) next
      
      ## eps is set to the median 5-nearest-neighbour distance, the usual
      ## data-driven default, and minPts to the dimension plus one.
      k_nn <- min(5L, max(1L, N - 1L))
      nn_d <- try(dbscan::kNNdist(y_train, k = k_nn), silent = TRUE)
      if (inherits(nn_d, "try-error")) {
        if (verbose) message("  dbscan: kNNdist() failed; skipping.")
        next
      }
      eps <- median(nn_d)
      if (!is.finite(eps) || eps <= 0) {
        if (verbose) message("  dbscan: invalid eps; skipping.")
        next
      }
      
      p_dim  <- ncol(y_train)
      minPts <- max(5L, p_dim + 1L)
      
      fit <- try(dbscan::dbscan(y_train, eps = eps, minPts = minPts),
                 silent = TRUE)
      if (inherits(fit, "try-error")) {
        if (verbose) message("  dbscan: dbscan() failed; skipping.")
        next
      }
      
      ## dbscan labels outliers 0; fold them into one extra cluster so the
      ## result is a partition, then collapse to at most M blocks.
      z <- as.integer(fit$cluster)
      if (all(z == 0L)) {
        if (verbose) message("  dbscan: all points noise; skipping.")
        next
      }
      if (any(z == 0L)) {
        z[z == 0L] <- max(z) + 1L
      }
      z <- relabel_partition(z)
      if (max(z) > M) {
        if (verbose) {
          message("  dbscan produced ", max(z),
                  " clusters; collapsing to at most M = ", M, ".")
        }
        z[z > M] <- M
        z <- relabel_partition(z)
      }
      cand_list[[paste0("dbscan_K", max(z))]] <- z
    }   else if (alg == "spectral") {
      if (!requireNamespace("kernlab", quietly = TRUE)) {
        warning("Skipping 'spectral': package 'kernlab' not installed.")
        next
      }
      for (K in K_seq) {
        if (K > N) next
        fit <- try(
          if (!is.null(spectral_nystrom_sample)) {
            kernlab::specc(y_train, centers = K, kernel = "rbfdot",
                           kpar = "automatic", nystrom.red = TRUE,
                           nystrom.sample = min(nrow(y_train) / 6,
                                                as.integer(spectral_nystrom_sample)))
          } else {
            kernlab::specc(y_train, centers = K, kernel = "rbfdot",
                           kpar = "automatic", nystrom.red = TRUE)
          },
          silent = TRUE
        )
        if (inherits(fit, "try-error")) {
          if (verbose) message("  spectral_K", K, " failed; skipping.")
          next
        }
        cand_list[[paste0("spectral_K", K)]] <- as.integer(fit)
      }
      
    } else {
      warning("Unknown algorithm '", alg, "'; skipping.")
    }
  }
  
  
  if (length(cand_list) == 0L) {
    if (verbose) message("get_ml_candidates: no candidates generated.")
    return(matrix(integer(0L), nrow = 0L, ncol = N))
  }
  
  ## Pool, canonicalize, and drop duplicates -- different algorithms often
  ## agree, and an identical partition must not be double-counted in the
  ## candidate weights.
  parts <- do.call(rbind, cand_list)
  rownames(parts) <- names(cand_list)
  
  if (relabel) {
    parts <- t(apply(parts, 1L, relabel_partition))
    rownames(parts) <- names(cand_list)
  }
  
  keys  <- apply(parts, 1L, paste, collapse = "_")
  parts <- parts[!duplicated(keys), , drop = FALSE]
  
  if (verbose) {
    message("Generated ", nrow(parts), " unique candidate partitions.")
  }
  
  parts
}



## ===========================================================================
##  Candidate-set diagnostics (not part of any method; off in production runs)
## ===========================================================================

## Print the size, block structure and, when the truth is supplied, the ARI of
## the first max_show candidate partitions.
print_candidate_summaries <- function(parts,
                                      M        = NULL,
                                      z_true   = NULL,
                                      max_show = 5) {
  if (!is.matrix(parts)) parts <- as.matrix(parts)
  
  n_cand <- nrow(parts)
  N      <- ncol(parts)
  if (is.null(M)) M <- max(parts)
  
  cat("\n---------------- ML candidate label sets ----------------\n")
  cat("Number of candidates:", n_cand, "\n")
  cat("Number of observations:", N, "\n\n")
  
  if (n_cand == 0L) {
    cat("No candidates to display.\n")
    cat("Head of label matrix (0 x", N, "):\n")
    print(parts)
    cat("----------------------------------------------------------\n")
    return(invisible(NULL))
  }
  
  idx_show <- seq_len(min(max_show, n_cand))
  
  for (i in idx_show) {
    z <- parts[i, ]
    
    name_i <- rownames(parts)[i]
    if (is.null(name_i) || nchar(name_i) == 0) {
      name_i <- paste0("Candidate ", i)
    }
    
    cat(">>> ", name_i, " (row ", i, ")\n", sep = "")
    
    K_eff <- length(unique(z[z > 0]))
    cl_sizes <- tabulate(z, nbins = M)
    singleton_count <- sum(cl_sizes == 1)
    
    cat("  K_eff          :", K_eff, "\n")
    cat("  cluster sizes  :", paste(cl_sizes, collapse = " "), "\n")
    cat("  singletons     :", singleton_count, "\n")
    
    if (!is.null(z_true)) {
      if (requireNamespace("mclust", quietly = TRUE)) {
        ari <- mclust::adjustedRandIndex(z_true, z)
        cat("  ARI w.r.t true :", sprintf("%.3f", ari), "\n")
      } else {
        cat("  ARI w.r.t true : <mclust not installed>\n")
      }
    }
    
    cat("\n")
  }
  
  cat("Head of label matrix (first", min(10, N), "obs, first",
      length(idx_show), "candidates):\n")
  print(parts[idx_show, 1:min(10, N), drop = FALSE])
  cat("----------------------------------------------------------\n")
  
  invisible(NULL)
}


## Scatter y against the second covariate, faceted by the true labels and the
## first max_show candidate partitions, to see what each generator proposes.
plot_truth_and_candidates <- function(y,
                                      x,
                                      parts,
                                      z_true   = NULL,
                                      max_show = 4) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package 'ggplot2' is required for plotting.")
    return(invisible(NULL))
  }
  library(ggplot2)
  
  if (!is.matrix(parts)) parts <- as.matrix(parts)
  
  N <- length(y)
  if (nrow(x) != N) {
    stop("x and y must have the same number of rows/length.")
  }
  
  x <- as.matrix(x)
  if (ncol(x) < 2L) {
    stop("x must have at least 2 columns; x[,2] is used on the x-axis.")
  }
  
  x2 <- x[, 2]
  
  n_cand   <- nrow(parts)
  idx_show <- seq_len(min(max_show, n_cand))
  
  df_list   <- list()
  lab_levels <- character(0)
  
  # Truth
  if (!is.null(z_true)) {
    if (length(z_true) != N) {
      stop("z_true must have length equal to length(y).")
    }
    df_truth <- data.frame(
      x2       = x2,
      y        = y,
      cluster  = factor(z_true),
      labeling = "Truth"
    )
    df_list[[length(df_list) + 1L]] <- df_truth
    lab_levels <- c(lab_levels, "Truth")
  }
  
  # Candidates
  for (i in idx_show) {
    z_i <- parts[i, ]
    name_i <- rownames(parts)[i]
    if (is.null(name_i) || nchar(name_i) == 0) {
      name_i <- paste0("Candidate ", i)
    } else {
      name_i <- paste0("Cand ", i, ": ", name_i)
    }
    
    df_i <- data.frame(
      x2       = x2,
      y        = y,
      cluster  = factor(z_i),
      labeling = name_i
    )
    df_list[[length(df_list) + 1L]] <- df_i
    lab_levels <- c(lab_levels, name_i)
  }
  
  df_long <- do.call(rbind, df_list)
  df_long$labeling <- factor(df_long$labeling, levels = lab_levels)
  
  p <- ggplot(df_long, aes(x = x2, y = y, color = cluster)) +
    geom_point(alpha = 0.7, size = 1.5) +
    facet_wrap(~ labeling, ncol = 2, scales = "free_y") +
    labs(
      x = "Second covariate (x[,2])",
      y = "Response y",
      color = "Cluster",
      title = "True labels vs ML candidate labelings"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  print(p)
  invisible(p)
}


## "MCMC": the unconstrained Gibbs sampler for the model of Equation (1), used
## as the comparison baseline.  Unlike the direct samplers this conditions on
## the current parameter values at every sweep, so successive draws are
## correlated and the chain must be burned in and thinned.
##
## Blocks: (1) pi | Z, (2) beta_k and sigma2_k | Z with beta_k drawn first and
## sigma2_k updated conditional on it, (3) Z | everything else.
##
## @param X  N x p covariate matrix; include an intercept column if wanted.
## @return list(beta [ntotal x M x p], sigma2, pi, z, kpost), where kpost is
##         the number of occupied components per retained draw.
full_mcmc_sampling_regression_dependent <- function(y, X, M=3, alpha=1,
                                                    omega=2, kappa=1, sig2_beta=100,
                                                    niter=1000, nburn=0, nthin=1){
  if(!is.matrix(X)) X <- as.matrix(X)
  N <- length(y)
  p <- ncol(X)
  stopifnot(nrow(X) == N)
  
  ntotal <- floor((niter - nburn) / nthin)
  
  rdirichlet1 <- function(alpha_vec){
    g <- rgamma(length(alpha_vec), shape = alpha_vec, rate = 1)
    as.numeric(g / sum(g))
  }
  rinvgamma_ <- function(n, shape, scale){  
    1 / rgamma(n, shape = shape, rate = scale)
  }
  
  # Initialize
  z    <- sample(1:M, N, replace = TRUE)
  beta <- matrix(0, nrow = M, ncol = p)  
  sig2 <- rep(1, M)
  w    <- rep(1/M, M)
  
  # Storage
  z_keep    <- matrix(NA_integer_, nrow = ntotal, ncol = N)
  beta_keep <- array(NA_real_, dim = c(ntotal, M, p))
  sig2_keep <- matrix(NA_real_, nrow = ntotal, ncol = M)
  pi_keep   <- matrix(NA_real_, nrow = ntotal, ncol = M)
  k_keep    <- matrix(NA_integer_, nrow = ntotal, ncol = 1)
  
  I_p <- diag(1, p)
  prior_prec <- (1 / sig2_beta) * I_p
  
  ii <- 1
  for(i in 1:niter){
    
    ## (1) Update weights w | z
    nk <- tabulate(z, nbins = M)
    w  <- rdirichlet1(alpha + nk)
    
    ## (2) Update component-specific parameters (beta_k, sig2_k)
    for(k in 1:M){
      idx_k <- which(z == k)
      n_k   <- length(idx_k)
      
      if(n_k > 0){
        Xk <- X[idx_k, , drop = FALSE]
        yk <- y[idx_k]
        
        # Posterior precision for beta_k (scaled by 1/sig2_k):
        # A_k = Xk'Xk + (1/sig2_beta) I_p
        XtX  <- crossprod(Xk)
        Xty  <- crossprod(Xk, yk)
      } else {
        # Empty cluster: X'X = 0, X'y = 0 
        XtX <- matrix(0, p, p)
        Xty <- rep(0, p)
      }
      
      A_k <- XtX + prior_prec  # p x p; free of sig2_k, so factor it once

      R <- chol(A_k)           # A_k = R'R
      mu_k <- backsolve(R, forwardsolve(t(R), Xty))   # A_k^{-1} X'y

      zstd <- rnorm(p)
      beta_k <- as.numeric(mu_k + sqrt(sig2[k]) * backsolve(R, zstd))
      beta[k, ] <- beta_k
      
      # Residual sum of squares for y | beta
      if(n_k > 0){
        res <- yk - as.numeric(Xk %*% beta_k)
        rss <- sum(res^2)
      } else {
        rss <- 0
      }
      ## The scale-dependent prior contributes beta'beta / sig2_beta to the
      ## inverse-gamma scale, and p to its shape.
      quad_prior <- sum(beta_k^2) / sig2_beta


      a_star <- omega + 0.5 * (n_k + p)
      b_star <- kappa + 0.5 * (rss + quad_prior)
      sig2[k] <- rinvgamma_(1, shape = a_star, scale = b_star)
    }
    
    ## (3) Update allocations z_i, with
    ##     log p_k = log w_k + log N(y_i | x_i' beta_k, sig2_k) + const.
    logpk <- matrix(NA_real_, nrow = M, ncol = N)
    for(k in 1:M){
      mu_i <- as.numeric(X %*% beta[k, ])
      logpk[k, ] <- dnorm(y, mean = mu_i, sd = sqrt(sig2[k]), log = TRUE) + log(w[k])
    }
    logpk <- sweep(logpk, 2, apply(logpk, 2, max), FUN = "-")   # stabilize
    pk    <- exp(logpk)
    pk    <- sweep(pk, 2, colSums(pk), FUN = "/")
    z <- apply(pk, 2, function(probs) sample.int(M, 1, prob = probs))

    krep <- sum(tabulate(z, nbins = M) > 0)   # occupied components
    
    ## (4) Save draws
    if((i > nburn) && (((i - nburn) %% nthin) == 0)){
      z_keep[ii, ]       <- z
      beta_keep[ii, , ]  <- beta
      sig2_keep[ii, ]    <- sig2
      pi_keep[ii, ]      <- w
      k_keep[ii, 1]      <- krep
      ii <- ii + 1
      if(ii > ntotal) break
    }
  }
  
  list(beta = beta_keep, sigma2 = sig2_keep, pi = pi_keep,
       z = z_keep, kpost = k_keep)
}

## Collapsed label chain: a Gibbs sampler for Equation (8), the conditional of
## Z_i given the other labels with pi, beta and sigma2 all integrated out.  This
## is the label half of "MC-MCMC"; the parameters are added afterwards by
## composition in direct_sampling_gibbs_regression().
##
## Naively, each of the N x M candidate reassignments per sweep would need its
## own p x p factorization.  Instead the per-cluster sufficient statistics
## X'X, X'y, y'y and the derived quantities A^{-1}, log|A| and
## Q = y'y - s_x'A^{-1}s_x are cached and updated by rank-one Sherman-Morrison
## corrections as an observation leaves or joins a cluster, which makes a sweep
## O(N M p^2) with no factorizations inside the inner loop.
##
## @return list(z [ntotal x N], kpost), labels only.
gibbs_labels_only_regression <- function(
    y, X, M = 25, alpha = 1, omega = 2, kappa = 1, sig2_beta = 100,
    niter = 2500, nburn = 500, nthin = 2, seed = NULL
){
  if (!is.null(seed)) set.seed(seed)
  if (!is.matrix(X)) X <- as.matrix(X)
  N <- length(y); p <- ncol(X); stopifnot(nrow(X) == N)
  ntotal <- floor((niter - nburn) / nthin); stopifnot(ntotal > 0L)
  if (length(alpha) == 1L) alpha <- rep(alpha, M)
  stopifnot(length(alpha) == M)
  
  ## Prior precision Lambda0 = I_p / sigma2_beta, and the quantities an empty
  ## cluster resets to.
  Lambda0 <- diag(1 / sig2_beta, p)
  invLambda0 <- diag(sig2_beta, p)
  logdet_Lambda0 <- p * log(1 / sig2_beta)  # = -p*log(sig2_beta)
  
  ## Sufficient statistics per cluster.
  nj <- integer(M)
  Sx <- array(0, dim = c(p, p, M))  # X'X
  sx <- matrix(0, nrow = p, ncol = M)  # X'y
  s2 <- numeric(M)                     # y'y
  
  ## Derived quantities maintained incrementally rather than recomputed.
  Ainv <- array(0, c(p, p, M))         # (X'X + Lambda0)^{-1}
  logdet_A <- rep(logdet_Lambda0, M)   # log|A|
  Ainv_sx <- matrix(0, p, M)           # A^{-1} s_x
  Qval <- numeric(M)                   # s2 - s_x' A^{-1} s_x
  
  ## Random initial labels, then build every cache from scratch once.
  z <- sample.int(M, N, replace = TRUE)
  for (j in 1:M) {
    idx <- (z == j)
    nj[j] <- sum(idx)
    if (nj[j] > 0L) {
      Xj <- X[idx, , drop = FALSE]
      yj <- y[idx]
      Sx[,,j] <- crossprod(Xj)
      sx[,j]  <- crossprod(Xj, yj)
      s2[j]   <- sum(yj * yj)
      Aj <- Sx[,,j] + Lambda0
      Rj <- chol(Aj)
      logdet_A[j] <- 2 * sum(log(diag(Rj)))
      Ainv[,,j] <- chol2inv(Rj)
    } else {
      Ainv[,,j] <- invLambda0
      logdet_A[j] <- logdet_Lambda0
    }
    Ainv_sx[,j] <- Ainv[,,j] %*% sx[,j]
    Qval[j] <- s2[j] - sum(sx[,j] * Ainv_sx[,j])
  }
  
  ## lgamma(omega + n/2) is needed for every candidate cluster size, so tabulate
  ## it over the only N + 1 values it can take.
  n_grid <- 0:N
  lg_w_n <- lgamma(omega + 0.5 * n_grid)
  lg_w <- function(n) lg_w_n[n + 1L]
  
  # Storage
  z_keep <- matrix(NA_integer_, nrow = ntotal, ncol = N)
  k_keep <- integer(ntotal)
  ii <- 1L
  Kcur <- sum(nj > 0L)
  
  ## Rank-one add.  The caller has already computed v = A^{-1}x, c = 1 + x'v
  ## and the updated t and Q while scoring the candidates, so committing the
  ## move costs no further solves.  Uses <<- to update the enclosing caches.
  add_update <- function(j, xi, yi, xi_yi, yi2, v_j, c_plus, t_new, Q_new) {
    nj[j] <<- nj[j] + 1L
    Sx[,,j] <<- Sx[,,j] + tcrossprod(xi)
    sx[,j]  <<- sx[,j]  + xi_yi
    s2[j]   <<- s2[j]   + yi2
    
    ## (A + xx')^{-1} = A^{-1} - (A^{-1}x)(A^{-1}x)' / (1 + x'A^{-1}x), and
    ## log|A + xx'| = log|A| + log(1 + x'A^{-1}x).
    Ainv_old <- Ainv[,,j]
    Ainv_j <- Ainv_old - (v_j %o% v_j) / c_plus
    Ainv[,,j] <<- Ainv_j
    
    logdet_A[j] <<- logdet_A[j] + log(c_plus)
    
    # cached from weight calc
    Ainv_sx[,j] <<- t_new
    Qval[j]     <<- Q_new
  }
  
  
  ## Rank-one remove, the inverse operation.  v and c_minus are computed from
  ## the cache with the observation still in the cluster.
  remove_update <- function(j, xi, yi, xi_yi, yi2) {
    v <- as.vector(Ainv[,,j] %*% xi)
    c_minus <- 1 - sum(xi * v)
    if (!is.finite(c_minus) || c_minus <= 1e-12) c_minus <- 1e-12
    
    nj[j] <<- nj[j] - 1L
    Sx[,,j] <<- Sx[,,j] - tcrossprod(xi)
    sx[,j]  <<- sx[,j]  - xi_yi
    s2[j]   <<- s2[j]   - yi2
    
    if (nj[j] == 0L) {
      ## An emptied cluster is reset to the prior exactly, rather than left as
      ## the accumulated result of many rank-one updates.
      Ainv[,,j] <<- invLambda0
      logdet_A[j] <<- logdet_Lambda0
      Ainv_sx[,j] <<- rep(0, p)
      Qval[j]     <<- 0
    } else {
      # A^{-1} <- A^{-1} + (v v') / (1 - x' A^{-1} x)
      Ainv_j <- Ainv[,,j] + (v %o% v) / c_minus
      Ainv[,,j] <<- Ainv_j
      logdet_A[j] <<- logdet_A[j] + log(c_minus)
      t_new <- Ainv_j %*% sx[,j]
      Ainv_sx[,j] <<- as.vector(t_new)
      Qval[j]     <<- s2[j] - sum(sx[,j] * t_new)
    }
  }
  
  ## Sweeps.  Observations are visited in random order each sweep.
  for (iter in 1:niter) {
    for (i in sample.int(N)) {
      j_old <- z[i]
      xi <- X[i, ]; yi <- y[i]; yi2 <- yi * yi
      xi_yi <- xi * yi
      
      ## Take observation i out of its current cluster.
      prev_n <- nj[j_old]
      remove_update(j_old, xi, yi, xi_yi, yi2)
      if (prev_n == 1L) Kcur <- Kcur - 1L
      
      ## Score every cluster: what log|A| and Q would become if i joined it.
      log_prior <- log(nj + alpha)
      n <- nj
      logdet_bef <- logdet_A
      Q_bef <- Qval
      
      logdet_aft <- numeric(M)
      Q_aft <- numeric(M)
      v_cache <- vector("list", M)
      t_cache <- vector("list", M)
      
      for (j in 1:M) {
        vj <- as.vector(Ainv[,,j] %*% xi)        # A^{-1} x
        v_cache[[j]] <- vj
        c <- 1 + sum(xi * vj)                    # 1 + x' A^{-1} x
        if (c <= 1e-12) c <- 1e-12
        logdet_aft[j] <- logdet_bef[j] + log(c)
        
        t1 <- Ainv_sx[,j] + yi * vj              # A^{-1} s_x + y A^{-1} x
        udot_t1 <- sum(xi * t1)
        tj <- t1 - vj * (udot_t1 / c)            # (A+xx')^{-1} (s_x + x y)
        t_cache[[j]] <- tj
        
        s_vec <- sx[,j] + xi_yi
        Q_aft[j] <- (s2[j] + yi2) - sum(s_vec * tj)
      }
      
      ## Equation (8): the ratio of collapsed marginals with and without
      ## observation i, times the Dirichlet-multinomial prior n_j + alpha.
      n_after <- n + 1L
      dm <- (lg_w(n_after) - lg_w(n)) +
        (omega + 0.5 * n)      * log(kappa + 0.5 * Q_bef) -
        (omega + 0.5 * n_after)* log(kappa + 0.5 * Q_aft) -
        0.5 * (logdet_aft - logdet_bef)
      logw <- log_prior + dm
      
      ## Sample the new label by inverse CDF on the stabilized weights.
      mmax <- max(logw)
      pj <- exp(logw - mmax)
      s <- sum(pj)
      j_new <- if (!is.finite(s) || s <= 0) which.max(logw) else {
        u <- runif(1) * s
        cs <- 0.0; pick <- 1L
        for (j in 1:M) { cs <- cs + pj[j]; if (u <= cs) { pick <- j; break } }
        pick
      }
      
      ## Commit, reusing the cached v, t and Q computed for j_new above.
      v_new <- v_cache[[j_new]]
      c_plus <- 1 + sum(xi * v_new); if (c_plus <= 1e-12) c_plus <- 1e-12
      t_new <- t_cache[[j_new]]
      Q_new <- Q_aft[j_new]
      add_update(j_new, xi, yi, xi_yi, yi2, v_new, c_plus, t_new, Q_new)
      if (nj[j_new] == 1L) Kcur <- Kcur + 1L
      z[i] <- j_new
    }
    
    if ((iter > nburn) && (((iter - nburn) %% nthin) == 0)) {
      z_keep[ii, ] <- z
      k_keep[ii] <- Kcur
      ii <- ii + 1L
      if (ii > ntotal) break
    }
  }
  
  list(z = z_keep, kpost = k_keep)
}


## "MC-MCMC": the method of composition applied to a collapsed label chain.
## Labels come from gibbs_labels_only_regression() (Equation 8); given each
## retained label draw, pi, sigma2 and beta are drawn exactly from Equations
## (5), (4) and (3), so only the label sequence is autocorrelated.
##
## @return list(beta [nsamp x M x p], sigma2, pi, Z, K, fitted.values).
direct_sampling_gibbs_regression <- function(
    y, X, M = 3, alpha = 1, sig2_beta = 100, omega = 2, kappa = 1, nsamp = 1000,
    niter = 3000, nburn = 1000, nthin = 2
){
  if (!is.matrix(X)) X <- as.matrix(X)
  N <- length(y); p <- ncol(X)
  stopifnot(nrow(X) == N)
  
  ## The label chain must supply at least nsamp draws for the composition step
  ## below; checked here so a mismatch is reported rather than surfacing as a
  ## subscript error the driver would record as a failed replicate.
  ntotal <- floor((niter - nburn) / nthin)
  if (ntotal < nsamp)
    stop("direct_sampling_gibbs_regression: label chain yields ", ntotal,
         " draws but nsamp = ", nsamp,
         ".  Increase niter or reduce nsamp/nthin.")

  gibbs_fit <- gibbs_labels_only_regression(
    y, X,
    M = M, alpha = alpha, omega = omega, kappa = kappa, sig2_beta = sig2_beta,
    niter = niter, nburn = nburn, nthin = nthin, seed = NULL
  )

  ## storage
  betas <- array(NA_real_, dim = c(nsamp, M, p))
  sigma2s <- matrix(NA_real_, nrow = nsamp, ncol = M)
  pis <- matrix(NA_real_, nrow = nsamp, ncol = M)
  Zs <- matrix(NA_integer_, nrow = nsamp, ncol = N)
  Ys <- matrix(NA_real_, nrow = nsamp, ncol = N)
  K <- integer(nsamp)
  
  rdirichlet1 <- function(alpha_vec) {
    g <- rgamma(length(alpha_vec), shape = alpha_vec, rate = 1)
    as.numeric(g / sum(g))
  }
  
  for (i in 1:nsamp){
    Z <- gibbs_fit$z[i, ]
    Zs[i, ] <- Z
    K[i] <- length(unique(Z))
    
    ## Sufficient statistics for this allocation.
    nk <- tabulate(Z, nbins = M)
    Sx <- array(0, dim = c(p, p, M))
    sx <- matrix(0, nrow = p, ncol = M)
    s2 <- numeric(M)
    
    for (k in 1:M){
      idx <- which(Z == k)
      if (length(idx) > 0L){
        Xk <- X[idx, , drop = FALSE]
        yk <- y[idx]
        Sx[,,k] <- crossprod(Xk)
        sx[,k] <- crossprod(Xk, yk)
        s2[k] <- sum(yk * yk)
      }
    }
    
    ## pi | Z ~ Dirichlet(alpha + n), Equation (5).
    pis[i, ] <- rdirichlet1(alpha + nk)

    ## sigma2_k then beta_k, Equations (4) and (3).  Empty components have
    ## A = I_p / sigma2_beta and s_x = 0, so both revert to their priors.
    for (k in 1:M){
      A <- Sx[,,k] + diag(1 / sig2_beta, p)
      R <- chol(A)
      Ainv_sx <- backsolve(R, forwardsolve(t(R), sx[,k]))
      mu_beta <- Ainv_sx
      
      a_star <- omega + 0.5 * nk[k]
      b_star <- kappa + 0.5 * (s2[k] - sum(sx[,k] * Ainv_sx))
      sigma2s[i, k] <- 1 / rgamma(1, shape = a_star, rate = b_star)
      
      zstd <- rnorm(p)
      beta_k <- as.numeric(mu_beta + sqrt(sigma2s[i, k]) * backsolve(R, zstd))
      betas[i, k, ] <- beta_k
    }
    
    ## Fitted values under this draw's own allocation.
    yh <- numeric(N)
    for (ii in 1:N) {
      yh[ii] <- sum(X[ii, ] * betas[i, Z[ii], ])
    }
    Ys[i, ] <- yh
  }
  
  list(beta = betas, sigma2 = sigma2s, pi = pis, Z = Zs, K = K, fitted.values = Ys)
}
