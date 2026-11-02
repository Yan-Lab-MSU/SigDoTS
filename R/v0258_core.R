
# SigDoTS core functions: the single-domain null model of a test (from the research code
# v0.2.5.8).
# ---------------------------------------------------------------------------
# At a node of the tree, the null model for the node's n spots x p genes is matrix normal,
#   X ~ MN(1 mu', Sigma_s, Sigma_g),   Sigma_s = eta * K(ell) + (1 - eta) * I  (+ nugget),
# with K the Matern 3/2 correlation of the distances between spots.
#
# Spatial parameters (estimate_node_spatial_params()). eta and ell are estimated from pairs of
# spots in the same leaf (candidate domain) of the subtree under the node. The correlogram of
# these pairs, with each leaf's mean removed and pooled over genes, is fitted by weighted
# least squares to its exact expectation under the model for leaf-centred data
# (centred_correlogram_expectation()). ell is searched on a log grid, in units of the median
# nearest-neighbour distance, and refined continuously between the best grid value's two
# neighbours; eta is optimised at every ell. Pairs that straddle two leaves are never used,
# because they carry the between-leaf mean differences that the test is looking for.
#
# Gene covariance (whitened_rows_without_group_means(), spiked_cov_est_from_whitened()). The
# data are whitened with Sigma_s^-1/2, the mean of every leaf is removed by generalised least
# squares (only the n - K whitened rows orthogonal to the whitened leaf indicators are kept),
# and Sigma_g is estimated from them as a spiked covariance: a few leading eigenvectors plus
# equal variance in all other directions.
#
# Spatial backends (build_exact_spatial_backend(), build_vecchia_spatial_backend()). The exact
# backend eigendecomposes K once (eta * K + (1 - eta) * I has K's eigenvectors, with
# eigenvalues eta * lambda + (1 - eta)); the Vecchia backend uses GpGp's (variance, range,
# nugget) parameterisation with variance = eta and nugget = ((1 - eta) + nugget) / eta, which
# gives the same matrix.
# ---------------------------------------------------------------------------

# -------------------------------
# Utility: Cluster Index (CI)
# -------------------------------
cluster_index <- function(Z, labels) {
  # Z: n x p feature matrix on which clustering was performed (e.g., PCA dims)
  # labels: integer or factor vector length n
  stopifnot(nrow(Z) == length(labels))
  labels <- as.factor(labels)
  n <- nrow(Z)
  mu <- colMeans(Z)
  tss <- sum(rowSums( (Z - matrix(mu, n, ncol(Z), byrow = TRUE))^2 ))
  wss <- 0
  for (lev in levels(labels)) {
    idx <- which(labels == lev)
    if (length(idx) > 0) {
      mu_k <- colMeans(Z[idx, , drop = FALSE])
      wss <- wss + sum(rowSums( (Z[idx, , drop = FALSE] - matrix(mu_k, length(idx), ncol(Z), byrow = TRUE))^2 ))
    }
  }
  ci <- wss / tss
  return(ci)
}

matern32_kernel <- function(coords, ell) {
  # Matérn 3/2 with unit marginal variance
  D <- as.matrix(dist(coords))
  sqrt3 <- sqrt(3)
  R <- (1 + sqrt3 * D / ell) * exp(-sqrt3 * D / ell)
  diag(R) <- 1
  return(R)
}

# ---- Spatial covariance backends -------------------------------------------
# v0.2.5: exact dense backend when feasible; Vecchia backend otherwise.

estimate_exact_peak_gib <- function(n, dense_matrix_copies = 5) {
  n <- as.double(n)
  bytes_dist <- 8 * n * (n - 1) / 2
  bytes_dense <- as.double(dense_matrix_copies) * 8 * n * n
  (bytes_dist + bytes_dense) / 1024^3
}

is_known_dense_failure_warning <- function(message) {
  grepl(
    paste(
      c(
        "integer overflow",
        "cannot allocate",
        "allocation",
        "size \\* size",
        "must be a finite number",
        "long vectors not supported"
      ),
      collapse = "|"
    ),
    message,
    ignore.case = TRUE
  )
}

new_spatial_backend_failure <- function(message) {
  structure(
    list(message = message, call = NULL),
    class = c("spatial_backend_failure", "error", "condition")
  )
}

build_exact_spatial_backend <- function(coords, ell, nugget = 0,
                                        eta = 1) {
  # Package: Matern 3/2 only (v0.2.5.8 also offered an rbf kernel, never used).
  coords <- as.matrix(coords)
  kernel <- "matern32"
  stopifnot(
    ncol(coords) >= 2L,
    all(is.finite(coords)),
    is.finite(ell),
    ell > 0,
    is.finite(nugget),
    nugget >= 0,
    is.finite(eta),
    eta > 0,
    eta <= 1
  )

  K <- matern32_kernel(coords, ell)
  # Sigma_s = eta * K + (1 - eta) * I, the same form the generator uses. The kernel
  # has a unit diagonal, so scaling by eta and adding (1 - eta) back to the diagonal
  # leaves the diagonal at exactly 1 and multiplies every off-diagonal by eta.
  # eta = 1 leaves K untouched, which is why the default reproduces v0.2.5.2.
  if (eta != 1) {
    K <- eta * K
    diag(K) <- diag(K) + (1 - eta)
  }
  diag(K) <- diag(K) + nugget
  ee <- eigen(K, symmetric = TRUE)
  tol <- max(dim(K)) * .Machine$double.eps * max(abs(ee$values), 1)
  if (any(!is.finite(ee$values)) || min(ee$values) < -tol) {
    stop("Exact spatial covariance is not numerically positive semidefinite.")
  }
  lambda <- pmax(ee$values, tol)
  structure(
    list(
      type = "dense_exact",
      n = nrow(coords),
      ell = ell,
      eta = eta,
      nugget = nugget,
      kernel = kernel,
      U = ee$vectors,
      lambda = lambda,
      selection = list(
        requested = "exact",
        selected = "dense_exact",
        fallback = FALSE,
        reason = NULL,
        estimated_exact_peak_gib = estimate_exact_peak_gib(nrow(coords))
      )
    ),
    class = "spatial_backend"
  )
}

build_vecchia_spatial_backend <- function(coords, ell, nugget = 0,
                                          n_neighbors = 60L,
                                          eta = 1,
                                          order_seed = NULL) {
  # The session's thread settings are left as they are. GpGp is an import of the package.
  coords <- as.matrix(coords)
  kernel <- "matern32"
  stopifnot(
    ncol(coords) >= 2L,
    all(is.finite(coords)),
    is.finite(ell),
    ell > 0,
    is.finite(nugget),
    nugget >= 0,
    is.finite(eta),
    eta > 0,
    eta <= 1,
    length(n_neighbors) == 1L,
    is.finite(n_neighbors),
    n_neighbors >= 1L
  )

  n <- nrow(coords)
  if (n < 2L) stop("Vecchia backend requires at least two spatial locations.")
  n_neighbors <- min(as.integer(n_neighbors), n - 1L)

  # GpGp::order_maxmin() and GpGp::find_ordered_nn() both draw random numbers (each jitters
  # the coordinates with stats::rnorm(); order_maxmin also samples). In v0.2.5.8 they drew
  # them from the global stream, so every test after the first depended on how earlier
  # replicates had been batched over cores. The package runs both under a seed of the test's
  # own (`order_seed`) and leaves the caller's stream untouched. order_seed = NULL keeps
  # v0.2.5.8's behaviour; only the comparisons with v0.2.5.8 use it.
  vecchia_order <- function() {
    ord <- GpGp::order_maxmin(coords)
    coords_ord <- coords[ord, , drop = FALSE]
    list(ord = ord, coords_ord = coords_ord,
         NNarray = GpGp::find_ordered_nn(coords_ord, m = n_neighbors))
  }
  vo <- if (is.null(order_seed)) vecchia_order() else .with_local_seed(order_seed, vecchia_order())
  ord <- vo$ord
  inv_ord <- order(ord)
  coords_ord <- vo$coords_ord
  NNarray <- vo$NNarray
  # GpGp's matern15_isotropic is variance * Matern(range) with variance * nugget
  # added to the diagonal, so variance = eta and nugget = ((1 - eta) + nugget) / eta
  # gives off-diagonal eta * K and diagonal eta + (1 - eta) + nugget = 1 + nugget,
  # matching the exact backend above (verified to 2.8e-17). At eta = 1 this is variance 1
  # with the nugget unchanged.
  covparms <- c(variance = eta, range = ell / sqrt(3),
                nugget = ((1 - eta) + nugget) / eta)
  Linv <- GpGp::vecchia_Linv(
    covparms = covparms,
    covfun_name = "matern15_isotropic",
    locs = coords_ord,
    NNarray = NNarray
  )
  if (any(!is.finite(Linv))) {
    stop("GpGp returned non-finite inverse-Cholesky entries.")
  }

  structure(
    list(
      type = "vecchia",
      n = n,
      ell = ell,
      eta = eta,
      nugget = nugget,
      kernel = kernel,
      n_neighbors = n_neighbors,
      order = ord,
      inverse_order = inv_ord,
      coords_ordered = coords_ord,
      NNarray = NNarray,
      Linv = Linv,
      covparms = covparms,
      covfun_name = "matern15_isotropic",
      GpGp_version = as.character(utils::packageVersion("GpGp")),
      selection = list(
        requested = "vecchia",
        selected = "vecchia",
        fallback = FALSE,
        reason = NULL,
        estimated_exact_peak_gib = estimate_exact_peak_gib(n)
      )
    ),
    class = "spatial_backend"
  )
}

build_spatial_backend_auto <- function(coords, ell, nugget = 0,
                                       max_exact_peak_gib = 16,
                                       vecchia_neighbors = 60L,
                                       force_method = c("vecchia", "exact"),
                                       eta = 1,
                                       order_seed = NULL) {
  # The caller chooses the backend (sigdots(spatial_method = )); there is no automatic choice
  # and no fallback from one backend to the other. The exact backend estimates its memory
  # first and stops when the estimate exceeds the limit.
  coords <- as.matrix(coords)
  force_method <- match.arg(force_method)
  n <- nrow(coords)
  estimated_gib <- estimate_exact_peak_gib(n)

  if (force_method == "exact") {
    preflight_reasons <- character(0)
    # Only the memory limit applies, which the caller sets; v0.2.5.8's additional limit of
    # floor(sqrt(.Machine$integer.max)) = 46,340 spots is not used.
    if (!is.finite(estimated_gib) || estimated_gib > max_exact_peak_gib) {
      preflight_reasons <- c(
        preflight_reasons,
        sprintf("estimated exact peak memory %.2f GiB exceeds limit %.2f GiB",
                estimated_gib, max_exact_peak_gib)
      )
    }
    if (length(preflight_reasons) > 0L) {
      # The message names the options of sigdots().
      stop(sprintf(paste0(
        "spatial_method = \"exact\" would need about %.2f GiB of memory at a test on %d spots, more than ",
        "max_exact_memory_gib = %g. Use spatial_method = \"vecchia\", or, if the machine has the memory, ",
        "raise it with control = sigdots_control(max_exact_memory_gib = ...)."),
        estimated_gib, n, max_exact_peak_gib), call. = FALSE)
    }
    exact_result <- tryCatch(
      withCallingHandlers(
        build_exact_spatial_backend(
          coords = coords,
          ell = ell,
          nugget = nugget,
          eta = eta
        ),
        warning = function(w) {
          msg <- conditionMessage(w)
          if (is_known_dense_failure_warning(msg)) {
            stop(new_spatial_backend_failure(msg))
          }
        }
      ),
      error = function(e) e
    )
    if (inherits(exact_result, "error")) {
      stop("Exact backend failed and fallback was disabled: ",
           conditionMessage(exact_result))
    }
    exact_result$selection <- list(
      requested = force_method,
      selected = "dense_exact",
      fallback = FALSE,
      reason = NULL,
      estimated_exact_peak_gib = estimated_gib
    )
    return(exact_result)
  }

  vecchia_result <- build_vecchia_spatial_backend(
    coords = coords,
    ell = ell,
    nugget = nugget,
    n_neighbors = vecchia_neighbors,
    eta = eta,
    order_seed = order_seed
  )
  vecchia_result$selection <- list(
    requested = force_method,
    selected = "vecchia",
    fallback = FALSE,
    reason = NULL,
    estimated_exact_peak_gib = estimated_gib
  )
  vecchia_result
}

apply_by_column <- function(X, FUN) {
  X <- as.matrix(X)
  out <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X),
                dimnames = dimnames(X))
  for (j in seq_len(ncol(X))) {
    out[, j] <- FUN(X[, j])
  }
  out
}

spatial_whiten_dense <- function(backend, X) {
  stopifnot(nrow(X) == backend$n)
  UtX <- crossprod(backend$U, X)
  scaled <- sweep(UtX, MARGIN = 1L, STATS = sqrt(backend$lambda), FUN = "/")
  backend$U %*% scaled
}

spatial_correlate_dense <- function(backend, E) {
  stopifnot(nrow(E) == backend$n)
  scaled <- sweep(crossprod(backend$U, E),
                  MARGIN = 1L, STATS = sqrt(backend$lambda), FUN = "*")
  backend$U %*% scaled
}

spatial_whiten_vecchia <- function(backend, X) {
  stopifnot(nrow(X) == backend$n)
  X_ord <- X[backend$order, , drop = FALSE]
  Y_ord <- apply_by_column(
    X_ord,
    FUN = function(x) {
      GpGp::Linv_mult(Linv = backend$Linv, z = x, NNarray = backend$NNarray)
    }
  )
  Y_ord[backend$inverse_order, , drop = FALSE]
}

spatial_correlate_vecchia <- function(backend, E) {
  stopifnot(nrow(E) == backend$n)
  E_ord <- E[backend$order, , drop = FALSE]
  Z_ord <- apply_by_column(
    E_ord,
    FUN = function(z) {
      GpGp::fast_Gp_sim_Linv(
        Linv = backend$Linv,
        NNarray = backend$NNarray,
        z = z
      )
    }
  )
  Z_ord[backend$inverse_order, , drop = FALSE]
}

spatial_whiten <- function(backend, X) {
  switch(
    backend$type,
    dense_exact = spatial_whiten_dense(backend, X),
    vecchia = spatial_whiten_vecchia(backend, X),
    stop("Unknown spatial backend type: ", backend$type)
  )
}

spatial_correlate <- function(backend, E) {
  switch(
    backend$type,
    dense_exact = spatial_correlate_dense(backend, E),
    vecchia = spatial_correlate_vecchia(backend, E),
    stop("Unknown spatial backend type: ", backend$type)
  )
}

# Gene-side (feature) covariance from the whitened node data: SigClust-style spiked model.
# (v0.2.5.8 kept two older, commented-out estimators here; removed.)
spiked_cov_est_from_whitened <- function(
  Xw,
  mu,
  r_cap = 80,
  ridge = 1e-6,
  max_iter = 50,
  tol = 1e-6,
  eps = 1e-12,
  return_whitened = FALSE
) {
  Xw <- as.matrix(Xw)
  n <- nrow(Xw); p <- ncol(Xw)
  if (n < 2 || p < 1) stop("Bad whitened matrix dimensions.")
  if (length(mu) != p) stop("mu length must match ncol(Xw).")

  sv <- svd(Xw, nu = 0)
  lam_feat <- (sv$d^2) / n
  trS <- sum(lam_feat)
  gamma <- p / n
  sigma2 <- max(trS / p, eps)
  m <- 0L

  for (it in seq_len(max_iter)) {
    mp_edge <- sigma2 * (1 + sqrt(gamma))^2
    m_new <- sum(lam_feat >= mp_edge)
    m_new <- min(m_new, r_cap, length(lam_feat))
    sum_spikes <- if (m_new > 0) sum(lam_feat[seq_len(m_new)]) else 0
    denom <- p - m_new - (p * m_new) / n
    if (!is.finite(denom) || denom <= eps) denom <- max(p - m_new, 1)
    sigma2_new <- max((trS - sum_spikes) / denom, eps)
    if (m_new == m && abs(sigma2_new - sigma2) / max(sigma2, eps) < tol) {
      m <- m_new; sigma2 <- sigma2_new
      break
    }
    m <- m_new; sigma2 <- sigma2_new
  }

  if (m == 0L) {
    out <- list(
      U_r = matrix(0, p, 0),
      d_r = numeric(0),
      sigma2 = sigma2,
      r = 0L,
      mu = mu,
      ridge = ridge,
      gamma = gamma
    )
    if (return_whitened) out$Xw <- Xw
    return(out)
  }

  # Preserve SigClust-DE's debiased total eigenvalue, then convert to excess variance.
  lambda_total <- pmax(
    lam_feat[seq_len(m)] - sigma2 * (p / n),
    0
  )
  d_r <- pmax(lambda_total - sigma2, 0)
  U_r <- sv$v[, seq_len(m), drop = FALSE]
  out <- list(
    U_r = U_r,
    d_r = d_r,
    sigma2 = sigma2,
    r = m,
    mu = mu,
    ridge = ridge,
    gamma = gamma
  )
  if (return_whitened) out$Xw <- Xw
  out
}


# v0.2.5.8: the whitened rows the gene side is estimated from when
# `feature_centring = "leaves_gls"`. `X_whitened` is Sigma_s^-1/2 (X - 1 xbar') for the
# node, `group` the leaf of each of its spots. With B the leaf indicators and
# W = Sigma_s^-1/2 B, the rows returned are Q_perp' X_whitened, Q_perp an orthonormal basis
# of the complement of span(W), computed through the Householder reflectors of qr(W)
# (qr.qty()) without forming Q. The indicators sum to one, so the intercept goes too.
whitened_rows_without_group_means <- function(spatial_backend, X_whitened, group) {
  X_whitened <- as.matrix(X_whitened)
  n <- nrow(X_whitened)
  if (length(group) != n || anyNA(group))
    stop("leaves_gls: `group` must give a leaf for every spot of the node.")
  g <- as.integer(factor(group))
  K <- max(g)
  if (K >= n) stop("leaves_gls: as many groups as spots.")
  B <- matrix(0, n, K)
  B[cbind(seq_len(n), g)] <- 1
  Wb <- spatial_whiten(spatial_backend, B)
  qw <- qr(Wb)
  if (qw$rank < K) stop("leaves_gls: the whitened leaf indicators are rank deficient.")
  rows <- qr.qty(qw, X_whitened)[(K + 1L):n, , drop = FALSE]
  list(rows = rows, n_groups = K, n_rows = nrow(rows),
    removed_energy_share = 1 - sum(rows^2) / sum(X_whitened^2),
    group_sizes = tabulate(g))
}

simulate_matrix_normal_spiked_backend <- function(n, p, mu, spatial_backend,
                                                  U_r, d_r, sigma2) {
  if (n != spatial_backend$n) stop("n does not match spatial backend.")
  W0 <- matrix(rnorm(n * p), n, p)
  X0 <- spatial_correlate(spatial_backend, W0) * sqrt(sigma2)

  if (length(d_r) > 0) {
    Wr <- matrix(rnorm(n * length(d_r)), n, length(d_r))
    spatial_spikes <- spatial_correlate(spatial_backend, Wr)
    spatial_spikes <- sweep(
      spatial_spikes,
      MARGIN = 2L,
      STATS = sqrt(d_r),
      FUN = "*"
    )
    Xr <- spatial_spikes %*% t(U_r)
  } else {
    Xr <- matrix(0, n, p)
  }
  M <- matrix(mu, n, p, byrow = TRUE)
  M + X0 + Xr
}


# ===========================================================================
# v0.2.5.4: estimating (ell, eta) from the node's own data, using only pairs
# that lie inside the same candidate child.
# ===========================================================================

# Matern 3/2 correlation with a spatial fraction. Distances in the same units as
# `ell`. r(0) = 1 by definition; for d > 0 the value is eta * k(d), so eta < 1
# leaves a jump of height 1 - eta at the origin. That jump is the nugget, and it is
# the one feature no choice of `ell` can reproduce.
matern32_corr <- function(d, ell, eta = 1) {
  a <- sqrt(3) * d / ell
  eta * (1 + a) * exp(-a)
}

# Run `expr` with its own RNG stream and put the caller's stream back afterwards.
# v0.2.5.2's select_lengthscale_by_morans_i() sampled features straight off the
# global stream, which silently shifted every draw that followed it inside a node.
# Nothing here is allowed to do that.
.with_local_seed <- function(seed, expr) {
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }
  set.seed(seed)
  force(expr)
}

# Distance bins and a sample of pairs in each, drawn ONLY from pairs whose two
# spots share a group label. `group` is an integer vector over the node's spots,
# normally the candidate split into two children; pairs that cross it are the ones
# carrying the between-child mean difference and are never enumerated.
#
# Bin edges are multiples of the median nearest-neighbour distance, so the binning
# is scale free and does not assume a lattice. Weighting is by the number of pairs
# a bin actually contains: inside a child there is no between-domain contrast left,
# so the long-range pairs are clean and are exactly where the variance-inflation
# factor accumulates. That is the difference between this and simply truncating the
# correlogram at a fixed distance, which throws the clean long-range pairs away too.
build_within_pair_bins <- function(coords, group, pairs_per_bin = 4000L,
                                   min_pairs_per_bin = 200L, seed = 20260921L) {
  coords <- as.matrix(coords)
  n <- nrow(coords)
  stopifnot(length(group) == n)
  nn <- RANN::nn2(coords, k = min(2L, n))
  h <- stats::median(nn$nn.dists[, ncol(nn$nn.dists)])
  if (!is.finite(h) || h <= 0) stop("Degenerate coordinates: nearest-neighbour distance is not positive.")

  keep <- !is.na(group)
  idx <- which(keep)
  by_group <- split(idx, group[idx])
  i1 <- integer(0); i2 <- integer(0)
  for (g in by_group) {
    m <- length(g)
    if (m < 2L) next
    i1 <- c(i1, rep.int(g[-m], times = (m - 1L):1L))
    i2 <- c(i2, g[sequence.default((m - 1L):1L, from = 2L:m)])
  }
  if (!length(i1)) stop("No within-group pairs are available at this node.")
  d <- sqrt((coords[i1, 1L] - coords[i2, 1L])^2 +
            (coords[i1, 2L] - coords[i2, 2L])^2)

  edges <- c(seq(0.5, 30.5, by = 1), 35.5, 42.5, 55.5, 75) * h
  edges <- edges[edges <= max(d) + h]
  if (length(edges) < 5L) edges <- seq(0.5 * h, max(d) + 1e-12, length.out = 6L)
  b <- findInterval(d, edges, rightmost.closed = TRUE)
  nb <- length(edges) - 1L

  bins <- vector("list", nb)
  counts <- integer(nb); centres <- rep(NA_real_, nb)
  .with_local_seed(seed, {
    for (k in seq_len(nb)) {
      hit <- which(b == k)
      counts[k] <- length(hit)
      if (counts[k] < min_pairs_per_bin) next
      take <- if (counts[k] > pairs_per_bin) sample(hit, pairs_per_bin) else hit
      bins[[k]] <- list(i1 = i1[take], i2 = i2[take])
      centres[k] <- mean(d[take])
    }
  })
  usable <- which(!vapply(bins, is.null, logical(1)))
  list(h = h, edges = edges, bins = bins, usable = usable,
       centre = centres, n_pairs = counts,
       total_within_pairs = length(i1),
       max_within_distance = max(d))
}

# Per-feature correlation at each usable bin, summarised across features by the
# median. A median rather than a mean because a minority of features can carry a
# mean shift the split has not captured; it costs nothing and bounds their effect.
#
# Centring and scaling are done WITHIN each group, not over the whole node. This is
# not a detail. Restricting the pairs while still dividing by the node-wide variance
# leaves the between-child mean difference in the denominator, which deflates every
# bin by the ratio of within-group to total variance and the fit reads that as a
# smaller eta. A test of the code caught exactly this: eta came back at 0.735, 0.361,
# 0.380 and 0.179 against truths of 1.0, 0.5, 0.5 and 0.25, consistently low, while
# ell was unaffected. Both the numerator and the denominator have to see the same
# restriction.
# v0.2.5.5 additionally returns `pooled`: the mean of the centred pair products
# over features AND pairs, divided by the mean within-group variance over features.
#
# This matters, and it is not a cosmetic addition. `median` and `mean` are
# cross-feature summaries of PER-FEATURE ratios, that is estimates of E[num/den].
# The analytic centring expectation derived above is E[num]/E[den]. Those two
# differ whenever the per-feature denominator is noisy, which is exactly when the
# field is smooth relative to the group: measured on a 20x20 grid with 6000
# features, the derivation matches `pooled` to 0.007-0.009 at every setting, while
# its gap to `median` grows 0.042 -> 0.108 -> 0.173 as the effective degrees of
# freedom per feature fall 9.7 -> 3.5 -> 2.0. Fitting the exact expectation to
# `median` would be fitting an exact model to the wrong target, so the analytic
# path uses `pooled` and the intercept path keeps `median` unchanged.
#
# The cost is that `pooled` is not robust to a minority of features carrying a
# mean shift, which is what `median` was there for. Route (c) removes the
# between-child contrast by excluding pairs rather than by trimming features, so
# the spatial exclusion carries that job; whether it carries it well enough is a
# question for the experiment, and the ALT-against-NULL absorption comparison
# measures precisely it.
within_group_correlogram <- function(X, pair_bins, group, feature_chunk = 250L) {
  X <- as.matrix(X)
  p <- ncol(X)
  grp <- split(seq_len(nrow(X)), group)
  grp <- grp[vapply(grp, length, integer(1)) >= 2L]
  nb <- length(pair_bins$bins)
  R <- matrix(NA_real_, p, nb)
  num_sum <- numeric(nb)                 # sum over features of the per-bin mean product
  den_sum <- 0
  # v0.2.5.6: the per-gene numerators and denominators, kept separately. A gene-
  # weighted ratio of sums needs them apart; R = Num / Den cannot be inverted back
  # exactly in floating point. `pooled` below keeps its own summation path so it
  # stays bit-identical to v0.2.5.5.
  Num <- matrix(NA_real_, p, nb)
  Den <- numeric(p)
  for (start in seq.int(1L, p, by = feature_chunk)) {
    cols <- start:min(start + feature_chunk - 1L, p)
    Xc <- X[, cols, drop = FALSE]
    ss <- numeric(length(cols)); nn <- 0L
    for (ix in grp) {
      Xc[ix, ] <- sweep(Xc[ix, , drop = FALSE], 2L,
                        colMeans(Xc[ix, , drop = FALSE]), "-")
      ss <- ss + colSums(Xc[ix, , drop = FALSE]^2)
      nn <- nn + length(ix)
    }
    s2c <- ss / nn                       # pooled within-group variance
    den_sum <- den_sum + sum(s2c)
    Den[cols] <- s2c
    for (k in pair_bins$usable) {
      pr <- pair_bins$bins[[k]]
      cm <- colMeans(Xc[pr$i1, , drop = FALSE] * Xc[pr$i2, , drop = FALSE])
      R[cols, k] <- cm / s2c
      Num[cols, k] <- cm
      num_sum[k] <- num_sum[k] + sum(cm)
    }
    rm(Xc)
  }
  pooled <- rep(NA_real_, nb)
  pooled[pair_bins$usable] <- (num_sum[pair_bins$usable] / p) / (den_sum / p)
  list(per_feature = R,
       median = apply(R, 2L, stats::median, na.rm = TRUE),
       mean = apply(R, 2L, mean, na.rm = TRUE),
       pooled = pooled,
       per_feature_num = Num, per_feature_den = Den)
}

# ===========================================================================
# v0.2.5.5: the ANALYTIC within-group centring correction.
#
# WHY. within_group_correlogram() centres each feature by its own mean INSIDE
# each group, so what it measures is not r(d) but the correlogram of a
# group-centred field. Writing R_ij = eta*k_ij for i != j and R_ii = 1, and with
# m_g = |G_g|, Rbar_{i,g} = mean_{l in G_g} R_il, Rbarbar_g = mean_{k,l in G_g} R_kl:
#
#   E[a_i a_j] / sigma^2 = R_ij - Rbar_{i,g} - Rbar_{j,g} + Rbarbar_g
#   E[s^2]     / sigma^2 = 1 - sum_g (m_g/n_v) * Rbarbar_g
#
# Substituting R and collecting the eta terms, with
#   S_{i,g} = sum_{l in G_g, l != i} k_il        (own-group row sum of the kernel)
#   T_g     = sum_{k != l in G_g} k_kl
#   A_ij    = k_ij - (S_{i,g} + S_{j,g})/m_g + T_g/m_g^2
#   B       = sum_g (m_g/n_v) * T_g/m_g^2
# the expected correlogram in bin b is
#
#   E[r_hat(b)] = ( eta*Abar_b - mean_b(1/m_g) ) / ( (1 - G/n_v) - eta*B )
#
# where Abar_b and mean_b(1/m_g) average over that bin's sampled pairs. This is
# the double-centring identity for a single group, with the row means taken inside
# each spot's OWN group. Checked against simulation.
#
# WHY THE FREE INTERCEPT IS NOT ENOUGH. The bias -Rbar_i - Rbar_j + Rbarbar is not
# a constant: Rbar_{i,g} depends on where spot i sits in its group, and long
# within-group pairs are disproportionately edge-to-edge, so the bias grows with
# distance. Measured on a 20x20 grid at ell = 10, eta = 1, the gap between
# eta*k(d) and the true group-centred expectation runs 0.074 at d = 1 to 0.700 at
# d = 8 -- a factor of 9.5 across the fitted range, and the truth goes NEGATIVE
# (-0.103) where eta*k(d) is +0.597. A single intercept fitted to that leaves a
# structured residual, which the fit reads as faster decay: ell biased low, v
# under-estimated, and the test over-splits. At ell = 3 the same gap runs 0.046 to
# 0.202 and flattens, which is why the intercept is adequate when ell << L_group.
#
# NOTE ON eta. The model above is a RATIO in eta, not linear in it, so the
# closed-form weighted regression of fit_matern_correlogram() no longer applies.
# eta is still a single bounded scalar, so each ell on the grid costs one
# one-dimensional optimisation.
# ===========================================================================

# Everything in the correction that does not depend on ell, tabulated once per
# node. Distances are tallied by their distinct values, which is exact on a
# regular lattice and is what makes each candidate ell a handful of dot products
# instead of an O(n_v^2) sweep. On an irregular geometry the distinct-value count
# can blow up, so it falls back to fine binning and SAYS SO in `mode` rather than
# doing it silently.
build_centring_basis <- function(coords, group, pair_bins,
                                 max_distinct = 4000L, bin_fraction = 0.05) {
  coords <- as.matrix(coords)
  n_v <- nrow(coords)
  keep <- !is.na(group)
  groups <- sort(unique(group[keep]))
  G <- length(groups)
  h <- pair_bins$h

  # One shared distance grid for every tally, so all of them index the same vector.
  d_all <- as.vector(stats::dist(coords))
  mode <- if (length(unique(round(d_all, 10L))) <= max_distinct) "exact" else "binned"
  keyof <- if (mode == "exact")
    function(d) round(d, 10L) else function(d) round(d / (bin_fraction * h)) * (bin_fraction * h)
  levels_d <- sort(unique(keyof(d_all)))
  idx_of <- function(d) match(keyof(d), levels_d)
  nd <- length(levels_d)

  m <- vapply(groups, function(g) sum(group == g, na.rm = TRUE), integer(1))
  names(m) <- as.character(groups)

  # T_g: the own-group kernel sum over ordered distinct pairs, as a distance tally.
  # Also C_g: per-spot own-group distance tallies, needed for S_{i,g}.
  tg <- matrix(0, G, nd)
  Clist <- vector("list", G)
  for (gi in seq_len(G)) {
    ix <- which(group == groups[gi])
    if (length(ix) < 2L) { Clist[[gi]] <- matrix(0, length(ix), nd); next }
    sub <- coords[ix, , drop = FALSE]
    C <- matrix(0, length(ix), nd)
    for (start in seq.int(1L, length(ix), by = 250L)) {
      rows <- start:min(start + 249L, length(ix))
      dd <- sqrt(outer(sub[rows, 1L], sub[, 1L], "-")^2 +
                 outer(sub[rows, 2L], sub[, 2L], "-")^2)
      dd[cbind(seq_along(rows), rows)] <- NA_real_        # drop l == i
      k <- idx_of(as.vector(dd))
      ok <- !is.na(k)
      rr <- rep(seq_along(rows), times = ncol(dd))[ok]
      # tabulate() rather than `C[cbind(i, j)] <- C[cbind(i, j)] + 1`: matrix
      # indexing assigns once per index PAIR, so with repeated pairs the last
      # write wins and the counts come out wrong. A spot has many own-group
      # members at the same distance, so the pairs repeat constantly. Sanity
      # check 7's (A)-vs-(B) comparison caught this; the error grew with ell,
      # up to 0.73 at ell = 10.
      add <- tabulate((rr - 1L) * nd + k[ok], nbins = length(rows) * nd)
      C[rows, ] <- C[rows, ] +
        matrix(add, nrow = length(rows), ncol = nd, byrow = TRUE)
    }
    Clist[[gi]] <- C
    tg[gi, ] <- colSums(C)          # sum over i of (own-group members at distance d)
  }

  # Per (bin, group): the pair-distance tally and the tally that produces
  # sum over the bin's pairs of (S_i + S_j).
  nb <- length(pair_bins$bins)
  pd <- array(0, c(nb, G, nd))      # counts of sampled pairs at each distance
  W  <- array(0, c(nb, G, nd))      # sum over pairs of (C_i + C_j)
  npair <- matrix(0L, nb, G)
  for (b in pair_bins$usable) {
    pr <- pair_bins$bins[[b]]
    gi_of <- match(group[pr$i1], groups)
    dd <- sqrt((coords[pr$i1, 1L] - coords[pr$i2, 1L])^2 +
               (coords[pr$i1, 2L] - coords[pr$i2, 2L])^2)
    ki <- idx_of(dd)
    for (gi in seq_len(G)) {
      sel <- which(gi_of == gi)
      if (!length(sel)) next
      npair[b, gi] <- length(sel)
      pd[b, gi, ] <- tabulate(ki[sel], nbins = nd)
      ix <- which(group == groups[gi])
      local1 <- match(pr$i1[sel], ix); local2 <- match(pr$i2[sel], ix)
      mult <- tabulate(c(local1, local2), nbins = length(ix))
      W[b, gi, ] <- as.vector(crossprod(Clist[[gi]], mult))
    }
  }
  list(mode = mode, levels_d = levels_d, h = h, n_v = n_v, G = G,
       groups = groups, m = m, tg = tg, pd = pd, W = W, npair = npair,
       usable = pair_bins$usable, n_distinct = nd)
}

# E[r_hat(b)] for the group-centred field, for one (ell, eta). Pure arithmetic on
# the tallies: no distance matrix is touched here, so the ell grid is cheap.
centred_correlogram_expectation <- function(basis, ell, eta, bins = basis$usable) {
  kv <- matern32_corr(basis$levels_d, ell, 1)
  Tg <- as.vector(basis$tg %*% kv)                      # length G
  B <- sum((basis$m / basis$n_v) * Tg / basis$m^2)
  den <- (1 - basis$G / basis$n_v) - eta * B
  out <- rep(NA_real_, length(bins))
  for (bi in seq_along(bins)) {
    b <- bins[bi]
    np <- sum(basis$npair[b, ])
    if (!np) next
    num <- 0
    for (gi in seq_len(basis$G)) {
      ng <- basis$npair[b, gi]
      if (!ng) next
      mg <- basis$m[gi]
      sum_k <- sum(basis$pd[b, gi, ] * kv)              # sum over pairs of k_ij
      sum_S <- sum(basis$W[b, gi, ] * kv)               # sum over pairs of (S_i + S_j)
      num <- num + eta * (sum_k - sum_S / mg + ng * Tg[gi] / mg^2) - ng / mg
    }
    out[bi] <- (num / np) / den
  }
  list(expected = out, denominator = den, T_g = Tg, B = B)
}

# Fit r_hat to the group-centred expectation. No free intercept: the offset the
# intercept used to absorb is now part of the model and depends on (ell, eta).
# For each ell on the grid, eta is found by one bounded one-dimensional search,
# because the expectation is a ratio in eta rather than linear in it.
fit_matern_correlogram_centred <- function(r_hat, w, basis, ell_grid,
                                           eta_bounds = c(1e-2, 1),
                                           bins = basis$usable) {
  ok <- is.finite(r_hat) & is.finite(w) & w > 0
  none <- list(ell = NA_real_, eta = NA_real_, offset = NA_real_, sse = NA_real_,
               n_bins = sum(ok), converged = FALSE)
  if (sum(ok) < 4L) return(none)
  r <- r_hat[ok]; ww <- w[ok]; bb <- bins[ok]
  best <- none; best$sse <- Inf
  for (l in ell_grid) {
    f <- function(eta) {
      e <- centred_correlogram_expectation(basis, l, eta, bb)$expected
      if (anyNA(e)) return(Inf)
      sum(ww * (r - e)^2)
    }
    # The two endpoints are evaluated EXPLICITLY alongside the interior search.
    # stats::optimize() cannot return an endpoint -- golden-section converges to an
    # interior point, so at the lower bound it returns 0.0100005 and at the upper
    # 0.999999. An exact `eta <= lower + 1e-12` test then never fires, and
    # eta_at_bound / eta_at_upper_bound would be silently dead on this path: a
    # degenerate fit would be reported as a healthy one.
    # Evaluating the endpoints makes eta exactly equal to the bound when the bound
    # wins, so the flags downstream stay exact comparisons.
    o <- try(stats::optimize(f, interval = eta_bounds, tol = 1e-8), silent = TRUE)
    cand <- list(eta_bounds[1], eta_bounds[2])
    if (!inherits(o, "try-error")) cand <- c(cand, list(o$minimum))
    for (eta in cand) {
      sse <- f(eta)
      if (is.finite(sse) && sse < best$sse)
        best <- list(ell = l, eta = eta, offset = NA_real_, sse = sse,
                     n_bins = sum(ok), converged = TRUE)
    }
  }
  if (!isTRUE(best$converged)) return(none)
  best
}

# v0.2.5.7: continuous refinement of the grid fit (`ell_search = "refine"`).
# The profile sum of squares S(ell) = min over eta of sum_b w_b (r_b - E_b(ell, eta))^2
# is minimised by golden-section search in log(ell) on the bracket formed by the best
# grid value's two neighbours (at either end of the grid, that end and its one
# neighbour). Eta is found at each ell exactly as fit_matern_correlogram_centred()
# finds it: a bounded optimise plus both endpoints, the first minimum in the order
# lower bound, upper bound, interior. That search is repeated here rather than
# factored out of the grid function, so the default path stays byte-for-byte
# v0.2.5.6. The refined fit replaces the grid fit only if its sum of squares is
# strictly lower; the grid fit is always kept as `ell_grid_best` / `sse_grid`.
refine_matern_correlogram_centred <- function(r_hat, w, basis, ell_grid, grid_fit,
                                              eta_bounds = c(1e-2, 1),
                                              bins = basis$usable, tol = 1e-6) {
  ok <- is.finite(r_hat) & is.finite(w) & w > 0
  r <- r_hat[ok]; ww <- w[ok]; bb <- bins[ok]
  best_eta <- function(l) {
    f <- function(eta) {
      e <- centred_correlogram_expectation(basis, l, eta, bb)$expected
      if (anyNA(e)) return(Inf)
      sum(ww * (r - e)^2)
    }
    o <- try(stats::optimize(f, interval = eta_bounds, tol = 1e-8), silent = TRUE)
    cand <- c(eta_bounds[1], eta_bounds[2])
    if (!inherits(o, "try-error")) cand <- c(cand, o$minimum)
    sse <- vapply(cand, f, numeric(1))
    k <- which.min(sse)
    list(eta = cand[k], sse = sse[k])
  }
  i <- which.min(abs(ell_grid - grid_fit$ell))
  lo <- ell_grid[max(1L, i - 1L)]
  hi <- ell_grid[min(length(ell_grid), i + 1L)]
  out <- grid_fit
  out$ell_grid_best <- grid_fit$ell
  out$sse_grid <- grid_fit$sse
  out$refine_improved <- FALSE
  out$refine_bracket <- c(lo, hi)
  o <- try(stats::optimize(function(z) best_eta(exp(z))$sse, interval = log(c(lo, hi)),
                           tol = tol), silent = TRUE)
  if (inherits(o, "try-error")) {
    out$refine_reason <- conditionMessage(attr(o, "condition"))
    return(out)
  }
  l <- exp(o$minimum)
  be <- best_eta(l)
  if (is.finite(be$sse) && be$sse < grid_fit$sse) {
    out$ell <- l
    out$eta <- be$eta
    out$sse <- be$sse
    out$refine_improved <- TRUE
  }
  out
}

# ===========================================================================
# The LEAF partition.
#
# Instead of grouping the node's spots by its two candidate children, group them by
# every LEAF of the candidate subtree under the node: pairs are taken only within a
# leaf, and each leaf is centred by its own mean. Centring within K groups is
# projecting out K group indicators, and the analytic expectation (above) is exact
# for ANY partition, so the correction carries over unchanged -- this is the
# method-of-moments counterpart of REML with a categorical trend.
#
# Why: a candidate child may still hold an unresolved sub-domain, whose mean
# differences would otherwise enter the correlogram; the leaves are the finest
# structure the method ever considers. Where both children are themselves leaves,
# the leaf partition IS the two-child partition and nothing changes.
#
# Known risk, stated before any test: the analytic expectation treats the
# partition as FIXED, but the leaves are clusters chosen on the same data. Centring
# within data-chosen groups can remove more structure than the formula accounts
# for, which would bias ell low. The two-child split carries the same effect, which
# was small in simulations; with more data-chosen groups it may be larger.
# ===========================================================================
leaves_under_node <- function(observed_tree, node_id) {
  nt <- observed_tree$node_table
  out <- integer(0)
  stack <- node_id
  while (length(stack)) {
    v <- stack[1]; stack <- stack[-1]
    row <- match(v, nt$node_id)
    if (is.na(row)) stop(sprintf("node %s is not in the observed tree", v))
    if (isTRUE(nt$is_leaf[row])) { out <- c(out, v); next }
    # Children go on the FRONT, left before right: a depth-first walk, so the leaves
    # come out in the tree's own left-to-right order. An earlier draft appended them
    # at the back and then sorted by node id; see leaf_group_for_node() for why the
    # order is not cosmetic.
    stack <- c(nt$left_child_id[row], nt$right_child_id[row], stack)
  }
  out
}

# Leaf labels 1..K over the node's own spots, in the tree's left-to-right order.
#
# The ORDER matters, although mathematically only the partition does:
# build_within_pair_bins() subsamples pairs group by group in label order, so the
# same partition under different labels draws different pairs. Left-to-right order
# makes label 1 the left child and label 2 the right one wherever both children are
# leaves -- exactly the labels of the two-child partition -- so there the two rules
# are identical bit for bit for ANY tree. The earlier draft labelled leaves by
# sorted node id, which gives the same labels only when the left child's id is the
# smaller, which is not a property of trees in general.
#
# Three ways a tree can be corrupt, and each is checked EXPLICITLY. An earlier draft
# relied on `g[match(...)] <- k` to fail on a leaf spot outside the node, and a
# comment said it would; it does not -- with a scalar replacement R silently drops
# NA subscripts.
leaf_group_for_node <- function(observed_tree, node_id, idx_v) {
  lv <- leaves_under_node(observed_tree, node_id)
  g <- rep(NA_integer_, length(idx_v))
  for (k in seq_along(lv)) {
    pos <- match(observed_tree$spot_indices[[lv[k]]], idx_v)
    if (anyNA(pos))
      stop(sprintf("node %s: leaf %s has spots outside the node", node_id, lv[k]))
    if (any(!is.na(g[pos])))
      stop(sprintf("node %s: leaf %s overlaps another leaf", node_id, lv[k]))
    g[pos] <- k
  }
  if (anyNA(g)) stop(sprintf("node %s: its leaves do not cover all of its spots", node_id))
  g
}

# v0.2.5.5 adds `centring`:
#   "intercept" (DEFAULT) the v0.2.5.4 behaviour -- fit eta*k(d) + c with a free
#               intercept. Kept as the default so nothing changes unless asked.
#   "analytic"  fit the exact group-centred expectation instead, with no free
#               intercept. See build_centring_basis() above for the derivation.
estimate_node_spatial_params <- function(X_v, coords_v, group,
                                         ell_grid_mult = exp(seq(log(0.3), log(40),
                                                                 length.out = 160L)),
                                         eta_bounds = c(1e-2, 1),
                                         pairs_per_bin = 4000L,
                                         min_pairs_per_bin = 200L,
                                         seed = 20260921L,
                                         verbose = FALSE) {
  # Package: the within-group correlogram fit as in the paper (v0.2.5.8 with centring = "analytic",
  # gene_weighting = "none", ell_search = "refine"). The free-intercept fit (v0.2.5.4),
  # Tukey gene weighting (v0.2.5.6) and the grid-only search were removed.
  pb <- tryCatch(build_within_pair_bins(coords_v, group, pairs_per_bin,
                                        min_pairs_per_bin, seed),
                 error = function(e) e)
  if (inherits(pb, "error")) {
    return(list(converged = FALSE, reason = conditionMessage(pb)))
  }
  if (length(pb$usable) < 4L) {
    return(list(converged = FALSE,
                reason = sprintf("only %d usable distance bins (need 4)", length(pb$usable)),
                pair_bins = pb))
  }
  cg <- within_group_correlogram(X_v, pb, group)
  u <- pb$usable
  # ell is searched on a grid expressed in units of the nearest-neighbour spacing,
  # so the same grid is meaningful on any geometry.
  basis <- tryCatch(build_centring_basis(coords_v, group, pb),
                    error = function(e) e)
  if (inherits(basis, "error"))
    return(list(converged = FALSE,
                reason = sprintf("centring basis failed: %s", conditionMessage(basis)),
                pair_bins = pb, correlogram = cg))
  # `pooled`, not `median`: see within_group_correlogram() above. The analytic
  # expectation is E[num]/E[den] and only `pooled` estimates that.
  fit <- fit_matern_correlogram_centred(cg$pooled[u], pb$n_pairs[u], basis,
                                        ell_grid_mult * pb$h, eta_bounds, u)
  if (isTRUE(fit$converged))
    fit <- refine_matern_correlogram_centred(cg$pooled[u], pb$n_pairs[u], basis,
                                             ell_grid_mult * pb$h, fit, eta_bounds, u)
  if (!isTRUE(fit$converged)) {
    return(list(converged = FALSE,
                reason = fit$reason %||% "correlogram fit did not converge",
                pair_bins = pb, correlogram = cg))
  }
  # Guards. eta near zero makes ell unidentifiable while v = 1 + eta * S(ell) still
  # multiplies by ell^2, so an unconstrained fit can return a lengthscale longer
  # than the tissue. Both bounds are recorded rather than applied silently.
  ell_cap <- 2 * pb$max_within_distance
  hit_cap <- fit$ell > ell_cap
  ell_used <- min(fit$ell, ell_cap)
  # v0.2.5.7: the grid's own best ell, whatever the continuous search did afterwards.
  grid_best <- fit$ell_grid_best %||% fit$ell
  if (verbose) {
    message(sprintf(
      "    [spatial] within-child fit (analytic centring): ell=%.6g (%.2f spacings), eta=%.4f, offset=n/a (in the model), bins=%d, pairs=%d%s",
      ell_used, ell_used / pb$h, fit$eta,
      fit$n_bins, pb$total_within_pairs,
      if (hit_cap) "  [ell capped at 2x node extent]" else ""))
  }
  # Package: the expected correlogram at the ell and eta the null uses, per usable bin, for
  # plot(type = "correlogram"). Only recorded; nothing downstream reads it.
  r_fitted <- centred_correlogram_expectation(basis, ell_used, fit$eta, u)$expected
  list(converged = TRUE, ell = ell_used, eta = fit$eta,
       r_fitted = r_fitted,
       centring_mode = basis$mode,
       centring_n_distinct = basis$n_distinct,
       ell_unconstrained = fit$ell, ell_cap = ell_cap, ell_capped = hit_cap,
       eta_at_bound = fit$eta <= eta_bounds[1] + 1e-12,
       eta_at_upper_bound = fit$eta >= eta_bounds[2] - 1e-12,
       n_groups = length(unique(group[!is.na(group)])),
       sse = fit$sse, n_bins = fit$n_bins,
       h = pb$h, ell_in_spacings = ell_used / pb$h,
       total_within_pairs = pb$total_within_pairs,
       max_within_distance = pb$max_within_distance,
       bin_centre = pb$centre[u], bin_pairs = pb$n_pairs[u],
       r_median = cg$median[u], r_mean = cg$mean[u], r_pooled = cg$pooled[u],
       ell_grid_best = grid_best,
       sse_grid = fit$sse_grid %||% fit$sse,
       refine_improved = fit$refine_improved %||% NA,
       ell_at_grid_lower = grid_best <= min(ell_grid_mult * pb$h) * (1 + 1e-12),
       ell_at_grid_upper = grid_best >= max(ell_grid_mult * pb$h) * (1 - 1e-12))
}
