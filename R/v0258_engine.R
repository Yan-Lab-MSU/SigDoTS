# spaSigClust v0.2.5.8HC -- Multi-cluster SHC-style wrapper
# Observed tree: BayesSpace(K0) labels as leaves -> retained-component
# adjacency when connected, otherwise all-component adjacency -> constrained HC
# Node test: node-local null -> null local tree with K=Kv -> root-split CI
# nolint start: object_usage_linter

.bootstrap_node_budget <- function(
  alpha_star,
  B_min,
  B_max_configured,
  bootstrap_B_max_cap
) {
  alpha_star <- as.numeric(alpha_star)
  B_min <- max(1L, as.integer(B_min))
  B_max_configured <- max(B_min, as.integer(B_max_configured))
  bootstrap_B_max_cap <- as.numeric(bootstrap_B_max_cap)

  if (length(alpha_star) != 1L || !is.finite(alpha_star) ||
      alpha_star <= 0 || alpha_star >= 1) {
    stop("alpha_star must be a finite scalar in (0, 1).")
  }
  if (length(bootstrap_B_max_cap) != 1L ||
      !is.finite(bootstrap_B_max_cap) ||
      bootstrap_B_max_cap < B_max_configured) {
    stop("bootstrap_B_max_cap must be finite and >= the configured B_max.")
  }

  required <- max(1, ceiling(1 / alpha_star) - 1)
  uncapped <- max(as.numeric(B_max_configured), required)
  effective <- min(uncapped, floor(bootstrap_B_max_cap))
  resolution_sufficient <- 1 / (effective + 1) <= alpha_star

  list(
    configured = B_max_configured,
    required = required,
    effective = as.integer(effective),
    cap = as.integer(floor(bootstrap_B_max_cap)),
    resolution_sufficient = resolution_sufficient,
    minimum_attainable_p = 1 / (effective + 1)
  )
}

.normalize_observed_clustering <- function(observed_clustering, n, K0) {
  is_result_list <- is.list(observed_clustering) &&
    !is.data.frame(observed_clustering)
  if (is_result_list) {
    if (is.null(observed_clustering$labels)) {
      stop("observed_clustering list must contain `labels`.")
    }
    labels <- observed_clustering$labels
    method <- observed_clustering$method %||% "precomputed"
    diagnostics <- observed_clustering$diagnostics %||% list()
  } else {
    labels <- observed_clustering
    method <- "precomputed"
    diagnostics <- list()
  }

  if ((!is.atomic(labels) && !is.factor(labels)) || !is.null(dim(labels))) {
    stop("observed_clustering labels must be an atomic vector or factor.")
  }
  if (length(labels) != n) {
    stop(sprintf(
      "observed_clustering must provide one label per X row (expected %d; got %d).",
      n,
      length(labels)
    ))
  }
  if (anyNA(labels)) {
    stop("observed_clustering labels cannot contain missing values.")
  }
  label_keys <- as.character(labels)
  if (any(!nzchar(trimws(label_keys)))) {
    stop("observed_clustering labels cannot contain empty values.")
  }

  unique_keys <- unique(label_keys)
  if (length(unique_keys) != as.integer(K0)) {
    stop(sprintf(
      "observed_clustering must contain exactly K0=%d distinct labels; got %d.",
      as.integer(K0),
      length(unique_keys)
    ))
  }
  if (length(method) != 1L || is.na(method) ||
      !nzchar(trimws(as.character(method)))) {
    stop("observed_clustering$method must be one non-empty value when supplied.")
  }
  if (!is.list(diagnostics)) {
    stop("observed_clustering$diagnostics must be a list when supplied.")
  }

  diagnostics$precomputed <- TRUE
  integer_labels <- suppressWarnings(as.integer(label_keys))
  preserves_existing_ids <- all(!is.na(integer_labels)) &&
    all(as.character(integer_labels) == label_keys) &&
    identical(sort(unique(integer_labels)), seq_len(as.integer(K0)))
  normalized_labels <- if (preserves_existing_ids) {
    integer_labels
  } else {
    as.integer(factor(label_keys, levels = sort(unique_keys)))
  }
  list(
    labels = normalized_labels,
    method = method,
    diagnostics = diagnostics,
    precomputed = TRUE
  )
}

# Compute CI from an explicit binary partition on node-local data.
# v0.2.3.2: as in v0.2.3.1, CI is computed directly on X_node (no per-node PCA),
# so the observed CI and null-replicate CIs share the same feature space; the
# caller is responsible for passing X in the desired representation (e.g.,
# scaled HVGs or a fixed global PCA embedding).
compute_ci_from_split <- function(X_node, left_idx_local, right_idx_local) {
  X_node <- as.matrix(X_node)
  n <- nrow(X_node)
  labs <- integer(n)
  labs[left_idx_local] <- 1L
  labs[right_idx_local] <- 2L
  if (any(labs == 0L)) stop("Split does not partition all node spots.")
  cluster_index(X_node, labs)
}

# Build label-level tree from labels and a provided spot graph object.
build_contiguous_tree_from_labels <- function(
  X,
  labels,
  spot_graph,
  ci_pca_dims = 30L,
  n_comp_min = 0.05,
  K_expected = NULL
) {
  build_label_level_tree_from_labels(
    X = X,
    labels = labels,
    spot_graph = spot_graph,
    ci_pca_dims = ci_pca_dims,
    n_comp_min = n_comp_min,
    K_expected = K_expected
  )
}

# Build observed global contiguous tree.
build_observed_global_tree <- function(
  X,
  coords,
  K0,
  k_nn = 6L,
  ci_pca_dims = 30L,
  n_comp_min = 0.05,
  seed = 1L,
  cluster_method = "kmeans",
  bayesspace_args = list(),
  observed_clustering = NULL
) {
  spot_graph <- .connect_pieces(build_knn_graph(coords = coords, k_nn = k_nn), coords)
  bs <- observed_clustering %||% run_clustering(
    X = X,
    coords = coords,
    q = K0,
    cluster_method = cluster_method,
    seed = seed,
    pca_dims = ci_pca_dims,
    bayesspace_args = bayesspace_args
  )
  tree <- build_contiguous_tree_from_labels(
    X = X,
    labels = bs$labels,
    spot_graph = spot_graph,
    ci_pca_dims = ci_pca_dims,
    n_comp_min = n_comp_min,
    K_expected = K0
  )
  if (as.integer(tree$node_table$leaf_count[tree$root_id]) != as.integer(K0)) {
    stop("Observed tree root leaf_count does not match K0.")
  }
  tree$spot_graph <- spot_graph
  tree$bayesspace <- bs
  tree$n_comp_min <- n_comp_min
  tree
}

# v0.2.5.6: every read of `null_args` uses `[[` rather than `$`. On a list, `$`
# PARTIALLY matches names, so `null_args$ell` returned the value of `ell_method`
# whenever a caller set `ell_method` without an explicit `ell` -- and `null_args$eta`
# returned `eta_min` the same way, which is silent. v0.2.5.4's own end-to-end check
# passed only because it happened to pass `ell = NULL` explicitly. `[[` matches
# exactly, so an absent name is NULL and falls through to its default.
# Fit node-local spatial single-domain null.
#
# `child_group` and `feature_group` are integer vectors over the node's spots giving the
# leaf (initial group) of each spot; the engine passes the same vector for both. The spatial fit
# estimates (ell, eta) from pairs of spots in the same leaf, so the between-leaf mean
# differences the node is testing for never enter the correlogram; the gene side is
# estimated after the leaf means are removed by GLS (core:
# whitened_rows_without_group_means()).
fit_node_spatial_null <- function(X_v, coords_v, null_args = list(), verbose = FALSE,
                                  child_group = NULL, feature_group = NULL,
                                  order_seed = NULL) {
  # Package: the spatial parameters always come from the within-group correlogram fit (v0.2.5.8's
  # ell_method = "within_child", centring = "analytic", gene_weighting = "none",
  # ell_search = "refine"), on the Matern 3/2 kernel, and the gene side always uses
  # feature_centring = "leaves_gls". The other options of v0.2.5.8 (fixed fifth-neighbour
  # rule, Moran's I, supplied ell/eta, rbf kernel, grand-mean gene side, z-scoring,
  # "auto" backend, fallback) are not part of the package.
  r_cap <- as.integer(null_args[["r_cap"]] %||% 80L)
  ridge <- as.numeric(null_args[["ridge"]] %||% 1e-6)
  spatial_method <- match.arg(
    null_args[["spatial_method"]] %||% "vecchia",
    choices = c("exact", "vecchia")
  )
  max_exact_peak_gib <- as.numeric(null_args[["max_exact_peak_gib"]] %||% 16)
  vecchia_neighbors <- as.integer(null_args[["vecchia_neighbors"]] %||% 60L)
  nugget <- as.numeric(null_args[["nugget"]] %||% ridge)
  if (is.null(child_group)) {
    stop("The spatial fit needs `child_group`; the caller did not supply the leaves under the node.")
  }
  if (is.null(feature_group)) {
    stop("The gene side needs `feature_group`; the caller did not supply the leaves under the node.")
  }

  ell_range <- as.numeric(null_args[["ell_range"]] %||% c(0.3, 40))
  ell_candidates <- as.integer(null_args[["ell_candidates"]] %||% 160L)
  spatial_estimate <- estimate_node_spatial_params(
    X_v = X_v,
    coords_v = coords_v,
    group = child_group,
    ell_grid_mult = exp(seq(log(ell_range[1]), log(ell_range[2]), length.out = ell_candidates)),
    eta_bounds = c(as.numeric(null_args[["eta_min"]] %||% 1e-2), 1),
    pairs_per_bin = as.integer(null_args[["correlogram_pairs_per_bin"]] %||% 4000L),
    min_pairs_per_bin = as.integer(null_args[["correlogram_min_pairs"]] %||% 200L),
    verbose = verbose
  )
  if (!isTRUE(spatial_estimate$converged)) {
    # Not an error of the run: the caller (node_pvalue_multicluster) catches this class, leaves
    # the node untested and reports it.
    pb <- spatial_estimate$pair_bins
    stop(structure(class = c("sigdots_spatial_fit_failure", "error", "condition"), list(
      message = sprintf("The spatial parameters could not be estimated at this node (%s).",
        spatial_estimate$reason %||% "unknown"),
      call = NULL,
      reason = spatial_estimate$reason %||% "unknown",
      usable_bins = if (is.null(pb)) NA_integer_ else length(pb$usable),
      largest_bin_pairs = if (is.null(pb) || !length(pb$n_pairs)) NA_integer_ else max(pb$n_pairs),
      within_group_pairs = if (is.null(pb)) NA_integer_ else pb$total_within_pairs
    )))
  }
  ell <- spatial_estimate$ell
  eta <- spatial_estimate$eta

  spatial_backend <- build_spatial_backend_auto(
    coords = coords_v,
    ell = ell,
    nugget = nugget,
    max_exact_peak_gib = max_exact_peak_gib,
    vecchia_neighbors = vecchia_neighbors,
    force_method = spatial_method,
    eta = eta,
    order_seed = order_seed
  )
  mu <- colMeans(X_v)
  X_centered <- sweep(X_v, 2L, mu, "-")
  X_whitened <- spatial_whiten(spatial_backend, X_centered)
  gr <- whitened_rows_without_group_means(spatial_backend, X_whitened, feature_group)
  X_whitened <- gr$rows
  gene_side <- list(feature_centring = "leaves_gls", n_groups = gr$n_groups,
    n_rows = gr$n_rows, removed_energy_share = gr$removed_energy_share,
    group_sizes = gr$group_sizes)
  fc <- spiked_cov_est_from_whitened(
    Xw = X_whitened,
    mu = mu,
    r_cap = r_cap,
    ridge = ridge
  )
  list(
    spatial = spatial_backend,
    feature = fc,
    mu = fc$mu,
    ell_used = ell,
    eta_used = eta,
    spatial_estimate = spatial_estimate,
    estimate_used = TRUE,
    gene_side = gene_side
  )
}

simulate_node_null <- function(null_fit, n, p) {
  simulate_matrix_normal_spiked_backend(
    n = n,
    p = p,
    mu = null_fit$mu,
    spatial_backend = null_fit$spatial,
    U_r = null_fit$feature$U_r,
    d_r = null_fit$feature$d_r,
    sigma2 = null_fit$feature$sigma2
  )
}

.summarize_spatial_backend <- function(spatial_backend) {
  selection <- spatial_backend$selection %||% list()
  list(
    requested_method = selection$requested %||% NA_character_,
    selected_method = selection$selected %||% spatial_backend$type %||%
      NA_character_,
    backend_type = spatial_backend$type %||% NA_character_,
    fallback_used = isTRUE(selection$fallback),
    fallback_reason = selection$reason %||% NA_character_,
    estimated_exact_peak_gib = selection$estimated_exact_peak_gib %||%
      NA_real_,
    vecchia_neighbors = as.integer(
      spatial_backend$n_neighbors %||% NA_integer_
    ),
    GpGp_version = spatial_backend$GpGp_version %||% NA_character_,
    # v0.2.5.4: the spatial parameters the null actually used, so a saved run can
    # be audited without reconstructing the fit. eta is 1 under "fixed_5nn".
    ell = spatial_backend$ell %||% NA_real_,
    eta = spatial_backend$eta %||% NA_real_,
    nugget = spatial_backend$nugget %||% NA_real_
  )
}

.summarize_tree_adjacency <- function(tree) {
  adjacency_diagnostics <- tree$adjacency_diagnostics %||% list()
  list(
    adjacency_source = tree$adjacency_source %||% NA_character_,
    component_filter_used_for_tree =
      tree$component_filter_used_for_tree %||% NA,
    component_filter_reason = tree$component_filter_reason %||% NA_character_,
    component_filter_message = tree$component_filter_message %||% NA_character_,
    retained_edge_count = as.integer(
      adjacency_diagnostics$retained_edge_count %||% NA_integer_
    ),
    retained_connected = adjacency_diagnostics$retained_connected %||% NA,
    unfiltered_edge_count = as.integer(
      adjacency_diagnostics$unfiltered_edge_count %||% NA_integer_
    ),
    unfiltered_connected = adjacency_diagnostics$unfiltered_connected %||% NA
  )
}

# Build local null subtree for one node replicate.
build_local_null_subtree <- function(
  X_v_rep,
  coords_v,
  spot_graph_v,
  Kv,
  seed,
  ci_pca_dims = 30L,
  n_comp_min = 0.05,
  cluster_method = "kmeans",
  bayesspace_args = list()
) {
  bs <- run_clustering(
    X = X_v_rep,
    coords = coords_v,
    q = Kv,
    cluster_method = cluster_method,
    seed = seed,
    pca_dims = ci_pca_dims,
    bayesspace_args = bayesspace_args
  )
  tree <- build_contiguous_tree_from_labels(
    X = X_v_rep,
    labels = bs$labels,
    spot_graph = spot_graph_v,
    ci_pca_dims = ci_pca_dims,
    n_comp_min = n_comp_min,
    K_expected = Kv
  )
  tree$bayesspace <- bs
  tree$n_comp_min <- n_comp_min
  tree
}

compute_node_ci_observed <- function(observed_tree, node_id, X) {
  split <- extract_observed_node_split(observed_tree, node_id)
  node_spots <- observed_tree$spot_indices[[node_id]]
  left_local <- match(split$left_spots, node_spots)
  right_local <- match(split$right_spots, node_spots)
  if (any(is.na(left_local)) || any(is.na(right_local))) {
    stop("Observed split mapping failed for node-local indexing.")
  }
  X_v <- X[node_spots, , drop = FALSE]
  compute_ci_from_split(X_v, left_local, right_local)
}

compute_node_ci_null <- function(X_v_rep, null_tree) {
  split <- extract_root_split_from_tree(null_tree)
  compute_ci_from_split(
    X_node = X_v_rep,
    left_idx_local = split$left_spots,
    right_idx_local = split$right_spots
  )
}

# CI = WSS/TSS, so smaller values indicate stronger clustering. The normal
# approximation must match the bootstrap lower tail P(CI_null <= CI_obs).
.normal_lower_tail_p <- function(x, mu, sigma) {
  if (!is.finite(sigma) || sigma <= 0) return(NA_real_)
  stats::pnorm(x, mean = mu, sd = sigma, lower.tail = TRUE)
}

# Node-level Monte Carlo p-value in multi-cluster setting.
# v0.2.3.2: SHC-style sequential null-replicate stopping. After each batch of
# valid replicates we update b_valid and c_b = #{CI_null <= CI_obs}, then check
#   p_lower = (1 + c_b) / (1 + B_max)
#   p_upper = (1 + c_b + (B_max - b_valid)) / (1 + B_max)
# Once b_valid >= B_min, if p_lower > alpha_star we certify (no remaining
# replicate can make the result significant); if p_upper <= alpha_star we
# reject (no remaining replicate can make it nonsignificant).
node_pvalue_multicluster <- function(
  observed_tree,
  node_id,
  X,
  coords,
  B_min = 20L,
  B_max = 120L,
  bootstrap_B_max_cap = 5000L,
  alpha = 0.05,
  n_total = NULL,
  ci_pca_dims = 30L,
  min_node_size = 20L,
  min_child_size = 10L,
  seed = 1L,
  n_comp_min = 0.05,
  cluster_method = "kmeans",
  bayesspace_args = list(),
  null_args = list(),
  verbose = FALSE,
  n_cores = 1L,
  pvalue_method = c("bootstrap", "normal_adaptive"),
  normal_n_initial = 10L,
  normal_n_refine = 50L,
  normal_band_multiplier = 10
) {
  # pvalue_method: "bootstrap" (the Monte Carlo p-value) or "normal_adaptive" (the normal
  # approximation of v0.2.5.8: normal_n_initial replicates, then up to normal_n_refine when
  # the first p-value lies within a factor normal_band_multiplier of the level).
  pvalue_method <- match.arg(pvalue_method)
  B_min_int <- as.integer(B_min)
  B_max_int <- as.integer(B_max)
  if (B_min_int < 1L) B_min_int <- 1L
  if (B_max_int < B_min_int) B_max_int <- B_min_int
  B_max_configured_int <- B_max_int
  bootstrap_budget <- NULL
  normal_n_initial_int <- max(2L, as.integer(normal_n_initial))
  normal_n_refine_int <- max(normal_n_initial_int, as.integer(normal_n_refine))
  normal_band_multiplier <- as.numeric(normal_band_multiplier)
  if (!is.finite(normal_band_multiplier) || normal_band_multiplier <= 1) {
    stop("normal_band_multiplier must be finite and > 1.")
  }
  if (identical(pvalue_method, "normal_adaptive")) {
    B_max_int <- max(B_max_int, normal_n_refine_int)
  }

  node <- observed_tree$node_table[node_id, , drop = FALSE]
  if (isTRUE(node$is_leaf)) {
    return(list(tested = FALSE, reason_skipped = "leaf_node"))
  }

  idx_v <- observed_tree$spot_indices[[node_id]]
  n_v <- length(idx_v)
  if (n_v < as.integer(min_node_size)) {
    return(list(tested = FALSE, reason_skipped = "below_min_node_size", n_v = n_v))
  }

  obs_split <- extract_observed_node_split(observed_tree, node_id)
  left_n <- length(obs_split$left_spots)
  right_n <- length(obs_split$right_spots)
  if (left_n < as.integer(min_child_size) || right_n < as.integer(min_child_size)) {
    return(list(
      tested = FALSE,
      reason_skipped = "below_min_child_size_observed_split",
      n_v = n_v
    ))
  }

  if (is.null(n_total)) n_total <- nrow(X)
  alpha_star <- as.numeric(alpha) * ((n_v - 1) / (as.numeric(n_total) - 1))

  Kv <- as.integer(observed_tree$node_table$leaf_count[node_id])
  if (Kv < 2L) {
    return(list(tested = FALSE, reason_skipped = "Kv_less_than_2", n_v = n_v, Kv = Kv))
  }

  if (identical(pvalue_method, "bootstrap")) {
    bootstrap_budget <- .bootstrap_node_budget(
      alpha_star = alpha_star,
      B_min = B_min_int,
      B_max_configured = B_max_configured_int,
      bootstrap_B_max_cap = bootstrap_B_max_cap
    )
    B_max_int <- bootstrap_budget$effective
    if (!isTRUE(bootstrap_budget$resolution_sufficient)) {
      if (verbose) {
        message(sprintf(
          paste0(
            "    [Node %d] Bootstrap skipped: required B_max=%s exceeds ",
            "cap=%d (alpha*=%.6g, minimum attainable p=%.6g)."
          ),
          node_id, format(bootstrap_budget$required, scientific = FALSE),
          bootstrap_budget$cap, alpha_star,
          bootstrap_budget$minimum_attainable_p
        ))
      }
      return(list(
        tested = FALSE,
        reason_skipped = "bootstrap_resolution_exceeds_cap",
        n_v = n_v,
        Kv = Kv,
        alpha_star = alpha_star,
        pvalue_method = pvalue_method,
        B_min = B_min_int,
        B_max = B_max_int,
        B_max_configured = bootstrap_budget$configured,
        B_max_required = bootstrap_budget$required,
        B_max_effective = bootstrap_budget$effective,
        bootstrap_B_max_cap = bootstrap_budget$cap,
        bootstrap_resolution_sufficient = FALSE,
        minimum_attainable_p = bootstrap_budget$minimum_attainable_p
      ))
    }
  }

  X_v <- X[idx_v, , drop = FALSE]
  coords_v <- coords[idx_v, , drop = FALSE]
  spot_graph_v <- subset_spot_graph(observed_tree$spot_graph, idx_v)

  CI_obs <- compute_node_ci_observed(observed_tree, node_id, X)

  # Package: the spot graph restricted to this part may fall into separate pieces (a layer cut
  # by a fold; spots whose neighbours all belong to other groups). A simulated dataset whose
  # clustering gives a piece a group of its own then has no connected tree, and v0.2.5.8
  # stopped (seen 2026-10-29 on DLPFC section 151673 with BANKSY groups). The pieces are joined
  # by their closest spots (.connect_pieces()); a connected spot graph gives a connected group
  # graph for every clustering. Nothing changes when the graph is already connected.
  pieces <- sort(igraph::components(spot_graph_v$graph)$csize, decreasing = TRUE)
  if (length(pieces) > 1L) {
    if (verbose) message(sprintf("    [Node %d] its spots form %d separate pieces (%s); joined by their closest spots.",
                                 node_id, length(pieces), paste(pieces, collapse = ", ")))
    spot_graph_v <- .connect_pieces(spot_graph_v, coords_v)
  }

  if (verbose) message(sprintf("    [Node %d] Fitting spatial null (kernel=matern32)...", node_id))
  # The leaves of the candidate subtree under the node group both the spatial estimate
  # (v0.2.5.8's centring_partition = "leaves") and the gene side (feature_centring =
  # "leaves_gls"); v0.2.5.8 built the same vector twice.
  leaf_group_v <- leaf_group_for_node(observed_tree, node_id, idx_v)
  # The option SigDoTS.v0258_compatible = TRUE restores two v0.2.5.8 behaviours that made
  # results depend on the number of cores (the Vecchia ordering drawn from the global random
  # stream; the stopping rule checked only after each batch). It exists for the package's
  # bit-for-bit checks against v0.2.5.8 only.
  v0258_compatible <- isTRUE(getOption("SigDoTS.v0258_compatible", FALSE))
  # The Vecchia ordering's seed: this test's own, distinct from its replicates' seeds
  # (seed + node_id * 100000 + b, b >= 1).
  order_seed <- if (v0258_compatible) NULL else as.integer(seed) + as.integer(node_id) * 100000L
  null_fit <- tryCatch(fit_node_spatial_null(
    X_v = X_v,
    coords_v = coords_v,
    null_args = null_args,
    verbose = verbose,
    child_group = leaf_group_v,
    feature_group = leaf_group_v,
    order_seed = order_seed
  ), sigdots_spatial_fit_failure = function(e) e)
  if (inherits(null_fit, "sigdots_spatial_fit_failure")) {
    if (verbose) message(sprintf("    [Node %d] %s", node_id, conditionMessage(null_fit)))
    return(list(
      tested = FALSE,
      reason_skipped = "spatial_fit_failed",
      n_v = n_v,
      Kv = Kv,
      CI_obs = CI_obs,
      alpha_star = alpha_star,
      spatial_fit_reason = null_fit$reason,
      spatial_fit_usable_bins = null_fit$usable_bins,
      spatial_fit_largest_bin_pairs = null_fit$largest_bin_pairs,
      spatial_fit_within_group_pairs = null_fit$within_group_pairs,
      min_pairs_per_bin = as.integer(null_args[["correlogram_min_pairs"]] %||% 200L)
    ))
  }
  spatial_backend_diagnostics <- .summarize_spatial_backend(null_fit$spatial)

  n_cores_use <- max(1L, min(as.integer(n_cores), B_max_int))
  batch_size <- n_cores_use
  piece_sizes <- if (length(pieces) > 1L) pieces else NULL
  if (verbose) {
    if (identical(pvalue_method, "bootstrap")) {
      message(sprintf(
        paste0(
          "    [Node %d] Null fitted. Sequential MC: B_min=%d, ",
          "B_max=%d (configured=%d, required=%s, cap=%d), ",
          "alpha*=%.6g, batch=%d"
        ),
        node_id, B_min_int, B_max_int, bootstrap_budget$configured,
        format(bootstrap_budget$required, scientific = FALSE),
        bootstrap_budget$cap, alpha_star, batch_size))
    } else {
      message(sprintf(
        "    [Node %d] Null fitted. Adaptive normal MC: n=%d/%d, alpha*=%.4f, band=x%.1f, batch=%d",
        node_id, normal_n_initial_int, normal_n_refine_int, alpha_star,
        normal_band_multiplier, batch_size))
    }
  }

  p_v <- ncol(X_v)
  rep_seeds <- as.integer(seed) + as.integer(node_id) * 100000L + seq_len(B_max_int)
  min_child_int <- as.integer(min_child_size)

  mc_one_rep <- function(b) {
    set.seed(rep_seeds[b])
    tryCatch({
      Xb <- simulate_node_null(null_fit = null_fit, n = n_v, p = p_v)

      t_null <- build_local_null_subtree(
        X_v_rep = Xb, coords_v = coords_v, spot_graph_v = spot_graph_v,
        Kv = Kv, seed = rep_seeds[b],
        ci_pca_dims = ci_pca_dims, n_comp_min = n_comp_min,
        cluster_method = cluster_method,
        bayesspace_args = bayesspace_args
      )

      rs <- extract_root_split_from_tree(t_null)
      if (length(rs$left_spots) < min_child_int ||
          length(rs$right_spots) < min_child_int) {
        stop("null_root_split_below_min_child_size")
      }

      ci_val <- compute_node_ci_null(X_v_rep = Xb, null_tree = t_null)

      list(
        ci = ci_val,
        tree_adjacency = .summarize_tree_adjacency(t_null)
      )
    }, error = function(e) e)
  }

  run_batch <- function(batch_idx) {
    if (n_cores_use > 1L && length(batch_idx) > 1L) {
      parallel::mclapply(batch_idx, mc_one_rep, mc.cores = n_cores_use,
                         mc.preschedule = TRUE)
    } else {
      lapply(batch_idx, mc_one_rep)
    }
  }

  CI_null <- numeric(0)
  null_adjacency_sources <- character(0)
  null_component_filter_reasons <- character(0)
  failures <- 0L
  early_stop_reason <- NA_character_
  mc_t0 <- proc.time()[3]
  next_idx <- 1L
  p_lower_now <- NA_real_
  p_upper_now <- NA_real_

  append_valid_rep <- function(res, rep_b) {
    CI_null <<- c(CI_null, as.numeric(res$ci))
    tree_adjacency <- res$tree_adjacency %||% list()
    null_adjacency_sources <<- c(
      null_adjacency_sources,
      tree_adjacency$adjacency_source %||% "unknown"
    )
    null_component_filter_reasons <<- c(
      null_component_filter_reasons,
      tree_adjacency$component_filter_reason %||% NA_character_
    )
    invisible(NULL)
  }

  summarize_null_tree_adjacency <- function() {
    n_valid <- length(null_adjacency_sources)
    n_retained <- sum(
      null_adjacency_sources == "retained_components_after_filtering"
    )
    n_unfiltered <- sum(
      null_adjacency_sources == "all_components_without_filtering"
    )
    reason_values <- null_component_filter_reasons[
      !is.na(null_component_filter_reasons) &
        nzchar(null_component_filter_reasons)
    ]
    list(
      n_valid_replicates = n_valid,
      retained_components_after_filtering_count = n_retained,
      all_components_without_filtering_count = n_unfiltered,
      unknown_adjacency_source_count =
        n_valid - n_retained - n_unfiltered,
      all_components_without_filtering_fraction = if (n_valid > 0L) {
        n_unfiltered / n_valid
      } else {
        NA_real_
      },
      component_filter_reason_counts = as.list(table(reason_values))
    )
  }

  node_diagnostics <- function(...) {
    diagnostics <- list(...)
    diagnostics$null_tree_adjacency <- summarize_null_tree_adjacency()
    diagnostics
  }

  # Package (2026-11-01): the normal approximation of v0.2.5.8 (wrapper lines 885-1035). A
  # replicate that fails stops the run, as in the Monte Carlo branch below (strict mode is the
  # only mode). The replicates used are always numbers 1..normal_n_initial, then up to
  # normal_n_refine, so the result does not depend on n_cores.
  if (identical(pvalue_method, "normal_adaptive")) {
    run_until_valid <- function(target_valid) {
      while (next_idx <= B_max_int && length(CI_null) < target_valid) {
        remaining_valid <- target_valid - length(CI_null)
        batch_end <- min(next_idx + min(batch_size, remaining_valid) - 1L, B_max_int)
        batch_idx <- next_idx:batch_end

        if (verbose) {
          elapsed_mc <- proc.time()[3] - mc_t0
          message(sprintf(
            "    [Node %d] Normal MC batch %d-%d/%d  (target_valid=%d, b_valid=%d, fail=%d, cores=%d, elapsed=%.0fs)",
            node_id, next_idx, batch_end, B_max_int, target_valid,
            length(CI_null), failures, length(batch_idx), elapsed_mc))
        }

        batch_results <- run_batch(batch_idx)
        for (i in seq_along(batch_results)) {
          res <- batch_results[[i]]
          if (is.null(res) || inherits(res, "error") || inherits(res, "try-error")) {
            failures <<- failures + 1L
            msg <- if (inherits(res, "error")) res$message else "MC worker failed"
            stop(msg)
          } else if (
            is.list(res) && "ci" %in% names(res) &&
              length(res$ci) == 1L && is.finite(res$ci)
          ) {
            append_valid_rep(res, batch_idx[i])
          } else {
            failures <<- failures + 1L
            stop("MC worker returned a malformed result.")
          }
        }

        next_idx <<- batch_end + 1L
      }
    }

    run_until_valid(normal_n_initial_int)
    if (length(CI_null) < normal_n_initial_int) {
      return(list(
        tested = FALSE,
        reason_skipped = "insufficient_valid_normal_replicates",
        n_v = n_v,
        Kv = Kv,
        CI_obs = CI_obs,
        alpha_star = alpha_star,
        pvalue_method = pvalue_method,
        n_valid_reps = length(CI_null),
        n_failed_reps = failures,
        piece_sizes = piece_sizes,
        diagnostics = node_diagnostics(
          null_fit = null_fit,
          spatial_backend = spatial_backend_diagnostics,
          failures = failures
        )
      ))
    }

    mu_initial <- mean(CI_null[seq_len(normal_n_initial_int)])
    sigma_initial <- stats::sd(CI_null[seq_len(normal_n_initial_int)])
    p_initial <- .normal_lower_tail_p(CI_obs, mu_initial, sigma_initial)
    if (!is.finite(p_initial)) {
      return(list(
        tested = FALSE,
        reason_skipped = "invalid_normal_approximation",
        n_v = n_v,
        Kv = Kv,
        CI_obs = CI_obs,
        alpha_star = alpha_star,
        pvalue_method = pvalue_method,
        n_valid_reps = length(CI_null),
        n_failed_reps = failures,
        piece_sizes = piece_sizes,
        diagnostics = node_diagnostics(
          null_fit = null_fit,
          spatial_backend = spatial_backend_diagnostics,
          normal_mu_initial = mu_initial,
          normal_sigma_initial = sigma_initial
        )
      ))
    }

    band_low <- alpha_star / normal_band_multiplier
    band_high <- alpha_star * normal_band_multiplier
    normal_used_refine <- p_initial >= band_low && p_initial <= band_high
    if (normal_used_refine) run_until_valid(normal_n_refine_int)

    b_valid <- length(CI_null)
    normal_mu <- mean(CI_null)
    normal_sigma <- stats::sd(CI_null)
    pval <- .normal_lower_tail_p(CI_obs, normal_mu, normal_sigma)
    if (!is.finite(pval)) {
      return(list(
        tested = FALSE,
        reason_skipped = "invalid_normal_approximation",
        n_v = n_v,
        Kv = Kv,
        CI_obs = CI_obs,
        alpha_star = alpha_star,
        pvalue_method = pvalue_method,
        n_valid_reps = b_valid,
        n_failed_reps = failures,
        piece_sizes = piece_sizes,
        diagnostics = node_diagnostics(
          null_fit = null_fit,
          spatial_backend = spatial_backend_diagnostics,
          normal_mu = normal_mu,
          normal_sigma = normal_sigma
        )
      ))
    }

    # v0.2.5.8 wrote "normal_refine_to_50" / "normal_stop_after_10" whatever the two numbers
    # were; the package writes the numbers used (the same strings at the defaults).
    early_stop_reason <- if (normal_used_refine) {
      sprintf("normal_refine_to_%d", normal_n_refine_int)
    } else {
      sprintf("normal_stop_after_%d", normal_n_initial_int)
    }
    return(list(
      tested = TRUE,
      reason_skipped = NA_character_,
      n_v = n_v,
      Kv = Kv,
      CI_obs = CI_obs,
      CI_null = CI_null,
      raw_p = pval,
      p_lower = NA_real_,
      p_upper = NA_real_,
      b_valid = b_valid,
      early_stop_reason = early_stop_reason,
      piece_sizes = piece_sizes,
      significant_decision = pval <= alpha_star,
      alpha_star = alpha_star,
      n_valid_reps = b_valid,
      n_failed_reps = failures,
      B_min = B_min_int,
      B_max = B_max_int,
      pvalue_method = pvalue_method,
      normal_n_reps = b_valid,
      normal_used_refine = normal_used_refine,
      normal_mu = normal_mu,
      normal_sigma = normal_sigma,
      normal_p_initial = p_initial,
      diagnostics = node_diagnostics(
        null_fit = null_fit,
        spatial_backend = spatial_backend_diagnostics,
        normal_mu_initial = mu_initial,
        normal_sigma_initial = sigma_initial,
        normal_band = c(lower = band_low, upper = band_high)
      )
    ))
  }

  # A null replicate that fails stops the run (v0.2.5.8's strict mode, now the only mode):
  # dropping it would change the p-value without saying so.
  while (next_idx <= B_max_int) {
    batch_end <- min(next_idx + batch_size - 1L, B_max_int)
    batch_idx <- next_idx:batch_end

    if (verbose) {
      elapsed_mc <- proc.time()[3] - mc_t0
      done_so_far <- next_idx - 1L
      rate <- if (done_so_far > 0L) elapsed_mc / done_so_far else NA_real_
      eta <- if (!is.na(rate)) rate * (B_max_int - done_so_far) else NA_real_
      message(sprintf(
        "    [Node %d] MC batch %d-%d/%d  (b_valid=%d, fail=%d, cores=%d, elapsed=%.0fs, ETA<=%.0fs)",
        node_id, next_idx, batch_end, B_max_int,
        length(CI_null), failures, length(batch_idx),
        elapsed_mc, ifelse(is.na(eta), 0, eta)))
    }

    batch_results <- run_batch(batch_idx)
    # The package checks the stopping rule after every replicate, in replicate order, and
    # ignores the rest of the batch once it stops; so the replicates used, b_valid and the
    # p-value are those of a one-core run, whatever n_cores is. (Decisions never depended on
    # the batch size: the rule stops only when the outcome is certain.) v0.2.5.8 checked once
    # per batch; v0258_compatible keeps that.
    for (i in seq_along(batch_results)) {
      res <- batch_results[[i]]
      if (is.null(res) || inherits(res, "error") || inherits(res, "try-error")) {
        failures <- failures + 1L
        msg <- if (inherits(res, "error")) res$message else "MC worker failed"
        stop(msg)
      } else if (
        is.list(res) && "ci" %in% names(res) &&
          length(res$ci) == 1L && is.finite(res$ci)
      ) {
        append_valid_rep(res, batch_idx[i])
      } else {
        failures <- failures + 1L
        stop("MC worker returned a malformed result.")
      }
      if (!v0258_compatible && i < length(batch_results) && length(CI_null) >= B_min_int) {
        b_now <- length(CI_null)
        c_now <- sum(CI_null <= CI_obs)
        if ((1 + c_now) / (1 + B_max_int) > alpha_star ||
            (1 + c_now + (B_max_int - b_now)) / (1 + B_max_int) <= alpha_star) {
          break
        }
      }
    }

    next_idx <- batch_end + 1L
    b_valid <- length(CI_null)
    c_b <- sum(CI_null <= CI_obs)
    p_lower_now <- (1 + c_b) / (1 + B_max_int)
    p_upper_now <- (1 + c_b + (B_max_int - b_valid)) / (1 + B_max_int)

    if (verbose) {
      message(sprintf(
        "    [Node %d]   -> b_valid=%d, c_b=%d, p in [%.4f, %.4f], alpha*=%.4f",
        node_id, b_valid, c_b, p_lower_now, p_upper_now, alpha_star))
    }

    if (b_valid >= B_min_int) {
      if (p_lower_now > alpha_star) {
        early_stop_reason <- "certify_p_lower_gt_alpha"
        if (verbose) message(sprintf(
          "    [Node %d]   EARLY STOP: certify (p_lower=%.4f > alpha*=%.4f)",
          node_id, p_lower_now, alpha_star))
        break
      }
      if (p_upper_now <= alpha_star) {
        early_stop_reason <- "reject_p_upper_le_alpha"
        if (verbose) message(sprintf(
          "    [Node %d]   EARLY STOP: reject (p_upper=%.4f <= alpha*=%.4f)",
          node_id, p_upper_now, alpha_star))
        break
      }
    }
  }

  b_valid <- length(CI_null)
  c_b <- sum(CI_null <= CI_obs)
  pval <- (1 + c_b) / (1 + b_valid)

  if (is.na(early_stop_reason)) {
    significant_decision <- pval <= alpha_star
  } else if (identical(early_stop_reason, "reject_p_upper_le_alpha")) {
    significant_decision <- TRUE
  } else {
    significant_decision <- FALSE
  }

  list(
    tested = TRUE,
    reason_skipped = NA_character_,
    n_v = n_v,
    Kv = Kv,
    CI_obs = CI_obs,
    CI_null = CI_null,
    raw_p = pval,
    p_lower = p_lower_now,
    p_upper = p_upper_now,
    b_valid = b_valid,
    early_stop_reason = early_stop_reason,
    piece_sizes = piece_sizes,
    significant_decision = significant_decision,
    alpha_star = alpha_star,
    n_valid_reps = b_valid,
    n_failed_reps = failures,
    B_min = B_min_int,
    B_max = B_max_int,
    B_max_configured = bootstrap_budget$configured,
    B_max_required = bootstrap_budget$required,
    B_max_effective = bootstrap_budget$effective,
    bootstrap_B_max_cap = bootstrap_budget$cap,
    bootstrap_resolution_sufficient =
      bootstrap_budget$resolution_sufficient,
    minimum_attainable_p = bootstrap_budget$minimum_attainable_p,
    pvalue_method = pvalue_method,
    normal_n_reps = NA_integer_,
    normal_used_refine = NA,
    normal_mu = NA_real_,
    normal_sigma = NA_real_,
    normal_p_initial = NA_real_,
    diagnostics = node_diagnostics(
      null_fit = null_fit,
      spatial_backend = spatial_backend_diagnostics
    )
  )
}

traverse_tree_topdown <- function(
  observed_tree,
  X,
  coords,
  B_min = 20L,
  B_max = 120L,
  bootstrap_B_max_cap = 5000L,
  alpha = 0.05,
  ci_pca_dims = 30L,
  min_node_size = 20L,
  min_child_size = 10L,
  seed = 1L,
  n_comp_min = 0.05,
  cluster_method = "kmeans",
  bayesspace_args = list(),
  null_args = list(),
  verbose = TRUE,
  n_cores = 1L,
  progress = FALSE,
  pvalue_method = c("bootstrap", "normal_adaptive"),
  normal_n_initial = 10L,
  normal_n_refine = 50L,
  normal_band_multiplier = 10
) {
  pvalue_method <- match.arg(pvalue_method)
  B_min_int <- as.integer(B_min)
  B_max_int <- as.integer(B_max)
  n_total <- nrow(X)
  # The Monte Carlo budget (B_max raised to what the level needs) means nothing for the normal
  # approximation; v0.2.5.8 reported NA then (wrapper lines 1265 and 1399).
  na_budget <- list(configured = NA_integer_, required = NA_real_, effective = NA_integer_,
                    cap = NA_integer_, resolution_sufficient = NA, minimum_attainable_p = NA_real_)
  node_budget <- function(alpha_star) {
    if (!identical(pvalue_method, "bootstrap")) return(na_budget)
    .bootstrap_node_budget(alpha_star = alpha_star, B_min = B_min_int,
                           B_max_configured = B_max_int, bootstrap_B_max_cap = bootstrap_B_max_cap)
  }
  internal_nodes <- get_internal_nodes_topdown(observed_tree)
  n_internal <- length(internal_nodes)
  tested_map <- setNames(vector("list", n_internal), as.character(internal_nodes))
  visited <- integer(0)
  certified_nodes <- integer(0)
  untested_nodes <- integer(0)
  visit_counter <- 0L
  traverse_t0 <- proc.time()[3]

  visit <- function(node_id) {
    node <- observed_tree$node_table[node_id, , drop = FALSE]
    if (isTRUE(node$is_leaf)) {
      if (verbose) message(sprintf("  [traverse] Node %d is a leaf -> certified", node_id))
      certified_nodes <<- c(certified_nodes, node_id)
      return(invisible(NULL))
    }

    visit_counter <<- visit_counter + 1L
    visited <<- c(visited, node_id)
    n_v <- length(observed_tree$spot_indices[[node_id]])
    Kv_display <- observed_tree$node_table$leaf_count[node_id]
    alpha_star <- alpha * ((n_v - 1) / (n_total - 1))
    display_budget <- node_budget(alpha_star)

    if (verbose) {
      elapsed_trav <- proc.time()[3] - traverse_t0
      message(sprintf(
        "\n  [traverse] === Node %d  (%d/%d internal nodes visited, %.0fs elapsed) ===",
        node_id, visit_counter, n_internal, elapsed_trav))
      if (identical(pvalue_method, "bootstrap")) {
        message(sprintf(
          paste0(
            "  [traverse]     n_v=%d spots, Kv=%d leaves, ",
            "alpha*=%.6g, B in [%d, %d] ",
            "(configured=%d, required=%s, cap=%d)"
          ),
          n_v, Kv_display, alpha_star, B_min_int,
          display_budget$effective, display_budget$configured,
          format(display_budget$required, scientific = FALSE),
          display_budget$cap
        ))
      } else {
        message(sprintf(
          "  [traverse]     n_v=%d spots, Kv=%d leaves, alpha*=%.6g, normal approximation (%d, up to %d replicates)",
          n_v, Kv_display, alpha_star, as.integer(normal_n_initial), as.integer(normal_n_refine)))
      }
    }

    node_t0 <- proc.time()[3]
    res <- node_pvalue_multicluster(
      observed_tree = observed_tree,
      node_id = node_id,
      X = X,
      coords = coords,
      B_min = B_min_int,
      B_max = B_max_int,
      bootstrap_B_max_cap = bootstrap_B_max_cap,
      alpha = alpha,
      n_total = n_total,
      ci_pca_dims = ci_pca_dims,
      min_node_size = min_node_size,
      min_child_size = min_child_size,
      seed = seed,
      n_comp_min = n_comp_min,
      cluster_method = cluster_method,
      bayesspace_args = bayesspace_args,
      null_args = null_args,
      verbose = verbose,
      n_cores = n_cores,
      pvalue_method = pvalue_method,
      normal_n_initial = normal_n_initial,
      normal_n_refine = normal_n_refine,
      normal_band_multiplier = normal_band_multiplier
    )
    node_elapsed <- proc.time()[3] - node_t0
    res$alpha_star <- res$alpha_star %||% alpha_star
    if (isTRUE(res$tested)) {
      if (!is.null(res$significant_decision)) {
        res$significant <- isTRUE(res$significant_decision)
      } else {
        res$significant <- !is.null(res$raw_p) && is.finite(res$raw_p) &&
          (res$raw_p <= alpha_star)
      }
    } else {
      res$significant <- FALSE
    }
    tested_map[[as.character(node_id)]] <<- res
    if (progress && !verbose) {
      outcome <- if (isTRUE(res$tested)) {
        sprintf("%s (%d simulated datasets)", if (isTRUE(res$significant)) "split" else "no split",
                as.integer(res$b_valid))
      } else if (identical(res$reason_skipped, "spatial_fit_failed")) {
        "NOT TESTED: the spatial parameters could not be estimated"
      } else sprintf("not tested (%s)", res$reason_skipped %||% "unknown")
      message(sprintf("  test at node %d (%d spots, %d candidate domains): %s [%.0f s]",
                      node_id, n_v, as.integer(Kv_display), outcome, node_elapsed))
    }

    if (verbose) {
      if (!isTRUE(res$tested)) {
        message(sprintf(
          "  [traverse] Node %d SKIPPED (%s) [%.1fs]",
          node_id, res$reason_skipped %||% "unknown", node_elapsed))
      } else {
        es_tag <- if (is.na(res$early_stop_reason %||% NA_character_)) "no_early_stop"
                  else res$early_stop_reason
        message(sprintf(
          "  [traverse] Node %d DONE: p=%.4f (b_valid=%d/%d, [%.4f,%.4f], stop=%s) -> %s [%.1fs]",
          node_id,
          res$raw_p %||% NA_real_,
          as.integer(res$b_valid %||% 0L),
          as.integer(res$B_max_effective %||% res$B_max %||% B_max_int),
          res$p_lower %||% NA_real_,
          res$p_upper %||% NA_real_,
          es_tag,
          if (isTRUE(res$significant)) "SIGNIFICANT (split)" else "NOT significant (certify)",
          node_elapsed))
      }
    }

    if (isTRUE(res$significant)) {
      visit(node$left_child_id)
      visit(node$right_child_id)
    } else if (identical(res$reason_skipped, "spatial_fit_failed")) {
      # Not a domain: the test could not be done, so neither this node nor anything below it
      # was tested.
      untested_nodes <<- c(untested_nodes, node_id)
    } else {
      certified_nodes <<- c(certified_nodes, node_id)
    }
    invisible(NULL)
  }

  visit(observed_tree$root_id)

  if (verbose) {
    total_elapsed <- proc.time()[3] - traverse_t0
    message(sprintf(
      "\n  [traverse] All node tests complete: %d visited, %d certified [%.1fs total]",
      length(visited), length(certified_nodes), total_elapsed))
  }

  rows <- vector("list", length(internal_nodes))
  for (i in seq_along(internal_nodes)) {
    nid <- internal_nodes[i]
    node <- observed_tree$node_table[nid, , drop = FALSE]
    n_v <- length(observed_tree$spot_indices[[nid]])
    alpha_star <- alpha * ((n_v - 1) / (n_total - 1))
    row_budget <- node_budget(alpha_star)
    res <- tested_map[[as.character(nid)]]
    if (is.null(res)) {
      rows[[i]] <- data.frame(
        node_id = nid,
        parent_id = node$parent_id,
        left_child_id = node$left_child_id,
        right_child_id = node$right_child_id,
        n_v = n_v,
        Kv = observed_tree$node_table$leaf_count[nid],
        raw_p = NA_real_,
        p_lower = NA_real_,
        p_upper = NA_real_,
        b_valid = 0L,
        early_stop_reason = NA_character_,
        alpha_star = alpha_star,
        significant = FALSE,
        tested = FALSE,
        reason_skipped = "ancestor_not_significant",
        n_valid_reps = 0L,
        n_failed_reps = 0L,
        B_max_configured = row_budget$configured,
        B_max_required = row_budget$required,
        B_max_effective = row_budget$effective,
        bootstrap_B_max_cap = row_budget$cap,
        bootstrap_resolution_sufficient =
          row_budget$resolution_sufficient,
        minimum_attainable_p = row_budget$minimum_attainable_p,
        pvalue_method = pvalue_method,
        normal_n_reps = NA_integer_,
        normal_used_refine = NA,
        normal_mu = NA_real_,
        normal_sigma = NA_real_,
        normal_p_initial = NA_real_,
        spatial_backend = NA_character_,
        spatial_backend_reason = NA_character_,
        stringsAsFactors = FALSE
      )
    } else {
      backend_info <- res$diagnostics$spatial_backend %||% list()
      rows[[i]] <- data.frame(
        node_id = nid,
        parent_id = node$parent_id,
        left_child_id = node$left_child_id,
        right_child_id = node$right_child_id,
        n_v = res$n_v %||% n_v,
        Kv = res$Kv %||% observed_tree$node_table$leaf_count[nid],
        raw_p = res$raw_p %||% NA_real_,
        p_lower = res$p_lower %||% NA_real_,
        p_upper = res$p_upper %||% NA_real_,
        b_valid = as.integer(res$b_valid %||% 0L),
        early_stop_reason = res$early_stop_reason %||% NA_character_,
        alpha_star = res$alpha_star %||% alpha_star,
        significant = isTRUE(res$significant),
        tested = isTRUE(res$tested),
        reason_skipped = res$reason_skipped %||% NA_character_,
        n_valid_reps = as.integer(res$n_valid_reps %||% 0L),
        n_failed_reps = as.integer(res$n_failed_reps %||% 0L),
        B_max_configured =
          res$B_max_configured %||% row_budget$configured,
        B_max_required = res$B_max_required %||% row_budget$required,
        B_max_effective = res$B_max_effective %||% row_budget$effective,
        bootstrap_B_max_cap =
          res$bootstrap_B_max_cap %||% row_budget$cap,
        bootstrap_resolution_sufficient =
          res$bootstrap_resolution_sufficient %||%
            row_budget$resolution_sufficient,
        minimum_attainable_p =
          res$minimum_attainable_p %||% row_budget$minimum_attainable_p,
        pvalue_method = res$pvalue_method %||% pvalue_method,
        normal_n_reps = as.integer(res$normal_n_reps %||% NA_integer_),
        normal_used_refine = res$normal_used_refine %||% NA,
        normal_mu = res$normal_mu %||% NA_real_,
        normal_sigma = res$normal_sigma %||% NA_real_,
        normal_p_initial = res$normal_p_initial %||% NA_real_,
        spatial_backend = as.character(
          backend_info$selected_method %||% backend_info$backend_type %||%
            NA_character_
        ),
        spatial_backend_reason = as.character(
          backend_info$fallback_reason %||% NA_character_
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  node_table <- do.call(rbind, rows)
  node_table <- node_table[order(node_table$node_id), , drop = FALSE]

  certified_nodes <- sort(unique(certified_nodes))
  certified_regions <- lapply(certified_nodes, function(nid) {
    list(node_id = nid, spot_indices = observed_tree$spot_indices[[nid]])
  })
  # Package: no node is certified when the spatial fit fails at the root (v0.2.5.8 stopped the run then);
  # paste0() of a zero-length vector returns "node_", so name only when non-empty.
  names(certified_regions) <- if (length(certified_nodes)) paste0("node_", certified_nodes) else character(0)
  untested_nodes <- sort(unique(untested_nodes))
  untested_regions <- lapply(untested_nodes, function(nid) {
    list(node_id = nid, spot_indices = observed_tree$spot_indices[[nid]])
  })
  # paste0() of a zero-length vector returns "node_", not character(0); name only when non-empty.
  names(untested_regions) <- if (length(untested_nodes)) paste0("node_", untested_nodes) else character(0)

  list(
    node_table = node_table,
    certified_nodes = certified_nodes,
    certified_regions = certified_regions,
    untested_nodes = untested_nodes,
    untested_regions = untested_regions,
    tested_node_ids = sort(unique(visited)),
    node_results = tested_map
  )
}

# The SigDoTS engine (internal; the user-facing function is sigdots()).
# 1) builds the observed label-level tree: initial clustering into K0 groups (k-means on the
#    first ci_pca_dims principal components by default), retained-component adjacency,
#    constrained Ward merging of group centroids in the same principal components;
# 2) tests its internal nodes top-down; at each node it fits the node-local null (the within-group correlogram fit for
#    eta and ell, leaves_gls for the gene side), simulates null replicates, clusters each
#    into Kv groups with `cluster_method`, and compares the root split's cluster index (on
#    the node's full expression matrix) with the observed one, sequentially;
# 3) the per-node level is alpha * (n_v - 1) / (n - 1).
# Arguments as in v0.2.5.8, less the options that are not part of the package: B,
# time_profile, profile, script_dir, core_script_path.
# `null_args` keeps r_cap, ridge, nugget, spatial_method ("exact" or "vecchia"),
# max_exact_peak_gib, vecchia_neighbors, eta_min, correlogram_pairs_per_bin and
# correlogram_min_pairs.
# Returns list(observed_tree, node_table, certified_regions, node_results, settings,
# diagnostics).
spa_sigtest_multicluster <- function(
  X,
  coords,
  K0,
  B_min = 20L,
  B_max = 120L,
  bootstrap_B_max_cap = 5000L,
  alpha = 0.05,
  k_nn = 6L,
  n_comp_min = 0.05,
  ci_pca_dims = 30L,
  min_node_size = 20L,
  min_child_size = 10L,
  seed = 1L,
  cluster_method = "kmeans",
  bayesspace_args = list(),
  null_args = list(),
  verbose = TRUE,
  n_cores = 1L,
  observed_clustering = NULL,
  initial_seed = seed,
  progress = FALSE,
  pvalue_method = c("bootstrap", "normal_adaptive"),
  normal_n_initial = 10L,
  normal_n_refine = 50L,
  normal_band_multiplier = 10
) {
  pvalue_method <- match.arg(pvalue_method)
  B_min <- as.integer(B_min)
  B_max <- as.integer(B_max)
  if (B_min < 1L) B_min <- 1L
  if (B_max < B_min) B_max <- B_min
  bootstrap_B_max_cap <- as.numeric(bootstrap_B_max_cap)
  if (identical(pvalue_method, "bootstrap") &&
      (length(bootstrap_B_max_cap) != 1L ||
       !is.finite(bootstrap_B_max_cap) ||
       bootstrap_B_max_cap < B_max)) {
    stop("bootstrap_B_max_cap must be finite and >= B_max.")
  }
  bootstrap_B_max_cap <- as.integer(floor(bootstrap_B_max_cap))

  X <- as.matrix(X)
  coords <- as.matrix(coords)
  n <- nrow(X)
  if (nrow(coords) != n) stop("coords and X must have the same number of rows.")
  if (K0 < 2L) stop("K0 must be >= 2.")
  if (!is.finite(n_comp_min) || n_comp_min <= 0 || n_comp_min > 1) {
    stop("n_comp_min must be in (0, 1].")
  }
  observed_clustering_precomputed <- !is.null(observed_clustering)
  observed_clustering_result <- if (observed_clustering_precomputed) {
    .normalize_observed_clustering(
      observed_clustering = observed_clustering,
      n = n,
      K0 = K0
    )
  } else {
    NULL
  }
  if (!is.null(seed)) set.seed(seed)

  n_cores <- max(1L, as.integer(n_cores))
  cluster_method <- match.arg(cluster_method, choices = c("bayesspace", "kmeans", "banksy"))
  if (verbose) {
    message(sprintf(
      paste0(
        "[spaSigClust-HC] n=%d, p=%d, K0=%d, baseline B in [%d, %d], node-specific cap=%d, ",
        "cluster=%s, n_cores=%d"
      ),
      n, ncol(X), as.integer(K0), B_min, B_max, bootstrap_B_max_cap,
      cluster_method, n_cores
    ))
  }

  t0 <- proc.time()[3]
  obs_tree <- build_observed_global_tree(
    X = X,
    coords = coords,
    K0 = as.integer(K0),
    k_nn = as.integer(k_nn),
    ci_pca_dims = as.integer(ci_pca_dims),
    n_comp_min = as.numeric(n_comp_min),
    seed = as.integer(initial_seed),
    cluster_method = cluster_method,
    bayesspace_args = bayesspace_args,
    observed_clustering = observed_clustering_result
  )

  if (verbose) {
    message(sprintf(
      "[spaSigClust-HC] Observed tree built: %d nodes, %d leaves",
      nrow(obs_tree$node_table), obs_tree$n_leaves
    ))
    if (!isTRUE(obs_tree$component_filter_used_for_tree)) {
      message(sprintf(
        "[spaSigClust-HC] Observed adjacency: %s",
        obs_tree$component_filter_message
      ))
    }
  }
  if (progress && !verbose) {
    message(sprintf("SigDoTS: %d spots, %d genes; %d candidate domains; testing from the top of the tree",
                    n, ncol(X), as.integer(K0)))
  }

  trav <- traverse_tree_topdown(
    observed_tree = obs_tree,
    X = X,
    coords = coords,
    B_min = B_min,
    B_max = B_max,
    bootstrap_B_max_cap = bootstrap_B_max_cap,
    alpha = alpha,
    ci_pca_dims = as.integer(ci_pca_dims),
    min_node_size = as.integer(min_node_size),
    min_child_size = as.integer(min_child_size),
    seed = as.integer(seed),
    n_comp_min = as.numeric(n_comp_min),
    cluster_method = cluster_method,
    bayesspace_args = bayesspace_args,
    null_args = null_args,
    verbose = verbose,
    n_cores = n_cores,
    progress = progress,
    pvalue_method = pvalue_method,
    normal_n_initial = as.integer(normal_n_initial),
    normal_n_refine = as.integer(normal_n_refine),
    normal_band_multiplier = as.numeric(normal_band_multiplier)
  )
  elapsed <- proc.time()[3] - t0

  diagnostics <- list(
    elapsed_seconds = elapsed,
    tested_nodes = trav$tested_node_ids,
    certified_node_ids = trav$certified_nodes,
    untested_node_ids = trav$untested_nodes,
    bayesspace_method_observed = obs_tree$bayesspace$method %||% "unknown",
    observed_leaf_count = as.integer(obs_tree$node_table$leaf_count[obs_tree$root_id]),
    observed_clustering_precomputed = observed_clustering_precomputed,
    observed_tree_adjacency = .summarize_tree_adjacency(obs_tree)
  )
  settings <- list(
    K0 = as.integer(K0),
    B_min = B_min,
    B_max = B_max,
    bootstrap_B_max_cap = bootstrap_B_max_cap,
    bootstrap_B_max_policy = paste(
      "max(configured B_max, ceiling(1 / alpha_star) - 1),",
      "capped at bootstrap_B_max_cap"
    ),
    alpha = alpha,
    k_nn = as.integer(k_nn),
    n_comp_min = as.numeric(n_comp_min),
    ci_pca_dims = as.integer(ci_pca_dims),
    min_node_size = as.integer(min_node_size),
    min_child_size = as.integer(min_child_size),
    seed = as.integer(seed),
    initial_seed = as.integer(initial_seed),
    cluster_method = cluster_method,
    tree_adjacency_policy = paste(
      "retained components after filtering when connected;",
      "otherwise all components without filtering"
    ),
    n_cores = n_cores,
    pvalue_method = pvalue_method,
    normal_n_initial = as.integer(normal_n_initial),
    normal_n_refine = as.integer(normal_n_refine),
    normal_band_multiplier = as.numeric(normal_band_multiplier),
    observed_clustering_precomputed = observed_clustering_precomputed,
    bayesspace_args = bayesspace_args,
    null_args = null_args
  )

  list(
    observed_tree = obs_tree,
    node_table = trav$node_table,
    certified_regions = trav$certified_regions,
    node_results = trav$node_results,
    untested_regions = trav$untested_regions,
    settings = settings,
    diagnostics = diagnostics
  )
}
# nolint end
