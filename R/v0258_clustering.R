# spaSigClust v0.2.5.8HC -- Clustering Utilities
# Thin adapters for deterministic clustering calls on spot-level matrices.
# Supported backends: "bayesspace", "kmeans" (more can be added).

# Encode a regular square lattice once, then reuse the same levels for every
# node-local subset. This preserves gaps instead of re-ranking each subset.
configure_bayesspace_reference_grid <- function(
  bayesspace_args = list(),
  coords,
  platform = "ST"
) {
  coords <- as.matrix(coords)
  if (ncol(coords) < 2L || any(!is.finite(coords[, 1:2]))) {
    stop("coords must contain two finite coordinate columns.")
  }
  bayesspace_args$platform <- platform
  bayesspace_args$array_grid <- list(
    x = sort(unique(as.numeric(coords[, 1]))),
    y = sort(unique(as.numeric(coords[, 2])))
  )
  bayesspace_args$require_spatial_neighbors <- TRUE
  bayesspace_args
}

.match_reference_levels <- function(values, reference, axis_name) {
  reference <- as.numeric(reference)
  if (length(reference) < 1L || any(!is.finite(reference))) {
    stop("Invalid BayesSpace reference levels for ", axis_name, ".")
  }
  matched <- match(as.numeric(values), reference)
  if (anyNA(matched)) {
    stop(
      "Some ", axis_name,
      " coordinates are absent from the global BayesSpace reference grid."
    )
  }
  as.integer(matched)
}

.bayesspace_array_coordinates <- function(coords, array_grid = NULL) {
  coords <- as.matrix(coords)
  if (is.null(array_grid)) {
    # Backward-compatible path for genuine platform array coordinates.
    return(list(
      array_row = as.numeric(coords[, 1]),
      array_col = as.numeric(coords[, 2]),
      coordinate_mode = "native"
    ))
  }
  if (!all(c("x", "y") %in% names(array_grid))) {
    stop("BayesSpace array_grid must contain global `x` and `y` levels.")
  }
  list(
    array_row = .match_reference_levels(coords[, 2], array_grid$y, "y"),
    array_col = .match_reference_levels(coords[, 1], array_grid$x, "x"),
    coordinate_mode = "reference_grid"
  )
}

.bayesspace_neighbor_offsets <- function(platform) {
  if (identical(platform, "Visium")) {
    return(cbind(
      col = c(-2, 2, -1, 1, -1, 1),
      row = c(0, 0, -1, -1, 1, 1)
    ))
  }
  if (platform %in% c("ST", "VisiumHD")) {
    return(cbind(
      col = c(0, 1, 0, -1),
      row = c(-1, 0, 1, 0)
    ))
  }
  stop("Unsupported BayesSpace platform: ", platform)
}

.bayesspace_neighbor_diagnostics <- function(array_row, array_col, platform) {
  coordinate_key <- function(col, row) paste(col, row, sep = ":")
  keys <- coordinate_key(array_col, array_row)
  if (anyDuplicated(keys)) {
    stop("BayesSpace array coordinates must be unique.")
  }
  offsets <- .bayesspace_neighbor_offsets(platform)
  neighbor_counts <- integer(length(keys))
  for (i in seq_len(nrow(offsets))) {
    candidate <- coordinate_key(
      array_col + offsets[i, "col"],
      array_row + offsets[i, "row"]
    )
    neighbor_counts <- neighbor_counts + as.integer(candidate %in% keys)
  }
  n_with_neighbors <- sum(neighbor_counts > 0L)
  list(
    platform = platform,
    n_spots = length(keys),
    n_with_neighbors = n_with_neighbors,
    neighbor_fraction = n_with_neighbors / length(keys),
    n_undirected_edges = sum(neighbor_counts) / 2,
    min_neighbors = min(neighbor_counts),
    median_neighbors = stats::median(neighbor_counts),
    max_neighbors = max(neighbor_counts)
  )
}

# Build a minimal SCE object with coordinate metadata and PCA embedding.
make_sce_for_bayesspace <- function(
  X,
  coords,
  pca_dims = 15L,
  platform = "Visium",
  array_grid = NULL,
  require_spatial_neighbors = FALSE
) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment is required for BayesSpace integration.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment is required for BayesSpace integration.")
  }
  X <- as.matrix(X)
  coords <- as.matrix(coords)
  n <- nrow(X)
  p <- ncol(X)
  if (n < 2L || p < 2L) stop("X must have at least 2 spots and 2 features.")
  if (nrow(coords) != n) stop("coords and X must have matching number of rows.")

  gene_names <- paste0("feature_", seq_len(p))
  spot_names <- paste0("spot_", seq_len(n))
  d <- min(as.integer(pca_dims), p, n - 1L)
  if (d < 1L) stop("Invalid PCA dimension.")

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = t(X))
  )
  rownames(sce) <- gene_names
  colnames(sce) <- spot_names

  array_coords <- .bayesspace_array_coordinates(
    coords = coords,
    array_grid = array_grid
  )
  neighbor_diagnostics <- .bayesspace_neighbor_diagnostics(
    array_row = array_coords$array_row,
    array_col = array_coords$array_col,
    platform = platform
  )
  if (isTRUE(require_spatial_neighbors) &&
      neighbor_diagnostics$n_with_neighbors == 0L) {
    stop(
      "BayesSpace array coordinates produce zero spatial neighbors. ",
      "Check `platform` and `array_grid`."
    )
  }

  cd <- SummarizedExperiment::colData(sce)
  cd$row       <- array_coords$array_row
  cd$col       <- array_coords$array_col
  cd$imagerow  <- as.numeric(coords[, 1])
  cd$imagecol  <- as.numeric(coords[, 2])
  cd$array_row <- array_coords$array_row
  cd$array_col <- array_coords$array_col
  SummarizedExperiment::colData(sce) <- cd

  pr <- stats::prcomp(X, center = TRUE, scale. = FALSE, rank. = d)
  pca_mat <- pr$x[, seq_len(d), drop = FALSE]
  SingleCellExperiment::reducedDim(sce, "PCA") <- pca_mat
  attr(sce, "pca_dims") <- d
  attr(sce, "bayesspace_coordinate_mode") <- array_coords$coordinate_mode
  attr(sce, "bayesspace_neighbor_diagnostics") <- neighbor_diagnostics
  sce
}

# Internal deterministic fallback if BayesSpace is unavailable or fails.
.fallback_kmeans_labels <- function(X, q, pca_dims, seed) {
  d <- min(as.integer(pca_dims), ncol(X), nrow(X) - 1L)
  pr <- stats::prcomp(X, center = TRUE, scale. = FALSE, rank. = d)
  Z <- pr$x[, seq_len(d), drop = FALSE]
  if (!is.null(seed)) set.seed(seed)
  km <- stats::kmeans(Z, centers = as.integer(q), nstart = 10, iter.max = 100)
  as.integer(km$cluster)
}

# Run BayesSpace with fixed q and deterministic seed behavior.
# A missing package or a failed BayesSpace run stops with an error; no other method is substituted.
# Returns stable integer labels with optional diagnostics.
run_bayesspace_labels <- function(
  X,
  coords,
  q,
  seed = 1L,
  pca_dims = 15L,
  bayesspace_args = list()
) {
  if (q < 2L) stop("BayesSpace target q must be >= 2.")

  have_bs <- requireNamespace("BayesSpace", quietly = TRUE)
  have_sce <- requireNamespace("SingleCellExperiment", quietly = TRUE)

  if (!have_bs || !have_sce) {
    stop("BayesSpace and SingleCellExperiment are required.")
  }

  args <- bayesspace_args
  platform <- args$platform %||% "Visium"
  init_method <- args$init.method %||% "mclust"
  model <- args$model %||% "t"
  nrep <- args$nrep %||% 5000L
  burn_in <- args$burn.in %||% 100L
  gamma <- args$gamma %||% 3
  save_chain <- isTRUE(args$save.chain)

  sce <- make_sce_for_bayesspace(
    X = X,
    coords = coords,
    pca_dims = pca_dims,
    platform = platform,
    array_grid = args$array_grid %||% NULL,
    require_spatial_neighbors = isTRUE(args$require_spatial_neighbors)
  )
  d <- attr(sce, "pca_dims")
  coordinate_mode <- attr(sce, "bayesspace_coordinate_mode")
  neighbor_diagnostics <- attr(sce, "bayesspace_neighbor_diagnostics")

  if (!is.null(seed)) set.seed(seed)
  out <- tryCatch(
    {
      BayesSpace::spatialCluster(
        sce,
        q = as.integer(q),
        d = as.integer(d),
        platform = platform,
        init.method = init_method,
        model = model,
        nrep = as.integer(nrep),
        burn.in = as.integer(burn_in),
        gamma = gamma,
        save.chain = save_chain
      )
    },
    error = function(e) e
  )

  if (inherits(out, "error")) {
    stop(sprintf("BayesSpace::spatialCluster failed: %s", out$message))
  }

  if (is.null(SummarizedExperiment::colData(out)$spatial.cluster)) {
    stop("BayesSpace output did not contain `spatial.cluster` labels.")
  }

  list(
    labels = as.integer(SummarizedExperiment::colData(out)$spatial.cluster),
    method = "bayesspace",
    diagnostics = list(
      platform = platform,
      init_method = init_method,
      model = model,
      coordinate_mode = coordinate_mode,
      neighbors = neighbor_diagnostics
    )
  )
}

# ============================================================
# BANKSY clustering backend
# ============================================================
run_banksy_labels <- function(
  X,
  coords,
  q,
  seed = 1L,
  pca_dims = 15L,
  banksy_args = list()
) {
  if (q < 2L) stop("BANKSY target q must be >= 2.")

  have_banksy <- requireNamespace("Banksy", quietly = TRUE)
  have_sce <- requireNamespace("SingleCellExperiment", quietly = TRUE)

  if (!have_banksy || !have_sce) {
    stop("Banksy and SingleCellExperiment are required.")
  }

  lambda   <- banksy_args$lambda   %||% 0.2
  k_geom   <- banksy_args$k_geom   %||% 15L
  use_agf  <- isTRUE(banksy_args$use_agf)
  npcs     <- as.integer(banksy_args$npcs %||% min(pca_dims, 20L))

  X <- as.matrix(X)
  coords <- as.matrix(coords)
  n <- nrow(X); p <- ncol(X)
  spot_names <- paste0("spot_", seq_len(n))
  gene_names <- if (!is.null(colnames(X))) colnames(X) else paste0("feature_", seq_len(p))
  rownames(X) <- spot_names

  se <- SingleCellExperiment::SingleCellExperiment(
    assays = list(normcounts = t(X))
  )
  colnames(se) <- spot_names
  rownames(se) <- gene_names
  cd <- SummarizedExperiment::colData(se)
  cd$x <- as.numeric(coords[, 1])
  cd$y <- as.numeric(coords[, 2])
  SummarizedExperiment::colData(se) <- cd

  # Package: Banksy prints a line per call ("Using seed=..."), once per simulated dataset;
  # silenced, so that sigdots()'s own progress lines stay readable.
  out <- tryCatch(suppressMessages({
    se <- Banksy::computeBanksy(se, assay_name = "normcounts",
                                 coord_names = c("x", "y"),
                                 k_geom = as.integer(k_geom),
                                 compute_agf = use_agf, verbose = FALSE)
    se <- Banksy::runBanksyPCA(se, use_agf = use_agf, lambda = lambda,
                                npcs = npcs, assay_name = "normcounts",
                                seed = seed)
    se <- Banksy::clusterBanksy(se, use_agf = use_agf, lambda = lambda,
                                 use_pcs = TRUE, npcs = npcs,
                                 assay_name = "normcounts",
                                 algo = "kmeans",
                                 kmeans.centers = as.integer(q),
                                 seed = seed)
    cnames <- grep("^clust_", colnames(SummarizedExperiment::colData(se)), value = TRUE)
    if (length(cnames) == 0L) stop("No cluster column found in BANKSY output.")
    as.integer(SummarizedExperiment::colData(se)[[cnames[1]]])
  }), error = function(e) e)

  if (inherits(out, "error")) {
    stop(sprintf("Banksy clustering failed: %s", out$message))
  }

  list(
    labels = out,
    method = "banksy",
    diagnostics = list(lambda = lambda, k_geom = k_geom,
                       use_agf = use_agf, npcs = npcs)
  )
}

# ============================================================
# Unified clustering dispatcher
# ============================================================
# cluster_method: "kmeans" (the package's default), "bayesspace", or "banksy".
# Returns list(labels, method, diagnostics) regardless of backend.
run_clustering <- function(
  X,
  coords,
  q,
  cluster_method = "kmeans",
  seed = 1L,
  pca_dims = 15L,
  bayesspace_args = list()
) {
  cluster_method <- match.arg(cluster_method, choices = c("bayesspace", "kmeans", "banksy"))

  if (cluster_method == "kmeans") {
    labs <- .fallback_kmeans_labels(X = X, q = q, pca_dims = pca_dims, seed = seed)
    return(list(
      labels = labs,
      method = "kmeans",
      diagnostics = list(nstart = 10L, pca_dims = pca_dims)
    ))
  }

  if (cluster_method == "banksy") {
    return(run_banksy_labels(
      X = X, coords = coords, q = q, seed = seed,
      pca_dims = pca_dims,
      banksy_args = bayesspace_args
    ))
  }

  run_bayesspace_labels(
    X = X, coords = coords, q = q, seed = seed,
    pca_dims = pca_dims, bayesspace_args = bayesspace_args
  )
}
