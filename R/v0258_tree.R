# spaSigClust v0.2.5.8HC -- Tree Utilities (label-level)
# Multi-cluster SHC-style utilities:
# kNN spot graph, retained-component label adjacency, label-level constrained HC,
# and tree bookkeeping helpers.

# Build an undirected symmetric kNN graph on spatial coordinates.
build_knn_graph <- function(coords, k_nn = 6L) {
  coords <- as.matrix(coords)
  n <- nrow(coords)
  if (n < 2L) stop("Need at least 2 spots to build a kNN graph.")
  if (ncol(coords) < 2L) stop("coords must have at least 2 columns.")
  k_eff <- min(as.integer(k_nn), n - 1L)
  if (k_eff < 1L) stop("k_nn must be >= 1.")

  nn <- RANN::nn2(coords, k = k_eff + 1L)
  idx <- nn$nn.idx[, -1, drop = FALSE]
  i <- rep(seq_len(n), each = k_eff)
  j <- as.vector(t(idx))
  adj <- Matrix::sparseMatrix(i = i, j = j, x = 1L, dims = c(n, n))
  adj <- (adj + Matrix::t(adj)) > 0
  Matrix::diag(adj) <- FALSE

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  list(adj = adj, graph = g, n = n, k_nn = k_eff)
}

# Build an induced subgraph object from selected global spot indices.
subset_spot_graph <- function(spot_graph, spot_indices) {
  idx <- sort(unique(as.integer(spot_indices)))
  if (length(idx) < 2L) stop("subset_spot_graph requires at least 2 spots.")
  adj_sub <- spot_graph$adj[idx, idx, drop = FALSE]
  g_sub <- igraph::graph_from_adjacency_matrix(adj_sub, mode = "undirected", diag = FALSE)
  list(
    adj = adj_sub,
    graph = g_sub,
    spot_indices_global = idx
  )
}

# Connected components per BayesSpace label and retained-component filtering.
split_labels_into_cc <- function(labels, spot_graph, n_comp_min = 0.05) {
  labels <- as.integer(labels)
  n <- length(labels)
  if (is.null(spot_graph$adj)) stop("spot_graph must contain sparse adjacency matrix `adj`.")
  if (n != nrow(spot_graph$adj)) stop("labels length must match spot graph size.")
  if (!is.finite(n_comp_min) || n_comp_min <= 0 || n_comp_min > 1) {
    stop("n_comp_min must be in (0, 1].")
  }

  label_values <- sort(unique(labels))
  component_members <- list()
  retained_members_by_label <- setNames(vector("list", length(label_values)), as.character(label_values))
  retained_spot_mask <- rep(FALSE, n)
  rows <- list()
  cc_global_id <- 0L

  for (lab in label_values) {
    lab_idx <- which(labels == lab)
    n_lab <- length(lab_idx)
    if (n_lab == 0L) next

    if (n_lab == 1L) {
      comp <- 1L
      cc_global_id <- cc_global_id + 1L
      comp_size <- 1L
      prop <- 1
      retained <- prop >= n_comp_min
      component_members[[cc_global_id]] <- lab_idx
      rows[[length(rows) + 1L]] <- data.frame(
        label = lab,
        cc_id_within_label = comp,
        cc_global_id = cc_global_id,
        size = comp_size,
        prop_within_label = prop,
        retained = retained,
        stringsAsFactors = FALSE
      )
      if (retained) {
        retained_members_by_label[[as.character(lab)]] <- list(lab_idx)
        retained_spot_mask[lab_idx] <- TRUE
      }
      next
    }

    adj_sub <- spot_graph$adj[lab_idx, lab_idx, drop = FALSE]
    g_sub <- igraph::graph_from_adjacency_matrix(adj_sub, mode = "undirected", diag = FALSE)
    comp_membership <- igraph::components(g_sub)$membership
    cc_ids <- sort(unique(comp_membership))
    retained_list <- list()

    for (cc in cc_ids) {
      local_idx <- which(comp_membership == cc)
      global_idx <- lab_idx[local_idx]
      comp_size <- length(global_idx)
      prop <- comp_size / n_lab
      retained <- prop >= n_comp_min
      cc_global_id <- cc_global_id + 1L

      component_members[[cc_global_id]] <- global_idx
      rows[[length(rows) + 1L]] <- data.frame(
        label = lab,
        cc_id_within_label = as.integer(cc),
        cc_global_id = cc_global_id,
        size = comp_size,
        prop_within_label = prop,
        retained = retained,
        stringsAsFactors = FALSE
      )
      if (retained) {
        retained_list[[length(retained_list) + 1L]] <- global_idx
        retained_spot_mask[global_idx] <- TRUE
      }
    }
    retained_members_by_label[[as.character(lab)]] <- retained_list
  }

  component_table <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()
  if (nrow(component_table) > 0L) {
    component_table$label <- as.integer(component_table$label)
    component_table$cc_id_within_label <- as.integer(component_table$cc_id_within_label)
    component_table$cc_global_id <- as.integer(component_table$cc_global_id)
    component_table$size <- as.integer(component_table$size)
    component_table$prop_within_label <- as.numeric(component_table$prop_within_label)
    component_table$retained <- as.logical(component_table$retained)
  }

  list(
    label_values = label_values,
    component_table = component_table,
    component_members = component_members,
    retained_members_by_label = retained_members_by_label,
    retained_spot_mask = retained_spot_mask
  )
}

# Label-wise spot membership using BayesSpace labels as leaves.
build_label_membership <- function(labels, K_expected = NULL) {
  labels <- as.integer(labels)
  label_values <- sort(unique(labels))
  if (!is.null(K_expected)) {
    exp_vals <- seq_len(as.integer(K_expected))
    if (length(label_values) != length(exp_vals) || any(label_values != exp_vals)) {
      stop(
        sprintf(
          "BayesSpace labels must be exactly 1..K_expected (expected %d labels; got: %s).",
          as.integer(K_expected),
          paste(label_values, collapse = ",")
        )
      )
    }
  }

  members <- lapply(label_values, function(lab) which(labels == lab))
  names(members) <- as.character(label_values)
  list(label_values = label_values, label_members = members)
}

# Build label-label adjacency using retained components.
build_label_adjacency_from_retained <- function(labels, spot_graph, cc_info) {
  labels <- as.integer(labels)
  label_values <- cc_info$label_values
  K <- length(label_values)
  label_to_pos <- setNames(seq_len(K), as.character(label_values))
  retained_mask <- cc_info$retained_spot_mask

  if (!any(retained_mask)) {
    # Empty K x K logical sparse matrix: K vertices, no edges.
    # `Matrix::Matrix(FALSE, K, K, sparse = TRUE)` would return an ldiMatrix
    # (all-FALSE is trivially diagonal), which igraph 1.5.1 cannot convert --
    # graph.adjacency.sparse() aborts with "argument is of length zero".
    # sparseMatrix() yields lgCMatrix, matching this function's non-empty
    # return class, so is_label_adjacency_connected() correctly reports FALSE
    # for K >= 2 and the caller falls back to the unfiltered adjacency.
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = logical(0), dims = c(K, K)
    ))
  }

  e <- Matrix::summary(Matrix::triu(spot_graph$adj, k = 1))
  keep <- retained_mask[e$i] & retained_mask[e$j] & (labels[e$i] != labels[e$j])
  if (!any(keep)) {
    # Empty K x K logical sparse matrix: K vertices, no edges.
    # `Matrix::Matrix(FALSE, K, K, sparse = TRUE)` would return an ldiMatrix
    # (all-FALSE is trivially diagonal), which igraph 1.5.1 cannot convert --
    # graph.adjacency.sparse() aborts with "argument is of length zero".
    # sparseMatrix() yields lgCMatrix, matching this function's non-empty
    # return class, so is_label_adjacency_connected() correctly reports FALSE
    # for K >= 2 and the caller falls back to the unfiltered adjacency.
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = logical(0), dims = c(K, K)
    ))
  }

  li <- label_to_pos[as.character(labels[e$i[keep]])]
  lj <- label_to_pos[as.character(labels[e$j[keep]])]
  adj <- Matrix::sparseMatrix(i = li, j = lj, x = 1L, dims = c(K, K))
  adj <- (adj + Matrix::t(adj)) > 0
  Matrix::diag(adj) <- FALSE
  adj
}

# Build label-label adjacency from every cross-label spatial contact, without
# applying the connected-component size filter.
build_label_adjacency_unfiltered <- function(
  labels,
  spot_graph,
  label_values = sort(unique(as.integer(labels)))
) {
  labels <- as.integer(labels)
  label_values <- as.integer(label_values)
  K <- length(label_values)
  if (length(labels) != nrow(spot_graph$adj)) {
    stop("labels length must match spot graph size.")
  }
  label_to_pos <- setNames(seq_len(K), as.character(label_values))
  e <- Matrix::summary(Matrix::triu(spot_graph$adj, k = 1))
  keep <- labels[e$i] != labels[e$j]
  if (!any(keep)) {
    # Empty K x K logical sparse matrix: K vertices, no edges.
    # `Matrix::Matrix(FALSE, K, K, sparse = TRUE)` would return an ldiMatrix
    # (all-FALSE is trivially diagonal), which igraph 1.5.1 cannot convert --
    # graph.adjacency.sparse() aborts with "argument is of length zero".
    # sparseMatrix() yields lgCMatrix, matching this function's non-empty
    # return class, so is_label_adjacency_connected() correctly reports FALSE
    # for K >= 2 and the caller falls back to the unfiltered adjacency.
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = logical(0), dims = c(K, K)
    ))
  }

  li <- label_to_pos[as.character(labels[e$i[keep]])]
  lj <- label_to_pos[as.character(labels[e$j[keep]])]
  if (anyNA(li) || anyNA(lj)) {
    stop("Encountered labels outside label_values while building adjacency.")
  }
  adj <- Matrix::sparseMatrix(i = li, j = lj, x = 1L, dims = c(K, K))
  adj <- (adj + Matrix::t(adj)) > 0
  Matrix::diag(adj) <- FALSE
  adj
}

is_label_adjacency_connected <- function(label_adjacency) {
  K <- nrow(label_adjacency)
  if (K != ncol(label_adjacency)) {
    stop("label_adjacency must be square.")
  }
  if (K <= 1L) return(TRUE)
  graph <- igraph::graph_from_adjacency_matrix(
    label_adjacency,
    mode = "undirected",
    diag = FALSE
  )
  igraph::is_connected(graph)
}

.label_adjacency_edge_count <- function(label_adjacency) {
  as.integer(sum(label_adjacency != 0) / 2)
}

# Compute group centroids in node-level PCA space (frozen representation).
compute_region_centroids_pca <- function(X, region_members, pca_dims = 30L) {
  X <- as.matrix(X)
  if (nrow(X) < 2L) stop("Need at least 2 rows in X.")
  d <- min(as.integer(pca_dims), ncol(X), nrow(X) - 1L)
  if (d < 1L) stop("Invalid PCA dimension; check node size and pca_dims.")

  pr <- stats::prcomp(X, center = TRUE, scale. = FALSE, rank. = d)
  Z <- pr$x[, seq_len(d), drop = FALSE]

  R <- length(region_members)
  cent <- matrix(NA_real_, nrow = R, ncol = d)
  for (r in seq_len(R)) {
    idx <- region_members[[r]]
    cent[r, ] <- colMeans(Z[idx, , drop = FALSE])
  }
  rownames(cent) <- paste0("region_", seq_len(R))

  list(centroids = cent, Z = Z, pca = pr)
}

# Internal helper: compute constrained Ward merge score between two groups.
.ward_delta <- function(size_a, size_b, centroid_a, centroid_b) {
  size_a <- as.numeric(size_a)
  size_b <- as.numeric(size_b)
  w <- (size_a * size_b) / (size_a + size_b)
  d2 <- sum((centroid_a - centroid_b)^2)
  w * d2
}

# Build a contiguous binary tree from initial groups and adjacency.
adjacency_constrained_ward_tree <- function(region_members, region_adj, region_centroids, leaf_values = NULL) {
  R <- length(region_members)
  if (R < 2L) stop("Need at least two leaves to build a binary tree.")
  if (nrow(region_adj) != R || ncol(region_adj) != R) stop("region_adj shape mismatch.")
  if (nrow(region_centroids) != R) stop("region_centroids row count must equal number of groups.")

  max_nodes <- 2L * R - 1L
  members <- vector("list", max_nodes)
  sizes <- integer(max_nodes)
  centroids <- matrix(NA_real_, nrow = max_nodes, ncol = ncol(region_centroids))
  parent_id <- rep(NA_integer_, max_nodes)
  left_id <- rep(NA_integer_, max_nodes)
  right_id <- rep(NA_integer_, max_nodes)
  merge_score <- rep(NA_real_, max_nodes)
  is_leaf <- rep(FALSE, max_nodes)
  labels_under_node <- vector("list", max_nodes)

  neighbors <- vector("list", max_nodes)
  leaf_values <- leaf_values %||% seq_len(R)

  for (r in seq_len(R)) {
    members[[r]] <- sort(unique(as.integer(region_members[[r]])))
    sizes[r] <- length(members[[r]])
    centroids[r, ] <- region_centroids[r, ]
    is_leaf[r] <- TRUE
    neighbors[[r]] <- which(region_adj[r, ])
    labels_under_node[[r]] <- as.integer(leaf_values[r])
  }

  active <- seq_len(R)
  next_id <- R + 1L

  find_best_pair <- function(active_ids, neighbors_local) {
    pairs_i <- integer(0)
    pairs_j <- integer(0)
    for (a in active_ids) {
      nei <- neighbors_local[[a]]
      nei <- nei[nei %in% active_ids]
      nei <- nei[nei > a]
      if (length(nei) > 0L) {
        pairs_i <- c(pairs_i, rep(a, length(nei)))
        pairs_j <- c(pairs_j, nei)
      }
    }
    if (length(pairs_i) == 0L) return(NULL)

    deltas <- mapply(function(a, b) {
      .ward_delta(
        size_a = sizes[a], size_b = sizes[b],
        centroid_a = centroids[a, ], centroid_b = centroids[b, ]
      )
    }, pairs_i, pairs_j)

    ord <- order(deltas, pmin(pairs_i, pairs_j), pmax(pairs_i, pairs_j))
    list(a = pairs_i[ord[1]], b = pairs_j[ord[1]], delta = deltas[ord[1]])
  }

  while (length(active) > 1L) {
    best <- find_best_pair(active, neighbors)
    if (is.null(best)) {
      stop(
        "No adjacent mergeable groups remain. ",
        "The supplied label-adjacency graph is disconnected."
      )
    }
    a <- best$a
    b <- best$b
    m <- next_id
    next_id <- next_id + 1L

    left_id[m] <- a
    right_id[m] <- b
    parent_id[a] <- m
    parent_id[b] <- m
    merge_score[m] <- best$delta

    members[[m]] <- sort(c(members[[a]], members[[b]]))
    sizes[m] <- sizes[a] + sizes[b]
    wa <- sizes[a] / sizes[m]
    wb <- sizes[b] / sizes[m]
    centroids[m, ] <- wa * centroids[a, ] + wb * centroids[b, ]
    labels_under_node[[m]] <- sort(unique(c(labels_under_node[[a]], labels_under_node[[b]])))

    nei_m <- sort(unique(c(neighbors[[a]], neighbors[[b]])))
    nei_m <- setdiff(nei_m, c(a, b))
    nei_m <- nei_m[nei_m %in% active]
    neighbors[[m]] <- nei_m

    for (u in nei_m) {
      neigh_u <- neighbors[[u]]
      neigh_u <- setdiff(neigh_u, c(a, b))
      neighbors[[u]] <- sort(unique(c(neigh_u, m)))
    }

    active <- setdiff(active, c(a, b))
    active <- c(active, m)
  }

  root <- active[1]
  n_nodes <- root
  node_table <- data.frame(
    node_id = seq_len(n_nodes),
    parent_id = parent_id[seq_len(n_nodes)],
    left_child_id = left_id[seq_len(n_nodes)],
    right_child_id = right_id[seq_len(n_nodes)],
    is_leaf = seq_len(n_nodes) <= R,
    n_spots = sizes[seq_len(n_nodes)],
    merge_score = merge_score[seq_len(n_nodes)],
    stringsAsFactors = FALSE
  )

  list(
    root_id = root,
    n_leaves = R,
    node_table = node_table,
    spot_indices = members[seq_len(n_nodes)],
    labels_under_node = labels_under_node[seq_len(n_nodes)],
    region_centroids = centroids[seq_len(n_nodes), , drop = FALSE],
    leaf_values = as.integer(leaf_values)
  )
}

# Label-level observed/local tree. Retained connected components define
# adjacency when they connect all labels; otherwise all components are used.
build_label_level_tree_from_labels <- function(
  X,
  labels,
  spot_graph,
  ci_pca_dims = 30L,
  n_comp_min = 0.05,
  K_expected = NULL
) {
  labels <- as.integer(labels)
  memb <- build_label_membership(labels = labels, K_expected = K_expected)
  cc_info <- split_labels_into_cc(labels = labels, spot_graph = spot_graph, n_comp_min = n_comp_min)

  cobj <- compute_region_centroids_pca(
    X = X,
    region_members = memb$label_members,
    pca_dims = ci_pca_dims
  )
  retained_adjacency <- build_label_adjacency_from_retained(
    labels = labels,
    spot_graph = spot_graph,
    cc_info = cc_info
  )
  unfiltered_adjacency <- build_label_adjacency_unfiltered(
    labels = labels,
    spot_graph = spot_graph,
    label_values = memb$label_values
  )
  retained_connected <- is_label_adjacency_connected(retained_adjacency)
  unfiltered_connected <- is_label_adjacency_connected(unfiltered_adjacency)

  if (retained_connected) {
    ladj <- retained_adjacency
    adjacency_source <- "retained_components_after_filtering"
    component_filter_used_for_tree <- TRUE
    component_filter_reason <- NA_character_
    component_filter_message <- paste(
      "Component-size filtering retained a connected label-adjacency graph;",
      "the hierarchical tree used retained components."
    )
  } else {
    if (!unfiltered_connected) {
      stop(paste(
        "No connected adjacency among cluster labels remained after",
        "component-size filtering, and the label-adjacency graph built from",
        "all spatial components without filtering is also disconnected.",
        "A spatially constrained complete hierarchy cannot be constructed."
      ))
    }
    ladj <- unfiltered_adjacency
    adjacency_source <- "all_components_without_filtering"
    component_filter_used_for_tree <- FALSE
    component_filter_reason <- "no_connected_adjacency_after_filtering"
    component_filter_message <- paste(
      "No connected adjacency among cluster labels remained after",
      "component-size filtering. The hierarchical tree was constructed using",
      "all spatial components without filtering."
    )
  }

  tree <- adjacency_constrained_ward_tree(
    region_members = memb$label_members,
    region_adj = ladj,
    region_centroids = cobj$centroids,
    leaf_values = memb$label_values
  )
  tree <- compute_leaf_counts(tree)
  tree$label_values <- memb$label_values
  tree$label_members <- memb$label_members
  tree$cc_info <- cc_info
  tree$label_adjacency <- ladj
  tree$retained_label_adjacency <- retained_adjacency
  tree$unfiltered_label_adjacency <- unfiltered_adjacency
  tree$adjacency_source <- adjacency_source
  tree$component_filter_used_for_tree <- component_filter_used_for_tree
  tree$component_filter_reason <- component_filter_reason
  tree$component_filter_message <- component_filter_message
  tree$adjacency_diagnostics <- list(
    n_labels = length(memb$label_values),
    retained_edge_count = .label_adjacency_edge_count(retained_adjacency),
    retained_connected = retained_connected,
    unfiltered_edge_count = .label_adjacency_edge_count(unfiltered_adjacency),
    unfiltered_connected = unfiltered_connected
  )
  tree
}

# Compute descendant leaf counts (Kv) for every node.
compute_leaf_counts <- function(tree) {
  node_table <- tree$node_table
  n_nodes <- nrow(node_table)
  leaf_count <- integer(n_nodes)

  recurse <- function(node_id) {
    row <- node_table[node_id, , drop = FALSE]
    if (isTRUE(row$is_leaf)) {
      leaf_count[node_id] <<- 1L
      return(1L)
    }
    lc <- recurse(row$left_child_id)
    rc <- recurse(row$right_child_id)
    leaf_count[node_id] <<- lc + rc
    lc + rc
  }

  recurse(tree$root_id)
  tree$node_table$leaf_count <- leaf_count
  tree
}

# Extract observed split for an internal node.
extract_observed_node_split <- function(tree, node_id) {
  node <- tree$node_table[node_id, , drop = FALSE]
  if (isTRUE(node$is_leaf)) stop("Cannot extract split from a leaf node.")
  left <- node$left_child_id
  right <- node$right_child_id
  list(
    left_spots = tree$spot_indices[[left]],
    right_spots = tree$spot_indices[[right]],
    left_child_id = left,
    right_child_id = right
  )
}

# Extract root split from a local tree.
extract_root_split_from_tree <- function(tree) {
  extract_observed_node_split(tree, tree$root_id)
}

# Internal nodes in root-to-leaf traversal order.
get_internal_nodes_topdown <- function(tree) {
  node_table <- tree$node_table
  q <- tree$root_id
  out <- integer(0)
  while (length(q) > 0L) {
    nid <- q[1]
    q <- q[-1]
    node <- node_table[nid, , drop = FALSE]
    if (isTRUE(node$is_leaf)) next
    out <- c(out, nid)
    q <- c(q, node$left_child_id, node$right_child_id)
  }
  out
}
