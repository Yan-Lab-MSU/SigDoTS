# Turns the engine's output into the object users read. Keeps everything print(), summary(),
# the accessors and the plots need, and drops the large per-test matrices (the exact
# backend's eigenvectors alone are n x n per test).
.sigdots_result <- function(fit, inp, settings, call, elapsed) {
  ot <- fit$observed_tree
  tr <- ot$node_table
  n <- nrow(inp$x)

  # Node ids in preorder (left before right), used to number the domains from left to right.
  preorder <- integer(0)
  walk <- function(id) {
    preorder <<- c(preorder, id)
    i <- match(id, tr$node_id)
    for (k in c(tr$left_child_id[i], tr$right_child_id[i])) if (!is.na(k)) walk(k)
  }
  walk(ot$root_id)
  cert <- vapply(fit$certified_regions, `[[`, integer(1), "node_id")
  domain_nodes <- preorder[preorder %in% cert]
  untested_nodes <- vapply(fit$untested_regions %||% list(), `[[`, integer(1), "node_id")
  domains <- rep(NA_integer_, n)
  for (k in seq_along(domain_nodes)) domains[ot$spot_indices[[match(domain_nodes[k], tr$node_id)]]] <- k
  untested_spot <- rep(NA_integer_, n)
  for (id in untested_nodes) untested_spot[ot$spot_indices[[match(id, tr$node_id)]]] <- id

  # The engine numbers the candidate domains 1..K0; with supplied labels, keep the user's names
  # for printing (character labels are numbered in sorted order, so "g10" may become 2).
  init <- as.integer(ot$bayesspace$labels)
  candidate_names <- if (isTRUE(settings$candidate_labels_supplied)) {
    lab <- as.character(inp$candidate_labels)
    vapply(seq_len(max(init)), function(k) lab[match(k, init)], character(1))
  } else as.character(seq_len(max(init)))

  tests <- .sigdots_tests_table(fit, preorder, candidate_names)
  null_statistics <- lapply(setNames(tests$node[tests$tested], tests$node[tests$tested]),
                            function(id) fit$node_results[[as.character(id)]]$CI_null)
  correlograms <- lapply(setNames(tests$node[tests$tested], tests$node[tests$tested]), function(id) {
    se <- fit$node_results[[as.character(id)]]$diagnostics$null_fit$spatial_estimate
    data.frame(distance = se$bin_centre / se$h, correlation = se$r_pooled,
               fitted = se$r_fitted %||% rep(NA_real_, length(se$r_pooled)), pairs = se$bin_pairs)
  })

  structure(list(
    domains = domains,
    domain_nodes = domain_nodes,
    untested_nodes = untested_nodes,
    untested_spot_node = untested_spot,
    candidate_domains = init,
    candidate_names = candidate_names,
    tree = data.frame(node = tr$node_id, parent = tr$parent_id, left = tr$left_child_id,
                      right = tr$right_child_id, is_candidate_domain = tr$is_leaf,
                      spots = tr$n_spots, candidate_domains = tr$leaf_count,
                      merge_score = tr$merge_score),
    root = ot$root_id,
    spot_indices = setNames(ot$spot_indices, tr$node_id),
    groups_under = setNames(ot$labels_under_node, tr$node_id),
    tests = tests,
    null_statistics = null_statistics,
    correlograms = correlograms,
    coords = inp$coords,
    spot_ids = inp$spot_ids,
    spot_pieces = ot$spot_graph$piece_sizes,
    n_spots = n,
    n_genes = ncol(inp$x),
    settings = settings,
    call = call,
    elapsed_sec = elapsed,
    version = as.character(utils::packageVersion("SigDoTS"))
  ), class = "sigdots")
}

# Why a reached node was not tested, in words (the engine's codes are kept in `skip_reason`).
.skip_reason_text <- c(
  spatial_fit_failed = "not tested: spatial fit",
  below_min_node_size = "not tested: too few spots",
  below_min_child_size_observed_split = "not tested: a side too small",
  bootstrap_resolution_exceeds_cap = "not tested: B_cap too small",
  Kv_less_than_2 = "not tested: one candidate domain",
  invalid_normal_approximation = "not tested: normal approximation impossible",
  insufficient_valid_normal_replicates = "not tested: too few simulated datasets")

# One row per node the top-down traversal reached, in preorder.
.sigdots_tests_table <- function(fit, preorder, candidate_names) {
  tr <- fit$observed_tree$node_table
  reached <- preorder[preorder %in% fit$diagnostics$tested_nodes]
  rows <- lapply(reached, function(id) {
    r <- fit$node_results[[as.character(id)]]
    i <- match(id, tr$node_id)
    se <- if (isTRUE(r$tested)) r$diagnostics$null_fit$spatial_estimate else NULL
    sp <- if (isTRUE(r$tested)) r$diagnostics$spatial_backend else NULL
    reason <- if (isTRUE(r$tested)) NA_character_ else r$reason_skipped %||% "unknown"
    decision <- if (isTRUE(r$tested)) {
      if (isTRUE(r$significant)) "split" else "no split"
    } else if (reason %in% names(.skip_reason_text)) .skip_reason_text[[reason]] else
      paste0("not tested: ", reason)
    ci_null <- r$CI_null
    data.frame(
      node = id,
      parent = tr$parent_id[i],
      spots = tr$n_spots[i],
      left_spots = tr$n_spots[match(tr$left_child_id[i], tr$node_id)],
      right_spots = tr$n_spots[match(tr$right_child_id[i], tr$node_id)],
      candidate_domains = paste(candidate_names[fit$observed_tree$labels_under_node[[i]]], collapse = ","),
      tested = isTRUE(r$tested),
      decision = decision,
      skip_reason = reason,
      statistic = r$CI_obs %||% NA_real_,
      null_median = if (length(ci_null)) stats::median(ci_null) else NA_real_,
      null_used = if (isTRUE(r$tested)) as.integer(r$b_valid) else NA_integer_,
      null_max = as.integer(r$B_max_effective %||% NA_integer_),
      clearer = if (length(ci_null)) as.integer(sum(ci_null <= r$CI_obs)) else NA_integer_,
      p_value = r$raw_p %||% NA_real_,
      level = r$alpha_star %||% NA_real_,
      stop_reason = r$early_stop_reason %||% NA_character_,
      p_initial = r$normal_p_initial %||% NA_real_,
      refined = r$normal_used_refine %||% NA,
      ell = se$ell_in_spacings %||% NA_real_,
      eta = se$eta %||% NA_real_,
      ell_at_lower_end = isTRUE(se$ell_at_grid_lower),
      ell_at_upper_end = isTRUE(se$ell_at_grid_upper),
      ell_capped = isTRUE(se$ell_capped),
      eta_at_lower_bound = isTRUE(se$eta_at_bound),
      eta_at_upper_bound = isTRUE(se$eta_at_upper_bound),
      backend = sp$backend_type %||% NA_character_,
      vecchia_neighbors = as.integer(sp$vecchia_neighbors %||% NA_integer_),
      pieces = if (length(r$piece_sizes)) paste(r$piece_sizes, collapse = ",") else NA_character_,
      spatial_fit_reason = r$spatial_fit_reason %||% NA_character_,
      spatial_fit_usable_bins = as.integer(r$spatial_fit_usable_bins %||% NA_integer_),
      min_pairs_per_bin = as.integer(r$min_pairs_per_bin %||% NA_integer_),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Warnings at the end of a run: parts whose spatial parameters could not be estimated, tests that
# would need more than B_cap simulated datasets, and spatial ranges at the upper end of the allowed
# range. Short: what happened, which nodes, what to do; the
# reasons are in vignette("faq").
.sigdots_warnings <- function(res) {
  t <- res$tests
  nofit <- t[t$decision == "not tested: spatial fit", , drop = FALSE]
  if (nrow(nofit)) {
    warning(sprintf(paste0(
      "The spatial parameters could not be estimated at %d part(s) of the tree (too few spots), so these parts ",
      "were not tested and their spots have no domain: %s. A smaller K0 makes the parts larger; see summary() and ",
      "vignette(\"faq\")."),
      nrow(nofit),
      paste(sprintf("node %d (%d spots, candidate domains %s)", nofit$node, nofit$spots, nofit$candidate_domains), collapse = "; ")),
      call. = FALSE)
  }
  capped <- t[t$skip_reason %in% "bootstrap_resolution_exceeds_cap", , drop = FALSE]
  if (nrow(capped)) {
    warning(sprintf(paste0(
      "At %d part(s) of the tree the test would need more than B_cap = %d simulated datasets: node(s) %s. ",
      "These parts were not tested and are reported as domains. To test them, raise B_cap in sigdots_control()."),
      nrow(capped), res$settings$control$B_cap, paste(capped$node, collapse = ", ")), call. = FALSE)
  }
  high <- t[t$tested & (t$ell_at_upper_end | t$ell_capped), , drop = FALSE]
  if (nrow(high)) {
    rng <- res$settings$control$ell_range
    warning(sprintf(paste0(
      "At %d test(s) the estimated spatial range ell is at the upper end of its allowed range (%g): node(s) %s. ",
      "Rerun with a wider range, for example control = sigdots_control(ell_range = c(%g, %g)); see vignette(\"faq\")."),
      nrow(high), rng[2], paste(high$node, collapse = ", "), rng[1], 2 * rng[2]), call. = FALSE)
  }
  invisible(NULL)
}

# The table summary() prints and returns: one row per part of the tree the testing reached.
.summary_tests <- function(fit) {
  t <- fit$tests
  data.frame(node = t$node, spots = t$spots, candidate_domains = t$candidate_domains, statistic = t$statistic,
             null_used = t$null_used, p_value = t$p_value, level = t$level, decision = t$decision,
             stringsAsFactors = FALSE)
}

#' The domain of every spot
#'
#' `domains(fit)` returns the result of [sigdots()] as a table with one row per spot, in the
#' order of the input data, so that it can be joined to any other per-spot information (the
#' true layers, cell types, a plot). For a Seurat or SpatialExperiment object, [add_domains()]
#' does this joining for you.
#'
#' @param fit A result of [sigdots()].
#' @return A data frame with one row per spot and three columns:
#'   * `spot`: the spot identifier (the row name of the expression matrix, or the row number
#'     when there are no row names);
#'   * `domain`: the number of the domain the spot belongs to, 1, 2, ... in the order of the
#'     split tree from left to right (as printed by `print(fit)`). `NA` for a spot in a part
#'     of the tree that could not be tested (see `summary(fit)` and the warning given at the
#'     end of the run);
#'   * `candidate_domain`: the candidate domain the spot was in before the tests (its number, or your
#'     own label when `candidate_labels` was given).
#' @details The domains are numbered by their position in the split tree, from left to right
#'   as `print(fit)` draws it, not by size. The spots of domain `k` are
#'   `domains(fit)$spot[domains(fit)$domain == k]`.
#' @seealso [summary.sigdots()] for the tests, [plot.sigdots()] for a map of the domains.
#' @examples
#' \donttest{
#' fit <- sigdots(sigdots_example$x, sigdots_example$coords, K0 = 5, verbose = FALSE)
#' head(domains(fit))
#' table(domains(fit)$domain)
#' table(domain = domains(fit)$domain, true_region = sigdots_example$region)
#' }
#' @export
domains <- function(fit) {
  stopifnot(inherits(fit, "sigdots"))
  init <- if (isTRUE(fit$settings$candidate_labels_supplied)) fit$candidate_names[fit$candidate_domains] else
    fit$candidate_domains
  data.frame(spot = fit$spot_ids, domain = fit$domains, candidate_domain = init, stringsAsFactors = FALSE)
}
