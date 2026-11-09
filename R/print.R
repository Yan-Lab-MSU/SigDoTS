# Characters for drawing the tree: box-drawing characters where the console can show them,
# plain ASCII otherwise. Written as escapes so the R code itself stays ASCII.
.tree_glyphs <- function(unicode = isTRUE(l10n_info()[["UTF-8"]])) {
  if (unicode) list(tee = "\u251c\u2500\u2500 ", last = "\u2514\u2500\u2500 ", pipe = "\u2502   ",
                    blank = "    ", mark = "[!]")
  else list(tee = "|-- ", last = "`-- ", pipe = "|   ", blank = "    ", mark = "[!]")
}

.candidates_text <- function(g) if (length(g) == 1L) sprintf("candidate domain %s", g) else
  sprintf("candidate domains %s", paste(g, collapse = ","))

# The split tree as lines of text. Above the domains every node shows its test's decision;
# inside a domain only the candidate domains and their sizes are listed.
.tree_lines <- function(fit, g = .tree_glyphs()) {
  tr <- fit$tree
  t <- fit$tests
  out <- character(0)
  rec <- function(id, prefix, last, root, inside) {
    i <- match(id, tr$node)
    groups <- fit$candidate_names[fit$groups_under[[as.character(id)]]]
    line <- if (inside) {
      sprintf("%s  %d spots", .candidates_text(groups), tr$spots[i])
    } else {
      k <- match(id, t$node)
      dom <- match(id, fit$domain_nodes)
      head <- sprintf("[%d] %d spots, %s", id, tr$spots[i], .candidates_text(groups))
      mid <- if (!is.na(k)) switch(t$decision[k],
          "split" = ": split",
          "no split" = ": no split",
          "not tested: spatial fit" = sprintf(": %s NOT TESTED (spatial parameters could not be estimated)", g$mark),
          "not tested: simulated splits too small" =
            sprintf(": %s NOT TESTED (too many simulated datasets with a side below min_child_size)", g$mark),
          sprintf(": %s", t$decision[k])) else ""
      tail <- if (!is.na(dom)) sprintf("   => DOMAIN %d", dom) else ""
      paste0(head, mid, tail)
    }
    out <<- c(out, if (root) line else paste0(prefix, if (last) g$last else g$tee, line))
    kids <- c(tr$left[i], tr$right[i])
    kids <- kids[!is.na(kids)]
    pre <- if (root) "" else paste0(prefix, if (last) g$blank else g$pipe)
    now_inside <- inside || id %in% fit$domain_nodes || id %in% fit$untested_nodes
    for (j in seq_along(kids)) rec(kids[j], pre, j == length(kids), FALSE, now_inside)
  }
  rec(fit$root, "", TRUE, TRUE, FALSE)
  out
}

.domain_line <- function(fit) {
  sizes <- tabulate(fit$domains[!is.na(fit$domains)], nbins = length(fit$domain_nodes))
  s <- if (length(sizes)) sprintf("Number of domains: %d   (spots per domain: %s)", length(fit$domain_nodes),
                                  paste(sizes, collapse = ", ")) else "Number of domains: 0"
  if (length(fit$untested_nodes))
    s <- c(s, sprintf("Not tested: %d part(s) of the tree with %d spots in all (marked [!] below)",
                      length(fit$untested_nodes), sum(is.na(fit$domains))))
  s
}

#' Print and summarise a SigDoTS result
#'
#' `print(fit)` (or typing `fit`) shows the number of domains and the split tree of a
#' [sigdots()] result. `summary(fit)` shows, in addition, the settings of the run and a table
#' with one row per test; the table is also returned, as `summary(fit)$tests`, for use in your
#' own code.
#'
#' @param x,object A result of [sigdots()] (for `print.summary.sigdots()`, a result of
#'   `summary()`).
#' @param ... Not used.
#'
#' @details
#' The printed summary has four parts: the settings (header lines), the number of domains, the
#' split tree, and the tests table, followed, where they apply, by the tests in which simulated
#' datasets were replaced, the parts that could not be tested, the parts whose spots fall into
#' separate pieces, and a note when the candidate domains were supplied by you.
#'
#' @section The header lines:
#' * `Data`: the number of spots and genes, and, for an object, where the expression came from.
#' * `Start`: how the candidate domains were made (`cluster_method`, or "labels supplied by the
#'   user") and how many there are.
#' * `Spatial`: the computation of the spatial covariance (Vecchia approximation with its number
#'   of neighbours, or exact).
#' * `Level`: `alpha`, the significance level for the whole tree.
#' * `P-values`: "Monte Carlo", or "normal approximation" with its numbers of simulated datasets
#'   (see `pvalue_method` in [sigdots()]).
#'
#' @section The split tree:
#' Every line that starts with a number in brackets, `[9]` say, is a part of the tree (a
#' *node*): the number of spots it holds and the candidate domains it is made of, then the
#' decision of the test at that node. SigDoTS tests from the top: the whole sample first, then
#' the two sides of every split. `split` means the two sides are significantly different and
#' each side is tested in turn; `no split` means they are not, and the part is a domain.
#' `=> DOMAIN k` marks a part that was not split further: it is the k-th domain, and the
#' lines under it list the candidate domains it contains with their numbers of spots. A single
#' candidate domain is never tested; when its parent is split, it is a domain by itself.
#' `[!] NOT TESTED` marks a part that could not be tested, because its spatial parameters could
#' not be estimated or because too many of its simulated datasets had a top split with a side
#' smaller than `min_child_size`; its spots get no domain. The tree is drawn with box-drawing
#' characters, or with `|--` and `` `-- `` where the console cannot show them.
#'
#' @section The tests table:
#' One row per part of the tree that the testing reached, from the top. The printed table shows
#' the tested parts; the returned table, `summary(fit)$tests`, has every reached part, with `NA`
#' numbers for those not tested. The columns:
#' * `node` (printed as `test`): the node number used in the tree;
#' * `spots`: the number of spots in that part;
#' * `candidate_domains`: the candidate domains in that part (in the returned table only);
#' * `statistic`: the split statistic of the data, between 0 and 1: the within-group sum of
#'   squares of the two sides of the split divided by the total sum of squares, over all genes.
#'   A smaller value means a clearer split;
#' * `null_used`: the number of simulated single-domain datasets the test used. A simulated
#'   dataset whose top split leaves fewer than `min_child_size` spots on a side is replaced by
#'   a new one and is not counted here; the summary lists the tests where this happened, with
#'   the number replaced;
#' * `p_value` (printed as `p`): the p-value of the test. With the Monte Carlo p-value it is
#'   (1 + the number of simulated datasets whose statistic was at most the observed one) /
#'   (1 + `null_used`), so with few simulated datasets it can take only a few values; the
#'   decision is what the test establishes. With the normal approximation it is the normal
#'   probability of a statistic at most the observed one, given the mean and standard deviation
#'   of the simulated statistics; it is printed in scientific notation when below 0.0001;
#' * `level`: the significance level of this test, `alpha * (spots - 1) / (n - 1)` for a sample
#'   of n spots, which keeps the probability of any false split in the whole tree at most
#'   `alpha`;
#' * `decision`: `split` when `p_value <= level`, `no split` otherwise, or why the part was not
#'   tested: `not tested: spatial fit` (the spatial parameters could not be estimated because
#'   the part has too few spots; its spots get no domain),
#'   `not tested: simulated splits too small` (the simulated datasets whose top split left
#'   fewer than `min_child_size` spots on a side outnumbered the valid ones; its spots get no
#'   domain), `not tested: too few spots` or
#'   `not tested: a side too small` (smaller than `min_node_size` or `min_child_size` of
#'   [sigdots_control()]; reported as a domain), `not tested: B_cap too small` (the level is so
#'   small that more than `B_cap` simulated datasets would be needed; reported as a domain).
#'
#' @section Separate pieces:
#' The tree needs the spots of each tested part to be connected through their nearest-neighbour
#' links (`k_nn` of [sigdots_control()]). Where they fall into separate pieces (a layer cut by a
#' fold of the tissue, or a few spots whose neighbours all belong to other candidate domains),
#' SigDoTS joins each smaller piece to the rest by a link between their closest spots, repeatedly
#' until all are connected, before it builds the trees of the data and of the simulated datasets.
#' Nothing changes for parts whose spots are connected. `summary()` lists the parts where this
#' happened, with the sizes of the pieces (for example "242, 49, 14 and 3 single spots").
#'
#' @return `print()` returns its argument invisibly. `summary()` returns a list of class
#'   `"summary.sigdots"` with two elements: `fit`, the result, and `tests`, the table described
#'   above (a data frame with the columns `node`, `spots`, `candidate_domains`, `statistic`,
#'   `null_used`, `p_value`, `level` and `decision`). Printing it shows the summary.
#' @seealso [domains()] for the domain of every spot, [plot.sigdots()] for the maps and the tree,
#'   [sigdots()] for the components of the result.
#' @examples
#' \donttest{
#' fit <- sigdots(sigdots_example$x, sigdots_example$coords, K0 = 5, verbose = FALSE)
#' fit
#' summary(fit)
#' summary(fit)$tests
#' }
#' @name summary.sigdots
NULL

#' @rdname summary.sigdots
#' @method print sigdots
#' @export
print.sigdots <- function(x, ...) {
  cat(sprintf("SigDoTS: %d spots, %d genes\n", x$n_spots, x$n_genes))
  writeLines(.domain_line(x))
  cat("\nSplit tree ([node] spots, candidate domains: decision):\n")
  writeLines(.tree_lines(x))
  cat("\nUse summary() for the tests.\n")
  invisible(x)
}

#' @rdname summary.sigdots
#' @method summary sigdots
#' @export
summary.sigdots <- function(object, ...) {
  structure(list(fit = object, tests = .summary_tests(object)), class = "summary.sigdots")
}

#' @rdname summary.sigdots
#' @method print summary.sigdots
#' @export
print.summary.sigdots <- function(x, ...) {
  fit <- x$fit
  s <- fit$settings
  t <- fit$tests
  cat(sprintf("SigDoTS %s\n", fit$version))
  cat(sprintf("Data:     %d spots x %d genes%s\n", fit$n_spots, fit$n_genes,
              if (!is.null(s$input) && s$input != "matrices") paste0(" (", s$input, ")") else ""))
  cat(sprintf("Start:    %s; %d candidate domains\n",
              if (isTRUE(s$candidate_labels_supplied)) "labels supplied by the user"
              else if (s$cluster_method == "kmeans") sprintf("k-means on %d principal components", s$n_pcs)
              else s$cluster_method, s$K0))
  cat(sprintf("Spatial:  %s\n", if (s$spatial_method == "exact") "exact computation"
              else sprintf("Vecchia approximation (%d neighbours)", s$vecchia_neighbors)))
  cat(sprintf("Level:    %g for the whole tree\n", s$alpha))
  cat(sprintf("P-values: %s\n\n", if (identical(s$pvalue_method, "normal_adaptive"))
    sprintf("normal approximation (%d simulated datasets per test, %d when p is near the level)",
            s$normal_n_initial, s$normal_n_refine) else "Monte Carlo"))
  writeLines(.domain_line(fit))
  cat("\nSplit tree:\n")
  writeLines(.tree_lines(fit))

  tested <- t[t$tested, , drop = FALSE]
  if (nrow(tested)) {
    cat("\nTests, from the top:\n")
    tab <- data.frame(test = sprintf("[%d]", tested$node), spots = tested$spots,
                      statistic = formatC(tested$statistic, format = "f", digits = 3),
                      null_used = tested$null_used,
                      p = .format_p(tested$p_value),
                      level = formatC(tested$level, format = "f", digits = 4),
                      decision = tested$decision, check.names = FALSE)
    print(tab, row.names = FALSE, right = TRUE)
  }
  repl <- tested[!is.na(tested$null_replaced) & tested$null_replaced > 0L, , drop = FALSE]
  if (nrow(repl)) {
    cat(sprintf("\nSimulated datasets replaced, because their top split left fewer than %d spots on a side (min_child_size):\n",
                s$control$min_child_size))
    writeLines(sprintf("  [%d] %d replaced, %d used", repl$node, repl$null_replaced, repl$null_used))
  }
  nofit <- t[t$decision == "not tested: spatial fit", , drop = FALSE]
  if (nrow(nofit)) {
    cat("\nNot tested, because the spatial parameters could not be estimated (too few spots); these spots have no domain:\n")
    writeLines(sprintf("  [%d] %d spots, candidate domains %s", nofit$node, nofit$spots, nofit$candidate_domains))
  }
  small <- t[t$decision == "not tested: simulated splits too small", , drop = FALSE]
  if (nrow(small)) {
    cat(sprintf(paste0("\nNot tested, because the simulated datasets whose top split left fewer than %d spots on a side ",
                       "(min_child_size) outnumbered the valid ones; these spots have no domain:\n"),
                s$control$min_child_size))
    writeLines(sprintf("  [%d] %d spots, candidate domains %s (%d simulated datasets replaced)", small$node, small$spots,
                       small$candidate_domains, small$null_replaced))
  }
  disc <- t[!is.na(t$pieces), , drop = FALSE]
  if (length(fit$spot_pieces) > 1L || nrow(disc)) {
    cat(sprintf("\nSpots in separate pieces on the %d-nearest-neighbour graph (joined by their closest spots for the tree):\n", s$control$k_nn))
    if (length(fit$spot_pieces) > 1L)
      cat(sprintf("  all spots: pieces of %s\n", .pieces_text(fit$spot_pieces)))
    if (nrow(disc)) writeLines(sprintf("  [%d] %d spots: pieces of %s", disc$node, disc$spots,
                                       vapply(strsplit(disc$pieces, ","), function(v) .pieces_text(as.integer(v)), character(1))))
  }
  if (isTRUE(s$candidate_labels_supplied))
    cat(sprintf("\nThe candidate domains were supplied; the simulated datasets were grouped with %s.\n", s$cluster_method))
  invisible(x)
}

# "242, 49, 14 and 3 single spots": the piece sizes, with single spots counted rather than listed.
.pieces_text <- function(sizes) {
  sizes <- sort(as.integer(sizes), decreasing = TRUE)
  big <- sizes[sizes > 1L]; n1 <- sum(sizes == 1L)
  paste0(paste(big, collapse = ", "),
         if (n1) sprintf("%s%d single spot%s", if (length(big)) " and " else "", n1, if (n1 > 1) "s" else "") else "")
}
