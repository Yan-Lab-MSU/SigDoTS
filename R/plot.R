# Plots of a sigdots() result (ggplot2).

# Categorical colours for domains, in a fixed order. Beyond the first three, some pairs are hard
# to tell apart for colour-blind readers, so the maps also label each domain and draw the
# boundaries between domains; identity never rests on colour alone. More than eight domains
# fall back to an HCL palette with the same labels and boundaries.
.domain_colours <- function(n) {
  base <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}
.untested_colour <- "#9a9a96"
.ink <- "#0b0b0b"
.ink_secondary <- "#52514e"

.map_theme <- function() {
  ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(colour = .ink, face = "bold"),
                   plot.subtitle = ggplot2::element_text(colour = .ink_secondary),
                   legend.title = ggplot2::element_text(colour = .ink_secondary),
                   legend.text = ggplot2::element_text(colour = .ink),
                   plot.margin = ggplot2::margin(6, 6, 6, 6))
}
.chart_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "#e6e5e1", linewidth = 0.3),
                   axis.text = ggplot2::element_text(colour = .ink_secondary),
                   axis.title = ggplot2::element_text(colour = .ink_secondary),
                   plot.title = ggplot2::element_text(colour = .ink, face = "bold"),
                   plot.subtitle = ggplot2::element_text(colour = .ink_secondary),
                   legend.position = "top", legend.title = ggplot2::element_blank())
}

# Default point size: smaller as the number of spots grows.
.point_size <- function(n) max(0.2, min(3, 45 / sqrt(n)))

# The spot nearest to the centre of each set of spots, so a label sits on the set itself.
.label_positions <- function(coords, labels) {
  ids <- sort(unique(labels[!is.na(labels)]))
  do.call(rbind, lapply(ids, function(k) {
    i <- which(labels == k)
    c0 <- colMeans(coords[i, , drop = FALSE])
    j <- i[which.min(rowSums((coords[i, , drop = FALSE] - rep(c0, each = length(i)))^2))]
    data.frame(x = coords[j, 1], y = coords[j, 2], label = as.character(k))
  }))
}

# Boundary between two sets of spots: for each pair of nearby spots with different labels, a
# short segment across the midpoint, perpendicular to the line joining them.
.boundary_segments <- function(coords, labels) {
  n <- nrow(coords)
  k <- min(7L, n)
  nn <- RANN::nn2(coords, k = k)
  h <- stats::median(nn$nn.dists[, 2])
  i <- rep(seq_len(n), k - 1L)
  j <- as.vector(nn$nn.idx[, -1, drop = FALSE])
  d <- as.vector(nn$nn.dists[, -1, drop = FALSE])
  # direct neighbours only: on a square grid the diagonal ones (1.41 spacings) would draw crosses
  keep <- i < j & d <= 1.2 * h & !is.na(labels[i]) & !is.na(labels[j]) & labels[i] != labels[j]
  i <- i[keep]; j <- j[keep]; d <- d[keep]
  if (!length(i)) return(data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0)))
  mx <- (coords[i, 1] + coords[j, 1]) / 2; my <- (coords[i, 2] + coords[j, 2]) / 2
  ux <- -(coords[j, 2] - coords[i, 2]) / d; uy <- (coords[j, 1] - coords[i, 1]) / d
  half <- 0.5 * d
  data.frame(x = mx - half * ux, y = my - half * uy, xend = mx + half * ux, yend = my + half * uy)
}

#' Plot a SigDoTS result
#'
#' `plot(fit, type = ...)` draws one of four plots of a [sigdots()] result: a map of the
#' domains, a map of the candidate domains, the split tree, or the spatial correlation
#' estimated at a test. Each plot is a ggplot2 object, so titles, colours and themes can be
#' changed by adding to it as usual, and it can be saved with `ggplot2::ggsave()`.
#'
#' @param x A result of [sigdots()].
#' @param type Which plot to draw, one of:
#'   * `"domains"` (the default): a map of the spots at their coordinates, each coloured by
#'     its domain, with the domain numbers and the boundaries between domains drawn in. Spots
#'     in a part of the tree that could not be tested are grey.
#'   * `"candidates"`: the same map coloured by candidate domain, with the boundaries of the final
#'     domains drawn in, so that you can see which candidate domains were merged.
#'   * `"tree"`: the split tree, as `print(x)` shows it, as a diagram: every tested part with
#'     its decision, p-value and level; the parts that became domains in the colours of
#'     `"domains"`; the candidate domains along the bottom.
#'   * `"correlogram"`: for one test (`test`), the spatial correlation of the noise, as
#'     SigDoTS estimated it for the single-domain model of that test. The points are the
#'     correlation between pairs of spots at each distance, after removing each candidate
#'     domain's mean and averaged over genes; the line is the correlation the fitted model gives
#'     (`eta` is the share of the noise variance that is spatially correlated, `ell` how far
#'     the correlation reaches). Distances are in units of the typical distance between
#'     neighbouring spots. The line following the points means the model describes the
#'     spatial correlation of the data.
#' @param test For `"correlogram"`: the node number of a tested part, as shown in brackets by
#'   `print(x)` and in the `node` column of `summary(x)$tests`. `NULL` (the default) means the
#'   first test, at the top of the tree. A part that was not tested has no correlogram and gives
#'   an error that lists the tested nodes.
#' @param point_size For the two maps: the size of the points, a positive number. `NULL` (the
#'   default) means `45 / sqrt(number of spots)`, kept between 0.2 and 3.
#' @param ... Not used.
#' @details
#' The maps draw each spot at its coordinates with equal scales on the two axes, and the y axis
#' pointing upwards. Image coordinates (such as the pixel rows of a Visium image) count
#' downwards, so a map may appear upside down compared with the tissue image; add
#' `ggplot2::scale_y_reverse()` to turn it. Domains are coloured in the order of their numbers
#' with a fixed palette of eight colours (more domains use `grDevices::hcl.colors()`); the
#' boundaries are drawn between neighbouring spots of different domains.
#'
#' In the correlogram, points close to the line mean that the single-domain model describes the
#' spatial correlation of the data at that test; points well above the line at short distances
#' mean the model has less spatial correlation than the data, which can make splits easier to keep.
#' @return A ggplot object; printing it draws the plot.
#' @seealso [summary.sigdots()] for the tests as a table, [domains()] for the domain of every
#'   spot.
#' @examples
#' \donttest{
#' fit <- sigdots(sigdots_example$x, sigdots_example$coords, K0 = 5, verbose = FALSE)
#' plot(fit)
#' plot(fit, type = "candidates")
#' plot(fit, type = "tree")
#' plot(fit, type = "correlogram")
#' plot(fit) + ggplot2::ggtitle("My title")
#' }
#' @method plot sigdots
#' @export
plot.sigdots <- function(x, type = c("domains", "candidates", "tree", "correlogram"),
                         test = NULL, point_size = NULL, ...) {
  type <- match.arg(type)
  switch(type,
         domains = .plot_domains(x, point_size),
         candidates = .plot_candidates(x, point_size),
         tree = .plot_tree(x),
         correlogram = .plot_correlogram(x, .which_test(x, test)))
}

.which_test <- function(fit, test) {
  t <- fit$tests
  done <- t$node[t$tested]
  if (!length(done)) stop("No test was carried out, so there is nothing to show.", call. = FALSE)
  if (is.null(test)) return(done[1])
  if (!test %in% t$node) stop(sprintf("Node %s was not reached by the testing. Tested nodes: %s.", test,
                                      paste(done, collapse = ", ")), call. = FALSE)
  if (!test %in% done) stop(sprintf("Node %s was not tested (%s).", test,
                                    t$decision[match(test, t$node)]), call. = FALSE)
  test
}

.plot_domains <- function(fit, point_size) {
  co <- fit$coords
  k <- length(fit$domain_nodes)
  lev <- c(if (k) as.character(seq_len(k)), if (anyNA(fit$domains)) "not tested")
  dom <- factor(ifelse(is.na(fit$domains), "not tested", fit$domains), levels = lev)
  pal <- c(stats::setNames(.domain_colours(k), if (k) as.character(seq_len(k))), "not tested" = .untested_colour)
  df <- data.frame(x = co[, 1], y = co[, 2], domain = dom)
  lab <- .label_positions(co, fit$domains)
  # boundaries also between two separate untested parts (negative node numbers keep them apart)
  seg <- .boundary_segments(co, ifelse(is.na(fit$domains), -fit$untested_spot_node, fit$domains))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(colour = domain), size = point_size %||% .point_size(nrow(df))) +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          colour = .ink, linewidth = 0.4, lineend = "round") +
    ggplot2::scale_colour_manual(values = pal[lev], drop = FALSE, name = "Domain") +
    ggplot2::coord_equal() + .map_theme() +
    ggplot2::labs(title = sprintf("SigDoTS: %d domain%s", k, if (k == 1) "" else "s"),
                  subtitle = sprintf("%d spots%s", fit$n_spots,
                                     if (anyNA(fit$domains)) sprintf("; %d not tested (grey)", sum(is.na(fit$domains))) else ""))
  if (!is.null(lab) && nrow(lab))
    p <- p + ggplot2::geom_label(data = lab, ggplot2::aes(x = x, y = y, label = label), colour = .ink,
                                 fill = "white", size = 3.5, label.size = 0.2, label.padding = ggplot2::unit(0.12, "lines"))
  p
}

.plot_candidates <- function(fit, point_size) {
  co <- fit$coords
  g <- fit$candidate_domains
  K <- length(fit$candidate_names)
  df <- data.frame(x = co[, 1], y = co[, 2], group = factor(fit$candidate_names[g], levels = fit$candidate_names))
  lab <- .label_positions(co, g)
  lab$label <- fit$candidate_names[as.integer(lab$label)]
  seg <- .boundary_segments(co, ifelse(is.na(fit$domains), -fit$untested_spot_node, fit$domains))
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(colour = group), size = point_size %||% .point_size(nrow(df))) +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          colour = .ink, linewidth = 0.5, lineend = "round") +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = x, y = y, label = label), colour = .ink, size = 3,
                       fontface = "bold") +
    ggplot2::scale_colour_manual(values = .domain_colours(K), name = "Candidate domain") +
    ggplot2::coord_equal() + .map_theme() +
    ggplot2::labs(title = sprintf("Candidate domains (%d) and final domains (%d)", K, length(fit$domain_nodes)),
                  subtitle = "Black lines: boundaries between final domains")
}

# Layout of the split tree: candidate domains along the bottom in the order of print(), each
# internal node centred above its two children, depth from the top.
.tree_layout <- function(fit) {
  tr <- fit$tree
  pos <- data.frame(node = tr$node, x = NA_real_, depth = NA_integer_)
  nxt <- 0
  rec <- function(id, depth) {
    i <- match(id, tr$node)
    kids <- c(tr$left[i], tr$right[i]); kids <- kids[!is.na(kids)]
    pos$depth[i] <<- depth
    if (!length(kids)) { nxt <<- nxt + 1; pos$x[i] <<- nxt; return(invisible()) }
    for (k in kids) rec(k, depth + 1L)
    pos$x[i] <<- mean(pos$x[match(kids, tr$node)])
  }
  rec(fit$root, 0L)
  pos
}

.plot_tree <- function(fit) {
  tr <- fit$tree
  t <- fit$tests
  pos <- .tree_layout(fit)
  maxd <- max(pos$depth)
  pos$y <- -pos$depth
  # domain of every node at or below a domain node; status of every node
  dom_of <- rep(NA_integer_, nrow(tr))
  for (k in seq_along(fit$domain_nodes)) {
    stack <- fit$domain_nodes[k]
    while (length(stack)) {
      id <- stack[1]; stack <- stack[-1]; i <- match(id, tr$node)
      dom_of[i] <- k
      stack <- c(stack, stats::na.omit(c(tr$left[i], tr$right[i])))
    }
  }
  untested_below <- rep(FALSE, nrow(tr))
  for (id in fit$untested_nodes) {
    stack <- id
    while (length(stack)) {
      i <- match(stack[1], tr$node); stack <- stack[-1]
      untested_below[i] <- TRUE
      stack <- c(stack, stats::na.omit(c(tr$left[i], tr$right[i])))
    }
  }
  status <- ifelse(tr$is_candidate_domain, "candidate domain", "not reached")
  k <- match(tr$node, t$node)
  status[!is.na(k)] <- ifelse(t$decision[k[!is.na(k)]] == "split", "split",
                       ifelse(t$decision[k[!is.na(k)]] == "no split", "no split", "not tested"))
  nodes <- data.frame(node = tr$node, x = pos$x, y = pos$y, status = status, domain = dom_of,
                      untested = untested_below)
  # elbow edges: parent's level down to each child
  e <- tr[!is.na(tr$parent), c("node", "parent")]
  pi <- match(e$parent, nodes$node); ci <- match(e$node, nodes$node)
  edges <- rbind(
    data.frame(x = nodes$x[pi], y = nodes$y[pi], xend = nodes$x[ci], yend = nodes$y[pi], domain = nodes$domain[ci],
               untested = nodes$untested[ci]),
    data.frame(x = nodes$x[ci], y = nodes$y[pi], xend = nodes$x[ci], yend = nodes$y[ci], domain = nodes$domain[ci],
               untested = nodes$untested[ci]))
  K <- length(fit$domain_nodes)
  pal <- .domain_colours(max(K, 1))
  edges$colour <- ifelse(edges$untested, .untested_colour,
                         ifelse(is.na(edges$domain), .ink_secondary, pal[pmax(1L, edges$domain)]))
  nodes$colour <- ifelse(nodes$untested, .untested_colour,
                         ifelse(is.na(nodes$domain), .ink, pal[pmax(1L, nodes$domain)]))
  nodes$shape <- ifelse(nodes$status == "candidate domain", 15, ifelse(nodes$status == "not tested", 1, 16))
  # labels
  tested <- t[t$tested, , drop = FALSE]
  lab_test <- data.frame(x = nodes$x[match(tested$node, nodes$node)], y = nodes$y[match(tested$node, nodes$node)],
                         label = sprintf("[%d] %d spots\n%s: p = %s, level %.4f%s", tested$node, tested$spots,
                                         tested$decision, .format_p(tested$p_value), tested$level,
                                         ifelse(is.na(match(tested$node, fit$domain_nodes)), "",
                                                sprintf("\n=> DOMAIN %d", match(tested$node, fit$domain_nodes)))))
  un <- t[!t$tested, , drop = FALSE]
  lab_un <- data.frame(x = nodes$x[match(un$node, nodes$node)], y = nodes$y[match(un$node, nodes$node)],
                       label = sprintf("[%d] %d spots\n%s%s", un$node, un$spots,
                                       ifelse(un$decision == "not tested: spatial fit", "[!] NOT TESTED (spatial fit)", un$decision),
                                       ifelse(is.na(match(un$node, fit$domain_nodes)), "",
                                              sprintf("\n=> DOMAIN %d", match(un$node, fit$domain_nodes)))))
  leaves <- nodes[nodes$status == "candidate domain", ]
  lab_leaf <- data.frame(x = leaves$x, y = leaves$y,
                         label = sprintf("%s\n%d", fit$candidate_names[vapply(fit$groups_under[as.character(leaves$node)], `[`, integer(1), 1)],
                                         tr$spots[match(leaves$node, tr$node)]))
  leaf_dom <- which(fit$domain_nodes %in% leaves$node)
  lab_dom <- if (length(leaf_dom)) data.frame(x = nodes$x[match(fit$domain_nodes[leaf_dom], nodes$node)],
                                             y = nodes$y[match(fit$domain_nodes[leaf_dom], nodes$node)],
                                             label = sprintf("DOMAIN %d", leaf_dom)) else NULL
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          colour = edges$colour, linewidth = 0.7, lineend = "round") +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y), colour = nodes$colour, shape = nodes$shape,
                        size = 3, stroke = 1) +
    ggplot2::geom_label(data = rbind(lab_test, lab_un), ggplot2::aes(x = x, y = y, label = label),
                        colour = .ink, fill = "white", size = 2.9, label.size = 0.15, vjust = 0,
                        nudge_y = 0.1, lineheight = 0.95) +
    ggplot2::geom_text(data = lab_leaf, ggplot2::aes(x = x, y = y, label = label), colour = .ink_secondary,
                       size = 2.7, vjust = 1.4, lineheight = 0.9) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0.6, 0.75))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.6)) +
    .map_theme() +
    ggplot2::labs(title = "SigDoTS split tree",
                  subtitle = "Bottom row: candidate domains (name, number of spots)")
  if (!is.null(lab_dom))
    p <- p + ggplot2::geom_text(data = lab_dom, ggplot2::aes(x = x, y = y, label = label), colour = .ink,
                                fontface = "bold", size = 3, vjust = -0.9)
  p
}

.plot_correlogram <- function(fit, test) {
  t <- fit$tests[match(test, fit$tests$node), ]
  cg <- fit$correlograms[[as.character(test)]]
  obs <- data.frame(distance = cg$distance, value = cg$correlation, what = "observed")
  fitd <- data.frame(distance = cg$distance, value = cg$fitted, what = "fitted")
  fitd <- fitd[!is.na(fitd$value), ]
  first <- fit$tests$node[fit$tests$tested][1]
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "#c9c8c3", linewidth = 0.4) +
    ggplot2::geom_line(data = fitd, ggplot2::aes(x = distance, y = value, colour = what), linewidth = 0.8) +
    ggplot2::geom_point(data = obs, ggplot2::aes(x = distance, y = value, colour = what), size = 2) +
    ggplot2::scale_colour_manual(values = c(observed = "#2a78d6", fitted = "#eb6834"),
                                 breaks = c("observed", "fitted")) +
    .chart_theme() +
    ggplot2::labs(x = "Distance", y = "Spatial correlation",
                  title = if (identical(test, first)) "Spatial correlation" else
                    sprintf("Spatial correlation, test at node %d (%d spots)", test, t$spots),
                  subtitle = sprintf("eta = %.3f, ell = %.2f", t$eta, t$ell))
}

utils::globalVariables(c("x", "y", "xend", "yend", "domain", "group", "label", "distance", "value", "what"))
