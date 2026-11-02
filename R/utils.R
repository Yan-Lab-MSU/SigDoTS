# Each of v0.2.5.8's four files defined this operator; the package defines it once.
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Joins the separate pieces of a spot graph (components of the k-NN graph, or of its restriction
# to a part of the tree) by adding, for every piece but the largest, an edge from its spot that
# is closest to a spot outside the piece to that spot; repeated until one piece remains. The
# graph is unchanged when it is already connected. Used by the engine so that every clustering
# of the spots has a connected group graph and the constrained tree can always be built.
.connect_pieces <- function(spot_graph, coords) {
  memb <- igraph::components(spot_graph$graph)$membership
  n_pieces <- max(memb)
  if (n_pieces == 1L) return(spot_graph)
  piece_sizes <- sort(as.integer(table(memb)), decreasing = TRUE)
  adj <- spot_graph$adj
  coords <- as.matrix(coords)
  repeat {
    pieces <- split(seq_along(memb), memb)
    if (length(pieces) == 1L) break
    big <- which.max(lengths(pieces))
    for (k in seq_along(pieces)) {
      if (k == big) next
      a <- pieces[[k]]; b <- setdiff(seq_along(memb), a)
      nn <- RANN::nn2(coords[b, , drop = FALSE], coords[a, , drop = FALSE], k = 1L)
      j <- which.min(nn$nn.dists[, 1]); i <- a[j]; jj <- b[nn$nn.idx[j, 1]]
      adj[i, jj] <- TRUE; adj[jj, i] <- TRUE
    }
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
    memb <- igraph::components(g)$membership
    if (max(memb) == 1L) break
  }
  spot_graph$adj <- adj
  spot_graph$graph <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  spot_graph$piece_sizes <- piece_sizes   # before joining; absent when the graph was connected
  spot_graph
}

# A p-value as printed: four decimals, or two significant digits in scientific notation when it is
# below 0.0001 (the normal approximation can give p-values of any size).
.format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p >= 1e-4, formatC(p, format = "f", digits = 4),
                                formatC(p, format = "e", digits = 1)))
}
