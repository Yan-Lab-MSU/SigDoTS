# A small dataset for the tests: spots on an nx x ny grid split into three vertical bands;
# the first q_signal genes are shifted by band, the noise is independent normal. Runs in
# seconds; it is not meant to resemble real data.
make_bands <- function(nx = 12, ny = 12, p = 40, q_signal = 15, shift = 1.5, seed = 1) {
  set.seed(seed)
  coords <- as.matrix(expand.grid(x = seq_len(nx), y = seq_len(ny)))
  band <- cut(coords[, 1], breaks = 3, labels = FALSE)
  x <- matrix(stats::rnorm(nrow(coords) * p), nrow(coords), p)
  x[, seq_len(q_signal)] <- x[, seq_len(q_signal)] + shift * c(-1, 0, 1)[band]
  rownames(x) <- sprintf("spot%03d", seq_len(nrow(x)))
  list(x = x, coords = coords, band = band)
}
