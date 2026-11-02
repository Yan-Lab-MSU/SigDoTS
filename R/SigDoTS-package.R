#' SigDoTS: significance testing for spatial domains
#'
#' Domain-detection methods divide a tissue section into spatial domains, but the number of
#' domains comes from a parameter, not from a test. SigDoTS tests which of the proposed domains
#' are supported by the data. Starting from more candidate domains than you expect (`K0` of
#' them), it arranges them in a tree by spatial adjacency and expression similarity and tests
#' the splits of the tree from the top down, each against a model in which the part being split
#' is a single domain, with the spatial correlation of the noise and the correlation between
#' genes estimated from the data. The parts of the tree that are not split further are the
#' domains; the probability that any of the kept splits is false is at most `alpha`.
#'
#' @section The functions:
#' * [sigdots()] runs the test; it takes an expression matrix and coordinates, a Seurat
#'   object or a SpatialExperiment object, and the number of candidate domains `K0` (or your own
#'   candidate domains).
#' * [sigdots_control()] holds the settings that rarely need changing.
#' * [sigdots_prepare()] makes the expression matrix from a raw counts matrix.
#' * `print()` and `summary()` ([summary.sigdots()]) show the split tree and the tests.
#' * [domains()] gives the domain of every spot as a table; [add_domains()] writes it into
#'   a Seurat or SpatialExperiment object.
#' * `plot()` ([plot.sigdots()]) draws the domains, the candidate domains, the tree and the
#'   spatial correlation.
#' * [sigdots_example] is a small simulated dataset to try the package on.
#'
#' The tutorials: `vignette("SigDoTS")` (getting started), `vignette("own-data")` (your own
#' data), `vignette("faq")` (errors, warnings and questions), `vignette("how-it-works")` (the
#' method step by step). The website, with a real-data application, is at
#' <https://yan-lab-msu.github.io/SigDoTS/>.
#'
#' @section How to cite:
#' Wang H, Cao F, Shen H, Liu Y, Yan G (2026). Statistical significance testing for spatial
#' domain detection in spatial transcriptomics. Manuscript.
#'
#' @keywords internal
#' @importFrom stats dist median prcomp rnorm setNames
#' @importFrom Matrix sparseMatrix
#' @importFrom fields rdist
"_PACKAGE"
