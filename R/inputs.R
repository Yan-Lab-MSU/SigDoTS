# Taking the objects users have: a Seurat object, a SpatialExperiment (or any
# SummarizedExperiment with coordinates given), or plain matrices. Everything is turned into
# the two matrices the engine uses, x (spots x genes) and coords (spots x 2), before
# .sigdots_check_inputs() sees them.

.need <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("To use %s, the package %s must be installed.", what, pkg), call. = FALSE)
}

# Spot coordinates from a data frame as Seurat's GetTissueCoordinates() returns it: columns
# x, y and cell (Seurat 5 FOVs), or imagecol, imagerow with the cells as row names (Visium
# images). Returned in the order of `cells`.
.coords_from_df <- function(tc, cells) {
  tc <- as.data.frame(tc)
  if ("cell" %in% names(tc)) { rownames(tc) <- tc$cell; tc$cell <- NULL }
  nm <- names(tc)
  pick <- if (all(c("x", "y") %in% nm)) c("x", "y") else
    if (all(c("imagecol", "imagerow") %in% nm)) c("imagecol", "imagerow") else
    nm[vapply(tc, is.numeric, logical(1))][1:2]
  if (anyNA(pick)) stop("Could not find two coordinate columns in the object's tissue coordinates.", call. = FALSE)
  missing <- setdiff(cells, rownames(tc))
  if (length(missing))
    stop(sprintf("%d of the %d cells have no tissue coordinates (first: %s). Give `coords` yourself.",
                 length(missing), length(cells), missing[1]), call. = FALSE)
  m <- as.matrix(tc[cells, pick])
  colnames(m) <- c("x", "y")
  m
}

# A column of the object's per-spot metadata, by name.
.meta_column <- function(x, name) {
  v <- if (inherits(x, "Seurat")) {
    if (!name %in% colnames(x[[]])) stop(sprintf("`candidate_labels`: there is no metadata column \"%s\" in the Seurat object.", name), call. = FALSE)
    x[[name, drop = TRUE]]
  } else {
    cd <- SummarizedExperiment::colData(x)
    if (!name %in% colnames(cd)) stop(sprintf("`candidate_labels`: there is no colData column \"%s\" in the object.", name), call. = FALSE)
    cd[[name]]
  }
  as.vector(v)
}

.sigdots_extract <- function(x, coords, genes, assay, candidate_labels) {
  stop_ <- function(...) stop(..., call. = FALSE)
  labels <- candidate_labels
  if (inherits(x, "Seurat")) {
    .need("SeuratObject", "a Seurat object as `x`")
    assay <- assay %||% SeuratObject::DefaultAssay(x)
    if (!assay %in% SeuratObject::Assays(x)) stop_(sprintf("The Seurat object has no assay \"%s\".", assay))
    if (!"data" %in% SeuratObject::Layers(x[[assay]]))
      stop_(sprintf("Assay \"%s\" has no normalised \"data\" layer. Run Seurat::NormalizeData() first.", assay))
    expr <- SeuratObject::LayerData(x, assay = assay, layer = "data")       # genes x cells
    if (is.null(genes)) {
      genes <- SeuratObject::VariableFeatures(x, assay = assay)
      if (!length(genes))
        stop_("No variable genes are set in the Seurat object. Run Seurat::FindVariableFeatures() first, or give `genes`.")
    }
    if (is.null(coords)) {
      if (!length(SeuratObject::Images(x)))
        stop_("The Seurat object has no image or FOV with the spot positions. Give `coords` (one row per cell, two columns).")
      coords <- .coords_from_df(SeuratObject::GetTissueCoordinates(x), colnames(expr))
    }
    if (is.character(labels) && length(labels) == 1L) labels <- .meta_column(x, labels)
    source <- sprintf("Seurat object, assay \"%s\", %d genes", assay, length(genes))
  } else if (inherits(x, "SummarizedExperiment")) {
    .need("SummarizedExperiment", "a SummarizedExperiment-based object as `x`")
    an <- SummarizedExperiment::assayNames(x)
    if (is.null(assay)) {
      if (!"logcounts" %in% an)
        stop_(sprintf("The object has no \"logcounts\" assay (it has: %s). Run scuttle::logNormCounts() first, or name the assay to use with `assay`.",
                      paste(an, collapse = ", ")))
      assay <- "logcounts"
    }
    if (!assay %in% an) stop_(sprintf("The object has no assay \"%s\" (it has: %s).", assay, paste(an, collapse = ", ")))
    expr <- SummarizedExperiment::assay(x, assay)
    if (is.null(colnames(expr))) colnames(expr) <- as.character(seq_len(ncol(expr)))
    if (is.null(genes)) {
      if (nrow(expr) > 5000)
        stop_(sprintf(paste0("The object has %d genes. Give `genes`, the genes to use (for example the 2,000 most variable: ",
                             "scran::getTopHVGs(scran::modelGeneVar(x), n = 2000))."), nrow(expr)))
      genes <- rownames(expr)
    }
    if (is.null(coords)) {
      if (!inherits(x, "SpatialExperiment"))
        stop_("Give `coords` (one row per spot, two columns): the object is not a SpatialExperiment, so it carries no spot positions.")
      .need("SpatialExperiment", "a SpatialExperiment object as `x`")
      coords <- SpatialExperiment::spatialCoords(x)
      rownames(coords) <- colnames(expr)
    }
    if (is.character(labels) && length(labels) == 1L) labels <- .meta_column(x, labels)
    source <- sprintf("%s, assay \"%s\", %d genes", class(x)[1], assay, length(genes))
  } else {
    if (is.null(coords)) stop_("`coords` is needed when `x` is a matrix: one row per spot, two columns.")
    if (!is.null(genes)) {
      if (inherits(x, "Matrix") || is.data.frame(x)) x <- as.matrix(x)
      missing <- setdiff(genes, colnames(x))
      if (length(missing)) stop_(sprintf("%d of the genes in `genes` are not columns of `x` (first: %s).", length(missing), missing[1]))
      x <- x[, genes, drop = FALSE]
    }
    return(list(x = x, coords = coords, candidate_labels = labels, source = "matrices"))
  }
  missing <- setdiff(genes, rownames(expr))
  if (length(missing)) stop_(sprintf("%d of the genes in `genes` are not in the object (first: %s).", length(missing), missing[1]))
  xm <- t(as.matrix(expr[genes, , drop = FALSE]))                           # spots x genes
  coords <- as.matrix(coords)
  if (!is.null(rownames(coords)) && !identical(rownames(coords), rownames(xm))) {
    missing <- setdiff(rownames(xm), rownames(coords))
    if (length(missing)) stop_(sprintf("%d spots of `x` have no row in `coords` (first: %s).", length(missing), missing[1]))
    coords <- coords[rownames(xm), , drop = FALSE]
  }
  list(x = xm, coords = coords, candidate_labels = labels, source = source)
}

#' Add the domains to a Seurat or SpatialExperiment object
#'
#' `add_domains()` writes the result of [sigdots()] into the per-spot metadata of the object
#' the data came from, so that the domains can be used like any other per-spot label (plotting,
#' subsetting, differential expression). Two columns are added: `sigdots_domain`, the domain
#' number of every spot (`NA` for a spot in a part of the tree that could not be tested), and
#' `sigdots_candidate_domain`, the candidate domain the spot was in.
#'
#' @param object A Seurat object, or a SpatialExperiment, SingleCellExperiment or
#'   SummarizedExperiment object. Its cells or spots must be exactly those `fit` was computed
#'   on, in the same order (which is the case when `fit` came from `sigdots(object, ...)`).
#' @param fit A result of [sigdots()].
#' @param prefix The first part of the two column names (`"sigdots"` by default, giving
#'   `sigdots_domain` and `sigdots_candidate_domain`). Use another prefix to keep the results of
#'   several runs side by side.
#' @details The spots are matched by their identifiers (the column names of the object and the
#'   spot identifiers of `fit`), so the object may hold its spots in a different order, but it
#'   must not hold spots that `fit` does not have; such spots stop the function with an error
#'   naming the first of them. Existing columns with the same names are replaced.
#' @return The object with the two columns added: in `object[[]]` (the metadata) for a Seurat
#'   object, in `colData(object)` for the others.
#' @seealso [domains()] for the same information as a plain table.
#' @examples
#' \donttest{
#' if (requireNamespace("SpatialExperiment", quietly = TRUE)) {
#'   spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(logcounts = t(sigdots_example$x)), spatialCoords = sigdots_example$coords)
#'   fit <- sigdots(spe, K0 = 5, verbose = FALSE)
#'   spe <- add_domains(spe, fit)
#'   table(SummarizedExperiment::colData(spe)$sigdots_domain)
#' }
#' }
#' @export
add_domains <- function(object, fit, prefix = "sigdots") {
  stopifnot(inherits(fit, "sigdots"))
  ids <- colnames(object)
  if (is.null(ids)) ids <- as.character(seq_len(ncol(object)))
  i <- match(ids, fit$spot_ids)
  if (anyNA(i))
    stop(sprintf("%d of the object's %d spots are not in the result (first: %s). The result must come from the same spots.",
                 sum(is.na(i)), length(ids), ids[which(is.na(i))[1]]), call. = FALSE)
  dom <- fit$domains[i]
  init <- if (isTRUE(fit$settings$candidate_labels_supplied)) fit$candidate_names[fit$candidate_domains[i]] else fit$candidate_domains[i]
  cols <- stats::setNames(list(dom, init), paste0(prefix, c("_domain", "_candidate_domain")))
  if (inherits(object, "Seurat")) {
    .need("SeuratObject", "a Seurat object")
    md <- as.data.frame(cols, row.names = ids, check.names = FALSE, stringsAsFactors = FALSE)
    return(SeuratObject::AddMetaData(object, metadata = md))
  }
  if (inherits(object, "SummarizedExperiment")) {
    .need("SummarizedExperiment", "this object")
    cd <- SummarizedExperiment::colData(object)
    for (nm in names(cols)) cd[[nm]] <- cols[[nm]]
    SummarizedExperiment::colData(object) <- cd
    return(object)
  }
  stop("`object` must be a Seurat object or a SummarizedExperiment-based object (SpatialExperiment, SingleCellExperiment).",
       call. = FALSE)
}

#' Expression matrix for sigdots() from a counts matrix
#'
#' `sigdots_prepare()` makes the expression matrix `x` that [sigdots()] takes from a matrix of
#' raw counts. It is the preparation used in the SigDoTS paper's simulations; any other
#' normalisation and gene selection (Seurat's, scran's, ...) can be used instead.
#'
#' @param counts The raw counts: one row per gene and one column per spot, as most tools
#'   return them, with the gene names as row names and the spot identifiers as column names.
#'   A matrix or a sparse `Matrix`.
#' @param n_genes The number of genes to keep (a whole number, 2000 by default). All genes
#'   are kept when there are fewer.
#' @param scale_factor Each spot's counts are scaled to this total before the log transform
#'   (10,000 by default).
#' @details Three steps:
#'   1. Normalisation: the counts of each spot are divided by the spot's total count and
#'      multiplied by `scale_factor`.
#'   2. Log transform: log(1 + value), as `Seurat::NormalizeData()` does by default.
#'   3. Gene selection: the variance of every gene's log-normalised values across the spots is
#'      computed, and the `n_genes` genes with the largest variance are kept, in decreasing order
#'      of variance.
#'
#'   The values are not scaled to variance 1. A sparse `dgCMatrix` stays sparse until the
#'   selected genes are turned into the dense result. Spots with a total count of 0 stop the
#'   function; remove them first.
#' @return A numeric matrix with one row per spot and one column per selected gene (the spot
#'   identifiers as row names, the gene names as column names), ready for [sigdots()].
#' @examples
#' counts <- t(round(expm1(sigdots_example$x) * 50))   # a toy counts matrix, genes x spots
#' x <- sigdots_prepare(counts, n_genes = 100)
#' dim(x)
#' @export
sigdots_prepare <- function(counts, n_genes = 2000, scale_factor = 1e4) {
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (!(is.matrix(counts) || inherits(counts, "Matrix"))) stop("`counts` must be a matrix: one row per gene, one column per spot.", call. = FALSE)
  if (is.null(rownames(counts))) stop("`counts` needs gene names as row names.", call. = FALSE)
  totals <- Matrix::colSums(counts)
  if (any(totals <= 0)) stop(sprintf("%d spot(s) have zero total counts; remove them first.", sum(totals <= 0)), call. = FALSE)
  lc <- if (inherits(counts, "dgCMatrix")) {
    m <- counts
    m@x <- log1p(m@x / rep.int(totals, diff(m@p)) * scale_factor)
    m
  } else log1p(sweep(as.matrix(counts), 2L, totals, "/") * scale_factor)
  n <- ncol(lc)
  mu <- Matrix::rowMeans(lc)
  v <- (Matrix::rowMeans(lc^2) - mu^2) * n / (n - 1)
  keep <- order(v, decreasing = TRUE)[seq_len(min(n_genes, nrow(lc)))]
  x <- t(as.matrix(lc[keep, , drop = FALSE]))
  if (is.null(rownames(x))) rownames(x) <- as.character(seq_len(nrow(x)))
  x
}
