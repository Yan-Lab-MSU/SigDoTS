# Inputs other than two matrices: a Seurat object, a SpatialExperiment object, labels named by
# a metadata column, counts prepared with sigdots_prepare(). The objects are built from
# a 20 x 20 corner of sigdots_example (400 spots, so that each run takes a few seconds); the
# results must equal the matrix path on the same genes.

sub <- local({
  co <- sigdots_example$coords
  keep <- co[, 1] <= sort(unique(co[, 1]))[20] & co[, 2] <= sort(unique(co[, 2]))[20]
  list(x = sigdots_example$x[keep, ], coords = sigdots_example$coords[keep, ])
})
make_counts <- function() Matrix::Matrix(t(round(expm1(sub$x) * 50)), sparse = TRUE)   # genes x spots, integers

test_that("a Seurat object is accepted: its variable genes, its tissue coordinates, add_domains()", {
  skip_if_not_installed("Seurat"); skip_if_not_installed("SeuratObject")
  counts <- make_counts()
  obj <- SeuratObject::CreateSeuratObject(counts = counts, assay = "Spatial")
  expect_error(sigdots(obj, K0 = 5, verbose = FALSE), "NormalizeData")
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  expect_error(sigdots(obj, K0 = 5, verbose = FALSE), "FindVariableFeatures")
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 120, verbose = FALSE)
  expect_error(sigdots(obj, K0 = 5, verbose = FALSE), "no image or FOV")
  cents <- SeuratObject::CreateCentroids(data.frame(x = sub$coords[, 1], y = sub$coords[, 2],
                                                    cell = colnames(counts)))
  obj[["slice1"]] <- SeuratObject::CreateFOV(cents, type = "centroids", assay = "Spatial")
  fit <- sigdots(obj, K0 = 5, B_min = 5, B_max = 5, verbose = FALSE)
  # the same computation through the matrix path, on the same genes
  vf <- SeuratObject::VariableFeatures(obj)
  x <- t(as.matrix(SeuratObject::LayerData(obj, layer = "data")))[, vf]
  fit_m <- sigdots(x, sub$coords, K0 = 5, B_min = 5, B_max = 5, verbose = FALSE)
  expect_identical(fit$domains, fit_m$domains)
  expect_identical(fit$tests, fit_m$tests)
  expect_identical(fit$spot_ids, colnames(counts))
  expect_match(fit$settings$input, "Seurat object")
  expect_output(print(summary(fit)), "Seurat object")
  # genes given by name; add_domains writes two metadata columns
  fit_g <- sigdots(obj, K0 = 5, genes = vf[1:50], B_min = 5, B_max = 5, verbose = FALSE)
  expect_equal(fit_g$n_genes, 50L)
  obj <- add_domains(obj, fit)
  expect_identical(unname(obj[["sigdots_domain", drop = TRUE]]), fit$domains)
  expect_identical(unname(obj[["sigdots_candidate_domain", drop = TRUE]]), fit$candidate_domains)
  # labels named by a metadata column
  obj$mygroups <- paste0("g", fit$candidate_domains)
  fit_l <- sigdots(obj, candidate_labels = "mygroups", B_min = 5, B_max = 5, verbose = FALSE)
  expect_true(fit_l$settings$candidate_labels_supplied)
  expect_identical(domains(fit_l)$candidate_domain, paste0("g", fit$candidate_domains))
  expect_error(sigdots(obj, candidate_labels = "nosuch", verbose = FALSE), "no metadata column")
})

test_that("a SpatialExperiment object is accepted, and add_domains() writes into colData", {
  skip_if_not_installed("SpatialExperiment"); skip_if_not_installed("SummarizedExperiment")
  spe <- SpatialExperiment::SpatialExperiment(assays = list(logcounts = t(sub$x)),
                                              spatialCoords = sub$coords)
  fit <- sigdots(spe, K0 = 5, B_min = 5, B_max = 5, verbose = FALSE)
  fit_m <- sigdots(sub$x, sub$coords, K0 = 5, B_min = 5, B_max = 5, verbose = FALSE)
  expect_identical(fit$domains, fit_m$domains)
  expect_identical(fit$tests, fit_m$tests)
  expect_match(fit$settings$input, "SpatialExperiment")
  spe <- add_domains(spe, fit)
  expect_identical(SummarizedExperiment::colData(spe)$sigdots_domain, fit$domains)
  # no logcounts: told to normalise; another assay by name; too many genes without `genes`
  spe2 <- SpatialExperiment::SpatialExperiment(assays = list(counts = make_counts()), spatialCoords = sub$coords)
  expect_error(sigdots(spe2, K0 = 5, verbose = FALSE), "logNormCounts")
  spe3 <- SpatialExperiment::SpatialExperiment(assays = list(mine = t(sub$x)), spatialCoords = sub$coords)
  fit3 <- sigdots(spe3, K0 = 5, assay = "mine", B_min = 5, B_max = 5, verbose = FALSE)
  expect_identical(fit3$domains, fit_m$domains)
  # a SingleCellExperiment without coordinates needs coords
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = t(sub$x)))
  expect_error(sigdots(sce, K0 = 5, verbose = FALSE), "not a SpatialExperiment")
  fit4 <- sigdots(sce, coords = sub$coords, K0 = 5, B_min = 5, B_max = 5, verbose = FALSE)
  expect_identical(fit4$domains, fit_m$domains)
})

test_that("sigdots_prepare() turns a counts matrix into x", {
  counts <- make_counts()
  x <- sigdots_prepare(counts, n_genes = 100)
  expect_equal(dim(x), c(400L, 100L))
  expect_identical(rownames(x), colnames(counts))
  # the same from a sparse matrix
  xd <- sigdots_prepare(as.matrix(counts), n_genes = 100)
  expect_equal(xd, x)
  # the 100 kept genes are the 100 with the largest variance of the log-normalised values
  lc <- log1p(sweep(as.matrix(counts), 2, Matrix::colSums(counts), "/") * 1e4)
  expect_setequal(colnames(x), names(sort(apply(lc, 1, var), decreasing = TRUE))[1:100])
  expect_error(sigdots_prepare(unname(as.matrix(counts))), "gene names")
  fit <- sigdots(x, sub$coords, K0 = 5, B_min = 5, B_max = 5, verbose = FALSE)
  expect_s3_class(fit, "sigdots")
})

test_that("`genes` subsets a matrix by column name", {
  d <- make_bands(nx = 20, ny = 20)
  colnames(d$x) <- paste0("g", seq_len(ncol(d$x)))
  fit <- sigdots(d$x, d$coords, K0 = 5, genes = paste0("g", 1:30), B_min = 5, B_max = 5, verbose = FALSE)
  expect_equal(fit$n_genes, 30L)
  expect_error(sigdots(d$x, d$coords, K0 = 5, genes = "nope", verbose = FALSE), "not columns of `x`")
  expect_error(sigdots(d$x, K0 = 5, verbose = FALSE), "`coords` is needed")
})
