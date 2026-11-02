# Add the domains to a Seurat or SpatialExperiment object

`add_domains()` writes the result of [`sigdots()`](sigdots.md) into the
per-spot metadata of the object the data came from, so that the domains
can be used like any other per-spot label (plotting, subsetting,
differential expression). Two columns are added: `sigdots_domain`, the
domain number of every spot (`NA` for a spot in a part of the tree that
could not be tested), and `sigdots_candidate_domain`, the candidate
domain the spot was in.

## Usage

``` r
add_domains(object, fit, prefix = "sigdots")
```

## Arguments

- object:

  A Seurat object, or a SpatialExperiment, SingleCellExperiment or
  SummarizedExperiment object. Its cells or spots must be exactly those
  `fit` was computed on, in the same order (which is the case when `fit`
  came from `sigdots(object, ...)`).

- fit:

  A result of [`sigdots()`](sigdots.md).

- prefix:

  The first part of the two column names (`"sigdots"` by default, giving
  `sigdots_domain` and `sigdots_candidate_domain`). Use another prefix
  to keep the results of several runs side by side.

## Value

The object with the two columns added: in `object[[]]` (the metadata)
for a Seurat object, in `colData(object)` for the others.

## Details

The spots are matched by their identifiers (the column names of the
object and the spot identifiers of `fit`), so the object may hold its
spots in a different order, but it must not hold spots that `fit` does
not have; such spots stop the function with an error naming the first of
them. Existing columns with the same names are replaced.

## See also

[`domains()`](domains.md) for the same information as a plain table.

## Examples

``` r
# \donttest{
if (requireNamespace("SpatialExperiment", quietly = TRUE)) {
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(logcounts = t(sigdots_example$x)), spatialCoords = sigdots_example$coords)
  fit <- sigdots(spe, K0 = 5, verbose = FALSE)
  spe <- add_domains(spe, fit)
  table(SummarizedExperiment::colData(spe)$sigdots_domain)
}
# }
```
