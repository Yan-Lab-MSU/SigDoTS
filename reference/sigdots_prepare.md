# Expression matrix for sigdots() from a counts matrix

`sigdots_prepare()` makes the expression matrix `x` that
[`sigdots()`](sigdots.md) takes from a matrix of raw counts. It is the
preparation used in the SigDoTS paper's simulations; any other
normalisation and gene selection (Seurat's, scran's, ...) can be used
instead.

## Usage

``` r
sigdots_prepare(counts, n_genes = 2000, scale_factor = 10000)
```

## Arguments

- counts:

  The raw counts: one row per gene and one column per spot, as most
  tools return them, with the gene names as row names and the spot
  identifiers as column names. A matrix or a sparse `Matrix`.

- n_genes:

  The number of genes to keep (a whole number, 2000 by default). All
  genes are kept when there are fewer.

- scale_factor:

  Each spot's counts are scaled to this total before the log transform
  (10,000 by default).

## Value

A numeric matrix with one row per spot and one column per selected gene
(the spot identifiers as row names, the gene names as column names),
ready for [`sigdots()`](sigdots.md).

## Details

Three steps:

1.  Normalisation: the counts of each spot are divided by the spot's
    total count and multiplied by `scale_factor`.

2.  Log transform: log(1 + value), as
    [`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html)
    does by default.

3.  Gene selection: the variance of every gene's log-normalised values
    across the spots is computed, and the `n_genes` genes with the
    largest variance are kept, in decreasing order of variance.

The values are not scaled to variance 1. A sparse `dgCMatrix` stays
sparse until the selected genes are turned into the dense result. Spots
with a total count of 0 stop the function; remove them first.

## Examples

``` r
counts <- t(round(expm1(sigdots_example$x) * 50))   # a toy counts matrix, genes x spots
x <- sigdots_prepare(counts, n_genes = 100)
dim(x)
#> [1] 900 100
```
