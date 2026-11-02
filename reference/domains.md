# The domain of every spot

`domains(fit)` returns the result of [`sigdots()`](sigdots.md) as a
table with one row per spot, in the order of the input data, so that it
can be joined to any other per-spot information (the true layers, cell
types, a plot). For a Seurat or SpatialExperiment object,
[`add_domains()`](add_domains.md) does this joining for you.

## Usage

``` r
domains(fit)
```

## Arguments

- fit:

  A result of [`sigdots()`](sigdots.md).

## Value

A data frame with one row per spot and three columns:

- `spot`: the spot identifier (the row name of the expression matrix, or
  the row number when there are no row names);

- `domain`: the number of the domain the spot belongs to, 1, 2, ... in
  the order of the split tree from left to right (as printed by
  `print(fit)`). `NA` for a spot in a part of the tree that could not be
  tested (see `summary(fit)` and the warning given at the end of the
  run);

- `candidate_domain`: the candidate domain the spot was in before the
  tests (its number, or your own label when `candidate_labels` was
  given).

## Details

The domains are numbered by their position in the split tree, from left
to right as `print(fit)` draws it, not by size. The spots of domain `k`
are `domains(fit)$spot[domains(fit)$domain == k]`.

## See also

[`summary.sigdots()`](summary.sigdots.md) for the tests,
[`plot.sigdots()`](plot.sigdots.md) for a map of the domains.

## Examples

``` r
# \donttest{
fit <- sigdots(sigdots_example$x, sigdots_example$coords, K0 = 5, verbose = FALSE)
head(domains(fit))
#>      spot domain candidate_domain
#> 1 spot001      3                1
#> 2 spot002      3                1
#> 3 spot003      3                1
#> 4 spot004      3                5
#> 5 spot005      3                5
#> 6 spot006      3                5
table(domains(fit)$domain)
#> 
#>   1   2   3 
#> 300 300 300 
table(domain = domains(fit)$domain, true_region = sigdots_example$region)
#>       true_region
#> domain   1   2   3
#>      1   0   0 300
#>      2   0 300   0
#>      3 300   0   0
# }
```
