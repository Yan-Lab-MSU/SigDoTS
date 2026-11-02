# A small simulated dataset with three spatial domains

900 spots on a 30 x 30 grid, 200 genes, three true domains of 300 spots
each (vertical bands). The first 40 genes ("marker" genes) have
different means in the three domains; the other 160 do not. The noise is
spatially correlated: half of its variance (`eta = 0.5`) is correlated
across spots, over a range of about 3 spot spacings (`ell = 3`, Matérn
3/2 correlation), and the genes are correlated with each other.
Expression values are log-normalised counts simulated from this model,
as in the simulation study of the SigDoTS paper.

## Usage

``` r
sigdots_example
```

## Format

A list with three elements:

- x:

  Numeric matrix, 900 spots x 200 genes, log-normalised expression. Row
  names are spot identifiers; column names mark the 40 marker genes.

- coords:

  Numeric matrix, 900 x 2, the spot coordinates (columns `x`, `y`).

- region:

  Integer vector of length 900, the true domain of each spot (1, 2 or
  3). SigDoTS does not use it; it is there to compare with the result.

## Source

Simulated for the package with the data generator of the SigDoTS paper's
simulation study: a 30 x 30 grid of spots divided into three vertical
bands of 10 columns; for each gene, a latent log-expression with a
band-specific mean for the 40 marker genes, spatially correlated noise
(`eta = 0.5`, `ell = 3` spot spacings) and correlation between genes;
counts drawn from a Poisson distribution given the latent values and a
spot-specific sequencing depth; then log(1 + 10,000 x count / total
count of the spot). Random seed 20261008.

## Details

`sigdots(sigdots_example$x, sigdots_example$coords, K0 = 5)` runs in
about 20 seconds on one core and recovers the three domains;
[`vignette("SigDoTS")`](../articles/SigDoTS.md) walks through it.

## Examples

``` r
str(sigdots_example, max.level = 1)
#> List of 3
#>  $ x     : num [1:900, 1:200] 4.92 4.93 4.02 3.07 4.05 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ coords: num [1:900, 1:2] 0 0.0345 0.069 0.1034 0.1379 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ region: int [1:900] 1 1 1 1 1 1 1 1 1 1 ...
table(sigdots_example$region)
#> 
#>   1   2   3 
#> 300 300 300 
```
