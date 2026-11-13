# Print and summarise a SigDoTS result

`print(fit)` (or typing `fit`) shows the number of domains and the split
tree of a [`sigdots()`](sigdots.md) result. `summary(fit)` shows, in
addition, the settings of the run and a table with one row per test; the
table is also returned, as `summary(fit)$tests`, for use in your own
code.

## Usage

``` r
# S3 method for class 'sigdots'
print(x, ...)

# S3 method for class 'sigdots'
summary(object, ...)

# S3 method for class 'summary.sigdots'
print(x, ...)
```

## Arguments

- x, object:

  A result of [`sigdots()`](sigdots.md) (for `print.summary.sigdots()`,
  a result of [`summary()`](https://rdrr.io/r/base/summary.html)).

- ...:

  Not used.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns its argument
invisibly. [`summary()`](https://rdrr.io/r/base/summary.html) returns a
list of class `"summary.sigdots"` with two elements: `fit`, the result,
and `tests`, the table described above (a data frame with the columns
`node`, `spots`, `candidate_domains`, `statistic`, `null_used`,
`p_value`, `level` and `decision`). Printing it shows the summary.

## Details

The printed summary has four parts: the settings (header lines), the
number of domains, the split tree, and the tests table, followed, where
they apply, by the tests in which simulated datasets were replaced, the
parts that could not be tested, the parts whose spots fall into separate
pieces, and a note when the candidate domains were supplied by you.

## The header lines

- `Data`: the number of spots and genes, and, for an object, where the
  expression came from.

- `Start`: how the candidate domains were made (`cluster_method`, or
  "labels supplied by the user") and how many there are.

- `Spatial`: the computation of the spatial covariance (Vecchia
  approximation with its number of neighbours, or exact).

- `Level`: `alpha`, the significance level for the whole tree.

- `P-values`: "Monte Carlo", or "normal approximation" with its numbers
  of simulated datasets (see `pvalue_method` in
  [`sigdots()`](sigdots.md)).

## The split tree

Every line that starts with a number in brackets, `[9]` say, is a part
of the tree (a *node*): the number of spots it holds and the candidate
domains it is made of, then the decision of the test at that node.
SigDoTS tests from the top: the whole sample first, then the two sides
of every split. `split` means the two sides are significantly different
and each side is tested in turn; `no split` means they are not, and the
part is a domain. `=> DOMAIN k` marks a part that was not split further:
it is the k-th domain, and the lines under it list the candidate domains
it contains with their numbers of spots. A single candidate domain is
never tested; when its parent is split, it is a domain by itself.
`[!] NOT TESTED` marks a part that could not be tested, because its
spatial parameters could not be estimated or because too many of its
simulated datasets had a top split with a side smaller than
`min_child_size`; its spots get no domain. The tree is drawn with
box-drawing characters, or with `|--` and `` `-- `` where the console
cannot show them.

## The tests table

One row per part of the tree that the testing reached, from the top. The
printed table shows the tested parts; the returned table,
`summary(fit)$tests`, has every reached part, with `NA` numbers for
those not tested. The columns:

- `node` (printed as `test`): the node number used in the tree;

- `spots`: the number of spots in that part;

- `candidate_domains`: the candidate domains in that part (in the
  returned table only);

- `statistic`: the split statistic of the data, between 0 and 1: the
  within-group sum of squares of the two sides of the split divided by
  the total sum of squares, over all genes. A smaller value means a
  clearer split;

- `null_used`: the number of simulated single-domain datasets the test
  used. A simulated dataset whose top split leaves fewer than
  `min_child_size` spots on a side is replaced by a new one and is not
  counted here; the summary lists the tests where this happened, with
  the number replaced;

- `p_value` (printed as `p`): the p-value of the test. With the Monte
  Carlo p-value it is (1 + the number of simulated datasets whose
  statistic was at most the observed one) / (1 + `null_used`), so with
  few simulated datasets it can take only a few values; the decision is
  what the test establishes. With the normal approximation it is the
  normal probability of a statistic at most the observed one, given the
  mean and standard deviation of the simulated statistics; it is printed
  in scientific notation when below 0.0001;

- `level`: the significance level of this test,
  `alpha * (spots - 1) / (n - 1)` for a sample of n spots, which keeps
  the probability of any false split in the whole tree at most `alpha`;

- `decision`: `split` when `p_value <= level`, `no split` otherwise, or
  why the part was not tested: `not tested: spatial fit` (the spatial
  parameters could not be estimated because the part has too few spots;
  its spots get no domain), `not tested: simulated splits too small`
  (the simulated datasets whose top split left fewer than
  `min_child_size` spots on a side outnumbered the valid ones; its spots
  get no domain), `not tested: too few spots` or
  `not tested: a side too small` (smaller than `min_node_size` or
  `min_child_size` of [`sigdots_control()`](sigdots_control.md);
  reported as a domain), `not tested: B_cap too small` (the level is so
  small that more than `B_cap` simulated datasets would be needed;
  reported as a domain).

## Separate pieces

The tree needs the spots of each tested part to be connected through
their nearest-neighbour links (`k_nn` of
[`sigdots_control()`](sigdots_control.md)). Where they fall into
separate pieces (a layer cut by a fold of the tissue, or a few spots
whose neighbours all belong to other candidate domains), SigDoTS joins
each smaller piece to the rest by a link between their closest spots,
repeatedly until all are connected, before it builds the trees of the
data and of the simulated datasets. Nothing changes for parts whose
spots are connected. [`summary()`](https://rdrr.io/r/base/summary.html)
lists the parts where this happened, with the sizes of the pieces (for
example "242, 49, 14 and 3 single spots").

## See also

[`domains()`](domains.md) for the domain of every spot,
[`plot.sigdots()`](plot.sigdots.md) for the maps and the tree,
[`sigdots()`](sigdots.md) for the components of the result.

## Examples

``` r
# \donttest{
fit <- sigdots(sigdots_example$x, sigdots_example$coords, K0 = 5, verbose = FALSE)
fit
#> SigDoTS: 900 spots, 200 genes
#> Number of domains: 3   (spots per domain: 300, 300, 300)
#> 
#> Split tree ([node] spots, candidate domains: decision):
#> [9] 900 spots, candidate domains 1,2,3,4,5: split
#> ├── [7] 300 spots, candidate domains 2,3: no split   => DOMAIN 1
#> │   ├── candidate domain 2  109 spots
#> │   └── candidate domain 3  191 spots
#> └── [8] 600 spots, candidate domains 1,4,5: split
#>     ├── [4] 300 spots, candidate domain 4   => DOMAIN 2
#>     └── [6] 300 spots, candidate domains 1,5: no split   => DOMAIN 3
#>         ├── candidate domain 1  103 spots
#>         └── candidate domain 5  197 spots
#> 
#> Use summary() for the tests.
summary(fit)
#> SigDoTS 0.2.5.10
#> Data:     900 spots x 200 genes
#> Start:    k-means on 15 principal components; 5 candidate domains
#> Spatial:  Vecchia approximation (60 neighbours)
#> Level:    0.05 for the whole tree
#> P-values: Monte Carlo
#> 
#> Number of domains: 3   (spots per domain: 300, 300, 300)
#> 
#> Split tree:
#> [9] 900 spots, candidate domains 1,2,3,4,5: split
#> ├── [7] 300 spots, candidate domains 2,3: no split   => DOMAIN 1
#> │   ├── candidate domain 2  109 spots
#> │   └── candidate domain 3  191 spots
#> └── [8] 600 spots, candidate domains 1,4,5: split
#>     ├── [4] 300 spots, candidate domain 4   => DOMAIN 2
#>     └── [6] 300 spots, candidate domains 1,5: no split   => DOMAIN 3
#>         ├── candidate domain 1  103 spots
#>         └── candidate domain 5  197 spots
#> 
#> Tests, from the top:
#>  test spots statistic null_used      p  level decision
#>   [9]   900     0.718        19 0.0500 0.0500    split
#>   [7]   300     0.935        10 0.2727 0.0166 no split
#>   [8]   600     0.903        30 0.0323 0.0333    split
#>   [6]   300     0.946        15 0.1250 0.0166 no split
summary(fit)$tests
#>   node spots candidate_domains statistic null_used    p_value      level
#> 1    9   900         1,2,3,4,5 0.7178312        19 0.05000000 0.05000000
#> 2    7   300               2,3 0.9345954        10 0.27272727 0.01662959
#> 3    8   600             1,4,5 0.9027189        30 0.03225806 0.03331479
#> 4    6   300               1,5 0.9464686        15 0.12500000 0.01662959
#>   decision
#> 1    split
#> 2 no split
#> 3    split
#> 4 no split
# }
```
