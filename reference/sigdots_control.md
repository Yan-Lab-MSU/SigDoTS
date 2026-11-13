# Settings of sigdots() that rarely need changing

`sigdots_control()` collects the settings of [`sigdots()`](sigdots.md)
that rarely need changing. Its defaults are the values used in the
SigDoTS paper. Pass its result as the `control` argument of
[`sigdots()`](sigdots.md), for example
`sigdots(x, coords, K0 = 10, control = sigdots_control(ell_range = c(0.3, 80)))`;
settings not named keep their defaults.

## Usage

``` r
sigdots_control(
  k_nn = 6,
  n_comp_min = 0.01,
  min_node_size = 20,
  min_child_size = 8,
  B_cap = 5000,
  eta_min = 0.01,
  ell_range = c(0.3, 40),
  ell_candidates = 160,
  pairs_per_bin = 4000,
  min_pairs_per_bin = 200,
  r_cap = 30,
  ridge = 1e-06,
  nugget = 1e-06,
  max_exact_memory_gib = 16,
  candidate_seed = NULL,
  normal_band = 10
)
```

## Arguments

- k_nn:

  The number of nearest spots each spot is linked to (a whole number, 6
  by default; on the hexagonal Visium array each spot has 6 nearest
  neighbours). Two candidate domains are neighbours, and can be merged
  into one part of the tree, when some of their spots are linked. The
  same links are used for every simulated dataset.

- n_comp_min:

  A number greater than 0 and at most 1 (0.01 by default). When the
  spots of a candidate domain fall into separate pieces on the
  nearest-neighbour links, pieces with fewer spots than this share of
  the candidate domain are ignored when deciding which candidate domains
  are neighbours, so that a few stray spots do not make distant
  candidate domains neighbours. If ignoring them would leave the
  candidate domains unconnected, all pieces are used.

- min_node_size:

  A whole number (20 by default). A part of the tree with fewer spots is
  not tested and is reported as a domain.

- min_child_size:

  A whole number (8 by default). A split that would leave fewer spots
  than this on either side is not tested, and the part is reported as a
  domain. A simulated dataset whose top split leaves fewer spots on a
  side is replaced by a new one, so that the test keeps its number of
  simulated datasets ([`summary()`](https://rdrr.io/r/base/summary.html)
  lists how many were replaced). If such simulated datasets outnumber
  the valid ones, the part is not tested, its spots get no domain, and
  [`sigdots()`](sigdots.md) warns.

- B_cap:

  The largest number of simulated datasets a single Monte Carlo test may
  use (a whole number, 5000 by default, at least `B_max`). A test on m
  of n spots needs about `(n - 1) / (alpha * (m - 1))` simulated
  datasets to keep a split, so a test on a small part of a large sample
  may need more than `B_cap`; that part is then reported as a domain
  without a test, and [`sigdots()`](sigdots.md) warns. Not used by
  `pvalue_method = "normal_adaptive"`.

- eta_min:

  The smallest value allowed for `eta`, the share of the noise variance
  that is spatially correlated (a number greater than 0 and smaller than
  1, 0.01 by default; the largest value is 1).

- ell_range:

  The range searched for `ell`, the distance over which the noise is
  spatially correlated: two increasing positive numbers, in spot
  spacings (the median distance from a spot to its nearest other spot in
  the part being tested), `c(0.3, 40)` by default. `ell` is also kept
  below twice the largest distance between two spots of the same
  candidate domain. When [`sigdots()`](sigdots.md) warns that `ell` is
  at the upper end of this range, the true range may be longer, which
  can make splits easier to keep; rerun with a wider range, for example
  `c(0.3, 80)`.

- ell_candidates:

  The number of values of `ell` tried in `ell_range`, evenly spaced on a
  log scale (a whole number from 2, 160 by default); the best of them is
  then refined between its two neighbouring values.

- pairs_per_bin, min_pairs_per_bin:

  `eta` and `ell` are estimated from pairs of spots in the same
  candidate domain, grouped by their distance into classes (one per spot
  spacing from 0.5 to 30.5 spacings, then 30.5-35.5, 35.5-42.5,
  42.5-55.5 and 55.5-75 spacings). At most `pairs_per_bin` pairs (a
  whole number, 4000 by default) are drawn at random from each class,
  and a class with fewer than `min_pairs_per_bin` pairs (200 by default)
  is not used. At least 4 usable classes are needed; a part with fewer
  is not tested, its spots get no domain, and [`sigdots()`](sigdots.md)
  warns.

- r_cap:

  The largest number of strong components kept in the correlation
  between genes of the simulated data (a whole number, 30 by default).
  The number actually kept is estimated from the data at every test.

- ridge, nugget:

  Small non-negative numbers added to the diagonals of the gene
  covariance and of the spatial covariance matrices, for numerical
  stability (both 1e-6 by default).

- max_exact_memory_gib:

  For `spatial_method = "exact"`: the largest memory, in GiB (1 GiB =
  2^30 bytes), a test may need (a positive number, 16 by default). If a
  test is estimated to need more, the run stops before it starts and
  suggests the Vecchia approximation or a larger value. 16 GiB is
  reached at about 19,800 spots.

- candidate_seed:

  The random seed used to make the candidate domains (k-means, BANKSY or
  BayesSpace), one number. `NULL` (the default) uses the `seed` of
  [`sigdots()`](sigdots.md). Setting it separately keeps the candidate
  domains fixed while `seed` changes the simulated datasets.

- normal_band:

  For `pvalue_method = "normal_adaptive"`: a number larger than 1 (10 by
  default). A test whose first p-value lies between its level divided by
  `normal_band` and its level multiplied by `normal_band` simulates more
  datasets, up to `normal_n_refine`, before deciding.

## Value

A list of class `"sigdots_control"` holding the settings, for the
`control` argument of [`sigdots()`](sigdots.md). Every value is checked
when the list is made; an invalid one stops with a message naming the
setting.

## Details

Most analyses need none of these changed. The ones you may need:
`ell_range` and `B_cap`, when [`sigdots()`](sigdots.md) warns about
them; `max_exact_memory_gib`, to use the exact computation on a large
sample; `candidate_seed`, to reproduce the candidate domains of an
earlier run. Values other than the defaults have not been evaluated.

## See also

[`sigdots()`](sigdots.md).

## Examples

``` r
sigdots_control()
#> $k_nn
#> [1] 6
#> 
#> $n_comp_min
#> [1] 0.01
#> 
#> $min_node_size
#> [1] 20
#> 
#> $min_child_size
#> [1] 8
#> 
#> $B_cap
#> [1] 5000
#> 
#> $eta_min
#> [1] 0.01
#> 
#> $ell_range
#> [1]  0.3 40.0
#> 
#> $ell_candidates
#> [1] 160
#> 
#> $pairs_per_bin
#> [1] 4000
#> 
#> $min_pairs_per_bin
#> [1] 200
#> 
#> $r_cap
#> [1] 30
#> 
#> $ridge
#> [1] 1e-06
#> 
#> $nugget
#> [1] 1e-06
#> 
#> $max_exact_memory_gib
#> [1] 16
#> 
#> $candidate_seed
#> NULL
#> 
#> $normal_band
#> [1] 10
#> 
#> attr(,"class")
#> [1] "sigdots_control"
sigdots_control(ell_range = c(0.3, 80), B_cap = 20000)
#> $k_nn
#> [1] 6
#> 
#> $n_comp_min
#> [1] 0.01
#> 
#> $min_node_size
#> [1] 20
#> 
#> $min_child_size
#> [1] 8
#> 
#> $B_cap
#> [1] 20000
#> 
#> $eta_min
#> [1] 0.01
#> 
#> $ell_range
#> [1]  0.3 80.0
#> 
#> $ell_candidates
#> [1] 160
#> 
#> $pairs_per_bin
#> [1] 4000
#> 
#> $min_pairs_per_bin
#> [1] 200
#> 
#> $r_cap
#> [1] 30
#> 
#> $ridge
#> [1] 1e-06
#> 
#> $nugget
#> [1] 1e-06
#> 
#> $max_exact_memory_gib
#> [1] 16
#> 
#> $candidate_seed
#> NULL
#> 
#> $normal_band
#> [1] 10
#> 
#> attr(,"class")
#> [1] "sigdots_control"
```
