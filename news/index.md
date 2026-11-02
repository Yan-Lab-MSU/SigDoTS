# Changelog

## SigDoTS 0.2.5.8

First version of the package.

- [`sigdots()`](../reference/sigdots.md) runs the test, from an
  expression matrix and coordinates, a Seurat object or a
  SpatialExperiment object;
  [`sigdots_control()`](../reference/sigdots_control.md) holds the
  settings that rarely change;
  [`sigdots_prepare()`](../reference/sigdots_prepare.md) makes the
  expression matrix from raw counts.
- [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`domains()`](../reference/domains.md),
  [`add_domains()`](../reference/add_domains.md) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (domains,
  candidate domains, tree, spatial correlation) read the result.
- The example dataset `sigdots_example`; four tutorials and a real-data
  application.
- Two ways to compute each test’s p-value: Monte Carlo (the default) and
  a normal approximation that needs fewer simulated datasets
  (`pvalue_method = "normal_adaptive"`).
- Results do not depend on the number of cores.
- Where the spatial parameters cannot be estimated (too few spots in a
  part of the tree), that part is left untested and marked, and the rest
  of the tree is still tested.
- Where the spots of a part of the tree fall into separate pieces on the
  nearest-neighbour graph, the pieces are joined by their closest spots
  before the tree is built.
