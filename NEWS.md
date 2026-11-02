# SigDoTS 0.2.5.8

First version of the package.

* `sigdots()` runs the test, from an expression matrix and coordinates, a Seurat object or a
  SpatialExperiment object; `sigdots_control()` holds the settings that rarely change;
  `sigdots_prepare()` makes the expression matrix from raw counts.
* `print()`, `summary()`, `domains()`, `add_domains()` and `plot()` (domains, candidate domains,
  tree, spatial correlation) read the result.
* The example dataset `sigdots_example`; four tutorials and a real-data application.
* Two ways to compute each test's p-value: Monte Carlo (the default) and a normal
  approximation that needs fewer simulated datasets (`pvalue_method = "normal_adaptive"`).
* Results do not depend on the number of cores.
* Where the spatial parameters cannot be estimated (too few spots in a part of the tree), that
  part is left untested and marked, and the rest of the tree is still tested.
* Where the spots of a part of the tree fall into separate pieces on the nearest-neighbour
  graph, the pieces are joined by their closest spots before the tree is built.
