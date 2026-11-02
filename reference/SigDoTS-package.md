# SigDoTS: significance testing for spatial domains

Domain-detection methods divide a tissue section into spatial domains,
but the number of domains comes from a parameter, not from a test.
SigDoTS tests which of the proposed domains are supported by the data.
Starting from more candidate domains than you expect (`K0` of them), it
arranges them in a tree by spatial adjacency and expression similarity
and tests the splits of the tree from the top down, each against a model
in which the part being split is a single domain, with the spatial
correlation of the noise and the correlation between genes estimated
from the data. The parts of the tree that are not split further are the
domains; the probability that any of the kept splits is false is at most
`alpha`.

## The functions

- [`sigdots()`](sigdots.md) runs the test; it takes an expression matrix
  and coordinates, a Seurat object or a SpatialExperiment object, and
  the number of candidate domains `K0` (or your own candidate domains).

- [`sigdots_control()`](sigdots_control.md) holds the settings that
  rarely need changing.

- [`sigdots_prepare()`](sigdots_prepare.md) makes the expression matrix
  from a raw counts matrix.

- [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html)
  ([`summary.sigdots()`](summary.sigdots.md)) show the split tree and
  the tests.

- [`domains()`](domains.md) gives the domain of every spot as a table;
  [`add_domains()`](add_domains.md) writes it into a Seurat or
  SpatialExperiment object.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  ([`plot.sigdots()`](plot.sigdots.md)) draws the domains, the candidate
  domains, the tree and the spatial correlation.

- [sigdots_example](sigdots_example.md) is a small simulated dataset to
  try the package on.

The tutorials: [`vignette("SigDoTS")`](../articles/SigDoTS.md) (getting
started), [`vignette("own-data")`](../articles/own-data.md) (your own
data), [`vignette("faq")`](../articles/faq.md) (errors, warnings and
questions), [`vignette("how-it-works")`](../articles/how-it-works.md)
(the method step by step). The website, with a real-data application, is
at <https://yan-lab-msu.github.io/SigDoTS/>.

## See also

Useful links:

- <https://github.com/Yan-Lab-MSU/SigDoTS>

- <https://yan-lab-msu.github.io/SigDoTS/>

- Report bugs at <https://github.com/Yan-Lab-MSU/SigDoTS/issues>

## Author

**Maintainer**: Huairui Wang <byeah8887@gmail.com>

Authors:

- Feiyu Cao

- Hui Shen

- Yufeng Liu

- Guanao Yan <gayan@msu.edu>
