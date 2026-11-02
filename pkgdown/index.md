<!-- The website's home page: README.md with its header image replaced by the title and the
     subtitle as text. Written when the site is built; edit README.Rmd, not this file. -->


<!-- README.md is generated from README.Rmd. Edit README.Rmd, then knit it. -->

# SigDoTS <img src="man/figures/logo.png" align="right" height="240" alt="SigDoTS logo" />

**Statistical significance testing for domain detection in spatial transcriptomics.**

<!-- badges: start -->
<!-- When the repository is public: the R-CMD-check badge from .github/workflows/R-CMD-check.yaml. -->
<!-- badges: end -->

Domain-detection methods (BANKSY, BayesSpace, Seurat/Scanpy clustering,
…) divide a tissue section into spatial domains, but the number of
domains comes from a parameter, not from a test: the same tissue can be
cut into 3, 5 or 8 domains depending on the setting, and some of the
cuts only separate variation that exists inside one domain. SigDoTS
tests the cuts. Starting from an over-segmentation (more candidate
domains than you expect), it arranges the candidate domains in a tree by
spatial adjacency and expression similarity and tests the splits from
the top of the tree down. Each test compares how clearly the two sides
of a split separate with how clearly the same clustering separates data
simulated from a *single-domain* model, which keeps the gene-gene
covariance and the spatial correlation of the noise estimated from the
data. A split is kept only when it is significant, with the family-wise
error rate controlled over the whole tree; the parts that are not split
further are the domains.

<img src="man/figures/workflow_steps.png" alt="The three steps of SigDoTS: initial candidate domains, a spatially constrained hierarchy, and tests from the root downward" width="100%" />

*The three steps of SigDoTS (Figure 1c of the paper): a domain-detection
pipeline gives `K0` initial candidate domains; spatially adjacent,
expression-similar candidates are merged into a hierarchy; the splits
are tested from the root downward, and a branch stops at the first split
that is not significant.*

## Key features

- **Tests, not thresholds.** Every retained split has passed a test
  against a single-domain null model that includes within-domain spatial
  correlation and gene-gene covariance; the family-wise error rate over
  all tests is controlled at `alpha` (0.05 by default).
- **Accounts for selection.** The same clustering procedure that
  produced the candidate domains is applied to every simulated dataset,
  so the “separation that clustering always finds” is part of the
  reference distribution.
- **Works with your objects.** Takes a Seurat object, a
  SpatialExperiment object, or plain matrices; writes the domains back
  with `add_domains()`.
- **Works with your pipeline.** Start from the built-in k-means, from
  BANKSY or BayesSpace, or from any labels you already have
  (`candidate_labels`).
- **Readable results.** `print()` draws the tested tree in the console;
  `summary()` lists every test with its statistic, p-value and level;
  `plot()` draws the domains, the candidate domains, the tree and the
  estimated spatial correlation.
- **Reproducible and parallel.** Same data, arguments and seed give the
  same result on any number of cores.
- **Scales.** A Vecchia approximation (default) keeps memory roughly
  linear in the number of spots; the exact computation is available for
  smaller samples.

## Installation

SigDoTS needs R 4.3.0 or later. It depends on the CRAN packages GpGp,
fields, igraph, Matrix, RANN and ggplot2, which are installed
automatically when SigDoTS is installed from GitHub.

``` r
# install.packages("remotes")
remotes::install_github("Yan-Lab-MSU/SigDoTS")
```

Or from a downloaded package file (the dependencies are then not
installed automatically):

``` r
install.packages(c("GpGp", "fields", "igraph", "RANN", "ggplot2"))
install.packages("path/to/SigDoTS_0.2.5.8.tar.gz", repos = NULL, type = "source")
```

Installation time, measured on 2026-10-16 on a Linux computer that
already had R 4.3.2 and the CRAN packages above except GpGp: 36 seconds
for GpGp (downloaded and compiled from CRAN) and 5 seconds for SigDoTS.
A computer without igraph, RANN and ggplot2 needs a few minutes more: in
clean Linux containers (2026-11-01; R 4.4.2 and R 4.6.1, empty package
library, binary packages from Posit Package Manager) the line above took
3 to 7 minutes, and SigDoTS 7 to 10 seconds. SigDoTS itself contains no
compiled code. GpGp does, but on Windows and macOS `install.packages()`
normally installs a ready-made binary of it from CRAN, so no compiler is
needed; only if R insists on compiling it from source are Rtools
(Windows) or the Xcode command-line tools (macOS) required.

On Linux, igraph needs the GLPK library of the operating system. If
`library(SigDoTS)` stops with “libglpk.so.40: cannot open shared object
file”, install the system libraries igraph uses and then igraph again,
for example on Ubuntu or Debian:

``` bash
sudo apt-get install libglpk-dev libgmp-dev libxml2-dev
```

(on Fedora: `sudo dnf install glpk-devel gmp-devel libxml2-devel`).

Tested on Windows, macOS and Linux: `R CMD check --as-cran` passes on
Windows (R 4.6.2), macOS (R 4.6.2, Apple silicon) and Ubuntu (R 4.5.3, R
4.6.2 and the development version of R) on GitHub Actions, and in clean
Ubuntu 24.04 containers with R 4.6.1; the package also runs with R
4.3.2, R 4.4.2 and R 4.5.1 on Linux.

Optional, only to make the candidate domains with BANKSY or BayesSpace:

``` r
install.packages("BiocManager")
BiocManager::install(c("Banksy", "SingleCellExperiment", "SummarizedExperiment"))   # or "BayesSpace"
```

## Quick start

The package includes a small simulated dataset: 900 spots on a 30 x 30
grid, 200 genes, three true domains of 300 spots each.

``` r
library(SigDoTS)
dim(sigdots_example$x)          # spots x genes
#> [1] 900 200
time <- system.time(
  fit <- sigdots(sigdots_example$x, sigdots_example$coords, K0 = 5)
)
#> SigDoTS: 900 spots, 200 genes; 5 candidate domains; testing from the top of the tree
#>   test at node 9 (900 spots, 5 candidate domains): split (19 simulated datasets) [10 s]
#>   test at node 7 (300 spots, 2 candidate domains): no split (10 simulated datasets) [2 s]
#>   test at node 8 (600 spots, 3 candidate domains): split (30 simulated datasets) [7 s]
#>   test at node 6 (300 spots, 2 candidate domains): no split (15 simulated datasets) [3 s]
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
```

Run time of this demo: 22 seconds on one core of an AMD EPYC 9654
(Linux). It also ran on a macOS laptop (time not measured).

The three domains match the three true regions:

``` r
table(domain = domains(fit)$domain, true_region = sigdots_example$region)
#>       true_region
#> domain   1   2   3
#>      1   0   0 300
#>      2   0 300   0
#>      3 300   0   0
```

``` r
plot(fit)
```

<img src="man/figures/README-domains-plot-1.png" width="70%" />

``` r
plot(fit, type = "tree")
```

<img src="man/figures/README-tree-plot-1.png" width="90%" />

`summary(fit)` adds a table of the tests (statistic, number of simulated
datasets, p-value, level, decision). The [Getting
started](https://yan-lab-msu.github.io/SigDoTS/articles/SigDoTS.html)
tutorial explains every line of this output.

## Input data format

`sigdots()` takes the data in the form you have it:

| You have                                                           | Call                                            |
|--------------------------------------------------------------------|-------------------------------------------------|
| a **Seurat** object (normalised, with variable genes set)          | `sigdots(obj, K0 = 10)`                         |
| a **SpatialExperiment** object (with a `logcounts` assay)          | `sigdots(spe, K0 = 10, genes = hvg)`            |
| a raw **counts matrix** (genes x spots) and a coordinate table     | `sigdots(sigdots_prepare(counts), xy, K0 = 10)` |
| a normalised **expression matrix** (spots x genes) and coordinates | `sigdots(x, xy, K0 = 10)`                       |

plus one number, `K0`, the number of candidate domains (more than the
number of domains you expect). Instead of `K0` you can give your own
over-segmentation with `candidate_labels`, a vector or the name of a
metadata column. `add_domains(obj, fit)` writes the domains back into a
Seurat or SpatialExperiment object. See [Using SigDoTS on your own
data](https://yan-lab-msu.github.io/SigDoTS/articles/own-data.html).

## Tutorials

- [Getting
  started](https://yan-lab-msu.github.io/SigDoTS/articles/SigDoTS.html):
  the demo above, line by line.
- [Real data
  application](https://yan-lab-msu.github.io/SigDoTS/articles/dlpfc.html):
  a public Visium section of human brain cortex with manually labelled
  layers, from download to the domains, including two regions of
  interest.
- [Using SigDoTS on your own
  data](https://yan-lab-msu.github.io/SigDoTS/articles/own-data.html):
  the five input forms (Seurat, SpatialExperiment, Scanpy file, counts,
  matrices), the choice of the candidate domains, exact or approximate
  spatial computation, run time and memory.
- [Troubleshooting and
  FAQ](https://yan-lab-msu.github.io/SigDoTS/articles/faq.html): every
  error and warning the package can give, and what to do.
- [How SigDoTS
  works](https://yan-lab-msu.github.io/SigDoTS/articles/how-it-works.html):
  the method step by step, without formulas.

All but the real-data application are also installed with the package:
`vignette("SigDoTS")`, `vignette("own-data")`, `vignette("faq")`,
`vignette("how-it-works")`.

## How long does it take?

Measured in the paper’s scalability experiment (simulated data with
2,000 genes and three domains, `K0 = 5`, exact spatial computation, 10
cores; mean of 2 datasets per size; peak memory of the whole run
including the 10 worker processes):

|  Spots | Run time | Peak memory |
|-------:|---------:|------------:|
|    400 |     14 s |     2.7 GiB |
|    900 |     32 s |     3.5 GiB |
|  1,600 |  2.0 min |     5.3 GiB |
|  2,500 |  4.3 min |     6.1 GiB |
|  5,041 |  9.2 min |     9.6 GiB |
| 10,000 |   18 min |    19.8 GiB |

With the default Vecchia approximation, 2,500 spots took a median of 4.0
minutes on three-domain data and 1.4 minutes on one-domain data (10
cores, 680 datasets). Beyond about 20,000 spots use the Vecchia
approximation: the exact computation at 50,176 spots took 12 hours and
206 GiB.

## How SigDoTS works, in brief

1.  **Candidate domains.** k-means on the first 15 principal components
    (or BayesSpace, BANKSY, or your own labels) divides the spots into
    `K0` candidate domains.
2.  **Tree.** Groups that touch in space are merged pairwise, the most
    similar first (Ward’s criterion on the candidate domains’ means in
    the first 15 principal components), into a binary tree.
3.  **Tests from the top.** At each split, the statistic is the
    within-group sum of squares divided by the total sum of squares,
    over all genes (smaller means a clearer split). A single-domain
    model is fitted to the spots of that part of the tree: a
    matrix-normal model whose spatial correlation (strength `eta`, range
    `ell`, Matérn 3/2) is estimated from pairs of spots within the same
    candidate domain, and whose gene covariance is estimated after
    removing each candidate domain’s mean. Datasets simulated from it
    are clustered in the same way, and the p-value is (1 + the number of
    them that split at least as clearly as the data) / (1 + the number
    of simulated datasets); a normal approximation that needs fewer
    simulated datasets is available
    (`pvalue_method = "normal_adaptive"`). A test on m of the n spots
    uses the level `alpha * (m - 1) / (n - 1)`.
4.  **Domains.** A part is split only when its test is significant;
    otherwise it is a domain.

## Reproducing the paper

The code of the simulation study and of the real-data analyses, and the
processed data, will be in a separate repository
(`Yan-Lab-MSU/SigDoTS_source`, to be published with the paper). The
DLPFC Visium data are available through the spatialLIBD package.

## Citation

Manuscript in preparation. Until it is published, please cite the
package: `citation("SigDoTS")`.

## Licence and contact

SigDoTS is released under the MIT licence ([LICENSE.md](LICENSE.md)).
Questions and bug reports: open an issue on GitHub, or write to the
maintainer (`maintainer("SigDoTS")`).
