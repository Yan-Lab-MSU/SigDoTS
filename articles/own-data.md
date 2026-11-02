# Using SigDoTS on your own data

This page covers what to prepare and which options to consider when you
run SigDoTS on a real dataset. Section 1 runs on small objects built
from the example data; the rest is shown without running, because it
needs your data. A complete run on a public Visium section (human
dorsolateral prefrontal cortex, 151673) is the [real data
application](dlpfc.md).

## 1. Input data: what `sigdots()` accepts

[`sigdots()`](../reference/sigdots.md) takes the data in whichever of
these forms you have:

| You have | Call | Where the pieces come from |
|----|----|----|
| a **Seurat** object (Visium, Xenium, Slide-seq, … loaded with [`Load10X_Spatial()`](https://satijalab.org/seurat/reference/Load10X_Spatial.html) or similar) | `sigdots(obj, K0 = 10)` | the normalised `"data"` layer of the default assay; the variable genes from [`FindVariableFeatures()`](https://satijalab.org/seurat/reference/FindVariableFeatures.html); the positions from the object’s image or FOV |
| a **SpatialExperiment** object (Bioconductor: [`read10xVisium()`](https://rdrr.io/pkg/SpatialExperiment/man/read10xVisium.html), spatialLIBD, STexampleData, …) | `sigdots(spe, K0 = 10, genes = hvg)` | the `"logcounts"` assay; the genes you name (all genes if there are at most 5,000); `spatialCoords(spe)` |
| a **SingleCellExperiment** without positions, or an AnnData file read with `zellkonverter::readH5AD()` | `sigdots(sce, coords = xy, K0 = 10)` | the `"logcounts"` assay; `coords` given by you |
| a **counts matrix** (genes x spots) and a coordinate table | `sigdots(sigdots_prepare(counts), xy, K0 = 10)` | [`sigdots_prepare()`](../reference/sigdots_prepare.md) log-normalises and keeps the 2,000 most variable genes |
| an **expression matrix** you normalised yourself | `sigdots(x, xy, K0 = 10)` | `x`: one row per spot, one column per gene; `xy`: one row per spot, two columns, same order |

Whatever the form, two things end up inside: an expression matrix with
one row per spot and one column per gene, normalised and
log-transformed, and the spots’ positions (any unit; SigDoTS measures
distances in units of the typical distance between neighbouring spots).
`candidate_labels` can be the name of a metadata column of the object
(for example a BANKSY or Seurat clustering) instead of a vector.

The examples below build small objects from the package’s example data,
so that they run here; with real data the objects come from the loading
functions named in the table.

All three examples start from the same toy counts matrix, made from the
example data:

``` r

library(SigDoTS)
counts <- Matrix::Matrix(t(round(expm1(sigdots_example$x) * 50)), sparse = TRUE)   # genes x spots
dim(counts)
#> [1] 200 900
```

### A Seurat object

``` r

library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 'SeuratObject' was built under R 4.6.0 but the current version is
#> 4.6.2; it is recomended that you reinstall 'SeuratObject' as the ABI
#> for R may have changed
#> 'SeuratObject' was built with package 'Matrix' 1.7.5 but the current
#> version is 1.7.6; it is recomended that you reinstall 'SeuratObject' as
#> the ABI for 'Matrix' may have changed
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
#> v5.6 adds significant improvements in speed and efficiency, particularly for large datasets.
#> For more information, see https://satijalab.org/seurat/articles/announcements.
# A toy Seurat object: the counts and the spot positions as a FOV. With real data:
# obj <- Load10X_Spatial("path/to/spaceranger/outs")
obj <- CreateSeuratObject(counts = counts, assay = "Spatial")
obj[["slice1"]] <- CreateFOV(CreateCentroids(data.frame(sigdots_example$coords, cell = colnames(counts))),
                             type = "centroids", assay = "Spatial")

obj <- NormalizeData(obj, verbose = FALSE)                       # the "data" layer SigDoTS uses
obj <- FindVariableFeatures(obj, nfeatures = 150, verbose = FALSE)   # the genes SigDoTS uses

fit <- sigdots(obj, K0 = 5, verbose = FALSE)
fit
#> SigDoTS: 900 spots, 150 genes
#> Number of domains: 3   (spots per domain: 300, 295, 305)
#> 
#> Split tree ([node] spots, candidate domains: decision):
#> [9] 900 spots, candidate domains 1,2,3,4,5: split
#> ├── [6] 300 spots, candidate domains 2,3: no split   => DOMAIN 1
#> │   ├── candidate domain 2  109 spots
#> │   └── candidate domain 3  191 spots
#> └── [8] 600 spots, candidate domains 1,4,5: split
#>     ├── [1] 295 spots, candidate domain 1   => DOMAIN 2
#>     └── [7] 305 spots, candidate domains 4,5: no split   => DOMAIN 3
#>         ├── candidate domain 4  154 spots
#>         └── candidate domain 5  151 spots
#> 
#> Use summary() for the tests.
obj <- add_domains(obj, fit)       # adds obj$sigdots_domain and obj$sigdots_candidate_domain
table(obj$sigdots_domain, sigdots_example$region)
#>    
#>       1   2   3
#>   1   0   0 300
#>   2 295   0   0
#>   3   5 300   0
```

With a Seurat object you can also name the assay (`assay = "SCT"`),
choose the genes (`genes = VariableFeatures(obj)[1:1000]`), give the
positions yourself (`coords = ...`), or start from a clustering stored
in the metadata (`candidate_labels = "seurat_clusters"`). Seurat objects
created with Seurat 4 (images of class `VisiumV1`) work the same way;
their coordinates are `imagecol` and `imagerow`, so `plot(fit)` may
appear flipped vertically relative to the tissue image (add
[`ggplot2::scale_y_reverse()`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)).

### A SpatialExperiment object

``` r

library(SpatialExperiment)
#> Loading required package: SingleCellExperiment
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> 
#> Attaching package: 'IRanges'
#> The following object is masked from 'package:sp':
#> 
#>     %over%
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
#> 
#> Attaching package: 'SummarizedExperiment'
#> The following object is masked from 'package:Seurat':
#> 
#>     Assays
#> The following object is masked from 'package:SeuratObject':
#> 
#>     Assays
library(scuttle)
library(scran)
# A toy SpatialExperiment: with real data, spe <- read10xVisium("path/to/spaceranger/outs")
spe <- SpatialExperiment(assays = list(counts = counts), spatialCoords = sigdots_example$coords)

spe <- logNormCounts(spe)                                        # the "logcounts" assay SigDoTS uses
#> Warning in .library_size_factors(assay(x, assay.type), ...): 'librarySizeFactors' is deprecated.
#> Use 'scrapper::centerSizeFactors' instead.
#> See help("Deprecated")
#> Warning in .local(x, ...): 'normalizeCounts' is deprecated.
#> Use 'scrapper::normalizeCounts' instead.
#> See help("Deprecated")
hvg <- getTopHVGs(modelGeneVar(spe), n = 150)                    # the genes SigDoTS uses
#> Warning in getTopHVGs(modelGeneVar(spe), n = 150): 'getTopHVGs' is deprecated.
#> Use 'scrapper::chooseHighlyVariableGenes' instead.
#> See help("Deprecated")
#> Warning in fitTrendVar(fm, fv, ...): 'fitTrendVar' is deprecated.
#> Use 'scrapper::fitVarianceTrend' instead.
#> See help("Deprecated")
#> Warning in combineBlocks(collected, method = method, equiweight = equiweight, : 'combineBlocks' is deprecated.
#> See help("Deprecated")

fit <- sigdots(spe, K0 = 5, genes = hvg, verbose = FALSE)
fit
#> SigDoTS: 900 spots, 94 genes
#> Number of domains: 3   (spots per domain: 300, 299, 301)
#> 
#> Split tree ([node] spots, candidate domains: decision):
#> [9] 900 spots, candidate domains 1,2,3,4,5: split
#> ├── [5] 300 spots, candidate domain 5   => DOMAIN 1
#> └── [8] 600 spots, candidate domains 1,2,3,4: split
#>     ├── [6] 299 spots, candidate domains 3,4: no split   => DOMAIN 2
#>     │   ├── candidate domain 3  160 spots
#>     │   └── candidate domain 4  139 spots
#>     └── [7] 301 spots, candidate domains 1,2: no split   => DOMAIN 3
#>         ├── candidate domain 1  144 spots
#>         └── candidate domain 2  157 spots
#> 
#> Use summary() for the tests.
spe <- add_domains(spe, fit)       # adds colData(spe)$sigdots_domain and $sigdots_candidate_domain
table(colData(spe)$sigdots_domain, sigdots_example$region)
#>    
#>       1   2   3
#>   1   0   0 300
#>   2 299   0   0
#>   3   1 300   0
```

For Visium data `spatialCoords(spe)` holds pixel positions
(`pxl_col_in_fullres`, `pxl_row_in_fullres`); any unit works. A
`SingleCellExperiment` without positions, for example from
`zellkonverter::readH5AD()` on a Scanpy file, needs `coords =`; Scanpy
keeps the positions in `obsm["spatial"]`, which becomes
`reducedDim(sce, "spatial")`.

### A counts matrix

``` r

# counts: genes x spots, raw; xy: spots x 2, rows in the same order as the columns of counts
x <- sigdots_prepare(counts, n_genes = 150)   # log(1 + counts per 10,000), 150 most variable genes
dim(x)
#> [1] 900 150
fit <- sigdots(x, sigdots_example$coords, K0 = 5, verbose = FALSE)
fit
#> SigDoTS: 900 spots, 150 genes
#> Number of domains: 3   (spots per domain: 300, 299, 301)
#> 
#> Split tree ([node] spots, candidate domains: decision):
#> [9] 900 spots, candidate domains 1,2,3,4,5: split
#> ├── [5] 300 spots, candidate domain 5   => DOMAIN 1
#> └── [8] 600 spots, candidate domains 1,2,3,4: split
#>     ├── [6] 299 spots, candidate domains 3,4: no split   => DOMAIN 2
#>     │   ├── candidate domain 3  171 spots
#>     │   └── candidate domain 4  128 spots
#>     └── [7] 301 spots, candidate domains 1,2: no split   => DOMAIN 3
#>         ├── candidate domain 1  154 spots
#>         └── candidate domain 2  147 spots
#> 
#> Use summary() for the tests.
```

### Notes on the genes

Any normalisation and gene selection can be used; the examples take the
most variable genes because that is what the paper did. The genes enter
two places: the principal components used for the candidate domains and
the tree, and the split statistic, which uses **all** the genes given.
Genes with larger variance weigh more in both; scale the columns first
if you want every gene to count equally. Give at least a few hundred
genes; the simulated data need a gene covariance to be estimated.

If an expression matrix has one row per gene, SigDoTS stops and tells
you to use `t(x)`.

## 2. The candidate domains

SigDoTS starts from an over-segmentation and decides which of its
candidate domains to merge. It assumes that **each candidate domain lies
inside one true domain**; it can merge candidate domains, but it cannot
split one.

**With the built-in clustering** give `K0`, the number of candidate
domains, larger than the number of domains you expect (twice as many is
a reasonable start):

``` r

fit <- sigdots(x, coords, K0 = 10)                                   # k-means (default)
fit <- sigdots(x, coords, K0 = 10, cluster_method = "banksy")        # needs Banksy
fit <- sigdots(x, coords, K0 = 10, cluster_method = "bayesspace")    # needs BayesSpace
```

BANKSY and BayesSpace are Bioconductor packages and are not installed
with SigDoTS; if they are missing, SigDoTS stops and prints the
installation command. Arguments for them go in `cluster_args`.

**With your own labels** (for example from Seurat, Scanpy or a BANKSY
run you have already tuned):

``` r

fit <- sigdots(x, coords, candidate_labels = my_labels, cluster_method = "kmeans")
```

The tree is built from your labels. The simulated datasets, however, are
still grouped with `cluster_method`, because SigDoTS must apply a
clustering procedure to them. The test is therefore calibrated for
labels made the way `cluster_method` makes them; when your labels come
from a very different method, keep this in mind when you read the
p-values.

**Size of the parts.** Every tested part of the tree needs enough spots
for the spatial parameters to be estimated (on a regular grid, about 80
or more). With many small candidate domains, the lower parts of the tree
may be too small; SigDoTS then leaves them untested and warns. A smaller
`K0` avoids this.

**A region of interest.** SigDoTS can be run on any subset of spots, for
example one anatomical layer, to ask whether that layer has
substructure. Subset `x` and `coords` to the spots of the region and run
as usual.

## 3. Exact or approximate spatial computation

The simulated datasets need the spots’ spatial correlation matrix, an n
x n matrix for a part of the tree with n spots.

- `spatial_method = "vecchia"` (default): the Vecchia approximation with
  `vecchia_neighbors = 60` neighbours. Memory and time grow roughly
  linearly with n.
- `spatial_method = "exact"`: the exact matrix and its
  eigendecomposition. Memory grows with n^2 (about 16 GiB at 20,000
  spots) and time faster still. SigDoTS stops before a test whose
  estimated memory exceeds `sigdots_control(max_exact_memory_gib = 16)`;
  raise it if your machine has more memory.

## 4. Run time, cores and memory

- `n_cores` runs the simulated datasets of a test in parallel (Linux and
  macOS; on Windows SigDoTS uses one core). The result is the same for
  any number of cores.
- The first test, on all spots, is the slowest; it uses 19 simulated
  datasets at `alpha = 0.05`.

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
cores, 680 datasets). At 50,176 spots the exact computation took 12
hours and 206 GiB; use Vecchia there.

Estimating the spatial parameters uses memory that grows with the square
of the number of spots in a part of the tree, whatever the spatial
method: measured for one estimate with 200 genes, 0.31 GiB at 2,500
spots, 2.05 GiB at 10,000, 7.42 GiB at 19,881 and 28.99 GiB at 40,000.
For very large samples, run SigDoTS on a computer with enough memory, or
on a region of interest.

## 5. Options you will rarely change

[`sigdots_control()`](../reference/sigdots_control.md) holds the
remaining settings; the defaults are those used in the paper. The most
likely to matter:

- `ell_range = c(0.3, 40)`: the range of `ell` (in spot spacings)
  searched. If SigDoTS warns that `ell` is at the upper end, rerun with
  a wider range, for example
  `control = sigdots_control(ell_range = c(0.3, 80))`.
- `min_node_size = 20`, `min_child_size = 8`: parts or sides smaller
  than this are not tested.
- `B_cap = 5000`: the largest number of simulated datasets for one test.
  A test on very few of very many spots has a tiny level and may need
  more; SigDoTS then warns.

Other arguments of [`sigdots()`](../reference/sigdots.md): `alpha`
(0.05), `B_min`, `B_max` (10 and 10), `n_pcs` (15), `seed` (1), and
`pvalue_method` (`"monte_carlo"`, or `"normal_adaptive"` for a normal
approximation of the p-value that needs fewer simulated datasets; see
[`?sigdots`](../reference/sigdots.md)).

## 6. Reading the result

`print(fit)`, `summary(fit)`, `domains(fit)` and `plot(fit, type = ...)`
are explained in [`vignette("SigDoTS")`](../articles/SigDoTS.md).
`add_domains(object, fit)` writes the domain of every spot into a Seurat
or SpatialExperiment object; `domains(fit)` returns one row per spot in
the order of the data for anything else.

For a figure of the domains on the tissue with your own styling, take
`domains(fit)$domain` and the coordinates and draw them with ggplot2, or
add layers to `plot(fit)`, which is a ggplot object.
