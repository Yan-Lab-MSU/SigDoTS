#' Significance test for spatial domains
#'
#' `sigdots()` finds the spatial domains of a spatial transcriptomics sample that are
#' supported by the data. It starts from an over-segmentation of the sample into `K0`
#' *candidate domains* (more than the number of domains you expect), arranges the candidate
#' domains in a tree by spatial adjacency and expression similarity, and tests the splits of
#' the tree from the top down, each against a model in which the part being split is a single
#' domain. A split is kept only when it is significant; the parts of the tree that are not
#' split further are the domains. The probability that the result contains any false split is
#' at most `alpha`. The section "How SigDoTS works" below describes each step.
#'
#' @param x The data, in one of three forms.
#'   * An **expression matrix**: one row per spot (or cell), one column per gene, with the
#'     spot identifiers as row names. The values must already be normalised, for example
#'     log(1 + 10,000 x count / total count of the spot), as `Seurat::NormalizeData()` or
#'     `scuttle::logNormCounts()` give; do not give raw counts ([sigdots_prepare()] makes the
#'     matrix from raw counts). Use the most variable genes, typically 1,000 to 3,000. A
#'     matrix, a data frame or a sparse `Matrix` is accepted. `coords` is then required.
#'   * A **Seurat** object. The expression comes from the normalised `"data"` layer of the
#'     assay given by `assay` (by default the object's default assay), restricted to the genes
#'     in `genes` (by default the variable genes set by `Seurat::FindVariableFeatures()`). The
#'     spot positions come from the object's image (`Seurat::GetTissueCoordinates()`) or FOV
#'     unless `coords` is given.
#'   * A **SpatialExperiment** object, or a SingleCellExperiment or SummarizedExperiment with
#'     `coords` given (for example a Scanpy `.h5ad` file read with `zellkonverter`). The
#'     expression comes from the `"logcounts"` assay (or the assay given by `assay`),
#'     restricted to the genes in `genes` (by default all genes when there are at most 5,000;
#'     otherwise `genes` is required). The positions come from
#'     `SpatialExperiment::spatialCoords()` unless `coords` is given.
#'
#'   The values are used as they are: genes are not scaled to variance 1, so a gene with a
#'   larger variance has more weight in the principal components and in the split statistic.
#' @param coords The positions of the spots: a matrix or data frame with one row per spot, in
#'   the same order as the rows of `x` (or, for an object, with the spot identifiers as row
#'   names), and two columns (x and y). Any unit can be used: SigDoTS measures distances in
#'   units of the typical distance between neighbouring spots (the median distance from a spot
#'   to its nearest other spot). Required for a matrix `x`; for an object it replaces the
#'   positions stored in the object.
#' @param K0 The number of candidate domains, a whole number from 2 up and smaller than the
#'   number of spots. SigDoTS merges candidate domains into domains and never splits a
#'   candidate domain, so `K0` must be larger than the number of domains you expect, and large
#'   enough that every candidate domain lies inside one true domain; see "Choosing K0" below.
#'   Not needed when `candidate_labels` is given (if given as well, it must equal the number of
#'   distinct labels).
#' @param genes A character vector of gene names: the genes to use. For a matrix `x` these must
#'   be column names of `x`; for an object, row names of the object. `NULL` (the default) means
#'   all columns of a matrix, the variable genes of a Seurat object, or all genes of a
#'   SpatialExperiment object that has at most 5,000 genes.
#' @param assay For an object only: the name of the assay to take the expression from, a
#'   character string. `NULL` (the default) means the Seurat object's default assay (its
#'   `"data"` layer), or the `"logcounts"` assay of a SpatialExperiment object.
#' @param cluster_method How the spots are divided into the `K0` candidate domains, and, in
#'   every test, how each simulated dataset is divided into candidate domains. One of
#'   * `"kmeans"` (the default): k-means clustering of the first `n_pcs` principal components
#'     of `x` (10 random starts, at most 100 iterations). Fast; it groups spots by expression
#'     only, so on real tissue its candidate domains can be scattered across the section.
#'   * `"banksy"`: the BANKSY method, which clusters each spot by its own expression together
#'     with the average expression of its neighbours and so gives spatially coherent candidate
#'     domains (k-means on the BANKSY principal components). Used in the SigDoTS paper for real
#'     tissue. Needs the Bioconductor package Banksy (R 4.4 or later).
#'   * `"bayesspace"`: the BayesSpace method, a Bayesian model with a spatial prior designed for
#'     the regular spot arrays of Visium and ST. Needs the Bioconductor package BayesSpace; it
#'     runs once for every simulated dataset, which makes it slow.
#' @param candidate_labels Your own candidate domains, instead of those of `cluster_method`: a
#'   vector with one label per spot (numbers, character strings or a factor, in the order of the
#'   spots, without missing values), or, for an object, the name of a metadata column that holds
#'   them (for example `"seurat_clusters"`). The tree is then built from these labels and `K0`
#'   is their number of distinct values. The simulated datasets are still divided with
#'   `cluster_method`, so the test is calibrated for labels made by that method: when your labels
#'   come from BANKSY or BayesSpace, set `cluster_method = "banksy"` or `"bayesspace"` as well.
#'   `NULL` (the default) means the candidate domains are made with `cluster_method`.
#' @param cluster_args A named list of settings passed to BANKSY or BayesSpace (not used for
#'   k-means). For `"banksy"`: `lambda` (the weight of the neighbourhood, 0.2 by default),
#'   `k_geom` (the number of neighbours, 15), `npcs` (the number of BANKSY principal components,
#'   `min(n_pcs, 20)`) and `use_agf` (`FALSE`). For `"bayesspace"`: `platform` (`"Visium"` or
#'   `"ST"`, `"Visium"` by default), `init.method` (`"mclust"`), `model` (`"t"`), `nrep`
#'   (5000 MCMC iterations), `burn.in` (100) and `gamma` (3). `list()` (the default) uses all
#'   these defaults.
#' @param n_pcs The number of principal components of `x` (a whole number, 15 by default) used
#'   by k-means to make the candidate domains, and to compute the mean of each candidate domain
#'   when the candidate domains are merged into the tree. The split statistic itself uses all
#'   genes.
#' @param spatial_method How the spatial covariance of the simulated data is computed:
#'   * `"vecchia"` (the default): the Vecchia approximation, in which each spot depends directly
#'     only on its `vecchia_neighbors` nearest spots. Memory and time grow in proportion to the
#'     number of spots, so it can be used for any sample size.
#'   * `"exact"`: the exact computation, which stores and decomposes a matrix with one row and
#'     one column per spot of the part being tested. Its memory grows with the square of the
#'     number of spots and its time with the cube; the run stops before it starts if a test would
#'     need more memory than `max_exact_memory_gib` of [sigdots_control()].
#' @param vecchia_neighbors For `spatial_method = "vecchia"`: the number of nearest spots each
#'   spot depends on (a whole number, 60 by default). More is closer to the exact computation
#'   and slower.
#' @param alpha The significance level for the whole tree, a number between 0 and 1 (0.05 by
#'   default): the probability that the result contains any false split is at most `alpha`. A
#'   test on m of the n spots uses the level `alpha * (m - 1) / (n - 1)`.
#' @param pvalue_method How the p-value of each test is computed from the simulated datasets:
#'   * `"monte_carlo"` (the default, used for the main results of the SigDoTS paper): the p-value is
#'     (1 + the number of simulated datasets that split at least as clearly as the data) /
#'     (1 + the number of simulated datasets). It makes no assumption about the distribution of
#'     the split statistic. The number of simulated datasets is set by `B_min` and `B_max`.
#'   * `"normal_adaptive"`: the p-value is the probability, under a normal distribution with the
#'     mean and standard deviation of the simulated statistics, of a statistic at most as large
#'     as the observed one. Each test first simulates `normal_n_initial` datasets; when the
#'     p-value from these lies between the test's level divided by `normal_band` and the level
#'     multiplied by `normal_band` (see [sigdots_control()]), it simulates more, up to
#'     `normal_n_refine` in all, and computes the p-value again from all of them. A test thus
#'     uses 10 or 50 simulated datasets with the defaults, whatever its level, whereas a
#'     Monte Carlo test needs about `1 / level` of them to keep a split; it is therefore
#'     much faster for the lower tests of a large tree. It relies on the simulated statistics
#'     having roughly a normal distribution, and on their mean and standard deviation, estimated
#'     from few datasets.
#' @param B_min,B_max For `pvalue_method = "monte_carlo"`: the smallest and the largest number
#'   of simulated datasets for one test (whole numbers with `1 <= B_min <= B_max`, both 10 by
#'   default).
#'   * A test simulates at least `B_min` datasets and stops as soon as its decision can no
#'     longer change, so a split that is clearly not significant stops early.
#'   * When `1 / (B_max + 1)` is larger than the test's level, the p-value could never reach
#'     the level; `B_max` is then raised for that test to `ceiling(1 / level) - 1` (19 at the
#'     level 0.05), up to `B_cap` of [sigdots_control()].
#'   * With the defaults, every test therefore uses the smallest number of simulated datasets
#'     that can establish its decision, as in the SigDoTS paper's simulations: a split is kept
#'     only when the data split more clearly than every simulated dataset.
#'   * Larger values give finer p-values and a better chance of keeping a true split that is
#'     only moderately clear; the run time grows in proportion. The error guarantee holds for
#'     any values.
#' @param normal_n_initial,normal_n_refine For `pvalue_method = "normal_adaptive"`: the number
#'   of simulated datasets every test starts with (a whole number from 2 up, 10 by default), and
#'   the number it uses in all when its first p-value is near the level (a whole number not
#'   smaller than `normal_n_initial`, 50 by default). The defaults are those described in the
#'   SigDoTS paper. Larger values make the estimated mean and standard deviation more precise.
#' @param n_cores The number of CPU cores used to simulate and cluster the datasets of each
#'   test in parallel (a whole number, 1 by default). The result does not depend on it. On
#'   Windows only 1 core is used, because the parallel computation relies on process forking
#'   (`parallel::mclapply()`).
#' @param seed The random seed, one number (1 by default). The same data, arguments and seed
#'   give the same result on any number of cores. Each simulated dataset has its own seed,
#'   derived from `seed`, the node and the dataset's number. The random-number state of your
#'   session is restored when `sigdots()` returns.
#' @param verbose How much progress to print: `TRUE` (the default) prints one line per test
#'   (its node, number of spots, decision, number of simulated datasets and time); `FALSE`
#'   prints nothing; `2` prints the detailed progress of every computation.
#' @param control The settings that rarely need changing, made with [sigdots_control()].
#'
#' @return An object of class `"sigdots"`. Read it with `print(fit)` (the split tree),
#'   `summary(fit)` (the settings and one row per test; see [summary.sigdots()]), [domains()]
#'   (the domain of every spot as a table), `plot(fit)` (maps, the tree and the spatial
#'   correlation; see [plot.sigdots()]) and [add_domains()] (writes the domains into a Seurat or
#'   SpatialExperiment object). For programming, the object is a list with these components:
#'   \describe{
#'     \item{`domains`}{Integer vector, one value per spot in the order of the input: the
#'       domain of the spot, numbered 1, 2, ... from left to right in the printed tree; `NA` for
#'       a spot in a part of the tree that could not be tested.}
#'     \item{`domain_nodes`}{The node numbers of the domains, in the order of their numbers.}
#'     \item{`untested_nodes`, `untested_spot_node`}{The node numbers of the parts that could
#'       not be tested (their spatial parameters could not be estimated, or too many of their
#'       simulated datasets had a top split with a side smaller than `min_child_size`), and,
#'       for every spot, the node number of the untested part it lies in (`NA` otherwise).}
#'     \item{`candidate_domains`, `candidate_names`}{The candidate domain of every spot,
#'       numbered 1 to `K0`, and the name of each candidate domain (your own label when
#'       `candidate_labels` was given).}
#'     \item{`tree`}{A data frame with one row per node of the tree: `node`, `parent`, `left`
#'       and `right` (node numbers; `NA` where there is none), `is_candidate_domain` (a bottom
#'       end of the tree), `spots`, `candidate_domains` (how many lie under the node) and
#'       `merge_score` (the increase in the within-group sum of squares caused by the merge that
#'       made the node).}
#'     \item{`root`}{The node number of the top of the tree.}
#'     \item{`spot_indices`, `groups_under`}{Lists named by node number: the row numbers of the
#'       spots under each node, and the candidate domains under each node.}
#'     \item{`tests`}{A data frame with one row per node the testing reached, from the top, with
#'       the columns of `summary(fit)$tests` and more: `parent`, `left_spots`, `right_spots`,
#'       `tested`, `skip_reason`, `null_median` (median statistic of the simulated datasets),
#'       `null_max` (the largest number of simulated datasets the test could have used),
#'       `null_replaced` (simulated datasets replaced because their top split left fewer than
#'       `min_child_size` spots on a side; they are not counted in `null_used`),
#'       `clearer` (simulated datasets that split at least as clearly as the data),
#'       `stop_reason`, `p_initial` and `refined` (normal approximation: the p-value from the
#'       first datasets, and whether more were simulated), `ell` and `eta` (the estimated
#'       spatial parameters), five flags for estimates at the ends of their ranges
#'       (`ell_at_lower_end`, `ell_at_upper_end`, `ell_capped`, `eta_at_lower_bound`,
#'       `eta_at_upper_bound`), `backend` and `vecchia_neighbors` (the spatial computation),
#'       `pieces` (sizes of separate pieces, see [summary.sigdots()]), and, for a part whose
#'       spatial parameters could not be estimated, `spatial_fit_reason`,
#'       `spatial_fit_usable_bins` and `min_pairs_per_bin`.}
#'     \item{`null_statistics`}{A list named by node number: the split statistics of the
#'       simulated datasets of each test, in the order they were simulated.}
#'     \item{`correlograms`}{A list named by node number: for each test, a data frame with the
#'       spatial correlation at each distance (`distance` in spot spacings, `correlation`
#'       observed, `fitted` by the model, `pairs` used), as drawn by
#'       `plot(fit, type = "correlogram")`.}
#'     \item{`coords`, `spot_ids`, `n_spots`, `n_genes`}{The positions and identifiers of the
#'       spots, and the size of the data.}
#'     \item{`spot_pieces`}{The sizes of the separate pieces of the spots' nearest-neighbour
#'       graph (one number when all spots are connected).}
#'     \item{`settings`}{A list of the arguments used, including `control`.}
#'     \item{`call`, `elapsed_sec`, `version`}{The call, the run time in seconds and the
#'       version of SigDoTS.}
#'   }
#'
#' @section How SigDoTS works:
#' 1. **Candidate domains.** The spots are divided into `K0` candidate domains with
#'    `cluster_method`, or taken from `candidate_labels`. SigDoTS never splits a candidate
#'    domain; it decides which candidate domains belong together. It therefore assumes that each
#'    candidate domain lies inside one true domain.
#' 2. **The tree.** Each spot is linked to its `k_nn` nearest spots (see [sigdots_control()]);
#'    two candidate domains are neighbours when some of their spots are linked. Of all pairs of
#'    neighbouring groups, the pair whose merge increases the within-group sum of squares the
#'    least (Ward's criterion, computed on the means of the groups in the first `n_pcs`
#'    principal components) is merged, and this is repeated until one group is left. The result
#'    is a binary tree whose top is the whole sample and whose bottom ends are the candidate
#'    domains; every merge in it joins two groups that touch in space.
#' 3. **The statistic.** The splits are tested from the top. The statistic of a split is the
#'    within-group sum of squares of its two sides divided by the total sum of squares of the
#'    part, summed over all genes: a number between 0 and 1, smaller for a clearer split.
#' 4. **The reference.** To judge whether the observed statistic is small, SigDoTS simulates
#'    datasets for the spots of the part from a model in which the part is a single domain (step
#'    5). Each simulated dataset is divided into as many candidate domains as the part has, with
#'    the same `cluster_method`; its tree is built in the same way, and the statistic of its top
#'    split is computed. Because the simulated data go through the same clustering and the same
#'    tree building as the real data, the comparison accounts for the separation that clustering
#'    finds even in data with only one domain.
#' 5. **The single-domain model.** In the model, every spot of the part has the same mean
#'    expression, and the noise is correlated both between genes and between nearby spots.
#'    * The spatial correlation is described by two numbers estimated from the part's own data:
#'      `eta`, the share of the noise variance that is spatially correlated (from `eta_min` to
#'      1), and `ell`, the distance over which the correlation reaches (Matérn 3/2 correlation
#'      function, in spot spacings).
#'    * They are estimated from pairs of spots in the same candidate domain, so that differences
#'      between candidate domains are not taken for spatial correlation: after the mean of each
#'      candidate domain is removed, the correlation between such pairs at each distance is
#'      matched to the correlation the model predicts (`plot(fit, type = "correlogram")` shows
#'      both).
#'    * The correlation between genes is then estimated from the data, with the candidate-domain
#'      means and the fitted spatial correlation removed, as a few strong components (at most
#'      `r_cap`) plus independent noise of equal variance.
#'    * Simulating from the model needs a factor of the spatial covariance matrix of the part's
#'      spots, computed exactly or with the Vecchia approximation (`spatial_method`).
#' 6. **The decision.** The p-value (`pvalue_method`) is compared with the level of the test,
#'    `alpha * (m - 1) / (n - 1)` for a part of m of the n spots. If it is at most the level,
#'    the split is kept and both sides are tested in turn; otherwise the part is a domain. With
#'    these levels, the probability that any kept split is false (the family-wise error rate) is
#'    at most `alpha` (the scheme of Kimes et al., 2017).
#' 7. **Parts that are not tested.** A part with fewer than `min_node_size` spots, a part whose
#'    split would leave fewer than `min_child_size` spots on a side, and a single candidate
#'    domain are domains without a test. A part where the spatial parameters cannot be estimated
#'    (too few pairs of spots in the same candidate domain) is left untested together with
#'    everything below it, and its spots get no domain (`NA`); `sigdots()` warns. Where the
#'    spots of a part fall into separate pieces on the nearest-neighbour graph, the pieces are
#'    joined by their closest spots before the trees are built (see [summary.sigdots()]).
#' 8. **Simulated datasets with a too-small split.** The data's split is tested only when both
#'    of its sides have at least `min_child_size` spots, so the same is asked of every
#'    simulated dataset: one whose top split leaves fewer spots on a side is replaced by the
#'    next one in seed order and is not counted, and the test still uses its number of
#'    simulated datasets (`summary()` lists how many were replaced). If the replaced ones
#'    outnumber the valid ones, the part is left untested together with everything below it,
#'    and its spots get no domain (`NA`); `sigdots()` warns.
#'
#' @section Choosing K0:
#' `K0` must be larger than the number of domains you expect, and fine enough that every
#' candidate domain lies inside one true domain, because a candidate domain that straddles two
#' true domains can never be split. About twice the expected number of domains is a reasonable
#' start. A very large `K0` makes the lower parts of the tree small: the spatial parameters need
#' about 80 or more spots in a part (on a regular grid of spots), and small parts have small
#' levels, so a Monte Carlo test there needs many simulated datasets to keep a split.
#'
#' @section Run time and memory:
#' Measured on simulated data with 2,000 genes and three domains (`K0 = 5`, exact computation,
#' 10 cores; peak memory of the whole run): 900 spots, 32 seconds and 3.5 GiB; 2,500 spots,
#' 4.3 minutes and 6.1 GiB; 10,000 spots, 18 minutes and 19.8 GiB. With the default Vecchia
#' approximation, 2,500 spots took a median of 4.0 minutes (three domains) and 1.4 minutes (one
#' domain) on 10 cores. Beyond about 20,000 spots use the Vecchia approximation: the exact
#' computation at 50,176 spots took 12 hours and 206 GiB. The estimation of the spatial
#' parameters computes the distances between all pairs of spots of the part being tested, so
#' its memory grows with the square of the number of spots whatever `spatial_method` is: one
#' estimation with 200 genes on a part of 10,000 spots used 2.1 GiB, of 19,881 spots 7.4 GiB,
#' and of 40,000 spots 29 GiB.
#'
#' @section Warnings at the end of a run:
#' * *The spatial parameters could not be estimated at ... part(s) of the tree*: a part has
#'   too few spots for the estimation; it is marked `[!] NOT TESTED` and its spots get no
#'   domain (`NA`). A smaller `K0` makes the parts larger.
#' * *The simulated datasets whose top split left fewer than min_child_size spots on a side
#'   outnumbered the valid ones*: the part is marked `[!] NOT TESTED` and its spots get no
#'   domain (`NA`). This happens where the spots of a part lie in many small separate pieces;
#'   see `vignette("faq")`.
#' * *The estimated spatial range ell is at the upper end of its allowed range*: the true
#'   range may be longer, which can make splits easier to keep; rerun with a wider `ell_range` in
#'   [sigdots_control()].
#' * *The test would need more than B_cap simulated datasets*: raise `B_cap` in
#'   [sigdots_control()], or accept that part as a domain.
#'
#' Before computing anything, `sigdots()` checks its input and stops with a message that names
#' the argument at fault (for example, when `x` has as many columns as `coords` has rows, it
#' suggests `t(x)`). A simulated dataset that fails also stops the run, rather than being
#' dropped (one whose top split is too small has not failed: it is replaced, as described
#' above). Every message is explained in `vignette("faq")`.
#'
#' @references
#' Wang H, Cao F, Shen H, Liu Y, Yan G (2026). Statistical significance testing for spatial
#' domain detection in spatial transcriptomics. Manuscript.
#'
#' Kimes PK, Liu Y, Hayes DN, Marron JS (2017). Statistical significance for hierarchical
#' clustering. *Biometrics* 73(3), 811-821.
#'
#' Vecchia AV (1988). Estimation and model identification for continuous spatial processes.
#' *Journal of the Royal Statistical Society Series B* 50(2), 297-312.
#'
#' Katzfuss M, Guinness J (2021). A general framework for Vecchia approximations of Gaussian
#' processes. *Statistical Science* 36(1), 124-141.
#'
#' Singhal V, Chou N, Lee J, et al. (2024). BANKSY unifies cell typing and tissue domain
#' segmentation for scalable spatial omics data analysis. *Nature Genetics* 56(3), 431-441.
#'
#' Zhao E, Stone MR, Ren X, et al. (2021). Spatial transcriptomics at subspot resolution with
#' BayesSpace. *Nature Biotechnology* 39(11), 1375-1384.
#'
#' @seealso [sigdots_control()] for the other settings, [summary.sigdots()] for reading the
#'   result, `vignette("SigDoTS")` for a worked example, `vignette("own-data")` for Seurat and
#'   SpatialExperiment objects and the choice of settings, `vignette("how-it-works")` for the
#'   method with figures.
#' @examples
#' \donttest{
#' # about 20 seconds on one core
#' fit <- sigdots(sigdots_example$x, sigdots_example$coords, K0 = 5)
#' fit
#' summary(fit)
#' table(domain = domains(fit)$domain, true_region = sigdots_example$region)
#' plot(fit)
#' plot(fit, type = "tree")
#'
#' # more simulated datasets per test (finer p-values), on 2 cores
#' fit2 <- sigdots(sigdots_example$x, sigdots_example$coords, K0 = 5,
#'                 B_min = 50, B_max = 50, n_cores = 2)
#'
#' # the normal approximation of the p-value
#' fit3 <- sigdots(sigdots_example$x, sigdots_example$coords, K0 = 5,
#'                 pvalue_method = "normal_adaptive")
#' summary(fit3)
#' }
#' @export
sigdots <- function(x, coords = NULL, K0 = NULL,
                    genes = NULL, assay = NULL,
                    cluster_method = c("kmeans", "bayesspace", "banksy"),
                    candidate_labels = NULL,
                    cluster_args = list(),
                    n_pcs = 15,
                    spatial_method = c("vecchia", "exact"),
                    vecchia_neighbors = 60,
                    alpha = 0.05,
                    pvalue_method = c("monte_carlo", "normal_adaptive"),
                    B_min = 10, B_max = 10,
                    normal_n_initial = 10, normal_n_refine = 50,
                    n_cores = 1,
                    seed = 1,
                    verbose = TRUE,
                    control = sigdots_control()) {
  call <- match.call()
  cluster_method <- match.arg(cluster_method)
  spatial_method <- match.arg(spatial_method)
  pvalue_method <- match.arg(pvalue_method)
  if (!inherits(control, "sigdots_control")) stop("`control` must come from sigdots_control().")

  # The user's random-number state is put back when sigdots() returns.
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  }, add = TRUE)

  ext <- .sigdots_extract(x, coords, genes, assay, candidate_labels)
  inp <- .sigdots_check_inputs(ext$x, ext$coords, K0, ext$candidate_labels, cluster_method, n_pcs,
                               vecchia_neighbors, alpha, B_min, B_max, n_cores, seed,
                               normal_n_initial, normal_n_refine)
  candidate_labels <- ext$candidate_labels
  n_cores <- inp$n_cores
  if (.Platform$OS.type == "windows" && n_cores > 1L) {
    message("On Windows SigDoTS runs on 1 core (parallel::mclapply() needs process forking). ",
            "The result does not depend on the number of cores.")
    n_cores <- 1L
  }

  null_args <- list(
    r_cap = control$r_cap, ridge = control$ridge, nugget = control$nugget,
    spatial_method = spatial_method, max_exact_peak_gib = control$max_exact_memory_gib,
    vecchia_neighbors = as.integer(vecchia_neighbors), eta_min = control$eta_min,
    correlogram_pairs_per_bin = control$pairs_per_bin,
    correlogram_min_pairs = control$min_pairs_per_bin,
    ell_range = control$ell_range, ell_candidates = control$ell_candidates)
  t0 <- proc.time()[["elapsed"]]
  fit <- spa_sigtest_multicluster(
    X = inp$x, coords = inp$coords, K0 = inp$K0,
    B_min = as.integer(B_min), B_max = as.integer(B_max),
    bootstrap_B_max_cap = control$B_cap, alpha = alpha,
    k_nn = control$k_nn, n_comp_min = control$n_comp_min,
    ci_pca_dims = as.integer(n_pcs), min_node_size = control$min_node_size,
    min_child_size = control$min_child_size, seed = as.integer(seed),
    cluster_method = cluster_method, bayesspace_args = cluster_args,
    null_args = null_args, verbose = isTRUE(verbose >= 2), n_cores = n_cores,
    observed_clustering = inp$candidate_labels,
    initial_seed = as.integer(control$candidate_seed %||% seed),
    progress = isTRUE(verbose >= 1) && !isTRUE(verbose >= 2),
    pvalue_method = if (pvalue_method == "monte_carlo") "bootstrap" else "normal_adaptive",
    normal_n_initial = as.integer(normal_n_initial), normal_n_refine = as.integer(normal_n_refine),
    normal_band_multiplier = control$normal_band)
  elapsed <- proc.time()[["elapsed"]] - t0

  settings <- list(K0 = inp$K0, cluster_method = cluster_method, input = ext$source,
    candidate_labels_supplied = !is.null(candidate_labels), n_pcs = as.integer(n_pcs),
    spatial_method = spatial_method, vecchia_neighbors = as.integer(vecchia_neighbors),
    alpha = alpha, pvalue_method = pvalue_method, B_min = as.integer(B_min), B_max = as.integer(B_max),
    normal_n_initial = as.integer(normal_n_initial), normal_n_refine = as.integer(normal_n_refine),
    n_cores = n_cores, seed = as.integer(seed), control = unclass(control))
  res <- .sigdots_result(fit, inp, settings, call, elapsed)
  .sigdots_warnings(res)
  res
}

#' Settings of sigdots() that rarely need changing
#'
#' `sigdots_control()` collects the settings of [sigdots()] that rarely need changing. Its
#' defaults are the values used in the SigDoTS paper. Pass its result as the `control` argument
#' of [sigdots()], for example
#' `sigdots(x, coords, K0 = 10, control = sigdots_control(ell_range = c(0.3, 80)))`; settings
#' not named keep their defaults.
#'
#' @param k_nn The number of nearest spots each spot is linked to (a whole number, 6 by
#'   default; on the hexagonal Visium array each spot has 6 nearest neighbours). Two candidate
#'   domains are neighbours, and can be merged into one part of the tree, when some of their
#'   spots are linked. The same links are used for every simulated dataset.
#' @param n_comp_min A number greater than 0 and at most 1 (0.01 by default). When the spots of
#'   a candidate domain fall into separate pieces on the nearest-neighbour links, pieces with
#'   fewer spots than this share of the candidate domain are ignored when deciding which
#'   candidate domains are neighbours, so that a few stray spots do not make distant candidate
#'   domains neighbours. If ignoring them would leave the candidate domains unconnected, all
#'   pieces are used.
#' @param min_node_size A whole number (20 by default). A part of the tree with fewer spots is
#'   not tested and is reported as a domain.
#' @param min_child_size A whole number (8 by default). A split that would leave fewer spots
#'   than this on either side is not tested, and the part is reported as a domain. A simulated
#'   dataset whose top split leaves fewer spots on a side is replaced by a new one, so that the
#'   test keeps its number of simulated datasets (`summary()` lists how many were replaced). If
#'   such simulated datasets outnumber the valid ones, the part is not tested, its spots get no
#'   domain, and [sigdots()] warns.
#' @param B_cap The largest number of simulated datasets a single Monte Carlo test may use (a
#'   whole number, 5000 by default, at least `B_max`). A test on m of n spots needs about
#'   `(n - 1) / (alpha * (m - 1))` simulated datasets to keep a split, so a test on a small part
#'   of a large sample may need more than `B_cap`; that part is then reported as a domain without
#'   a test, and [sigdots()] warns. Not used by `pvalue_method = "normal_adaptive"`.
#' @param eta_min The smallest value allowed for `eta`, the share of the noise variance that is
#'   spatially correlated (a number greater than 0 and smaller than 1, 0.01 by default; the
#'   largest value is 1).
#' @param ell_range The range searched for `ell`, the distance over which the noise is
#'   spatially correlated: two increasing positive numbers, in spot spacings (the median
#'   distance from a spot to its nearest other spot in the part being tested), `c(0.3, 40)` by
#'   default. `ell` is also kept below twice the largest distance between two spots of the same
#'   candidate domain. When [sigdots()] warns that `ell` is at the upper end of this range, the
#'   true range may be longer, which can make splits easier to keep; rerun with a wider range, for
#'   example `c(0.3, 80)`.
#' @param ell_candidates The number of values of `ell` tried in `ell_range`, evenly spaced on a
#'   log scale (a whole number from 2, 160 by default); the best of them is then refined between
#'   its two neighbouring values.
#' @param pairs_per_bin,min_pairs_per_bin `eta` and `ell` are estimated from pairs of spots in
#'   the same candidate domain, grouped by their distance into classes (one per spot spacing from
#'   0.5 to 30.5 spacings, then 30.5-35.5, 35.5-42.5, 42.5-55.5 and 55.5-75 spacings). At most
#'   `pairs_per_bin` pairs (a whole number, 4000 by default) are drawn at random from each
#'   class, and a class with fewer than `min_pairs_per_bin` pairs (200 by default) is not used.
#'   At least 4 usable classes are needed; a part with fewer is not tested, its spots get no
#'   domain, and [sigdots()] warns.
#' @param r_cap The largest number of strong components kept in the correlation between genes
#'   of the simulated data (a whole number, 30 by default). The number actually kept is
#'   estimated from the data at every test.
#' @param ridge,nugget Small non-negative numbers added to the diagonals of the gene
#'   covariance and of the spatial covariance matrices, for numerical stability (both 1e-6 by
#'   default).
#' @param max_exact_memory_gib For `spatial_method = "exact"`: the largest memory, in GiB
#'   (1 GiB = 2^30 bytes), a test may need (a positive number, 16 by default). If a test is
#'   estimated to need more, the run stops before it starts and suggests the Vecchia
#'   approximation or a larger value. 16 GiB is reached at about 19,800 spots.
#' @param candidate_seed The random seed used to make the candidate domains (k-means, BANKSY or
#'   BayesSpace), one number. `NULL` (the default) uses the `seed` of [sigdots()]. Setting it
#'   separately keeps the candidate domains fixed while `seed` changes the simulated datasets.
#' @param normal_band For `pvalue_method = "normal_adaptive"`: a number larger than 1 (10 by
#'   default). A test whose first p-value lies between its level divided by `normal_band` and its
#'   level multiplied by `normal_band` simulates more datasets, up to `normal_n_refine`, before
#'   deciding.
#'
#' @details Most analyses need none of these changed. The ones you may need: `ell_range` and
#'   `B_cap`, when [sigdots()] warns about them; `max_exact_memory_gib`, to use the exact
#'   computation on a large sample; `candidate_seed`, to reproduce the candidate domains of an
#'   earlier run. Values other than the defaults have not been evaluated.
#'
#' @return A list of class `"sigdots_control"` holding the settings, for the `control` argument
#'   of [sigdots()]. Every value is checked when the list is made; an invalid one stops with a
#'   message naming the setting.
#' @seealso [sigdots()].
#' @examples
#' sigdots_control()
#' sigdots_control(ell_range = c(0.3, 80), B_cap = 20000)
#' @export
sigdots_control <- function(k_nn = 6, n_comp_min = 0.01, min_node_size = 20,
                            min_child_size = 8, B_cap = 5000, eta_min = 0.01,
                            ell_range = c(0.3, 40), ell_candidates = 160,
                            pairs_per_bin = 4000, min_pairs_per_bin = 200,
                            r_cap = 30, ridge = 1e-6, nugget = 1e-6,
                            max_exact_memory_gib = 16, candidate_seed = NULL,
                            normal_band = 10) {
  chk <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
  pos_int <- function(v) length(v) == 1L && is.finite(v) && v >= 1 && v == round(v)
  chk(pos_int(k_nn), "k_nn must be a whole number >= 1.")
  chk(length(n_comp_min) == 1L && is.finite(n_comp_min) && n_comp_min > 0 && n_comp_min <= 1,
      "n_comp_min must be in (0, 1].")
  chk(pos_int(min_node_size), "min_node_size must be a whole number >= 1.")
  chk(pos_int(min_child_size), "min_child_size must be a whole number >= 1.")
  chk(pos_int(B_cap), "B_cap must be a whole number >= 1.")
  chk(length(eta_min) == 1L && is.finite(eta_min) && eta_min > 0 && eta_min < 1,
      "eta_min must be in (0, 1).")
  chk(length(ell_range) == 2L && all(is.finite(ell_range)) && ell_range[1] > 0 &&
        ell_range[2] > ell_range[1], "ell_range must be two increasing positive numbers.")
  chk(pos_int(ell_candidates) && ell_candidates >= 2, "ell_candidates must be a whole number >= 2.")
  chk(pos_int(pairs_per_bin), "pairs_per_bin must be a whole number >= 1.")
  chk(pos_int(min_pairs_per_bin), "min_pairs_per_bin must be a whole number >= 1.")
  chk(pos_int(r_cap), "r_cap must be a whole number >= 1.")
  chk(length(ridge) == 1L && is.finite(ridge) && ridge >= 0, "ridge must be >= 0.")
  chk(length(nugget) == 1L && is.finite(nugget) && nugget >= 0, "nugget must be >= 0.")
  chk(length(max_exact_memory_gib) == 1L && is.finite(max_exact_memory_gib) &&
        max_exact_memory_gib > 0, "max_exact_memory_gib must be > 0.")
  chk(is.null(candidate_seed) || (length(candidate_seed) == 1L && is.finite(candidate_seed)),
      "candidate_seed must be NULL or one number.")
  chk(length(normal_band) == 1L && is.finite(normal_band) && normal_band > 1,
      "normal_band must be one number larger than 1.")
  structure(list(k_nn = as.integer(k_nn), n_comp_min = n_comp_min,
                 min_node_size = as.integer(min_node_size),
                 min_child_size = as.integer(min_child_size), B_cap = as.integer(B_cap),
                 eta_min = eta_min, ell_range = as.numeric(ell_range),
                 ell_candidates = as.integer(ell_candidates),
                 pairs_per_bin = as.integer(pairs_per_bin),
                 min_pairs_per_bin = as.integer(min_pairs_per_bin), r_cap = as.integer(r_cap),
                 ridge = ridge, nugget = nugget, max_exact_memory_gib = max_exact_memory_gib,
                 candidate_seed = candidate_seed, normal_band = as.numeric(normal_band)),
            class = "sigdots_control")
}

# Checks every input before anything is computed, so that a mistake is reported at once and
# in terms of the user's arguments.
.sigdots_check_inputs <- function(x, coords, K0, candidate_labels, cluster_method, n_pcs,
                                  vecchia_neighbors, alpha, B_min, B_max, n_cores, seed,
                                  normal_n_initial = 10, normal_n_refine = 50) {
  stop_ <- function(...) stop(..., call. = FALSE)
  if (inherits(x, "Matrix") || is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) stop_("`x` must be a numeric matrix: one row per spot, one column per gene.")
  if (inherits(coords, "Matrix") || is.data.frame(coords)) coords <- as.matrix(coords)
  if (!is.matrix(coords) || !is.numeric(coords)) stop_("`coords` must be a numeric matrix with one row per spot and two columns.")
  if (ncol(coords) != 2L) stop_(sprintf("`coords` must have exactly two columns (it has %d).", ncol(coords)))
  if (nrow(x) != nrow(coords)) {
    if (ncol(x) == nrow(coords))
      stop_(sprintf(paste0("`x` has %d rows but `coords` has %d. `x` has %d columns, so it is probably ",
                           "one row per gene; use t(x) to give one row per spot."), nrow(x), nrow(coords), ncol(x)))
    stop_(sprintf("`x` has %d rows (spots) but `coords` has %d.", nrow(x), nrow(coords)))
  }
  n <- nrow(x)
  if (n < 4L) stop_("At least 4 spots are needed.")
  if (anyNA(x) || any(!is.finite(x))) stop_("`x` contains missing or infinite values.")
  if (anyNA(coords) || any(!is.finite(coords))) stop_("`coords` contains missing or infinite values.")
  nn <- RANN::nn2(coords, k = 2L)$nn.dists[, 2L]
  if (stats::median(nn) <= 0)
    stop_(paste0("At least half of the spots share their coordinates with another spot (the median distance ",
                 "from a spot to its nearest other spot is 0). SigDoTS needs distinct spot positions."))

  if (!is.null(candidate_labels)) {
    if (length(candidate_labels) != n)
      stop_(sprintf("`candidate_labels` has %d values; it needs one per spot (%d).", length(candidate_labels), n))
    if (anyNA(candidate_labels)) stop_("`candidate_labels` contains missing values.")
    k_lab <- length(unique(as.character(candidate_labels)))
    if (!is.null(K0) && as.integer(K0) != k_lab)
      stop_(sprintf("K0 = %d, but `candidate_labels` has %d distinct labels.", as.integer(K0), k_lab))
    K0 <- k_lab
  }
  if (is.null(K0)) stop_("Give K0, the number of candidate domains (or `candidate_labels`).")
  if (length(K0) != 1L || !is.finite(K0) || K0 < 2 || K0 != round(K0)) stop_("K0 must be a whole number >= 2.")
  if (K0 >= n) stop_(sprintf("K0 = %d must be smaller than the number of spots (%d).", as.integer(K0), n))
  pos_int <- function(v) length(v) == 1L && is.finite(v) && v >= 1 && v == round(v)
  if (!pos_int(n_pcs)) stop_("n_pcs must be a whole number >= 1.")
  if (!pos_int(vecchia_neighbors)) stop_("vecchia_neighbors must be a whole number >= 1.")
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) stop_("alpha must be in (0, 1).")
  if (!pos_int(B_min) || !pos_int(B_max) || B_max < B_min) stop_("B_min and B_max must be whole numbers with 1 <= B_min <= B_max.")
  if (!pos_int(normal_n_initial) || normal_n_initial < 2)
    stop_("normal_n_initial must be a whole number >= 2.")
  if (!pos_int(normal_n_refine) || normal_n_refine < normal_n_initial)
    stop_("normal_n_refine must be a whole number >= normal_n_initial.")
  if (!pos_int(n_cores)) stop_("n_cores must be a whole number >= 1.")
  if (length(seed) != 1L || !is.finite(seed)) stop_("seed must be one number.")

  needs <- switch(cluster_method,
    bayesspace = c("BayesSpace", "SingleCellExperiment", "SummarizedExperiment"),
    banksy = c("Banksy", "SingleCellExperiment", "SummarizedExperiment"),
    character(0))
  missing_pkgs <- needs[!vapply(needs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs))
    stop_(sprintf(paste0("cluster_method = \"%s\" needs the Bioconductor package(s) %s. Install with:\n",
                         "  install.packages(\"BiocManager\")\n  BiocManager::install(c(%s))"),
                  cluster_method, paste(missing_pkgs, collapse = ", "),
                  paste0("\"", missing_pkgs, "\"", collapse = ", ")))

  spot_ids <- rownames(x) %||% as.character(seq_len(n))
  list(x = x, coords = coords, K0 = as.integer(K0), candidate_labels = candidate_labels,
       n_cores = as.integer(n_cores), spot_ids = spot_ids)
}
