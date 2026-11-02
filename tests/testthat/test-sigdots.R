# The datasets come from make_bands() (helper-data.R). On the 20 x 20 grid every part of the tree
# is large enough for the spatial fit; on the 12 x 12 grid the two lower tests have 48 spots each, too few
# for the spatial fit (about 80 are needed on a grid), while the first test still runs.

test_that("a 20 x 20 three-band dataset runs end to end and the accessors agree", {
  d <- make_bands(nx = 20, ny = 20)
  fit <- sigdots(d$x, d$coords, K0 = 5, B_min = 5, B_max = 5, verbose = FALSE)
  expect_s3_class(fit, "sigdots")
  dm <- domains(fit)
  expect_equal(nrow(dm), nrow(d$x))
  expect_identical(dm$spot, rownames(d$x))
  expect_false(anyNA(dm$domain))
  expect_length(fit$untested_nodes, 0L)
  expect_setequal(unique(dm$domain), seq_along(fit$domain_nodes))
  expect_identical(names(dm), c("spot", "domain", "candidate_domain"))
  # the domains are exactly the spots under the domain nodes, in the tree's order
  for (k in seq_along(fit$domain_nodes))
    expect_setequal(fit$spot_ids[fit$spot_indices[[as.character(fit$domain_nodes[k])]]], dm$spot[dm$domain == k])
  st <- summary(fit)$tests
  expect_identical(names(st), c("node", "spots", "candidate_domains", "statistic", "null_used", "p_value", "level", "decision"))
  tt <- fit$tests
  expect_true(tt$tested[1])
  # the p-value is (1 + clearer) / (1 + null_used) at every test
  ok <- tt$tested
  expect_equal(st$p_value[ok], (1 + tt$clearer[ok]) / (1 + st$null_used[ok]))
  # each test's level is alpha * (m - 1) / (n - 1)
  expect_equal(st$level[ok], 0.05 * (st$spots[ok] - 1) / (nrow(d$x) - 1))
  # a node is split exactly when its decision says so, and the domains are the unsplit parts
  expect_identical(st$decision[ok] == "split", st$p_value[ok] <= st$level[ok])
  expect_output(print(fit), "Number of domains")
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("^SigDoTS 0\\.2\\.5\\.8", out)))
  expect_true(any(grepl("Tests, from the top:", out)))
  expect_false(any(grepl("eta|ell|clearer|null_median", out)))
  expect_identical(fit$version, as.character(packageVersion("SigDoTS")))
})

test_that("ell at the upper end of its allowed range is flagged and warned about", {
  # the example data were simulated with ell = 3 spot spacings; a range that ends at 1.5 cannot hold it
  keep <- sigdots_example$coords[, "y"] < 0.5   # the lower half of the grid: 15 rows of 30 spots (coordinates run from 0 to 1)
  expect_warning(
    fit <- sigdots(sigdots_example$x[keep, ], sigdots_example$coords[keep, ], K0 = 3, B_min = 3, B_max = 3,
                   verbose = FALSE, control = sigdots_control(ell_range = c(0.3, 1.5))),
    "upper end of its allowed range \\(1\\.5\\).*ell_range = c\\(0\\.3, 3\\)")
  expect_true(fit$tests$ell_at_upper_end[1])
  expect_lte(fit$tests$ell[1], 1.5 * (1 + 1e-8))
  # the default range holds it: no such warning, and the estimate is well inside
  expect_silent(fit2 <- sigdots(sigdots_example$x[keep, ], sigdots_example$coords[keep, ], K0 = 3, B_min = 3, B_max = 3,
                                verbose = FALSE))
  expect_false(fit2$tests$ell_at_upper_end[1])
  expect_gt(fit2$tests$ell[1], 1.5)
})

test_that("the result does not depend on the number of cores, and the user's seed is kept", {
  skip_on_os("windows")
  d <- make_bands()
  set.seed(123); before <- .Random.seed
  f1 <- suppressWarnings(sigdots(d$x, d$coords, K0 = 5, B_min = 5, B_max = 5, n_cores = 1, verbose = FALSE))
  expect_identical(.Random.seed, before)
  f2 <- suppressWarnings(sigdots(d$x, d$coords, K0 = 5, B_min = 5, B_max = 5, n_cores = 2, verbose = FALSE))
  expect_identical(f1$tests, f2$tests)
  expect_identical(f1$null_statistics, f2$null_statistics)
  expect_identical(f1$domains, f2$domains)
  # 3 cores as well, except under R CMD check --as-cran, which allows at most 2 (parallel stops with
  # "simultaneous processes spawned" when _R_CHECK_LIMIT_CORES_ is set).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_"))
  if (!nzchar(chk) || chk == "false") {
    f3 <- suppressWarnings(sigdots(d$x, d$coords, K0 = 5, B_min = 5, B_max = 5, n_cores = 3, verbose = FALSE))
    expect_identical(f1$tests, f3$tests)
    expect_identical(f1$null_statistics, f3$null_statistics)
    expect_identical(f1$domains, f3$domains)
  }
})

test_that("the normal approximation computes its p-values as defined, on any number of cores", {
  skip_on_os("windows")
  d <- make_bands(nx = 20, ny = 20)
  f1 <- suppressWarnings(sigdots(d$x, d$coords, K0 = 5, pvalue_method = "normal_adaptive", verbose = FALSE))
  f2 <- suppressWarnings(sigdots(d$x, d$coords, K0 = 5, pvalue_method = "normal_adaptive", n_cores = 2,
                                 verbose = FALSE))
  expect_identical(f1$tests, f2$tests)
  expect_identical(f1$null_statistics, f2$null_statistics)
  t <- f1$tests[f1$tests$tested, ]
  expect_gt(nrow(t), 0L)
  for (i in seq_len(nrow(t))) {
    nul <- f1$null_statistics[[as.character(t$node[i])]]
    # 10 simulated datasets, or 50 when the first p-value lies within a factor of 10 of the level
    near <- t$p_initial[i] >= t$level[i] / 10 && t$p_initial[i] <= 10 * t$level[i]
    expect_identical(t$refined[i], near)
    expect_identical(length(nul), if (near) 50L else 10L)
    expect_identical(t$null_used[i], length(nul))
    expect_equal(t$p_initial[i], stats::pnorm(t$statistic[i], mean(nul[1:10]), stats::sd(nul[1:10])))
    expect_equal(t$p_value[i], stats::pnorm(t$statistic[i], mean(nul), stats::sd(nul)))
    expect_identical(t$decision[i] == "split", t$p_value[i] <= t$level[i])
  }
  # the first 10 simulated datasets of a test are those of the Monte Carlo run with the same seed
  fm <- suppressWarnings(sigdots(d$x, d$coords, K0 = 5, verbose = FALSE))
  root <- as.character(f1$root)
  expect_identical(f1$null_statistics[[root]][1:10], fm$null_statistics[[root]][1:10])
  expect_identical(f1$settings$pvalue_method, "normal_adaptive")
  expect_output(print(summary(f1)),
                "normal approximation \\(10 simulated datasets per test, 50 when p is near the level\\)")
  expect_output(print(summary(fm)), "P-values: Monte Carlo")
  expect_error(sigdots(d$x, d$coords, K0 = 5, pvalue_method = "normal_adaptive", normal_n_initial = 1),
               "normal_n_initial")
  expect_error(sigdots(d$x, d$coords, K0 = 5, normal_n_initial = 20, normal_n_refine = 10), "normal_n_refine")
  expect_error(sigdots_control(normal_band = 1), "normal_band")
})

test_that("the exact backend runs and supplied labels are used as the candidate domains", {
  d <- make_bands(nx = 20, ny = 20)
  labs <- stats::kmeans(prcomp(d$x)$x[, 1:5], centers = 5, nstart = 5)$cluster
  fit <- sigdots(d$x, d$coords, candidate_labels = paste0("g", labs), spatial_method = "exact",
                 B_min = 5, B_max = 5, verbose = FALSE)
  expect_equal(fit$settings$K0, 5L)
  expect_true(fit$settings$candidate_labels_supplied)
  expect_identical(fit$tests$backend[1], "dense_exact")
  # the user's own labels come back, spot by spot
  expect_identical(domains(fit)$candidate_domain, paste0("g", labs))
  expect_output(print(summary(fit)), "labels supplied by the user")
})

test_that("where the spatial parameters cannot be estimated below the top, those parts are marked and the rest is tested", {
  d <- make_bands()   # 12 x 12
  expect_warning(fit <- sigdots(d$x, d$coords, K0 = 5, B_min = 5, B_max = 5, verbose = FALSE),
                 "could not be estimated")
  tt <- summary(fit)$tests
  expect_identical(tt$decision[1], "split")
  expect_true(any(tt$decision == "not tested: spatial fit"))
  expect_setequal(fit$untested_nodes, tt$node[tt$decision == "not tested: spatial fit"])
  # exactly the spots under the untested parts have no domain
  under <- unlist(lapply(fit$untested_nodes, function(id) fit$spot_ids[fit$spot_indices[[as.character(id)]]]))
  expect_setequal(domains(fit)$spot[is.na(domains(fit)$domain)], under)
  expect_output(print(fit), "[!] NOT TESTED", fixed = TRUE)
  expect_output(print(summary(fit)), "Not tested, because the spatial parameters could not be estimated")
})

test_that("where the spatial parameters cannot be estimated at the top, the run still returns, with no domain", {
  d <- make_bands(nx = 6, ny = 5, p = 20, q_signal = 10)   # 30 spots
  expect_warning(fit <- sigdots(d$x, d$coords, K0 = 3, B_min = 5, B_max = 5, verbose = FALSE),
                 "could not be estimated")
  tt <- summary(fit)$tests
  expect_identical(tt$decision[1], "not tested: spatial fit")
  expect_true(all(is.na(domains(fit)$domain)))
  expect_length(fit$domain_nodes, 0L)
  expect_output(print(fit), "NOT TESTED")
  expect_output(print(fit), "Number of domains: 0")
})

test_that("spots that form separate pieces are joined for the tree, and the run completes", {
  d <- make_bands(nx = 20, ny = 20)
  cut <- d$coords[, 1] %in% c(10, 11)              # remove two columns: two pieces of 180 spots
  fit <- sigdots(d$x[!cut, ], d$coords[!cut, ], K0 = 4, B_min = 5, B_max = 5, verbose = FALSE)
  expect_identical(summary(fit)$tests$decision[1], "split")
  expect_identical(fit$spot_pieces, c(180L, 180L))
  expect_false(anyNA(domains(fit)$domain))
  expect_output(print(summary(fit)), "all spots: pieces of 180, 180")
  # the joining itself: two pieces become one, by one added edge
  g <- SigDoTS:::build_knn_graph(d$coords[!cut, ], 6L)
  expect_equal(igraph::components(g$graph)$no, 2L)
  g2 <- SigDoTS:::.connect_pieces(g, d$coords[!cut, ])
  expect_equal(igraph::components(g2$graph)$no, 1L)
  expect_equal(sum(g2$adj) - sum(g$adj), 2L)
  expect_identical(SigDoTS:::.connect_pieces(g2, d$coords[!cut, ])$adj, g2$adj)
})

test_that("input mistakes are reported in terms of the user's arguments", {
  d <- make_bands()
  expect_error(sigdots(t(d$x), d$coords, K0 = 5, verbose = FALSE), "use t\\(x\\)")
  expect_error(sigdots(d$x, cbind(d$coords, 1), K0 = 5, verbose = FALSE), "exactly two columns")
  expect_error(sigdots(d$x, d$coords, verbose = FALSE), "Give K0")
  expect_error(sigdots(d$x, d$coords, K0 = 4, candidate_labels = rep(1:5, length.out = nrow(d$x)), verbose = FALSE),
               "distinct labels")
  dup <- d$coords; dup[seq(2, nrow(dup), 2), ] <- dup[seq(1, nrow(dup), 2), ]
  expect_error(sigdots(d$x, dup, K0 = 5, verbose = FALSE), "share their coordinates")
  expect_error(sigdots_control(ell_range = c(5, 1)), "increasing")
  # the exact backend stops, naming the options, when its memory estimate exceeds the limit
  expect_error(sigdots(d$x, d$coords, K0 = 5, spatial_method = "exact", verbose = FALSE,
                      control = sigdots_control(max_exact_memory_gib = 1e-6)), "max_exact_memory_gib")
})
