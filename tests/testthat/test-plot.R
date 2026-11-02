test_that("every plot type draws, with and without untested parts", {
  d <- make_bands()   # 12 x 12: two lower parts cannot be tested
  fit <- suppressWarnings(sigdots(d$x, d$coords, K0 = 5, B_min = 5, B_max = 5, verbose = FALSE))
  for (type in c("domains", "candidates", "tree", "correlogram")) {
    p <- plot(fit, type = type)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
  expect_error(plot(fit, type = "null"), "should be one of")
  # a test that was not done is refused with the list of tests that were
  skipped <- fit$tests$node[fit$tests$decision == "not tested: spatial fit"][1]
  expect_error(plot(fit, type = "correlogram", test = skipped), "was not tested")
  # the fitted correlogram is stored for every test that was done
  for (cg in fit$correlograms) expect_false(anyNA(cg$fitted))
  # the correlogram's labels are the plain ones
  p <- plot(fit, type = "correlogram")
  expect_identical(p$labels$x, "Distance")
  expect_identical(p$labels$y, "Spatial correlation")
  expect_identical(p$labels$title, "Spatial correlation")
  expect_null(p$labels$caption)
})
