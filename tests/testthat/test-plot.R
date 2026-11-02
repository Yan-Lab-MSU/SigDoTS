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

test_that("each number is written inside its set, away from the boundaries", {
  # a ring 5 spots wide around a 10 x 10 square: the ring's mean position is the square's centre,
  # so the ring's spot nearest to it lies on the boundary between the two
  co <- as.matrix(expand.grid(x = 1:20, y = 1:20))
  lab <- ifelse(pmax(abs(co[, 1] - 10.5), abs(co[, 2] - 10.5)) < 5, 2L, 1L)
  pos <- SigDoTS:::.label_positions(co, lab)
  expect_identical(pos$label, c("1", "2"))
  j <- match(paste(pos$x, pos$y), paste(co[, 1], co[, 2]))
  expect_identical(lab[j], 1:2)
  # distance from each number's spot to the nearest spot of the other set
  gap <- vapply(1:2, function(k) min(sqrt(colSums((t(co[lab != k, ]) - co[j[k], ])^2))), numeric(1))
  expect_gte(gap[1], 2)   # in the middle of the ring's width, not next to the square
  expect_gte(gap[2], 4)   # near the centre of the square
  # untested spots (NA) count as outside every set and get no number
  lab_na <- replace(lab, co[, 1] <= 2, NA)
  pos_na <- SigDoTS:::.label_positions(co, lab_na)
  expect_identical(pos_na$label, c("1", "2"))
  expect_true(all(pos_na$x > 2))
})
