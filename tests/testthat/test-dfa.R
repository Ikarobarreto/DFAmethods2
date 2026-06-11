test_that("dfa() returns the staged-C contract: $F (legacy F^2), $F2, $F_sqrt, $alpha", {
  set.seed(1)
  x <- cumsum(rnorm(2000))                       # random walk: alpha ~ 1.5
  fy <- dfa(x, np = 30, overlap = FALSE)
  expect_type(fy, "list")
  expect_named(fy, c("s", "F", "F2", "F_sqrt", "alpha"))
  # In the 1.x line $F is still the squared fluctuation (= $F2), preserved for
  # backward compatibility. $F_sqrt is the Peng-convention sqrt.
  expect_identical(fy$F, fy$F2)
  expect_equal(sqrt(fy$F2), fy$F_sqrt)
  # alpha equals the slope of log F_sqrt vs log s (Peng convention)
  pos <- fy$F_sqrt > 0
  slope_peng <- unname(coef(stats::lm(log(fy$F_sqrt[pos]) ~ log(fy$s[pos])))[[2]])
  expect_equal(fy$alpha, slope_peng)
  # Equivalently: alpha is half the slope of log F^2 vs log s
  slope_F2 <- unname(coef(stats::lm(log(fy$F2[pos]) ~ log(fy$s[pos])))[[2]])
  expect_equal(fy$alpha, 0.5 * slope_F2)
})

test_that("integrated random walk has DFA exponent near 1.5 (within sampling noise)", {
  set.seed(2)
  x <- cumsum(rnorm(4000))
  fy <- dfa(x, np = 40, overlap = FALSE)
  expect_true(abs(fy$alpha - 1.5) < 0.1)
})

test_that("stationary white noise has DFA exponent near 0.5", {
  set.seed(3)
  z <- rnorm(4000)
  fy <- dfa(z, np = 40, overlap = FALSE)
  expect_true(abs(fy$alpha - 0.5) < 0.1)
})

test_that("plotdfa() works on the new dfa() output", {
  set.seed(4)
  x <- cumsum(rnorm(800))
  expect_silent(suppressMessages(p <- plotdfa(dfa(x, np = 20, overlap = FALSE))))
  expect_s3_class(p, "ggplot")
})

test_that("rhodcca coefficient columns stay within [-1, 1]", {
  set.seed(42)
  x <- cumsum(rnorm(500))
  y <- x + cumsum(rnorm(500))
  out <- rhodcca(cbind(x, y))

  expect_s3_class(out, "data.frame")
  expect_true("s" %in% names(out))

  coef_cols <- setdiff(names(out), "s")
  rho_values <- unlist(out[coef_cols])
  rho_values <- rho_values[is.finite(rho_values)]
  expect_true(length(rho_values) > 0)
  expect_true(all(rho_values >= -1 - 1e-8 & rho_values <= 1 + 1e-8))
})
