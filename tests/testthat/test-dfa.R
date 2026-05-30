test_that("dfa returns a tibble of scales and fluctuations of matching length", {
  set.seed(1)
  x <- cumsum(rnorm(500))
  out <- dfa(x)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("s", "F"))
  expect_equal(length(out$s), length(out$F))
  expect_true(any(is.finite(out$s)))
})

test_that("rhodcca coefficient column stays within [-1, 1]", {
  set.seed(42)
  x <- cumsum(rnorm(500))
  y <- x + cumsum(rnorm(500))
  out <- rhodcca(cbind(x, y))

  expect_s3_class(out, "data.frame")
  expect_true("s" %in% names(out))

  # every column other than the scale column holds rho-DCCA coefficients
  coef_cols <- setdiff(names(out), "s")
  rho_values <- unlist(out[coef_cols])
  rho_values <- rho_values[is.finite(rho_values)]
  expect_true(length(rho_values) > 0)
  expect_true(all(rho_values >= -1 - 1e-8 & rho_values <= 1 + 1e-8))
})
