test_that("fracreg variance: 'inverse' equals 'marginal' for a single predictor", {
  set.seed(1)
  n <- 400
  d <- cbind(y = cumsum(rnorm(n)), x = cumsum(rnorm(n)))
  vi <- fracreg(d, dpo = 1, int = TRUE, np = 15, vcov = "inverse")$VDFA
  vm <- fracreg(d, dpo = 1, int = TRUE, np = 15, vcov = "marginal")$VDFA
  # the 1x1 inverse is exactly 1 / value, so the two estimators coincide
  expect_equal(vi, vm)
})

test_that("collinearity inflates the 'inverse' variance and the VIF", {
  set.seed(2)
  n  <- 400
  e1 <- rnorm(n); e2 <- rnorm(n)
  x1 <- cumsum(e1)
  x2 <- cumsum(0.9 * e1 + sqrt(1 - 0.81) * e2)   # correlated with x1
  y  <- x1 - 0.5 * x2 + cumsum(rnorm(n))
  fi <- fracreg(cbind(y, x1, x2), dpo = 1, int = TRUE, np = 15, vcov = "inverse")
  fm <- fracreg(cbind(y, x1, x2), dpo = 1, int = TRUE, np = 15, vcov = "marginal")

  expect_true(mean(fi$VIF, na.rm = TRUE) > 1)
  expect_true(stats::median(fi$VDFA / fm$VDFA, na.rm = TRUE) > 1)
  expect_true(all(fi$VIF   >= 1 - 1e-8, na.rm = TRUE))
  expect_true(all(fi$R2adj <= 1 + 1e-8, na.rm = TRUE))
})

test_that("fracreg.diag returns scale-wise diagnostics", {
  set.seed(3)
  n <- 400
  d <- data.frame(y = cumsum(rnorm(n)), x1 = cumsum(rnorm(n)),
                  x2 = cumsum(rnorm(n)))
  dg <- fracreg.diag(d, np = 15)
  expect_s3_class(dg, "data.frame")
  expect_true(all(c("s", "VIF_x1", "VIF_x2", "kappa", "R2adj") %in% names(dg)))
  expect_equal(nrow(dg), 15)
})

test_that("fracreg validates vcov and abs", {
  set.seed(4)
  d <- cbind(cumsum(rnorm(300)), cumsum(rnorm(300)))
  expect_error(fracreg(d, dpo = 1, int = TRUE, np = 15, vcov = "bogus"))
  expect_error(fracreg(d, dpo = 1, int = TRUE, np = 15, abs = "yes"), "abs")
})
