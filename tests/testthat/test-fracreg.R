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

test_that("vcov='HC' runs and changes only the variance", {
  set.seed(5)
  n  <- 500
  x1 <- cumsum(rnorm(n)); x2 <- cumsum(rnorm(n))
  y  <- x1 - 0.5 * x2 + cumsum(rnorm(n))
  fi <- fracreg(cbind(y, x1, x2), dpo = 1, int = TRUE, np = 12,
                overlap = FALSE, vcov = "inverse")
  fh <- fracreg(cbind(y, x1, x2), dpo = 1, int = TRUE, np = 12,
                overlap = FALSE, vcov = "HC")
  expect_equal(fh$BDFA, fi$BDFA)                      # HC leaves the coefficients unchanged
  expect_true(all(is.finite(fh$VDFA)))
  expect_true(all(fh$VDFA >= 0))
  expect_false(isTRUE(all.equal(fh$VDFA, fi$VDFA)))   # but the variance differs
})

test_that("fracreg.WB runs and i.i.d. weights reproduce the HC variance", {
  set.seed(1)
  n <- 500; np <- 10
  x1 <- cumsum(rnorm(n)); x2 <- cumsum(rnorm(n))
  d <- cbind(y = x1 - 0.5 * x2 + cumsum(rnorm(n)), x1, x2)

  wb <- fracreg.WB(d, B = 199, weights = "rademacher", np = np)
  expect_s3_class(wb, "data.frame")
  expect_equal(nrow(wb), np)
  expect_true(all(c("s", "beta_x1", "lower_x1", "upper_x1", "p_x1") %in% names(wb)))

  # the i.i.d. wild-bootstrap variance equals the HC sandwich variance (no HC1)
  pb <- .fracreg_perbox(d, dpo = 1, int = TRUE, np = np, overlap = FALSE)
  ps <- pb$perscale[[4]]; Ts <- ps$Ts
  set.seed(2)
  W  <- matrix(sample(c(-1, 1), Ts * 8000, replace = TRUE), Ts, 8000)
  BS <- ps$beta + ps$Fxx_inv %*% (crossprod(ps$scores, W) / Ts)
  var_wb <- apply(BS, 1, var)
  Om <- crossprod(ps$scores) / Ts
  var_hc <- diag((ps$Fxx_inv %*% Om %*% ps$Fxx_inv) / Ts)
  expect_equal(unname(var_wb), unname(var_hc), tolerance = 0.05)

  expect_error(fracreg.WB(d, weights = "bogus"))
})

test_that("fracreg validates vcov and abs", {
  set.seed(4)
  d <- cbind(cumsum(rnorm(300)), cumsum(rnorm(300)))
  expect_error(fracreg(d, dpo = 1, int = TRUE, np = 15, vcov = "bogus"))
  expect_error(fracreg(d, dpo = 1, int = TRUE, np = 15, abs = "yes"), "abs")
})
