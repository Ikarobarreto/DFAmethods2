# Well-posed designs throughout: stationary predictors with int = TRUE (DFA
# integrates the stationary input into the profile). Non-overlapping boxes are
# the sampling units of the analytic inference.

test_that("fracreg variance: 'inv' equals 'marginal' for a single predictor", {
  set.seed(1)
  n <- 800
  x <- rnorm(n)
  y <- 0.6 * x + rnorm(n)
  d <- cbind(y = y, x = x)
  vi <- fracreg(d, dpo = 1, int = TRUE, np = 15, overlap = FALSE, variance = "inv", min_boxes = 5)$VDFA
  vm <- fracreg(d, dpo = 1, int = TRUE, np = 15, overlap = FALSE, variance = "marginal", min_boxes = 5)$VDFA
  # the 1x1 inverse is exactly 1 / value, so the two estimators coincide
  expect_equal(vi, vm)
})

test_that("collinearity inflates the 'inv' variance and the VIF", {
  set.seed(2)
  n  <- 800
  e1 <- rnorm(n); e2 <- rnorm(n)
  x1 <- e1
  x2 <- 0.9 * e1 + sqrt(1 - 0.81) * e2          # stationary, correlated with x1
  y  <- x1 - 0.5 * x2 + rnorm(n)
  fi <- fracreg(cbind(y, x1, x2), dpo = 1, int = TRUE, np = 15,
                overlap = FALSE, variance = "inv", min_boxes = 5)
  fm <- fracreg(cbind(y, x1, x2), dpo = 1, int = TRUE, np = 15,
                overlap = FALSE, variance = "marginal", min_boxes = 5)

  expect_true(mean(fi$VIF, na.rm = TRUE) > 1)
  expect_true(stats::median(fi$VDFA / fm$VDFA, na.rm = TRUE) > 1)
  expect_true(all(fi$VIF   >= 1 - 1e-8, na.rm = TRUE))
  expect_true(all(fi$R2adj <= 1 + 1e-8, na.rm = TRUE))
})

test_that("inv_corrected scales the inverse variance by the memory factor", {
  set.seed(20)
  n  <- 1000
  x1 <- rnorm(n); x2 <- rnorm(n)
  y  <- x1 - 0.5 * x2 + rnorm(n)
  fi <- fracreg(cbind(y, x1, x2), dpo = 1, int = TRUE, np = 20,
                overlap = FALSE, variance = "inv", min_boxes = 5)
  fc <- fracreg(cbind(y, x1, x2), dpo = 1, int = TRUE, np = 20,
                overlap = FALSE, variance = "inv_corrected", min_boxes = 5)
  expect_equal(fc$VDFA, fi$VDFA * fc$kappa_factor[1])     # M4.2': inv x kappa(H)
  expect_true(fc$kappa_factor[1] > 0 && fc$kappa_factor[1] <= 1)
  expect_true(is.finite(fc$H_resid))
})

test_that("fracreg.diag returns scale-wise diagnostics", {
  set.seed(3)
  n <- 800
  d <- data.frame(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
  dg <- fracreg.diag(d, np = 15, overlap = FALSE)
  expect_s3_class(dg, "data.frame")
  expect_true(all(c("s", "VIF_x1", "VIF_x2", "kappa", "R2adj") %in% names(dg)))
  expect_equal(nrow(dg), 15)
})

test_that("variance = 'hc' runs and changes only the variance", {
  set.seed(5)
  n  <- 800
  x1 <- rnorm(n); x2 <- rnorm(n)
  y  <- x1 - 0.5 * x2 + rnorm(n)
  fi <- fracreg(cbind(y, x1, x2), dpo = 1, int = TRUE, np = 12,
                overlap = FALSE, variance = "inv", min_boxes = 5)
  fh <- fracreg(cbind(y, x1, x2), dpo = 1, int = TRUE, np = 12,
                overlap = FALSE, variance = "hc", min_boxes = 5)
  expect_equal(fh$BDFA, fi$BDFA)                      # HC leaves the coefficients unchanged
  expect_true(all(is.finite(fh$VDFA)))
  expect_true(all(fh$VDFA >= 0))
  expect_false(isTRUE(all.equal(fh$VDFA, fi$VDFA)))   # but the variance differs
})

test_that("overlap = TRUE returns point estimates only (variance = 'none')", {
  set.seed(21)
  n <- 600
  d <- cbind(y = rnorm(n), x = rnorm(n))
  expect_message(fo <- fracreg(d, dpo = 1, int = TRUE, np = 15, overlap = TRUE),
                 "overlap")
  expect_identical(fo$variance_method, "none")
  expect_true(all(is.na(fo$UCIB)))
})

test_that("the analytic variance uses the box count T_s = floor(N/s)", {
  # With the corrected denominator df = T_s - k, the standard error at a fixed
  # scale must SHRINK as N grows (more boxes), not stay tied to the box size.
  s_fixed <- 25
  se_at <- function(N) {
    set.seed(7)
    x <- rnorm(N); y <- 0.5 * x + rnorm(N)
    fr <- fracreg(cbind(y, x), dpo = 1, int = TRUE, np = 40, overlap = FALSE, min_boxes = 5)
    k  <- which.min(abs(fr$s - s_fixed))
    sqrt(fr$VDFA[1, 1, k])
  }
  expect_true(se_at(4000) < se_at(1000))
})

test_that("too few boxes (T_s <= k) yield NA limits and a warning", {
  set.seed(12)
  N <- 1000
  m <- as.data.frame(matrix(rnorm(7 * N), N, 7))
  names(m) <- c("y", paste0("x", 1:6))                # 6 predictors
  expect_warning(fracreg(m, dpo = 1, int = TRUE, np = 20, overlap = FALSE),
                 "T_s|floor")
  f <- suppressWarnings(fracreg(m, dpo = 1, int = TRUE, np = 20, overlap = FALSE))
  expect_true(any(is.na(f$VDFA)))                     # the largest scales are NA
})

test_that("a strong-memory error (H > 3/4) triggers the warning and reports H_resid", {
  set.seed(10)
  n <- 800
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  # supply H_eps directly so the regime is deterministic
  expect_warning(fracreg(d, dpo = 1, int = TRUE, np = 20, overlap = FALSE,
                         H_eps = 0.85, min_boxes = 5),
                 "Hermite|3/4|H_resid")
  fr <- suppressWarnings(
    fracreg(d, dpo = 1, int = TRUE, np = 20, overlap = FALSE,
            H_eps = 0.85, min_boxes = 5))
  expect_equal(fr$H_resid, 0.85)
  expect_equal(unname(fr$kappa_factor[1]), (2 * 0.85 + 1)^2 / 21)
  expect_true(all(c("alpha", "H_resid", "kappa_factor") %in% names(fr)))
  expect_length(fr$alpha, 2)
})

test_that("well-posed stationary inputs do not trip the strong-memory warning", {
  set.seed(11)
  n <- 1500
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  # narrow the scale range so all scales clear the min_boxes floor too
  expect_warning(fracreg(d, dpo = 1, int = TRUE, np = 12, overlap = FALSE,
                         min_boxes = 5),
                 regexp = NA)
})

test_that("min_boxes suppresses the standard errors of underpopulated scales", {
  set.seed(31)
  n <- 1000
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  expect_warning(fracreg(d, dpo = 1, int = TRUE, np = 30, min_boxes = 15),
                 "min_boxes = 15")
  fr <- suppressWarnings(
    fracreg(d, dpo = 1, int = TRUE, np = 30, min_boxes = 15))
  big <- fr$s > floor(n / 15)                          # T_s < 15 here
  expect_true(any(big))                                # there are suppressed scales
  expect_true(all(is.na(fr$UCIB[1, 1, big])))          # suppressed -> NA CI
  expect_true(all(is.finite(fr$BDFA[1, 1, big])))      # but the point estimate stays
})

test_that("N < 500 triggers a separate Likens-style warning", {
  set.seed(32)
  d <- data.frame(y = rnorm(300), x = rnorm(300))
  expect_warning(fracreg(d, dpo = 1, int = TRUE, np = 10, min_boxes = 5),
                 "Likens|< 500")
})

test_that("fracreg.WB drops the overlap argument by construction", {
  set.seed(33)
  d <- cbind(y = rnorm(600), x = rnorm(600))
  expect_error(suppressWarnings(fracreg.WB(d, B = 49, np = 12, overlap = FALSE)),
               "unused argument|overlap")
})

test_that("fracreg.WB min_boxes blanks the CI of small-T_s scales", {
  set.seed(34)
  d <- cbind(y = rnorm(800), x = rnorm(800))
  expect_warning(fracreg.WB(d, B = 49, np = 20, min_boxes = 15),
                 "min_boxes = 15")
  wb <- suppressWarnings(fracreg.WB(d, B = 49, np = 20, min_boxes = 15))
  big <- wb$s > floor(800 / 15)
  expect_true(any(big))
  expect_true(all(is.na(wb$lower_x[big])))             # CI suppressed at small T_s
  expect_true(all(is.finite(wb$beta_x[big])))          # point estimate still computed
})

test_that("fracreg.WB runs and i.i.d. weights reproduce the HC variance", {
  set.seed(1)
  n <- 500; np <- 10
  x1 <- cumsum(rnorm(n)); x2 <- cumsum(rnorm(n))
  d <- cbind(y = x1 - 0.5 * x2 + cumsum(rnorm(n)), x1, x2)

  wb <- suppressWarnings(fracreg.WB(d, B = 199, weights = "rademacher", np = np, min_boxes = 3))
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

test_that("fracreg validates variance and abs", {
  set.seed(4)
  d <- cbind(rnorm(300), rnorm(300))
  expect_error(fracreg(d, dpo = 1, int = TRUE, np = 15, variance = "bogus"))
  expect_error(fracreg(d, dpo = 1, int = TRUE, np = 15, abs = "yes"), "abs")
})
