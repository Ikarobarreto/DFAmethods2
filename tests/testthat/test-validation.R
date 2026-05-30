test_that("input validation stops on common user mistakes", {
  set.seed(3)
  n  <- 600
  x2 <- cumsum(rnorm(n))
  y  <- x2 + cumsum(rnorm(n))
  d2 <- data.frame(y = y, x2 = x2)

  # series too short / missing values / bad np
  expect_error(dfa(rnorm(20)), "too short")
  expect_error(dfa(c(y[-1], NA)), "missing")
  expect_error(dfa(y, np = 1), "np")

  # not enough columns
  expect_error(rhodcca(d2[, 1, drop = FALSE]), "at least 2")
  expect_error(rhodpcca(d2), "at least 3")

  # non-numeric column and bad B
  bad <- data.frame(y = y, z = rep(letters, length.out = n))
  expect_error(rhodcca(bad), "numeric")
  expect_error(fracreg.PStest(d2, B = 0), "B")
})

test_that("valid inputs still return the expected objects", {
  set.seed(4)
  n  <- 600
  x2 <- cumsum(rnorm(n))
  x3 <- cumsum(rnorm(n))
  y  <- 0.5 * x2 - 0.3 * x3 + cumsum(rnorm(n))
  d3 <- data.frame(y = y, x2 = x2, x3 = x3)

  expect_s3_class(rhodpcca(d3, np = 30), "data.frame")
  expect_s3_class(dmc2(d3, np = 30), "data.frame")
  expect_s3_class(effsizeDFA(d3, np = 30), "data.frame")
})
