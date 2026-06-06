# The per-box C primitives (rdfa_box / rdcca_box) expose the detrended
# (co)variance S(s, v) of every box; the mean over the T_s boxes of a scale must
# reproduce the corresponding F^2(s). They are the shared core for the HC,
# wild-bootstrap and surrogate machinery.

test_that("rdfa_box / rdcca_box per-box values average to F^2(s)", {
  set.seed(1)
  n <- 400; np <- 15; mx <- round(n / 5)
  x <- cumsum(rnorm(n)); y <- cumsum(rnorm(n))

  cfg <- as.integer(cbind(n, 2, 1, np, 1, 10, mx))   # NPTS,NFIT,IFLAG,NR,SW,MINBOX,MAXBOX
  a <- .C("rdfa_box", cfg, as.numeric(x), as.integer(numeric(np + 1)),
          numeric(np + 1), numeric(np * n), PACKAGE = "DFAmethods2")
  s  <- a[[3]][2:(np + 1)]
  f2 <- a[[4]][2:(np + 1)]
  box <- matrix(a[[5]], nrow = n)                     # column i = scale i, rows = boxes
  for (i in seq_len(np)) {
    Ts <- n - s[i] + 1                                # overlapping boxes
    expect_equal(mean(box[1:Ts, i]), f2[i])
  }

  cfgc <- as.integer(cbind(n, 2, 1, np, 1, 10, mx, 0))   # + ABSFLAG
  ac <- .C("rdcca_box", cfgc, as.numeric(x), as.numeric(y),
           as.integer(numeric(np + 1)), numeric(np + 1), numeric(np * n),
           PACKAGE = "DFAmethods2")
  sc  <- ac[[4]][2:(np + 1)]
  f2c <- ac[[5]][2:(np + 1)]
  boxc <- matrix(ac[[6]], nrow = n)
  for (i in seq_len(np)) {
    Ts <- n - sc[i] + 1
    expect_equal(mean(boxc[1:Ts, i]), f2c[i])
  }
})
