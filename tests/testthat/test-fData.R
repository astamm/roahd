test_that("`fData() correctly creates `fData` objects", {
  withr::local_seed(1234)
  N <- 100
  P <- 1000
  time_grid <- seq(0, 1, length.out = P)
  Cov <- exp_cov_function(time_grid, alpha = 0.3, beta = 0.4)
  Data <- generate_gauss_fdata(N, centerline = sin(2 * pi * time_grid), Cov = Cov)
  fD <- fData(time_grid, Data)

  save_png <- function(code, width = 400, height = 400) {
    path <- tempfile(fileext = ".png")
    png(path, width = width, height = height)
    on.exit(dev.off())
    code

    path
  }

  expect_snapshot_value(fData(time_grid, Data), style = "serialize")
  expect_snapshot_value(fData(time_grid, 1:P), style = "serialize")

  expect_snapshot_plot("plot_fData", plot(
    fD,
    xlab = 'time', ylab = 'values',
    main = 'A functional dataset'
  ))
})

test_that("statistical summaries work as expected on `fData` objects", {
  withr::local_seed(1234)
  N <- 100
  P <- 1000
  time_grid <- seq(0, 1, length.out = P)
  Cov <- exp_cov_function(time_grid, alpha = 0.3, beta = 0.4)
  Data <- generate_gauss_fdata(N, centerline = sin(2 * pi * time_grid), Cov = Cov)
  fD <- fData(time_grid, Data)

  expect_equal(
    as.numeric(mean(fD)$values),
    colMeans(fD$values)
  )
  expect_equal(
    as.numeric(median_fData(fD)$value),
    fD$values[which.max(MBD(fD$values)), ]
  )
  expect_equal(
    as.numeric(median_fData(fD, type = 'MBD', manage_ties = TRUE)$values),
    fD$values[which.max(MBD(fD$values, manage_ties = TRUE)), ]
  )
  expect_equal(
    as.numeric(median_fData(fD, type = 'MHRD')$value),
    fD$values[which.max(MHRD(fD$values)), ]
  )
  expect_equal(
    cov(fD$values),
    cov_fun(fD)$values
  )
  expect_equal(
    cov(fD$values, (fD + 1:P)$values),
    cov_fun(fD, fD + 1:P)$values
  )

  expect_error(cov_fun(1))
  expect_error(cov_fun(fD, fD[-1, ]))
  expect_error(cov_fun(fD, fD[, 1:10]))
})

test_that("statistical summaries work as expected on `mfData` objects", {
  withr::local_seed(1234)
  N <- 100
  P <- 100
  time_grid <- seq(0, 1, length.out = P)
  Cov <- exp_cov_function(time_grid, alpha = 0.3, beta = 0.4)
  Data <- generate_gauss_fdata(
    N,
    centerline = sin(2 * pi * time_grid),
    Cov = Cov
  )
  fD <- fData(time_grid, Data)
  mfD <- mfData(time_grid, list('comp1' = Data, 'comp2' = Data, 'comp3' = Data))
  mfD2 <- mfData(time_grid, list(Data, Data, Data))
  mfD3 <- mfData(time_grid, list(Data, Data))
  rhs <- lapply(1:(2 * mfD$L), function( i ) cov_fun(mfD$fDList[[1]]))
  names(rhs) <- c('comp1_1', 'comp1_2', 'comp1_3', 'comp2_2', 'comp2_3', 'comp3_3')

  expect_equal(cov_fun(mfD, mfD2), rhs)
  expect_equal(
    cov_fun(mfD2, fD),
    lapply(1:mfD$L, function(i) cov_fun(mfD$fDList[[1]]))
  )
  expect_equal(
    names(cov_fun(mfD)),
    c(
      'comp1_comp1', 'comp1_comp2', 'comp1_comp3',
      'comp2_comp2', 'comp2_comp3', 'comp3_comp3'
    )
  )
  expect_equal(
    names(cov_fun(mfD2)),
    c('1_1', '1_2', '1_3', '2_2', '2_3', '3_3')
  )

  expect_error(cov_fun(mfD, mfD3))
  expect_error(cov_fun(mfD, fD[, 1:10]))
})

test_that("Algebraic operations work as expected on `fData` objects", {
  fD <- fData(
    grid = seq(0, 1, length.out = 10),
    values = matrix(seq(1, 10), nrow = 21, ncol = 10, byrow = TRUE)
  )
  fD.2 <- fData(
    grid = seq(0, 1, length.out = 11),
    values = matrix(seq(1, 11), nrow = 21, ncol = 11, byrow = TRUE)
  )

  expect_equal(
    sum((fD + 1:10)$values -
          matrix(2 * seq(1, 10),  nrow = 21, ncol = 10, byrow = TRUE)),
    0
  )
  expect_equal(sum((fD - 1:10)$values), 0)
  expect_equal(
    sum((fD + array(1, dim = c(1, 10)))$values -
          matrix(seq(2, 11), nrow = 21, ncol = 10, byrow = TRUE)),
    0
  )
  expect_equal(
    sum((fD + fD - matrix(
      2 * seq(1, 10),
      nrow = 21, ncol = 10,
      byrow = TRUE
    ))$values),
    0
  )
  expect_equal(fD * 2, fD + fD)
  expect_equal((fD * 4) / 2, fD + fD)

  expect_error(fD + fD.2, regexp = 'Error.*')

})

test_that("Subsetting operations work as expected on `fData` objects", {
  fD <- fData(
    grid = seq(0, 1, length.out = 10),
    values = matrix(seq(1, 10), nrow = 21, ncol = 10, byrow = TRUE)
  )

  expect_identical(fD[1, ], fData(seq(0, 1, length.out = 10), 1:10))
  expect_identical(
    fD[, 1:2],
    fData(
      grid = seq(0, 1, length.out = 10)[1:2],
      values = matrix(1:2, nrow = 21, ncol = 2, byrow = TRUE)
    )
  )
  expect_identical(
    fD[1:2, 1:2, as_fData = FALSE],
    matrix(seq(1, 2), nrow = 2, ncol = 2, byrow = TRUE)
  )

  fD_logical_subset <- fD[, c(rep(FALSE, 5), rep(TRUE, 5))]
  expect_identical(
    fD_logical_subset$values,
    matrix(seq(6, 10), nrow = 21, ncol = 5, byrow = TRUE)
  )
  expect_identical(
    c(
      fD_logical_subset$t0,
      fD_logical_subset$tP,
      fD_logical_subset$h,
      fD_logical_subset$P
    ),
    c(5 * fD$h, 1, fD$h, 5)
  )

  fD_logical_subset <- fD[, c(rep(TRUE, 5), rep(FALSE, 5))]

  expect_identical(
    fD_logical_subset$values,
    matrix(seq(1, 5), nrow = 21, ncol = 5, byrow = TRUE)
  )
  expect_identical(
    c(
      fD_logical_subset$t0,
      fD_logical_subset$tP,
      fD_logical_subset$h,
      fD_logical_subset$P
    ),
    c(0, 4 * fD$h, fD$h, 5)
  )
})

test_that("Subsetting operations work as expected on `mfData` objects", {
  N <- 3
  P <- 20
  L <- 2
  grid <- seq(0, 1, length.out = P)
  values1 <- rep(1, P)
  values2 <- cos(2 * pi * grid)
  values3 <- sin(2 * pi * grid)
  values4 <- rep(2, P)
  values5 <- cos(4 * pi * grid)
  values6 <- sin(4 * pi * grid)
  mfD <- mfData(grid, list(
    matrix(c(values1, values2, values3), nrow = 3, ncol = P, byrow = TRUE),
    matrix(c(values4, values5, values6), nrow = 3, ncol = P, byrow = TRUE)
  ))

  expect_equal(mfD[1, ]$fDList[[1]]$values, t(as.matrix(values1)))
  expect_equal(mfD[1, ]$fDList[[2]]$values, t(as.matrix(values4)))
  expect_equal(
    mfD[1:2, ]$fDList[[1]]$values,
    matrix(c(values1, values2), nrow = 2, ncol = P, byrow = TRUE)
  )
  expect_equal(
    mfD[1:2, ]$fDList[[2]]$values,
    matrix(c(values4, values5), nrow = 2, ncol = P, byrow = TRUE)
  )
  expect_equal(
    mfD[, 1:5]$fDList[[1]]$values,
    matrix(
      c(values1[1:5], values2[1:5], values3[1:5]),
      nrow = 3, ncol = 5, byrow = TRUE
    )
  )
  expect_equal(
    mfD[, 1:5]$fDList[[2]]$values,
    matrix(
      c(values4[1:5], values5[1:5], values6[1:5]),
      nrow = 3, ncol = 5, byrow = TRUE
    )
  )
  expect_equal(
    append_mfData(mfD[1:2, ], mfD[3, ])$fDList[[1]]$values,
    mfD$fDList[[1]]$values
  )
  expect_equal(
    append_fData(mfD[1:2, ]$fDList[[1]], mfD[3, ]$fDList[[1]])$values,
    mfD$fDList[[1]]$values
  )

  expect_true(
    (mfD[1:2, ]$N == 2) & (mfD[1:2, ]$P == P) & (mfD[1:2, ]$L == 2 )
  )

  expect_error(mfD[ , c(1:5, 7:8)])
})

test_that("`mfData() correctly creates `mfData` objects", {
  withr::local_seed(1234)
  N <- 100
  P <- 1000
  time_grid <- seq(0, 1, length.out = P)
  Cov <- exp_cov_function(time_grid, alpha = 0.3, beta = 0.4)
  Data_1 <- generate_gauss_fdata(N, centerline = sin(2 * pi * time_grid), Cov = Cov)
  Data_2 <- generate_gauss_fdata(N, centerline = sin(2 * pi * time_grid), Cov = Cov)

  save_png <- function(code, width = 400, height = 400) {
    path <- tempfile(fileext = ".png")
    png(path, width = width, height = height)
    on.exit(dev.off())
    code

    path
  }

  expect_snapshot_value(mfData(time_grid, list(Data_1, Data_2)), style = "serialize")
  expect_snapshot_value(mfData(time_grid, list(1:P, 1:P)), style = "serialize")

  expect_snapshot_plot("plot_mfData", plot(
    mfData(time_grid, list(Data_1, Data_2)),
    xlab = 'time', ylab = list('values', 'values'),
    main = list('First Component', 'Second Component')
  ))

  mfD <- mfData(time_grid, list(Data_1, Data_2))

  expect_snapshot_plot("plot_mfData_sub", plot(mfD$fDList[[1]]))

  expect_snapshot_value(as.mfData(list(
    fData(time_grid, Data_1),
    fData(time_grid, Data_2)
  )), style = "serialize")

  expect_identical(
    toListOfValues(mfData(time_grid, list(Data_1, Data_2))),
    list(Data_1, Data_2)
  )
})

test_that("`unfold()` works as expected", {
  P <- 1000
  time_grid <- seq(0, 1, length.out = P)
  D <- matrix(c(
    sin(2 * pi * time_grid),
    cos(2 * pi * time_grid),
    sin(10 * pi * time_grid) * time_grid + 2
    ), ncol = P, nrow = 3, byrow = TRUE)
  fD <- fData(time_grid, D)
  fD_unfold <- unfold(fD)

  # Ana's implementation
  mon_func <- function(x) {
    x<- as.matrix(x)
    if (ncol(x) == 1) x <- t(x)
    P <- ncol(x)
    x_mon <- x
    diff_x <- t(abs(apply(x_mon, 1, diff)))
    for (j in 2:P)
      x_mon[, j] <- x_mon[, j - 1] + diff_x[, j - 1]
    x_mon
  }

  expect_true(all(apply(fD_unfold$values, 1, function(x) all(diff(x) >= 0))))
  expect_equal(fD_unfold$values, mon_func(fD$values))
})

test_that("`warp()` works as expected", {
  withr::local_seed(1234)
  N <- 30
  P <- 1001
  t0 <- 0
  t1 <- 1
  time_grid <- seq(t0, t1, length.out = P)
  means <- round(runif(N, t0 + (t1 - t0) / 8, t1 - (t1 - t0) / 8), 3)
  Data <- matrix(
    sapply(means, function(m) dnorm(time_grid, mean = m, sd = 0.05)),
    ncol = P, nrow = N, byrow = TRUE
  )
  fD <- fData(time_grid, Data)

  # Piecewise linear warpings
  template_warping <- function(m) c(
    time_grid[time_grid <= 0.5] * m / 0.5,
    (time_grid[time_grid > 0.5] - 0.5) * (1 - m) / 0.5 + m
  )
  warpings <- matrix(
    sapply(means, template_warping),
    ncol = P, nrow = N, byrow = TRUE
  )
  wfD <- fData(time_grid, warpings)
  fD_warped <- warp(fD, wfD)

  expect_true(all(
    maxima(fD_warped) - dnorm(0.5, 0.5, 0.05) <= .Machine$double.eps
  ))
  expect_true(all(
    maxima(fD_warped, which = TRUE)$grid - 0.5 <= .Machine$double.eps
  ))
})
