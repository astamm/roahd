# Functional boxplot for univariate functional data -----------------------

test_that("fbplot() works as expected for univariate data", {
  time_grid <- seq(0, 1, length.out = 1e2)
  D <- matrix(c(
    sin(2 * pi * time_grid) + 0,
    sin(2 * pi * time_grid) + 1,
    sin(2 * pi * time_grid) + 2,
    sin(2 * pi * time_grid) + 3,
    sin(2 * pi * time_grid) + 4,
    sin(2 * pi * time_grid) + 5,
    sin(2 * pi * time_grid) + 6,
    sin(2 * pi * time_grid) + 7,
    sin(2 * pi * time_grid) + 8,
    sin(2 * pi * time_grid) + 9,
    sin(2 * pi * time_grid) + 10,
    sin(2 * pi * time_grid) - 1,
    sin(2 * pi * time_grid) - 2,
    sin(2 * pi * time_grid) - 3,
    sin(2 * pi * time_grid) - 4,
    sin(2 * pi * time_grid) - 5,
    sin(2 * pi * time_grid) - 6,
    sin(2 * pi * time_grid) - 7,
    sin(2 * pi * time_grid) - 8,
    sin(2 * pi * time_grid) - 9,
    sin(2 * pi * time_grid) - 10
  ), nrow = 21, ncol = length(time_grid), byrow = TRUE)
  fD <- fData(time_grid, D)

  expect_snapshot_value(fbplot(fD, Fvalue = 10, display = FALSE))
  expect_snapshot_value(fbplot(
    fD,
    display = FALSE,
    xlab = 'time',
    ylab = 'value',
    main = 'My Functional Boxplot'
  ))
})

# Adjusted functional boxplot for univariate data -------------------------

test_that("`fbplot()` with adjustment works as expected for univariate data", {
  withr::local_seed(1234)
  time_grid <- seq(0, 1, length.out = 1e2)
  N <- 5e2
  Data <- generate_gauss_fdata(
    N,
    centerline = sin(2 * pi * time_grid),
    Cov = exp_cov_function(time_grid, alpha = 0.3, beta = 0.4)
  )
  fD <- fData(time_grid, Data)

  expect_warning(expect_warning(fbplot(
    fD,
    adjust = list(N_trials = 1, trial_size = N, foo = 'bar', baz = 'qux'),
    display = FALSE
  )))

  skip_on_cran()
  expect_snapshot_value(fbplot(
    fD,
    adjust = list(N_trials = 10, trial_size = N, VERBOSE = FALSE),
    display = FALSE,
    xlab = 'time', ylab = 'Values',
    main = 'My adjusted functional boxplot'
  ))

})

# For multivariate functional data ----------------------------------------

test_that("`fbplot()` works as expected for multivariate data", {
  time_grid <- seq(0, 1, length.out = 1e2)
  D <- matrix(c(
    sin(2 * pi * time_grid) - 10,
    sin(2 * pi * time_grid) - 9,
    sin(2 * pi * time_grid) - 8,
    sin(2 * pi * time_grid) - 7,
    sin(2 * pi * time_grid) - 6,
    sin(2 * pi * time_grid) - 5,
    sin(2 * pi * time_grid) - 4,
    sin(2 * pi * time_grid) - 3,
    sin(2 * pi * time_grid) - 2,
    sin(2 * pi * time_grid) - 1,
    sin(2 * pi * time_grid) + 0,
    sin(2 * pi * time_grid) + 1,
    sin(2 * pi * time_grid) + 2,
    sin(2 * pi * time_grid) + 3,
    sin(2 * pi * time_grid) + 4,
    sin(2 * pi * time_grid) + 5,
    sin(2 * pi * time_grid) + 6,
    sin(2 * pi * time_grid) + 7,
    sin(2 * pi * time_grid) + 8,
    sin(2 * pi * time_grid) + 9,
    sin(2 * pi * time_grid) + 10
  ), nrow = 21, ncol = length(time_grid), byrow = TRUE)
  mfD <- mfData(time_grid, list(D, D * abs(1:21 - 11) / 5))

  expect_snapshot_value(fbplot(mfD, Fvalue = 3, display = FALSE))
  expect_error(fbplot(mfD, adjust = list(N_trials = 2), display = FALSE))
})

test_that("`fbplot()` works for randomly generated multivariate data", {
  withr::local_seed(1234)
  P <- 1e2
  N <- 1e2
  L <- 3
  time_grid <- seq(0, 1, length.out = P)
  C1 <- exp_cov_function(time_grid, alpha = 0.3, beta = 0.4)
  C2 <- exp_cov_function(time_grid, alpha = 0.3, beta = 0.4)
  C3 <- exp_cov_function(time_grid, alpha = 0.3, beta = 0.4)
  Data <- generate_gauss_mfdata(
    N, L,
    centerline = matrix(sin(2 * pi * time_grid), nrow = 3, ncol = P, byrow = TRUE ),
    correlations = rep(0.5, 3),
    listCov = list(C1, C2, C3)
  )
  mfD <- mfData(time_grid, Data)

  expect_snapshot_value(fbplot(mfD, Fvalue = 2.5, display = FALSE))
})
