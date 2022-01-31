test_that("`multiMBD()` works as expected", {
  # Arrange
  withr::local_seed(1234)
  N <- 1e2
  P <- 1e3
  time_grid <- seq(0, 10, length.out = P)
  Cov <- exp_cov_function(time_grid, alpha = 0.2, beta = 0.8)
  Data_1 <- generate_gauss_fdata(N, centerline = rep(0, P), Cov = Cov)
  Data_2 <- generate_gauss_fdata(N, centerline = rep(0, P), Cov = Cov)

  # Act
  actual_wo_ties <- multiMBD(
    list(Data_1, Data_2),
    weights = 'uniform',
    manage_ties = FALSE
  ) - (MBD(Data_1, manage_ties = FALSE) +
         MBD(Data_2, manage_ties = FALSE)) / 2
  actual_with_ties <- multiMBD(
    list(Data_1, Data_2),
    weights = 'uniform',
    manage_ties = TRUE
  ) - (MBD(Data_1, manage_ties = TRUE) +
         MBD(Data_2, manage_ties = TRUE)) / 2
  actual_non_uniform_weights <- multiMBD(
    list(Data_1, Data_2),
    weights = c(1/3, 2/3),
    manage_ties = FALSE
  ) - (1/3 * MBD(Data_1, manage_ties = FALSE) +
         2/3 * MBD(Data_2, manage_ties = FALSE))

  # Assert
  expect_equal(actual_with_ties, rep(0, N))
  expect_equal(actual_wo_ties, rep(0, N))
  expect_equal(actual_non_uniform_weights, rep(0, N))
  expect_error(multiMBD(list(Data_1, Data_2), weights = c(1/2, 1)))
  expect_error(multiMBD(list(Data_1, Data_2), weights = 'unif'))
})

test_that("`multiBD()` works as expected", {
  # Arrange
  withr::local_seed(1234)
  N <- 1e2
  P <- 1e3
  time_grid <- seq(0, 10, length.out = P)
  Cov <- exp_cov_function(time_grid, alpha = 0.2, beta = 0.8)
  Data_1 <- generate_gauss_fdata(N, centerline = rep(0, P), Cov = Cov)
  Data_2 <- generate_gauss_fdata(N, centerline = rep(0, P), Cov = Cov)

  # Act
  actual_uniform_weights <- multiBD(
    list(Data_1, Data_2),
    weights = 'uniform'
  ) - (BD(Data_1) + BD(Data_2)) / 2
  actual_non_uniform_weights <- multiBD(
    list(Data_1, Data_2),
    weights = c(1/3, 2/3)
  ) - (1/3 * BD(Data_1) + 2/3 * BD(Data_2))

  # Assert
  expect_equal(actual_uniform_weights, rep(0, N))
  expect_equal(actual_non_uniform_weights, rep(0, N))
  expect_error(multiBD(list(Data_1, Data_2), weights = c(1/2, 1)))
  expect_error(multiBD(list(Data_1, Data_2), weights = 'unif'))
})
