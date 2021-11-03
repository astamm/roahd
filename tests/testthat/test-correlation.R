# maxima() & minima() -----------------------------------------------------

test_that("maxima() works for functional data when `which = TRUE`", {
  # Arrange
  P <- 1e4
  time_grid <- seq(0, 1, length.out = P)
  h <- time_grid[2] - time_grid[1]
  Data <- matrix(c(1 * time_grid,
                   2 *  time_grid,
                   3 * ( 0.5 - abs( time_grid - 0.5))),
                 nrow = 3, ncol = P, byrow = TRUE)
  fD <- fData(time_grid, Data)

  # Act
  actual <- maxima(fD, which = TRUE)

  # Assert
  expected_grid <- c(1, 1, 0 + 4999 * h)
  expected_value <- c(1, 2, 3 * ( 0.5 - abs( 0.5 - 4999 * h)))
  expect_equal(actual$grid, expected_grid)
  expect_equal(actual$value, expected_value)
})

test_that("minima() works for functional data when `which = TRUE`", {
  # Arrange
  P <- 1e4
  time_grid <- seq(0, 1, length.out = P)
  Data <- matrix(c(1 * time_grid,
                   2 *  time_grid,
                   3 * ( 0.5 - abs( time_grid - 0.5))),
                 nrow = 3, ncol = P, byrow = TRUE)
  fD <- fData(time_grid, Data)

  # Act
  actual <- minima(fD, which = TRUE)

  # Assert
  expected_grid <- rep(0, 3)
  expected_value <- rep(0, 3)
  expect_equal(actual$grid, expected_grid)
  expect_equal(actual$value, expected_value)
})

test_that("maxima() works for functional data when `which = FALSE`", {
  # Arrange
  P <- 1e4
  time_grid <- seq(0, 1, length.out = P)
  h <- time_grid[2] - time_grid[1]
  Data <- matrix(c(1 * time_grid,
                   2 *  time_grid,
                   3 * ( 0.5 - abs( time_grid - 0.5))),
                 nrow = 3, ncol = P, byrow = TRUE)
  fD <- fData(time_grid, Data)

  # Act
  actual <- maxima(fD, which = FALSE)

  # Assert
  expected <- c(1, 2, 3 * (0.5 - abs( 0.5 - 4999 * h)))
  expect_equal(actual, expected)
})

test_that("minima() works for functional data when `which = FALSE`", {
  # Arrange
  P <- 1e4
  time_grid <- seq(0, 1, length.out = P)
  Data <- matrix(c(1 * time_grid,
                   2 *  time_grid,
                   3 * ( 0.5 - abs( time_grid - 0.5))),
                 nrow = 3, ncol = P, byrow = TRUE)
  fD <- fData(time_grid, Data)

  # Act
  actual <- minima(fD, which = FALSE)

  # Assert
  expected <- rep(0, 3)
  expect_equal(actual, expected)
})

# area_under_curve() ------------------------------------------------------

test_that("area_under_curve() works for functional data", {
  # Arrange
  P <- 1e4
  time_grid <- seq(0, 1, length.out = P)
  fD_1 <- fData(time_grid,
                matrix(c(1 * time_grid,
                         2 *  time_grid,
                         3 * ( 0.5 - abs( time_grid - 0.5))),
                       nrow = 3, ncol = P, byrow = TRUE))
  fD_2 <- fData(time_grid,
                matrix(c(sin(2 * pi * time_grid),
                         cos(2 * pi * time_grid),
                         4 * time_grid * (1 - time_grid)),
                       nrow = 3, ncol = P, byrow = TRUE))

  # Act
  actual <- area_under_curve(fD_1)

  # Assert
  expected <- c(0.5, 1, 0.75)
  expect_equal(actual, expected)
  expect_true(all(c(
    area_under_curve(fD_2)[1:2],
    abs(area_under_curve(fD_2[3, ]) - 2/3)
  ) <= .Machine$double.eps^0.5))
})

# Ordering functions ------------------------------------------------------

test_that("max_ordered() works as expected", {
  # Arrange
  P <- 1e3
  time_grid <- seq(0, 1, length.out = P)
  h <- time_grid[2] - time_grid[1]

  Data_1 <- matrix(
    c(1 * time_grid, 2 *  time_grid),
    nrow = 2, ncol = P, byrow = TRUE
  )
  Data_2 <- matrix(
    3 * (0.5 - abs(time_grid - 0.5)),
    nrow = 1, byrow = TRUE
  )
  Data_3 <- rbind(Data_1, Data_1)

  fD_1 <- fData(time_grid, Data_1)
  fD_2 <- fData(time_grid, Data_2)
  fD_3 <- fData(time_grid, Data_3)

  # Act
  actual_max_1 <- max_ordered(fD_1, fD_2)
  actual_max_2 <- max_ordered(fD_2, fD_1)
  actual_max_3 <- max_ordered(fD_2, fD_3)
  actual_max_4 <- max_ordered(fD_3, fD_2)

  # Assert
  expected_max_1 <- c(TRUE, FALSE)
  expect_equal(actual_max_1, expected_max_1)

  expected_max_2 <- c(FALSE, TRUE)
  expect_equal(actual_max_2, expected_max_2)

  expect_error(max_ordered(fD_1, fD_3))
  expect_error(max_ordered(fD_3, fD_1))

  expected_max_3 <- c(FALSE, TRUE, FALSE, TRUE)
  expect_equal(actual_max_3, expected_max_3)

  expected_max_4 <- c(TRUE, FALSE, TRUE, FALSE)
  expect_equal(actual_max_4, expected_max_4)
})

test_that("area_ordered() works as expected", {
  # Arrange
  P <- 1e3
  time_grid <- seq(0, 1, length.out = P)
  h <- time_grid[2] - time_grid[1]

  Data_1 <- matrix(
    c(1 * time_grid, 2 *  time_grid),
    nrow = 2, ncol = P, byrow = TRUE
  )
  Data_2 <- matrix(
    3 * (0.5 - abs(time_grid - 0.5)),
    nrow = 1, byrow = TRUE
  )
  Data_3 <- rbind(Data_1, Data_1)

  fD_1 <- fData(time_grid, Data_1)
  fD_2 <- fData(time_grid, Data_2)
  fD_3 <- fData(time_grid, Data_3)

  # Act
  actual_area_1 <- area_ordered(fD_1, fD_2)
  actual_area_2 <- area_ordered(fD_2, fD_1)
  actual_area_3 <- area_ordered(fD_2, fD_3)
  actual_area_4 <- area_ordered(fD_3, fD_2)

  # Assert
  expected_area_1 <- c(TRUE, FALSE)
  expect_equal(actual_area_1, expected_area_1)

  expected_area_2 <- c(FALSE, TRUE)
  expect_equal(actual_area_2, expected_area_2)

  expect_error(area_ordered(fD_1, fD_3))
  expect_error(area_ordered(fD_3, fD_1))

  expected_area_3 <- c(FALSE, TRUE, FALSE, TRUE)
  expect_equal(actual_area_3, expected_area_3)

  expected_area_4 <- c(TRUE, FALSE, TRUE, FALSE)
  expect_equal(actual_area_4, expected_area_4)
})

# cor_kendall() & cor_spearman() ------------------------------------------

test_that("cor_kendall() and cor_spearman() work as expected", {
  # Arrange
  withr::local_seed(1234)
  N <- 2e2
  P <- 1e3
  time_grid <- seq(0, 1, length.out = P)
  Cov <- exp_cov_function(time_grid, alpha = 0.3, beta = 0.4)
  Data_1 <- generate_gauss_fdata(
    N,
    centerline = sin(2 * pi * time_grid),
    Cov = Cov
  )
  Data_2 <- generate_gauss_fdata(
    N,
    centerline = sin(2 * pi * time_grid),
    Cov = Cov
  )
  mfD <- mfData(time_grid, list(Data_1, Data_2))

  # Act
  actual_kendall_max  <- cor_kendall(mfD, ordering = 'max')
  actual_kendall_area <- cor_kendall(mfD, ordering = 'area')
  actual_spearman_mei <- cor_spearman(mfD, ordering = 'MEI')
  actual_spearman_mhi <- cor_spearman(mfD, ordering = 'MHI')

  # Assert
  expect_snapshot_value(actual_kendall_max,  style = "serialize")
  expect_snapshot_value(actual_kendall_area, style = "serialize")
  expect_snapshot_value(actual_spearman_mei, style = "serialize")
  expect_snapshot_value(actual_spearman_mhi, style = "serialize")
})

# Case studies from Dalia Valencia, Rosa Lillo, Juan Romo -----------------

test_that("cor_kendall() & cor_spearman() work on Case Study 1 from Dalia Valencia, Rosa Lillo, Juan Romo.", {
  withr::local_seed(1234)
  N <- 50
  P <- 50
  time_grid<- seq(0, 1, length.out = P)

  sigma_12 <- 0.8
  R <- matrix(c(1, sigma_12, sigma_12, 1), ncol = 2, nrow = 2)
  Z <- matrix(rnorm(N * 2, 0, 1), ncol = 2, nrow = N) %*% chol(R)
  X <- (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 1])^3 +
    (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 1])^2 +
    (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 1]) * 3
  Y <- (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 2])^2 +
    (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 2]) * 7 / 8 +
    - 10

  mfD <- mfData(time_grid, list(X, Y))

  expect_snapshot_value(cor_kendall(mfD, ordering = 'max'), style = "serialize")
  expect_snapshot_value(cor_kendall(mfD, ordering = 'area'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MEI'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MHI'), style = "serialize")
})

test_that("cor_kendall() & cor_spearman() work on Case Study 2 from Dalia Valencia, Rosa Lillo, Juan Romo.", {
  withr::local_seed(1234)
  N <- 50
  P <- 50
  time_grid<- seq(0, 1, length.out = P)

  sigma_12 <- - 0.7
  R <- matrix(c(1, sigma_12, sigma_12, 1), ncol = 2, nrow = 2)
  Z <- matrix(rnorm(N * 2, 0, 1), ncol = 2, nrow = N) %*% chol(R)
  X <- sin(matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 1])
  Y <- cos(matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 2])

  mfD <- mfData(time_grid, list(X, Y))

  expect_snapshot_value(cor_kendall(mfD, ordering = 'max'), style = "serialize")
  expect_snapshot_value(cor_kendall(mfD, ordering = 'area'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MEI'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MHI'), style = "serialize")
})

test_that("cor_kendall() & cor_spearman() work on Case Study 3 from Dalia Valencia, Rosa Lillo, Juan Romo.", {
  withr::local_seed(1234)
  N <- 50
  P <- 50
  time_grid<- seq(0, 1, length.out = P)

  sigma_12 <- 1
  Z <- rnorm(N, 0, 1)
  X <- (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z)^2
  Y <- (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z)^4

  mfD <- mfData(time_grid, list(X, Y))

  expect_equal(cor_kendall(mfD, ordering = 'max' ), 1)
  expect_equal(cor_kendall(mfD, ordering = 'area'), 1)

  expect_equal(cor_spearman(mfD, ordering = 'MEI'), 1)
  expect_equal(cor_spearman(mfD, ordering = 'MHI'), 1)
})

test_that("cor_kendall() & cor_spearman() work on Case Study 4 from Dalia Valencia, Rosa Lillo, Juan Romo.", {
  withr::local_seed(1234)
  N <- 50
  P <- 50
  time_grid<- seq(0, 1, length.out = P)

  sigma_12 <- 1
  Z <- rnorm(N, 0, 1)
  X <- (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z)^2 +
    (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z) * 7 +
    2
  Y <- ((matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z)^2 +
          (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z) * 7 +
          2)^3

  mfD <- mfData(time_grid, list(X, Y))

  expect_equal(cor_kendall(mfD, ordering = 'max' ), 1)
  expect_equal(cor_kendall(mfD, ordering = 'area'), 1)

  expect_equal(cor_spearman(mfD, ordering = 'MEI'), 1)
  expect_equal(cor_spearman(mfD, ordering = 'MHI'), 1)
})

test_that("cor_kendall() & cor_spearman() work on Case Study 5 from Dalia Valencia, Rosa Lillo, Juan Romo.", {
  withr::local_seed(1234)
  N <- 50
  P <- 50
  time_grid<- seq(0, 1, length.out = P)

  sigma_12 <- 1
  Z <- rnorm(N, 0, 1)
  X <- (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z)^2 +
    (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z) * 7 +
    2
  Y <- 1 - ((matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z)^2 +
              (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z) * 7 +
              2)^3

  mfD <- mfData(time_grid, list(X, Y))

  expect_equal(cor_kendall(mfD, ordering = 'max' ), -1)
  expect_equal(cor_kendall(mfD, ordering = 'area'), -1)

  expect_equal(cor_spearman(mfD, ordering = 'MEI'), -1)
  expect_equal(cor_spearman(mfD, ordering = 'MHI'), -1)
})

test_that("cor_kendall() & cor_spearman() work on Case Study 6 from Dalia Valencia, Rosa Lillo, Juan Romo.", {
  withr::local_seed(1234)
  N <- 50
  P <- 50
  time_grid<- seq(0, 1, length.out = P)

  sigma_12 <- 0.6
  R <- matrix(c(1, sigma_12, sigma_12, 1), ncol = 2, nrow = 2)
  Z <- matrix(rnorm(N * 2, 0, 1), ncol = 2, nrow = N) %*% chol(R)
  X <- exp(matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 1])
  Y <- (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 2])^3 +
    (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 2])^2 +
    (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 2]) * 3

  mfD <- mfData(time_grid, list(X, Y))

  expect_snapshot_value(cor_kendall(mfD, ordering = 'max'), style = "serialize")
  expect_snapshot_value(cor_kendall(mfD, ordering = 'area'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MEI'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MHI'), style = "serialize")
})

test_that("cor_kendall() & cor_spearman() work on Case Study 7 from Dalia Valencia, Rosa Lillo, Juan Romo.", {
  withr::local_seed(1234)
  N <- 50
  P <- 50
  time_grid<- seq(0, 1, length.out = P)

  sigma_12 <- -0.8
  R <- matrix(c(1, sigma_12, sigma_12, 1), ncol = 2, nrow = 2)
  Z <- matrix(rnorm(N * 2, 0, 1), ncol = 2, nrow = N) %*% chol(R)
  X <- exp(matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 1])^2
  Y <- cos((matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 2]))

  mfD <- mfData(time_grid, list(X, Y))

  expect_snapshot_value(cor_kendall(mfD, ordering = 'max'), style = "serialize")
  expect_snapshot_value(cor_kendall(mfD, ordering = 'area'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MEI'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MHI'), style = "serialize")
})

test_that("cor_kendall() & cor_spearman() work on Case Study 8 from Dalia Valencia, Rosa Lillo, Juan Romo.", {
  withr::local_seed(1234)
  N <- 50
  P <- 50
  time_grid<- seq(0, 1, length.out = P)

  sigma_12 <- 0.4
  R <- matrix(c(1, sigma_12, sigma_12, 1), ncol = 2, nrow = 2)
  Z <- matrix(rnorm(N * 2, 0, 1), ncol = 2, nrow = N) %*% chol(R)
  X <- sin(matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 1])
  Y <- (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 2])^2

  mfD <- mfData(time_grid, list(X, Y))

  expect_snapshot_value(cor_kendall(mfD, ordering = 'max'), style = "serialize")
  expect_snapshot_value(cor_kendall(mfD, ordering = 'area'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MEI'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MHI'), style = "serialize")
})

test_that("cor_kendall() & cor_spearman() work on Case Study 9 from Dalia Valencia, Rosa Lillo, Juan Romo.", {
  withr::local_seed(1234)
  N <- 50
  P <- 50
  time_grid<- seq(0, 1, length.out = P)

  sigma_12 <- 1
  Z <- rnorm(N, 0, 1)
  X <- (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z)^2 +
    (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z) * 9 +
    - 5
  Y <- cos(matrix(3 * time_grid, nrow = N, ncol = P, byrow = TRUE) + Z)

  mfD <- mfData(time_grid, list(X, Y))

  expect_snapshot_value(cor_kendall(mfD, ordering = 'max'), style = "serialize")
  expect_snapshot_value(cor_kendall(mfD, ordering = 'area'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MEI'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MHI'), style = "serialize")
})

test_that("cor_kendall() & cor_spearman() work on Case Study 10 from Dalia Valencia, Rosa Lillo, Juan Romo.", {
  withr::local_seed(1234)
  N <- 50
  P <- 50
  time_grid<- seq(0, 1, length.out = P)

  sigma_12 <- 0.9
  R <- matrix(c(1, sigma_12, sigma_12, 1), ncol = 2, nrow = 2)
  Z <- matrix(rnorm(N * 2, 0, 1), ncol = 2, nrow = N) %*% chol(R)
  X <- exp(matrix(time_grid, nrow = N, ncol = P, byrow = TRUE)^2 + Z[, 1])
  Y <- (matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 2])^2 +
    matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) * (-8) +
    (matrix(0, nrow = N, ncol = P, byrow = TRUE) + Z[, 2])

  mfD <- mfData(time_grid, list(X, Y))

  expect_snapshot_value(cor_kendall(mfD, ordering = 'max'), style = "serialize")
  expect_snapshot_value(cor_kendall(mfD, ordering = 'area'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MEI'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MHI'), style = "serialize")
})

test_that("cor_kendall() & cor_spearman() work on Case Study 11 from Dalia Valencia, Rosa Lillo, Juan Romo.", {
  withr::local_seed(1234)
  N <- 50
  P <- 50
  time_grid <- seq(0, 1, length.out = P)

  sigma_12 <- 0.
  R <- matrix(c(1, sigma_12, sigma_12, 1), ncol = 2, nrow = 2)
  Z <- matrix(rnorm(N * 2, 0, 1), ncol = 2, nrow = N) %*% chol(R)
  X <- exp(matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 1])
  Y <- sin(matrix(time_grid, nrow = N, ncol = P, byrow = TRUE) + Z[, 2])

  mfD <- mfData(time_grid, list(X, Y))

  expect_snapshot_value(cor_kendall(mfD, ordering = 'max'), style = "serialize")
  expect_snapshot_value(cor_kendall(mfD, ordering = 'area'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MEI'), style = "serialize")
  expect_snapshot_value(cor_spearman(mfD, ordering = 'MHI'), style = "serialize")
})
