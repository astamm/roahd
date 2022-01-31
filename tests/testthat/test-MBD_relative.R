test_that("`MBD_relative()` works for a single test function in row matrix form", {
  time_grid <- seq(0, 1, length.out = 1e2)
  Data_ref <- matrix(c(
     0 + sin(2 * pi * time_grid),
     1 + sin(2 * pi * time_grid),
    -1 + sin(2 * pi * time_grid)
  ), nrow = 3, ncol = length(time_grid), byrow = TRUE)

  Data_test_1 <- matrix(c(
    0.6 + sin(2 * pi * time_grid)
    ), nrow = 1, ncol = length(time_grid), byrow = TRUE)

  expect_equal(MBD_relative(Data_test_1, Data_ref), 2/3)
})

test_that("`MBD_relative()` works for a single test function in column matrix form", {
  time_grid <- seq(0, 1, length.out = 1e2)
  Data_ref <- matrix(c(
     0 + sin(2 * pi * time_grid),
     1 + sin(2 * pi * time_grid),
    -1 + sin(2 * pi * time_grid)
  ), nrow = 3, ncol = length(time_grid), byrow = TRUE)

  Data_test_2 <- matrix(c(
    0.6 + sin(2 * pi * time_grid)
    ), nrow = length(time_grid), ncol = 1, byrow = TRUE)

  expect_equal(MBD_relative(Data_test_2, Data_ref), 2/3)
})

test_that("`MBD_relative()` works for a single test function in vector form", {
  time_grid <- seq(0, 1, length.out = 1e2)
  Data_ref <- matrix(c(
     0 + sin(2 * pi * time_grid),
     1 + sin(2 * pi * time_grid),
    -1 + sin(2 * pi * time_grid)
  ), nrow = 3, ncol = length(time_grid), byrow = TRUE)

  Data_test_3 <- 0.6 + sin(2 * pi * time_grid)

  expect_equal(MBD_relative(Data_test_3, Data_ref), 2/3)
})

test_that("`MBD_relative()` works for a single test function in 1D array form", {
  time_grid <- seq(0, 1, length.out = 1e2)
  Data_ref <- matrix(c(
     0 + sin(2 * pi * time_grid),
     1 + sin(2 * pi * time_grid),
    -1 + sin(2 * pi * time_grid)
  ), nrow = 3, ncol = length(time_grid), byrow = TRUE)

  Data_test_4 <- array(0.6 + sin(2 * pi * time_grid), dim = length(time_grid))

  expect_equal(MBD_relative(Data_test_4, Data_ref), 2/3)
})

test_that("`MBD_relative()` works for a single test function in row-like 2D array form", {
  time_grid <- seq(0, 1, length.out = 1e2)
  Data_ref <- matrix(c(
     0 + sin(2 * pi * time_grid),
     1 + sin(2 * pi * time_grid),
    -1 + sin(2 * pi * time_grid)
  ), nrow = 3, ncol = length(time_grid), byrow = TRUE)

  Data_test_5 <- array(
    0.6 + sin(2 * pi * time_grid),
    dim = c(1, length(time_grid))
  )

  expect_equal(MBD_relative(Data_test_5, Data_ref), 2/3)
})

test_that("`MBD_relative()` works for a single test function in column-like 2D array form", {
  time_grid <- seq(0, 1, length.out = 1e2)
  Data_ref <- matrix(c(
     0 + sin(2 * pi * time_grid),
     1 + sin(2 * pi * time_grid),
    -1 + sin(2 * pi * time_grid)
  ), nrow = 3, ncol = length(time_grid), byrow = TRUE)

  Data_test_6 <- array(
    0.6 + sin(2 * pi * time_grid),
    dim = c(length(time_grid), 1)
  )

  expect_equal(MBD_relative(Data_test_6, Data_ref), 2/3)
})

test_that("`MBD_relative()` works for multiple test functions", {
  time_grid <- seq(0, 1, length.out = 1e2)
  Data_ref <- matrix(c(
     0 + sin(2 * pi * time_grid),
     1 + sin(2 * pi * time_grid),
    -1 + sin(2 * pi * time_grid)
  ), nrow = 3, ncol = length(time_grid), byrow = TRUE)

  Data_test_7 <- matrix(c(
     0.5 + sin(2 * pi * time_grid),
    -0.5 + sin(2 * pi * time_grid),
     1.1 + sin(2 * pi * time_grid)
  ), nrow = 3, ncol = length(time_grid), byrow = TRUE)

  expect_equal(MBD_relative(Data_test_7, Data_ref), c(2/3, 2/3, 0))
})
