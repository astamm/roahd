test_that("it handles properly a single test function in row matrix form.", {
  # Arrange
  time_grid <- seq(0, 1, length.out = 1e2)

  Data_ref <- matrix(c( 0  + sin( 2 * pi * time_grid ),
                        1  + sin( 2 * pi * time_grid ),
                        -1 + sin( 2 * pi * time_grid )
  ),
  nrow = 3, ncol = length( time_grid ), byrow = TRUE )

  Data_test <- matrix( c( 0.6 + sin( 2 * pi * time_grid ) ),
                        nrow = 1, ncol = length( time_grid ), byrow = TRUE )

  # Act
  actual <- BD_relative(Data_test, Data_ref)

  # Assert
  expected <- 2/3
  expect_equal(actual, expected)
})

test_that("it handles properly a single test function in column matrix form.", {
  # Arrange
  time_grid <- seq(0, 1, length.out = 1e2)

  Data_ref <- matrix(c( 0  + sin( 2 * pi * time_grid ),
                        1  + sin( 2 * pi * time_grid ),
                        -1 + sin( 2 * pi * time_grid )
  ),
  nrow = 3, ncol = length( time_grid ), byrow = TRUE )

  Data_test <- matrix( c( 0.6 + sin( 2 * pi * time_grid ) ),
                        nrow = length( time_grid ), ncol = 1, byrow = TRUE )

  # Act
  actual <- BD_relative(Data_test, Data_ref)

  # Assert
  expected <- 2/3
  expect_equal(actual, expected)
})

test_that("it handles properly a single test function in vector form.", {
  # Arrange
  time_grid <- seq(0, 1, length.out = 1e2)

  Data_ref <- matrix(c( 0  + sin( 2 * pi * time_grid ),
                        1  + sin( 2 * pi * time_grid ),
                        -1 + sin( 2 * pi * time_grid )
  ),
  nrow = 3, ncol = length( time_grid ), byrow = TRUE )

  Data_test <- 0.6 + sin( 2 * pi * time_grid )

  # Act
  actual <- BD_relative(Data_test, Data_ref)

  # Assert
  expected <- 2/3
  expect_equal(actual, expected)
})

test_that("it handles properly a single test function in 1D array form.", {
  # Arrange
  time_grid <- seq(0, 1, length.out = 1e2)

  Data_ref <- matrix(c( 0  + sin( 2 * pi * time_grid ),
                        1  + sin( 2 * pi * time_grid ),
                        -1 + sin( 2 * pi * time_grid )
  ),
  nrow = 3, ncol = length( time_grid ), byrow = TRUE )

  Data_test <- array( 0.6 + sin( 2 * pi * time_grid ), dim = length( time_grid ) )

  # Act
  actual <- BD_relative(Data_test, Data_ref)

  # Assert
  expected <- 2/3
  expect_equal(actual, expected)
})

test_that("it handles properly a single test function in row-like 2D array form.", {
  # Arrange
  time_grid <- seq(0, 1, length.out = 1e2)

  Data_ref <- matrix(c( 0  + sin( 2 * pi * time_grid ),
                        1  + sin( 2 * pi * time_grid ),
                        -1 + sin( 2 * pi * time_grid )
  ),
  nrow = 3, ncol = length( time_grid ), byrow = TRUE )

  Data_test <- array( 0.6 + sin( 2 * pi * time_grid ), dim = c( 1, length( time_grid ) ) )

  # Act
  actual <- BD_relative(Data_test, Data_ref)

  # Assert
  expected <- 2/3
  expect_equal(actual, expected)
})

test_that("it handles properly a single test function in column-like 2D array form.", {
  # Arrange
  time_grid <- seq(0, 1, length.out = 1e2)

  Data_ref <- matrix(c( 0  + sin( 2 * pi * time_grid ),
                        1  + sin( 2 * pi * time_grid ),
                        -1 + sin( 2 * pi * time_grid )
  ),
  nrow = 3, ncol = length( time_grid ), byrow = TRUE )

  Data_test <- array( 0.6 + sin( 2 * pi * time_grid ), dim = c( length( time_grid ), 1 ) )

  # Act
  actual <- BD_relative(Data_test, Data_ref)

  # Assert
  expected <- 2/3
  expect_equal(actual, expected)
})

test_that("it handles properly multiple test function.", {
  # Arrange
  time_grid <- seq(0, 1, length.out = 1e2)

  Data_ref <- matrix(c( 0  + sin( 2 * pi * time_grid ),
                        1  + sin( 2 * pi * time_grid ),
                        -1 + sin( 2 * pi * time_grid )
  ),
  nrow = 3, ncol = length( time_grid ), byrow = TRUE )

  Data_test <- matrix( c( 0.5  + sin( 2 * pi * time_grid ),
                           -0.5 + sin( 2 * pi * time_grid ),
                           1.1 + sin( 2 * pi * time_grid ) ),
                        nrow = 3, ncol = length( time_grid ), byrow = TRUE )

  # Act
  actual <- BD_relative(Data_test, Data_ref)

  # Assert
  expected <- c(2/3, 2/3, 0)
  expect_equal(actual, expected)
})
