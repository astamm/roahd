test_that("BD function works as expected.", {
  # Arrange
  time_grid <- seq(0, 1, length.out = 1e2)

  D <- matrix(data = c(
    1 + sin(2 * pi * time_grid),
    0 + sin(4 * pi * time_grid),
    1 - sin(pi * (time_grid - 0.2)),
    0.1 + cos(2 * pi * time_grid),
    0.5 + sin(3 * pi + time_grid),
    -2 + sin(pi * time_grid)),
    nrow = 6,
    ncol = length(time_grid),
    byrow = TRUE
  )

  # Act
  actual <- BD(D)

  # Assert
  expected <- c(1/3, 1/3, 1/3, 1/3, 14/30, 1/3)
  expect_equal(actual, expected)
})
