test_that("`MBD()` works as expected in presence of ties", {
  # Arrange
  D <- matrix(c(
    c(1, 0.5, 0.25, 0.1, 0.05),
    c(1, 0.75, 0.25, 0.2, 0.1),
    c(1, 0.7, 0.25, 0.25, 0.15),
    c(1, 0.9, 0.35, 0.3, 0.25),
    c(1, 0.6, 0.25, 0.2, 0.2) ,
    c(0.9, 0.8, 0.25, 0.1, 0.08),
    c(1, 0.4, 0.3, 0.2, 0.1),
    c(1, 0.4, 0.3, 0.2, 0.1)
  ), ncol = 5, nrow = 8, byrow = TRUE)

  # Act
  N <- nrow(D)
  P <- ncol(D)
  depths <- rep(0, N)
  for (i in 1:N)
  {
    for (j in 1:(N - 1))
    {
      for (k in (j + 1):N)
      {
        for (r in 1:P)
        {
          if ((D[j, r] - D[i, r]) * (D[k, r] - D[i, r]) <= 0)
            depths[i] <- depths[i] + 1
        }
      }
    }
  }
  depths <- depths / (N * (N - 1) / 2 * P)

  # Assert
  expect_equal(depths, MBD(D, manage_ties = TRUE))
})

test_that("`MBD()` works as expected without checking for ties", {
  N <- 3
  P <- 1e2
  time_grid <- seq(0, 1, length.out = P)
  Data <- matrix(c(
     0 + sin(2 * pi * time_grid),
     1 + sin(2 * pi * time_grid),
    -1 + sin(2 * pi * time_grid)
  ), nrow = 3, ncol = length(time_grid), byrow = TRUE)

  expect_equal(
    MBD(Data, manage_ties = TRUE),
    MBD(Data, manage_ties = FALSE)
  )
})
