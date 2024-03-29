test_that("hi() & mhi() work as expected", {
  N <- 20
  time_grid <- seq(0, 1, length.out = N * 1e2)

  Data <- matrix(0, nrow = N, ncol = length(time_grid))
  for (iObs in 1:N)
    Data[iObs, ] <- as.numeric(time_grid >= (iObs - 1) / N & time_grid < iObs / N)
  Data[N, length(time_grid)] <- 1

  expect_equal(HI(Data), rep(1 / N, N))
  expect_equal(MHI(Data), rep((N + (N - 1)^2) / N^2 , N))
})

test_that("hi() works on test by James Long (TAMU)", {
  yints <- c( 1.27, .927, 1/2, .217, 0)
  slopes <- c(-1, -1, 0, 1, 1)
  time_grid <- (0:100) / 100

  Data <- matrix(0, nrow = length(yints), ncol = length(time_grid))
  for (i in 1:length(yints))
    Data[i, ] <- yints[i] + time_grid * slopes[i]

  expect_equal(HI(Data), c(0.4, 0.2, 0.2, 0.4, 0.2))
})
