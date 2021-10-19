test_that("`depthgram()` works as expected for all supported input types", {
  # Arrange
  withr::local_seed(1234)
  N <- 100
  P <- 100
  grid <- seq(0, 1, length.out = P)
  Cov <- exp_cov_function(grid, alpha = 0.3, beta = 0.4)

  Data <- list()
  Data[[1]] <- generate_gauss_fdata(
    N,
    centerline = sin(2 * pi * grid),
    Cov = Cov
  )
  Data[[2]] <- generate_gauss_fdata(
    N,
    centerline = sin(2 * pi * grid),
    Cov = Cov
  )
  names <- paste0("id_", 1:nrow(Data[[1]]))
  fD <- fData(grid, Data[[1]])
  mfD <- mfData(grid, Data)

  # Act
  actual_list <- depthgram(Data, marginal_outliers = TRUE, ids = names)
  actual_fData <- depthgram(fD, marginal_outliers = TRUE, ids = names)
  actual_mfData <- depthgram(mfD, marginal_outliers = TRUE, ids = names)

  # Assert
  expect_snapshot(actual_list)
  expect_snapshot(actual_fData)
  expect_snapshot(actual_mfData)
})
