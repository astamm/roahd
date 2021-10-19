test_that("`BD()` handles equally arrays and fData objects", {
  expect_identical(
    BD(fD_restyling),
    BD(Data_restyling)
  )
})

test_that("`MBD(, manage_ties = TRUE)` handles equally arrays and fData objects", {
  expect_identical(
    MBD(fD_restyling, manage_ties = TRUE),
    MBD(Data_restyling, manage_ties = TRUE)
  )
})

test_that("`MBD(, manage_ties = FALSE)` handles equally arrays and fData objects", {
  expect_identical(
    MBD(fD_restyling, manage_ties = FALSE),
    MBD(Data_restyling, manage_ties = FALSE)
  )
})

test_that("`EI()` handles equally arrays and fData objects", {
  expect_identical(
    EI(fD_restyling),
    EI(Data_restyling)
  )
})

test_that("`MEI()` handles equally arrays and fData objects", {
  expect_identical(
    MEI(fD_restyling),
    MEI(Data_restyling)
  )
})

test_that("`HI()` handles equally arrays and fData objects", {
  expect_identical(
    HI(fD_restyling),
    HI(Data_restyling)
  )
})

test_that("`MHI()` handles equally arrays and fData objects", {
  expect_identical(
    MHI(fD_restyling),
    MHI(Data_restyling)
  )
})

test_that("`HRD()` handles equally arrays and fData objects", {
  expect_identical(
    HRD(fD_restyling),
    HRD(Data_restyling)
  )
})

test_that("`MHRD()` handles equally arrays and fData objects", {
  expect_identical(
    MHRD(fD_restyling),
    MHRD(Data_restyling)
  )
})

test_that("`BD_relative()` and `MBD_relative()` handle equally arrays and fData objects", {
  centerline_1 <- sin(2 * pi * grid)
  centerline_2 <- sin(2 * pi * grid) + 0.5
  Data_reference <- generate_gauss_fdata(N, centerline_1, Cov)
  Data_target <- generate_gauss_fdata(N, centerline_2, Cov)
  fD_target <- fData(grid, Data_target)
  fD_reference <- fData(grid, Data_reference)

  expect_identical(
    BD_relative(fD_target, fD_reference),
    BD_relative(Data_target, Data_reference)
  )
  expect_identical(
    MBD_relative(fD_target, fD_reference),
    MBD_relative(Data_target, Data_reference)
  )
})
