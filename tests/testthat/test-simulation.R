test_that("`generate_gauss_fdata()` works as expected using `Cov` argument", {
  withr::local_seed(1234)
  expect_snapshot_value(generate_gauss_fdata(N, centerline, Cov = C1))
})

test_that("`generate_gauss_fdata()` works as expected using `CholCov` argument", {
  withr::local_seed(1234)
  expect_snapshot_value(generate_gauss_fdata(N, centerline, CholCov = CholC1))
})

test_that("`generate_gauss_mfdata()` works as expected using `listCov` argument", {
  withr::local_seed(1234)
  expect_snapshot_value(generate_gauss_mfdata(
    N, L,
    centerlines,
    correlations = c(0.5, 0.5, 0.5),
    listCov = list(C1, C2, C3)
  ))
})

test_that("`generate_gauss_mfdata()` works as expected using `listCholCov` argument", {
  withr::local_seed(1234)
  expect_snapshot_value(generate_gauss_mfdata(
    N, L,
    centerlines,
    correlations = c(0.5, 0.5, 0.5),
    listCholCov = list(CholC1, CholC2, CholC3)
  ))
})

test_that("`generate_gauss_mfdata()` fails when dimensions mismatch (1/3)", {
  withr::local_seed(1234)
  expect_error(generate_gauss_mfdata(
    N, L,
    centerlines,
    correlations = c(0.5, 0.5, 0.5),
    listCholCov = list(CholC1[-1, ], CholC2[-1, ], CholC3[-1, ])
  ))
})

test_that("`generate_gauss_mfdata()` fails when dimensions mismatch (2/3)", {
  withr::local_seed(1234)
  expect_error(generate_gauss_mfdata(
    N, L,
    centerlines[-1, ],
    correlations = c(0.5, 0.5, 0.5),
    listCov = c(C1, C2, C3),
    listCholCov = list(CholC1[-1, ], CholC2[-1, ], CholC3[-1, ])
  ))
})

test_that("`generate_gauss_mfdata()` fails when dimensions mismatch (3/3)", {
  withr::local_seed(1234)
  expect_error(generate_gauss_mfdata(
    N, L,
    centerlines[, -1],
    correlations = c(0.5, 0.5),
    listCov = c(C1, C2, C3),
    listCholCov = list(CholC1, CholC2, CholC3)
  ))
})
