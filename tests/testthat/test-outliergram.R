test_that("`outliergram()` works as expected", {
  expect_snapshot(outliergram(fDout, display = FALSE))
})

test_that("`outliergram()` correctly identifies outliers with Fvalue = 1.5)", {
  expect_snapshot(outliergram(fDout, display = FALSE)$ID_outliers)
})

test_that("`outliergram()` correctly identifies outliers with Fvalue = 2.5)", {
  expect_snapshot(outliergram(fDout, Fvalue = 2.5, display = FALSE)$ID_outliers)
})

test_that("`outliergram()` correctly identifies outliers with auto-adjusted Fvalue", {
  skip_on_cran()
  expect_snapshot(
    outliergram(
      fDout,
      adjust = list(
        N_trials = 10,
        trial_size = 5 * fDout$N,
        TPR = 0.01,
        VERBOSE = FALSE
      ),
      display = FALSE)$ID_outliers
  )
})

test_that("`outliergram()` warns if unrecognized argument in `adjust` list", {
  expect_warning(expect_warning(
    outliergram(
      fDout,
      adjust = list(N_trials = 1, trial_size = fDout$N, foo = 'bar', baz = 'qux'),
      display = FALSE
    )
  ))
})
