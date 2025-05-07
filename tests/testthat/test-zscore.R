library(testthat)

# Unit tests
test_that("z_score.cal function works correctly", {
  # Test basic z-score calculation
  x <- c(1, 2, 3, 4, 5)
  expected_z_scores <- c(-1.2649111, -0.6324555,  0.0000000,  0.6324555,  1.2649111)
  actual_z_scores <- .z_score_cal(x)

  expect_equal(
    actual_z_scores,
    expected_z_scores,
    tolerance=1e-6
  )

  # Test with NA values
  x_with_na <- c(1, 2, NA, 4, 5)
  expected_na_z_scores <- c(-1.0954451, -0.5477226,  0.5477226,  1.0954451)

  actual_na_z_scores <- .z_score_cal(x_with_na, NA.rm=TRUE)

  # Modify test to handle NA carefully
  expect_true(
        all.equal(actual_na_z_scores[-3], expected_na_z_scores,
                  tolerance=1e-6)
  )

  # Test with constant vector
  constant_vec <- rep(5, 5)
  expect_equal(
    .z_score_cal(constant_vec),
    rep(0, 5)
  )

  # Test rescale01 option
  x <- c(1, 2, 3, 4, 5)
  rescaled <- .z_score_cal(x, rescale01=TRUE)
  expect_true(min(rescaled) >= 0)
  expect_true(max(rescaled) <= 1)

})
