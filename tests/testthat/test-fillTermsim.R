library(testthat)

testthat::test_that(".fillTermsim works as expected", {
    setClass("DummyTermSim", slots = c(termsim = "matrix"))
    m <- matrix(
        c(1, NA, 0.2,
          NA, 1, 0.3,
          0.4, 0.5, 1),
        nrow = 3, byrow = TRUE
    )
    rownames(m) <- colnames(m) <- c("A", "B", "C")
    x <- new("DummyTermSim", termsim = m)
    keep <- c("A", "C")
    out <- .fillTermsim(x, keep)
    testthat::expect_equal(dim(out), c(2, 2))
    testthat::expect_equal(rownames(out), keep)
    testthat::expect_equal(colnames(out), keep)
    testthat::expect_false(any(is.na(out)))
    testthat::expect_true(isSymmetric(out))
    testthat::expect_equal(as.vector(diag(out)), c(1, 1))
    testthat::expect_equal(out[1, 2], 0.2 + 0.4)
    testthat::expect_equal(out[2, 1], 0.2 + 0.4)
    expected <- matrix(c(1, 0.6, 0.6, 1), nrow = 2)
    rownames(expected) <- colnames(expected) <- keep
    testthat::expect_equal(out, expected)
})

testthat::test_that(".fillTermsim handles edge cases", {
    setClass("DummyTermSim2", slots = c(termsim = "matrix"))
    m <- matrix(1, nrow = 1)
    rownames(m) <- colnames(m) <- "A"
    x <- new("DummyTermSim2", termsim = m)
    out <- .fillTermsim(x, "A")
    testthat::expect_equal(dim(out), c(1, 1))
    testthat::expect_equal(as.vector(out), 1)

    m <- matrix(c(1, NA, NA, 1), nrow = 2)
    rownames(m) <- colnames(m) <- c("A", "B")
    x <- new("DummyTermSim2", termsim = m)
    out <- .fillTermsim(x, c("A", "B"))
    testthat::expect_equal(out[1, 2], 0)
    testthat::expect_equal(out[2, 1], 0)
})
