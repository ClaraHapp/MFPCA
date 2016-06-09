context("Testing FCP_TPA functionality")


test_that("normVec", {
 expect_equal(MFPCA:::normVec(c(1,0,0,0)), 1)
 
 x <- runif(4)
 expect_equal(MFPCA:::normVec(x), sqrt(sum(x^2)))
})