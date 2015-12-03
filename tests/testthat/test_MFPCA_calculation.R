context("Testing functions in MFCPA_calculation.R")

# rather plausibility checks
test_that("test .calcBasisIntegrals", {
  basis1D <- univExpansion("splines1D", scores = matrix(1:30, nrow = 3), argvals = seq(0, 1, 0.01), 
                           functions = NULL, params = list(k = 10, m = 2, bs = "ps"))
  calcBasis1D <- MFPCA:::.calcBasisIntegrals(basis1D@X, dimSupp = dimSupp(basis1D), argvals = basis1D@argvals)
  
  basis2D <- univExpansion("splines2D", scores = matrix(1:360, nrow = 3), argvals = list(seq(0, 1, 0.01), seq(-2, 2, 0.02)),
                                       functions = NULL, params = list(k = c(10,12), m = c(2,3), bs = "ps"))
  calcBasis2D <- MFPCA:::.calcBasisIntegrals(basis2D@X, dimSupp = dimSupp(basis2D), argvals = basis2D@argvals)
  
  # check 1D row sums
  expect_equal(rowSums(calcBasis1D), c(139.7564, 146.9259, 154.0953), tolerance = 1e-6)
  # check symmetry
  expect_equal(calcBasis1D, t(calcBasis1D))
  
  # check 2D row sums
  expect_equal(rowSums(calcBasis2D), c(3552908, 3572102, 3591295), tolerance = 1e-6)
  # check symmetry
  expect_equal(calcBasis2D, t(calcBasis2D))
})

# rather plausibility checks
test_that("test MFPCA main function", {
  # see also MFPCA examples
  set.seed(1)
  sim <-  simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
                          M = 5, eFunType = "Poly", eValType = "linear", N = 100)
  
  # check errors
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"))), 
                 "Function MFPCA_multidim: multivariate functional data object and univariate expansions must have the same length!")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"), list(type = "uFPCA")), bootstrap = TRUE), 
               "Specify number of bootstrap iterations.")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"), list(type = "uFPCA")), 
                     bootstrap = TRUE, nBootstrap = 10, bootstrapAlpha = -.1), 
               "Significance level for bootstrap confidence bands must be in (0,1).")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"), list(type = "uFPCA")), 
                     bootstrap = TRUE, nBootstrap = 10, bootstrapAlpha = 1.5), 
               "Significance level for bootstrap confidence bands must be in (0,1).")
  
  # check functionality
  expect_warning(uFPCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
                                                          list(type = "uFPCA"))), 
                 "Calculating a large percentage of principal components, approximation may not be appropriate.
            'approx.eigen' set to FALSE.")
  splines <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "splines1D", k = 10),
                                                            list(type = "splines1D", k = 10)))
  
  # values
  expect_identical(length(uFPCA$values), length(splines$values))
  expect_equal(sum(uFPCA$values), sum(uFPCA$values))
  expect_equal(uFPCA$values[1], 1.05175127)
  expect_equal(splines$values[1], 1.05266096)
  
  # functions
  expect_identical(nObs(uFPCA$functions), nObs(splines$functions))
  expect_equal(norm(uFPCA$functions), norm(splines$functions))
  expect_equal(norm(uFPCA$functions[[1]])[1], 0.57954971)
  expect_equal(norm(splines$functions[[1]])[1], 0.57956679)
})


test_that("test univariate decompositions 1D", {
  set.seed(1)
  f1 <- simFunData(seq(0,1,0.01), M = 10, eFunType = "Poly", eValType = "linear", N = 10)$simData

  spline1D <- MFPCA:::splineBasis1D(f1, bs = "ps", m = 3, k = 10)
  expect_identical(dim(spline1D$scores), c(10,10))
  expect_equal(mean(spline1D$scores),  17.22414491) 
  expect_identical(dim(spline1D$B), c(10,10))
  expect_equal(mean(spline1D$B), 0.0101531851)
  expect_false(spline1D$ortho)  
  expect_true(is.null(spline1D$functions))  
  expect_equal(spline1D$settings, list(bs = "ps", k = 10, m= c(3,3))) 
  
  spline1Dpen <- MFPCA:::splineBasis1Dpen(f1, bs = "ps", m = 3, k = 10)
  expect_identical(dim(spline1Dpen$scores), c(10,10))
  expect_equal(mean(spline1Dpen$scores),  17.19100716) 
  expect_identical(dim(spline1Dpen$B), c(10,10))
  expect_equal(mean(spline1Dpen$B), 0.0101531851)
  expect_false(spline1Dpen$ortho)  
  expect_true(is.null(spline1Dpen$functions))  
  expect_equal(spline1Dpen$settings, list(bs = "ps", k = 10, m= c(3,3))) 
  
  pca1D <- PACE(f1, pve = 0.95)
  expect_identical(pca1D$npc, 5)
  expect_identical(nObs(pca1D$fit), 10)
  expect_equal(mean(norm(pca1D$fit)), 4.41694367)
  expect_identical(dim(pca1D$scores), c(10,5))
  expect_equal(mean(pca1D$scores), -0.0107182127)
  expect_identical(nObs(pca1D$mu), 1)
  expect_equal(norm(pca1D$mu), 0.54801585)
  expect_identical(nObs(pca1D$functions), 5)
  expect_equal(norm(pca1D$functions), rep(1,5))
  expect_equal(sum(pca1D$values), 3.77553890)
  expect_equal(pca1D$values[1], 1.35690262)
  expect_equal(pca1D$sigma2, 0.0131059337)
})