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
  oldPar <- par(no.readonly = TRUE)
  
  set.seed(1)
  
  ### simulate data (one-dimensional domains)
  sim <-  simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
                          M = 5, eFunType = "Poly", eValType = "linear", N = 100)
  
  # MFPCA based on univariate FPCA
  uFPCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
                                                          list(type = "uFPCA")))
  
  # MFPCA based on univariate spline expansions
  splines <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "splines1D", k = 10),
                                                            list(type = "splines1D", k = 10)))
  
  # flip to make results more clear
  uFPCA$functions <- flipFuns(sim$trueFuns, uFPCA$functions)
  splines$functions <- flipFuns(sim$trueFuns, splines$functions)
})


test_that("test univariate decompositions 1D", {
  set.seed(1)
  f1 <- simFunData(seq(0,1,0.01), M = 10, eFunType = "Poly", eValType = "linear", N = 10)$simData

  spline1D <- splineBasis1D(f1, bs = "ps", m = 3, k = 10)
  expect_equal(dim(spline1D$scores), c(10,10))
  expect_equal(mean(spline1D$scores),  -6.00325721, tolerance=1e-7) 
  expect_equal(dim(spline1D$B), c(10,10))
  expect_equal(mean(spline1D$B), 0.010153185, tolerance=1e-7)
  expect_false(spline1D$ortho)  
  expect_true(is.null(spline1D$functions))  
  expect_equal(spline1D$settings, list(bs = "ps", k = 10, m= c(3,3))) 
  
  spline1Dpen <- splineBasis1Dpen(f1, bs = "ps", m = 3, k = 10)
  expect_equal(dim(spline1Dpen$scores), c(10,10))
  expect_equal(mean(spline1Dpen$scores),  -6.04131817, tolerance=1e-7) 
  expect_equal(dim(spline1Dpen$B), c(10,10))
  expect_equal(mean(spline1Dpen$B), 0.010153185, tolerance=1e-7)
  expect_false(spline1Dpen$ortho)  
  expect_true(is.null(spline1Dpen$functions))  
  expect_equal(spline1Dpen$settings, list(bs = "ps", k = 10, m= c(3,3))) 
  
  pca1D <- PACE(f1, pve = 0.95)
  expect_equal(pca1D$npc, 6)
  expect_equal(nObs(pca1D$fit), 10)
  expect_equal(mean(norm(pca1D$fit)), 4.518176, tolerance = 1e-7)
  expect_equal(dim(pca1D$scores), c(10,6))
  expect_equal(mean(pca1D$scores), -0.004868896, tolerance = 1e-7)
  expect_equal(nObs(pca1D$mu), 1)
  expect_equal(norm(pca1D$mu), 0.4576708, tolerance = 1e-7)
  expect_equal(nObs(pca1D$functions), 6)
  expect_equal(norm(pca1D$functions), rep(1,6))
  expect_equal(pca1D$values, c(1.69311466400647, 0.912435761461481, 0.742524381973502, 0.339119548289354, 0.128649288432991, 0.112615157127846))
  expect_equal(pca1D$sigma2, 0.03328641, tolerance = 1e-7)
})