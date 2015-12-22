context("Testing functions in univDecomp.R")

test_that("test univariate decompositions 1D", {
  set.seed(1)
  f1 <- simFunData(seq(0,1,0.01), M = 10, eFunType = "Poly", eValType = "linear", N = 10)$simData
  
  spline1D <- MFPCA:::splineBasis1D(f1, bs = "ps", m = 3, k = 10)
  expect_equal(dim(spline1D$scores), c(10,10))
  expect_equal(mean(spline1D$scores),  17.22414491) 
  expect_equal(dim(spline1D$B), c(10,10))
  expect_equal(mean(spline1D$B), 0.0101531851)
  expect_false(spline1D$ortho)  
  expect_null(spline1D$functions) 
  expect_equal(spline1D$settings, list(bs = "ps", k = 10, m= c(3,3))) 
  
  spline1Dpen <- MFPCA:::splineBasis1Dpen(f1, bs = "ps", m = 3, k = 10)
  expect_equal(dim(spline1Dpen$scores), c(10,10))
  expect_equal(mean(spline1Dpen$scores),  17.19100716) 
  expect_equal(dim(spline1Dpen$B), c(10,10))
  expect_equal(mean(spline1Dpen$B), 0.0101531851)
  expect_false(spline1Dpen$ortho)  
  expect_null(spline1Dpen$functions) 
  expect_equal(spline1Dpen$settings, list(bs = "ps", k = 10, m= c(3,3))) 
  
  fpca <- MFPCA:::fpcaBasis(f1, pve = 0.95)
  expect_equal(dim(fpca$scores), c(10,5))
  expect_equal(mean(fpca$scores), -0.0107182127) 
  expect_null(fpca$B) 
  expect_true(fpca$ortho)  
  expect_false(is.null(fpca$functions))  
  expect_equal(nObs(fpca$functions), 5)
  expect_equal(norm(fpca$functions), rep(1,5))
  
  # wrapper function
  decompSpline1D <- MFPCA:::univDecomp(type = "splines1D", data = f1, params = list(bs = "ps", m = 3, k = 10))
  expect_equal(decompSpline1D, spline1D)
  
  decompSpline1Dpen <- MFPCA:::univDecomp(type = "splines1Dpen", data = f1, params = list(bs = "ps", m = 3, k = 10))
  expect_equal(decompSpline1Dpen, spline1Dpen)
  
  decompFPCA1D <- MFPCA:::univDecomp(type = "uFPCA", data = f1, params = list(pve = 0.95))
  expect_equal(decompFPCA1D, fpca)
})

test_that("PACE function", {
  set.seed(1)
  f1 <- simFunData(seq(0,1,0.01), M = 10, eFunType = "Poly", eValType = "linear", N = 10)$simData
  
  # see also 1D decompositions, fpcaBasis
  pca1D <- PACE(f1, pve = 0.95)
  expect_equal(pca1D$npc, 5)
  expect_equal(nObs(pca1D$fit), 10)
  expect_equal(mean(norm(pca1D$fit)), 4.41694367)
  expect_equal(dim(pca1D$scores), c(10,5))
  expect_equal(mean(pca1D$scores), -0.0107182127)
  expect_equal(nObs(pca1D$mu), 1)
  expect_equal(norm(pca1D$mu), 0.54801585)
  expect_equal(nObs(pca1D$functions), 5)
  expect_equal(norm(pca1D$functions), rep(1,5))
  expect_equal(sum(pca1D$values), 3.77553890)
  expect_equal(pca1D$values[1], 1.35690262)
  expect_equal(pca1D$sigma2, 0.0131059337)
})

test_that("test univariate decompositions 2D", {
  set.seed(1)
  x1 <- seq(0,1,length.out=50)
  x2 <- seq(-1,1, length.out=75)
  f2 <- funData(argvals = list(x1, x2),
                X = aperm(replicate(10, outer(x1, cos(pi*x2))+matrix(rnorm(50*75, sd = 0.1), nrow = 50)), c(3,1,2)))
  
  spline2D <- MFPCA:::splineBasis2D(f2, bs = "ps", m = c(2,3), k = c(8,10))
  expect_equal(dim(spline2D$scores), c(10,80))
  expect_equal(mean(spline2D$scores),  -0.040453269) 
  expect_equal(dim(spline2D$B), c(80,80))
  expect_equal(sum(spline2D$B), 5.85240898)
  expect_false(spline2D$ortho)  
  expect_null(spline2D$functions)  
  expect_equal(spline2D$settings, list(bs = "ps", k = c(8,10), m =list(c(2,2), c(3,3)))) 
  
  spline2Dpen <- MFPCA:::splineBasis2Dpen(extractObs(f2,1:2), bs = "ps", m = c(2,3), k = c(8,10))
  expect_equal(dim(spline2Dpen$scores), c(2,80))
  expect_equal(mean(spline2Dpen$scores),  -0.0501778768) 
  expect_equal(dim(spline2Dpen$B), c(80,80))
  expect_equal(sum(spline2Dpen$B), 2.00353973)
  expect_false(spline2Dpen$ortho)  
  expect_null(spline2Dpen$functions)
  expect_equal(spline2Dpen$settings, list(bs = "ps", k = c(8,10), m =list(c(2,2), c(3,3)))) 
  
  # wrapper function
  decompSpline2D <- MFPCA:::univDecomp(type = "splines2D", data = f2, params = list(bs = "ps", m = c(2,3), k = c(8,10)))
  expect_equal(decompSpline2D, spline2D)
  
  decompSpline2Dpen <- MFPCA:::univDecomp(type = "splines2Dpen", data = extractObs(f2,1:2), params = list(bs = "ps", m = c(2,3), k = c(8,10)))
  expect_equal(decompSpline2Dpen, spline2Dpen)
})