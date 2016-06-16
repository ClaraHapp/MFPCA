context("Testing functions in univExpansion.R")

test_that("test expandBasis function", {
  expect_error(expandBasisFunction(scores = matrix(nrow = 2, ncol = 5), functions = funData(1:5, matrix(nrow = 3, ncol = 5))),
               "expandBasisFunction: number of scores for each observation and number of eigenfunctions does not match.")  
})

test_that("test univariate expansions 1D", {
  set.seed(1)
  scores <- sapply(1:5, function(x){rnorm(20, sd = exp(-x))})
  argvals <- list(seq(0, 1, 0.01))
  
  default1D <- MFPCA:::expandBasisFunction(scores = scores, argvals = argvals, functions = funData:::efPoly(argvals[[1]], M = 5))
  expect_equal(nObs(default1D), 20)
  expect_equal(nObsPoints(default1D), 101)
  expect_equal(mean(norm(default1D)),  0.12734517) 
  expect_equal(norm(default1D)[1], 0.0706192863)
  
  spline1D <- MFPCA:::splineFunction1D(scores = scores, argvals = argvals, bs = "ps", m = 3, k = 5)
  expect_equal(nObs(spline1D), 20)
  expect_equal(nObsPoints(spline1D), 101)
  expect_equal(mean(norm(spline1D)),  0.112358657) 
  expect_equal(norm(spline1D)[1], 0.0532589197)
  
  
  # wrapper function
  expandDefault1D <- MFPCA:::univExpansion(type = "default", scores = scores, argvals = argvals, funData:::efPoly(argvals[[1]], M = 5))
  expect_equal(expandDefault1D, default1D)
  
  expandSpline1D <- MFPCA:::univExpansion(type = "splines1D", scores = scores, argvals = argvals, functions = NULL, params = list(bs = "ps", m = 3, k = 5))
  expect_equal(expandSpline1D, spline1D)
  
  expandSpline1Dpen <- MFPCA:::univExpansion(type = "splines1Dpen", scores = scores, argvals = argvals, functions = NULL, params = list(bs = "ps", m = 3, k = 5))
  expect_equal(expandSpline1Dpen, spline1D) # spline1D, spline1Dpen have the same basis
  
  expandFPCA1D <- MFPCA:::univExpansion(type = "uFPCA", scores = scores, argvals = argvals, functions = funData:::efPoly(argvals[[1]], M = 5))
  expect_equal(expandFPCA1D, default1D)
})

test_that("test univariate expansions 2D", {
  set.seed(2)
  scores <- sapply(25:1, function(x){rnorm(20, sd = x/25)})
  argvals <- list(seq(0, 1, 0.01), seq(-1, 1, 0.02))
  
  default2D <- MFPCA:::expandBasisFunction(scores = scores, argvals = argvals, 
                                       functions = tensorProduct(funData:::efPoly(argvals[[1]], M = 5), funData:::efWiener(argvals[[2]], M = 5)))
  expect_equal(nObs(default2D), 20)
  expect_equal(nObsPoints(default2D), c(101,101))
  expect_equal(mean(norm(default2D)),  10.2855118) 
  expect_equal(norm(default2D)[1], 12.0950688)
  
  spline2D <- MFPCA:::splineFunction2D(scores = scores, argvals = argvals, bs = "ps", m = 3, k = 5)
  expect_equal(nObs(spline2D), 20)
  expect_equal(nObsPoints(spline2D), c(101,101))
  expect_equal(mean(norm(spline2D)),  2.84807347) 
  expect_equal(norm(spline2D)[1], 2.27016721)
  
  spline2Dpen <- MFPCA:::splineFunction2Dpen(scores = scores, argvals = argvals, bs = "ps", m = 3, k = 5)
  expect_equal(nObs(spline2Dpen), 20)
  expect_equal(nObsPoints(spline2Dpen), c(101,101))
  expect_equal(mean(norm(spline2Dpen)),  2.80048833) 
  expect_equal(norm(spline2Dpen)[1], 2.15059827)
  
  # wrapper function
  expandDefault2D <- MFPCA:::univExpansion(type = "default", scores = scores, argvals = argvals, 
                                           functions = tensorProduct(funData:::efPoly(argvals[[1]], M = 5), funData:::efWiener(argvals[[2]], M = 5)))
  expect_equal(expandDefault2D, default2D)
  
  expandUMPCA2D <- MFPCA:::univExpansion(type = "UMPCA", scores = scores, argvals = argvals, 
                                         functions = tensorProduct(funData:::efPoly(argvals[[1]], M = 5), funData:::efWiener(argvals[[2]], M = 5)))
  expect_equal(expandUMPCA2D, default2D)
  
  expandFCPTPA2D <- MFPCA:::univExpansion(type = "FCP_TPA", scores = scores, argvals = argvals, 
                                         functions = tensorProduct(funData:::efPoly(argvals[[1]], M = 5), funData:::efWiener(argvals[[2]], M = 5)))
  expect_equal(expandFCPTPA2D, default2D)
  
  expandSpline2D <- MFPCA:::univExpansion(type = "splines2D", scores = scores, argvals = argvals, functions = NULL, params = list(bs = "ps", m = 3, k = 5))
  expect_equal(expandSpline2D, spline2D)
  
  expandSpline2Dpen <- MFPCA:::univExpansion(type = "splines2Dpen", scores = scores, argvals = argvals, functions = NULL, params = list(bs = "ps", m = 3, k = 5))
  expect_equal(expandSpline2Dpen, spline2Dpen)
})


test_that("test univariate expansions 3D", {
  set.seed(3)
  scores <- sapply(60:1, function(x){rnorm(20, sd = x/25)})
  argvals <- list(seq(0, 1, 0.01), seq(-1, 1, 0.02), seq(-0.5, 0.5, 0.05))
  
  default3D <- MFPCA:::expandBasisFunction(scores = scores, argvals = argvals, 
                                       functions = tensorProduct(funData:::efPoly(argvals[[1]], M = 3),
                                                                 funData:::efWiener(argvals[[2]], M = 4),
                                                                 funData:::efFourier(argvals[[3]], M = 5)))
  expect_equal(nObs(default3D), 20)
  expect_equal(nObsPoints(default3D), c(101,101,21))
  expect_equal(mean(norm(default3D)),  117.481898) 
  expect_equal(norm(extractObs(default3D,1)), 118.67828)
  
  # wrapper function
  expandDefault3D <- MFPCA:::univExpansion(type = "default", scores = scores, argvals = argvals, 
                                           functions = tensorProduct(funData:::efPoly(argvals[[1]], M = 3),
                                                                     funData:::efWiener(argvals[[2]], M = 4),
                                                                     funData:::efFourier(argvals[[3]], M = 5)))
  expect_equal(expandDefault3D, default3D)
})
