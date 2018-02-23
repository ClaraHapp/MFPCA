context("Testing functions in univExpansion.R")

test_that("test expandBasis function", {
  expect_error(expandBasisFunction(scores = matrix(nrow = 2, ncol = 5), functions = funData(1:5, matrix(nrow = 3, ncol = 5))),
               "expandBasisFunction: number of scores for each observation and number of eigenfunctions does not match.")  
})

set.seed(1)
scores <- sapply(1:5, function(x){rnorm(20, sd = exp(-x))})
argvals <- list(seq(0, 1, 0.01))
functions <- funData:::efPoly(argvals[[1]], M = 5)

test_that("test univExpansion", {
  expect_error(MFPCA:::univExpansion(type = NULL, scores = scores, argvals = argvals, functions = functions),
              "Parameter 'type' is missing.")
  expect_error(MFPCA:::univExpansion(type = 5, scores = scores, argvals = argvals, functions = functions),
               "Parameter 'type' must be a character string. See ?univExpansion for details.", fixed = TRUE)
  expect_error(MFPCA:::univExpansion(type = "default", scores = NULL, argvals = argvals, functions = functions),
               "Parameter 'scores' is missing.")
  expect_error(MFPCA:::univExpansion(type = "default", scores = 1:5, argvals = argvals, functions = functions),
               "Parameter 'scores' must be passed as a matrix.")
  expect_error(MFPCA:::univExpansion(type = "default", scores = scores, argvals = NULL, functions = NULL),
               "Must pass 'argvals' if 'functions' is NULL.")
  expect_error(MFPCA:::univExpansion(type = "default", scores = scores, argvals = "Test", functions = NULL),
               "Parameter 'argvals' must be passed as a list.")
  expect_warning(MFPCA:::univExpansion(type = "default", scores = scores, argvals = argvals[[1]], functions = functions),
               "Parameter 'argvals' was passed as a vector and transformed to a list.")
  expect_error(MFPCA:::univExpansion(type = "default", scores = scores, argvals = argvals, functions = 5),
               "Parameter 'functions' must be a funData object.")
  expect_error(MFPCA:::univExpansion(type = "default", scores = scores[,1:4], argvals = argvals, functions = functions),
               "Number of scores per curve does not match the number of basis functions.")
  expect_error(MFPCA:::univExpansion(type = "default", scores = scores, argvals = argvals, functions = extractObs(functions, argvals = seq(0,0.5,0.01))),
               "The parameter 'argvals' does not match the argument values of 'functions'.")
  expect_error(MFPCA:::univExpansion(type = "default", scores = scores, argvals = argvals, functions = functions, params = c(a = 4)),
               "The parameter 'params' must be passed as a list.")
})

test_that("test univariate expansions 1D", {
  default1D <- MFPCA:::expandBasisFunction(scores = scores, argvals = argvals, functions = functions)
  expect_equal(nObs(default1D), 20)
  expect_equal(nObsPoints(default1D), 101)
  expect_equal(mean(norm(default1D)),  0.12735, tolerance = 1e-5)
  expect_equal(norm(default1D)[1], 0.07062, tolerance = 1e-5)
  
  spline1D <- MFPCA:::splineFunction1D(scores = scores, argvals = argvals, bs = "ps", m = 3, k = 5)
  expect_equal(nObs(spline1D), 20)
  expect_equal(nObsPoints(spline1D), 101)
  expect_equal(mean(norm(spline1D)),  0.11236, tolerance = 1e-5) 
  expect_equal(norm(spline1D)[1], 0.05326, tolerance = 1e-5)
  
  
  # wrapper function
  expandDefault1D <- MFPCA:::univExpansion(type = "default", scores = scores, argvals = argvals, functions = functions)
  expect_equal(expandDefault1D, default1D)
  
  expandSpline1D <- MFPCA:::univExpansion(type = "splines1D", scores = scores, argvals = argvals, functions = NULL, params = list(bs = "ps", m = 3, k = 5))
  expect_equal(expandSpline1D, spline1D)
  
  expandSpline1Dpen <- MFPCA:::univExpansion(type = "splines1Dpen", scores = scores, argvals = argvals, functions = NULL, params = list(bs = "ps", m = 3, k = 5))
  expect_equal(expandSpline1Dpen, spline1D) # spline1D, spline1Dpen have the same basis
  
  expandFPCA1D <- MFPCA:::univExpansion(type = "uFPCA", scores = scores, argvals = argvals, functions = functions)
  expect_equal(expandFPCA1D, default1D)
  
  expandGiven <- MFPCA:::univExpansion(type = "given", scores = scores, argvals = argvals, functions = funData:::efPoly(argvals[[1]], M = 5))
  expect_equal(expandGiven, default1D)
})

set.seed(2)
scores <- sapply(25:1, function(x){rnorm(20, sd = x/25)})
argvals <- list(seq(0, 1, 0.01), seq(-1, 1, 0.02))

test_that("test univariate expansions 2D", {
  default2D <- MFPCA:::expandBasisFunction(scores = scores, argvals = argvals, 
                                       functions = tensorProduct(funData:::efPoly(argvals[[1]], M = 5), funData:::efWiener(argvals[[2]], M = 5)))
  expect_equal(nObs(default2D), 20)
  expect_equal(nObsPoints(default2D), c(101,101))
  if(packageVersion("funData") != "1.0") # in 1.0, the tensor product was defined differently
  {
    expect_equal(mean(norm(default2D)),   10.27173, tolerance = 1e-5) 
    expect_equal(norm(default2D)[1], 12.08518, tolerance = 1e-5)
  }
  
  spline2D <- MFPCA:::splineFunction2D(scores = scores, argvals = argvals, bs = "ps", m = 3, k = 5)
  expect_equal(nObs(spline2D), 20)
  expect_equal(nObsPoints(spline2D), c(101,101))
  expect_equal(mean(norm(spline2D)),  2.84807, tolerance = 1e-5) 
  expect_equal(norm(spline2D)[1], 2.27017, tolerance = 1e-5)
  
  spline2Dpen <- MFPCA:::splineFunction2Dpen(scores = scores, argvals = argvals, bs = "ps", m = 3, k = 5)
  expect_equal(nObs(spline2Dpen), 20)
  expect_equal(nObsPoints(spline2Dpen), c(101,101))
  if(.Platform$endian == "big" | .Machine$sizeof.longdouble == 0)
    skip("Regression tests for spline2Dpen skipped on this architecture.")
  else
  {
    expect_equal(mean(norm(spline2Dpen)),  2.80049, tolerance = 1e-5) 
    expect_equal(norm(spline2Dpen)[1], 2.15060, tolerance = 1e-5)
  }
  
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
  if(packageVersion("funData") != "1.0") # in 1.0, the tensor product was defined differently
  {
    expect_equal(mean(norm(default3D)),  117.44730, tolerance = 1e-5) 
    expect_equal(norm(extractObs(default3D,1)), 118.64854, tolerance = 1e-5)
  }

  # wrapper function
  expandDefault3D <- MFPCA:::univExpansion(type = "default", scores = scores, argvals = argvals, 
                                           functions = tensorProduct(funData:::efPoly(argvals[[1]], M = 3),
                                                                     funData:::efWiener(argvals[[2]], M = 4),
                                                                     funData:::efFourier(argvals[[3]], M = 5)))
  expect_equal(expandDefault3D, default3D)
})

test_that("test univariate expansions 4D and higher", {
  set.seed(4)
  scores <- sapply(10:1, function(x){rnorm(20, sd = x/25)})
  argvals <- list(1:5,1:4,1:3,1:2)
  X <- array(runif(10*5*4*3*2), dim = c(10,5,4,3,2))
  
  default4D <- MFPCA:::expandBasisFunction(scores = scores, functions = funData(argvals, X))
  expect_equal(nObs(default4D), 20)
  expect_equal(nObsPoints(default4D), c(5,4,3,2))
  expect_equal(default4D@X[1,1,1,1,], c(0.29713, 0.57452), tol = 1e-5) # minimal check, as norm etc. are not implemented for 4D data
  
  # wrapper function
  expandDefault4D <- MFPCA:::univExpansion(type = "default", scores = scores, argvals = argvals,
                                           functions = funData(argvals, X))
  expect_equal(expandDefault4D, default4D)
})
