context("Testing functions in univDecomp.R")

test_that("test univariate decompositions 1D", {
  set.seed(1)
  f1 <- simFunData(seq(0,1,0.01), M = 10, eFunType = "Poly", eValType = "linear", N = 10)$simData
  
  # splines1D
  # check error
  expect_error(MFPCA:::splineBasis1D(tensorProduct(f1,f1), bs = "ps", m = 3, k = 10), 
               "splines1D is implemented for 1D functional data only.")
  # check functionality
  spline1D <- MFPCA:::splineBasis1D(f1, bs = "ps", m = 3, k = 10)
  expect_equal(dim(spline1D$scores), c(10,10))
  expect_equal(mean(spline1D$scores),  17.22414491) 
  expect_equal(dim(spline1D$B), c(10,10))
  expect_equal(mean(spline1D$B), 0.0101531851)
  expect_false(spline1D$ortho)  
  expect_null(spline1D$functions) 
  expect_equal(spline1D$settings, list(bs = "ps", k = 10, m= c(3,3))) 
  
  # splines1Dpen
  # check error
  expect_error(MFPCA:::splineBasis1Dpen(tensorProduct(f1,f1), bs = "ps", m = 3, k = 10), 
               "splines1Dpen is implemented for 1D functional data only.")
  # check functionality
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
  decompSpline1D <- MFPCA:::univDecomp(type = "splines1D", funDataObject = f1, bs = "ps", m = 3, k = 10)
  expect_equal(decompSpline1D, spline1D)
  
  decompSpline1Dpen <- MFPCA:::univDecomp(type = "splines1Dpen", funDataObject = f1, bs = "ps", m = 3, k = 10)
  expect_equal(decompSpline1Dpen, spline1Dpen)
  
  decompFPCA1D <- MFPCA:::univDecomp(type = "uFPCA", funDataObject = f1, pve = 0.95)
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

test_that("test UMPCA functionality", {
  A <- array(1:24, dim = c(3,4,2))
  
  # check ttv errors
  expect_error(MFPCA:::ttv(A, list(rep(1,3)), 2), "A and v have wrong dimensions!")
  expect_error(MFPCA:::ttv(A, list(rep(1,3), rep(2,3)), 1), "The parameters 'dim' and 'v' must have the same length!")
  
  # check ttv functionality 
  expect_equal(MFPCA:::ttv(A, list(rep(1,3)), 1), colSums(A,1))
  expect_equal(MFPCA:::ttv(A, list(rep(1,4)), 2), apply(A,3, rowSums))
  expect_equal(MFPCA:::ttv(A, list(rep(1,2)), 3), apply(A,2, rowSums))
  expect_equal(MFPCA:::ttv(A, list(rep(1,2), rep(1,3), rep(1,4)), c(3,1,2)), sum(A))
  
  # check internal function (double checked with matlab)
  expect_equal(MFPCA:::ttvCalculation(A, list(rep(1,3)), 1), colSums(A,1))
  expect_equal(MFPCA:::ttvCalculation(A, list(rep(1,4)), 2), apply(A,3, rowSums))
  expect_equal(MFPCA:::ttvCalculation(A, list(rep(1,2)), 3), apply(A,2, rowSums))
  
  # check maxeig
  expect_equal(as.numeric(MFPCA:::maxeig(diag(1:5))$lambda), 5)
  expect_equal(as.vector(MFPCA:::maxeig(diag(1:5))$x), c(0,0,0,0,1), tolerance = 1e-7)
  
  # see also 1D decompositions
  umpca2D <- MFPCA:::UMPCA(A, numP = 3)
  expect_equal(length(umpca2D$Us), 2)
  expect_equal(umpca2D$Us[[1]], matrix(1/sqrt(3), nrow = 3))
  expect_equal(umpca2D$Us[[2]], matrix(0.5, nrow = 4))
  expect_equal(dim(umpca2D$TXmean), c(3,4,1))
  expect_equal(umpca2D$TXmean[,,1], matrix(7:18, nrow = 3, ncol = 4))
  expect_equal(umpca2D$odrIdx, 1)
})

test_that("test univariate decompositions 2D", {
  set.seed(1)
  x1 <- seq(0,1,length.out=50)
  x2 <- seq(-1,1, length.out=75)
  f2 <- funData(argvals = list(x1, x2),
                X = aperm(replicate(10, outer(x1, cos(pi*x2))+matrix(rnorm(50*75, sd = 0.1), nrow = 50)), c(3,1,2)))
  
  # splines2D
  # check error
  expect_error(MFPCA:::splineBasis2D(funData(x1, rnorm(10)%o%x1), bs = "ps", m = 3, k = 10),
               "splines2D is implemented for 2D functional data (images) only.", fixed = TRUE)
  # check functionality
  spline2D <- MFPCA:::splineBasis2D(f2, bs = "ps", m = c(2,3), k = c(8,10))
  expect_equal(dim(spline2D$scores), c(10,80))
  expect_equal(mean(spline2D$scores),  -0.040453269) 
  expect_equal(dim(spline2D$B), c(80,80))
  expect_equal(sum(spline2D$B), 5.85240898)
  expect_false(spline2D$ortho)  
  expect_null(spline2D$functions)  
  expect_equal(spline2D$settings, list(bs = "ps", k = c(8,10), m =list(c(2,2), c(3,3)))) 
  
  # splines2Dpen
  # check error
  expect_error(MFPCA:::splineBasis2Dpen(funData(x1, rnorm(10)%o%x1), bs = "ps", m = 3, k = 10),
               "splines2Dpen is implemented for 2D functional data (images) only.", fixed = TRUE)
  # check functionality
  spline2Dpen <- MFPCA:::splineBasis2Dpen(extractObs(f2,1:2), bs = "ps", m = c(2,3), k = c(8,10))
  expect_equal(dim(spline2Dpen$scores), c(2,80))
  # expect_equal(mean(spline2Dpen$scores),  -0.0501778768) # different results for different mgcv versions...
  expect_equal(dim(spline2Dpen$B), c(80,80))
  expect_equal(sum(spline2Dpen$B), 2.00353973)
  expect_false(spline2Dpen$ortho)  
  expect_null(spline2Dpen$functions)
  expect_equal(spline2Dpen$settings, list(bs = "ps", k = c(8,10), m =list(c(2,2), c(3,3)))) 
  
  expect_warning(umpca2D <- MFPCA:::umpcaBasis(f2, npc = 4))
  expect_error(MFPCA:::umpcaBasis(funData(x1, t(sapply(1:5, function(x){x*x1}))), npc = 4), "UMPCA is implemented for (2D) image data only!", fixed = TRUE)
  expect_equal(dim(umpca2D$scores), c(10, 4))
  expect_equal(mean(umpca2D$scores),  0) 
  expect_equal(dim(umpca2D$B), c(4,4))
  expect_equal(sum(umpca2D$B), 0.00219103918)
  expect_equal(nObs(umpca2D$functions), 4)
  expect_equal(norm(umpca2D$functions)[1], 0.000523190474)
   
  set.seed(2)
  fcptpa2D <- MFPCA:::fcptpaBasis(f2, npc = 4, alphaRange = list(v = c(1e-4, 1e4), w = c(1e-4, 1e4)))
  expect_error(MFPCA:::fcptpaBasis(funData(x1, t(sapply(1:5, function(x){x*x1}))), npc = 4), "FCP_TPA is implemented for (2D) image data only!", fixed = TRUE)
  expect_equal(dim(fcptpa2D$scores), c(10, 4))
  expect_equal(mean(fcptpa2D$scores),  -6.34674399, tolerance = 1e-6) 
  expect_equal(dim(fcptpa2D$B), c(4,4))
  expect_equal(sum(fcptpa2D$B), 0.00213955703)
  expect_equal(nObs(fcptpa2D$functions), 4)
  expect_equal(norm(fcptpa2D$functions)[1], 0.000520573393)
  
  # wrapper function
  decompSpline2D <- MFPCA:::univDecomp(type = "splines2D", funDataObject = f2, bs = "ps", m = c(2,3), k = c(8,10))
  expect_equal(decompSpline2D, spline2D)
  
  decompSpline2Dpen <- MFPCA:::univDecomp(type = "splines2Dpen", funDataObject = extractObs(f2,1:2), bs = "ps", m = c(2,3), k = c(8,10))
  expect_equal(decompSpline2Dpen, spline2Dpen)
  
  expect_warning(decompUMPCA2D <- MFPCA:::univDecomp(type = "UMPCA", funDataObject = f2, npc = 4))
  expect_equal(decompUMPCA2D, umpca2D)
  
  set.seed(2)
  decompFCPTPA2D <- MFPCA:::univDecomp(type = "FCP_TPA", funDataObject = f2, npc = 4, alphaRange = list(v = c(1e-4, 1e4), w = c(1e-4, 1e4)))
  expect_equal(decompFCPTPA2D, fcptpa2D)
})