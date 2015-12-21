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
               'Significance level for bootstrap confidence bands must be in (0,1).')
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
  expect_equal(length(uFPCA$values), length(splines$values))
  expect_equal(sum(uFPCA$values), sum(uFPCA$values))
  expect_equal(uFPCA$values[1], 1.05175127)
  expect_equal(splines$values[1], 1.05266096)
  
  # functions
  expect_equal(nObs(uFPCA$functions), nObs(splines$functions))
  expect_equal(norm(uFPCA$functions), norm(splines$functions))
  expect_equal(norm(uFPCA$functions[[1]])[1], 0.57954971)
  expect_equal(norm(splines$functions[[1]])[1], 0.57956679)
})


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
  
  dct2D <- MFPCA:::dctBasis2D(f2, qThresh = 0.95)
  expect_equal(dim(dct2D$scores), c(10, 3750))
  expect_equal(length(dct2D$scores@i),  1880) 
  expect_equal(dct2D$scores@x[1],  -0.0196469613) 
  expect_equal(dim(dct2D$B), c(3750, 3750))
  expect_equal(dct2D$B@x[1], 0.202642367)
  expect_equal(var(diff(dct2D$B@x)), 0)
  expect_false(dct2D$ortho)  
  expect_null(dct2D$functions)
  
  # wrapper function
  decompSpline2D <- MFPCA:::univDecomp(type = "splines2D", data = f2, params = list(bs = "ps", m = c(2,3), k = c(8,10)))
  expect_equal(decompSpline2D, spline2D)
  
  decompSpline2Dpen <- MFPCA:::univDecomp(type = "splines2Dpen", data = extractObs(f2,1:2), params = list(bs = "ps", m = c(2,3), k = c(8,10)))
  expect_equal(decompSpline2Dpen, spline2Dpen)
  
  decompDCT2D<- MFPCA:::univDecomp(type = "DCT2D", data = f2, params = list(qThresh = 0.95))
  expect_equal(decompDCT2D, dct2D)
})


test_that("test univariate decompositions 3D", {
  set.seed(3)
  x1 <- seq(0, 1, length.out = 40)
  x2 <- seq(-1, 1, length.out = 30)
  x3 <- seq(0, 0.5, length.out = 20)
  f3 <- funData(argvals = list(x1, x2, x3), X = replicate(20, array(rnorm(10*40*30), dim = c(10, 40, 30))))
  
  dct3D <- MFPCA:::dctBasis3D(f3, qThresh = 0.95)
  expect_equal(dim(dct3D$scores), c(10, 23998))
  expect_equal(length(dct3D$scores@i),  12000) 
  expect_equal(dct3D$scores@x[1],  -0.071277532) 
  expect_equal(dim(dct3D$B), c(23998, 23998))
  expect_equal(dct3D$B@x[1], 0.032251534)
  expect_equal(var(diff(dct3D$B@x)), 0)
  expect_false(dct3D$ortho)  
  expect_null(dct3D$functions)
  
  # wrapper function
  decompDCT3D<- MFPCA:::univDecomp(type = "DCT3D", data = f3, params = list(qThresh = 0.95))
  expect_equal(decompDCT3D, dct3D)
})

test_that("Test fftw", {
  set.seed(4)
  image2D <- outer(seq(0,1,0.01), sin(seq(-1,1,0.02))) + matrix(rnorm(101^2, sd = 0.2), nrow= 101)
  image3D <- aperm(sapply(seq(0,10,0.5), 
                          function(x){x*outer(seq(0,1,0.01), sin(seq(-1,1,0.02))) + matrix(rnorm(101^2, sd = 0.2), nrow= 101)}, 
                          simplify = "array"), c(3,1,2))
  
  expect_error(MFPCA:::dct2D(image3D, qThresh = 0.9), "dct2D can handle only 2D images")
  expect_error(MFPCA:::dct3D(image2D, qThresh = 0.9), "dct3D can handle only 3D images")
  
  fftw2D <- MFPCA:::dct2D(image2D, qThresh = 0.95)
  expect_equal(length(fftw2D$ind), 510)
  expect_equal(fftw2D$ind[1:3], c(2,92,102))
  expect_equal(mean(fftw2D$val), -0.00049897171)
  expect_equal(fftw2D$val[1], 0.018282019)
  
  fftw3D <- MFPCA:::dct3D(image3D, qThresh = 0.95)
  expect_equal(length(fftw3D$ind), 10711)
  expect_equal(fftw3D$ind[1:3], c(102, 103, 110))
  expect_equal(mean(fftw3D$val), 0.00027522909)
  expect_equal(fftw3D$val[1], -0.32200299)
  
})
