context("Testing functions in univDecomp.R")

set.seed(1)
s1 <- simFunData(seq(0,1,0.01), M = 10, eFunType = "Poly", eValType = "linear", N = 10)
f1 <- s1$simData

set.seed(1)
x1 <- seq(0,1,length.out=50)
x2 <- seq(-1,1, length.out=75)
f2 <- funData(argvals = list(x1, x2),
              X = aperm(replicate(10, outer(x1, cos(pi*x2))+matrix(rnorm(50*75, sd = 0.1), nrow = 50)), c(3,1,2)))

test_that("test univDecomp", {
  expect_error(MFPCA:::univDecomp(type = NULL, funDataObject = f1), 
               "Parameter 'type' is missing.")   
  expect_error(MFPCA:::univDecomp(type = 5, funDataObject = f1), 
               "Parameter 'type' must be a character string. See ?univDecomp for details.", fixed = TRUE)   
  expect_error(MFPCA:::univDecomp(type = "Test", funDataObject = 5), 
               "Parameter 'funDataObject' must be a funData object.")   
          })

test_that("test univariate decompositions 1D", {
  # splines1D
  # check error
  expect_error(MFPCA:::splineBasis1D(tensorProduct(f1,f1), bs = "ps", m = 3, k = 10), 
               "splines1D is implemented for 1D functional data only.")
  # check functionality
  spline1D <- MFPCA:::splineBasis1D(f1, bs = "ps", m = 3, k = 10)
  expect_equal(dim(spline1D$scores), c(10,10))
  expect_equal(mean(spline1D$scores),  17.22414, tolerance = 1e-5) 
  expect_equal(dim(spline1D$B), c(10,10))
  expect_equal(mean(spline1D$B), 0.01015, tolerance = 1e-5)
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
  expect_equal(mean(spline1Dpen$scores),  17.19101, tolerance = 1e-5) 
  expect_equal(dim(spline1Dpen$B), c(10,10))
  expect_equal(mean(spline1Dpen$B), 0.01015, tolerance = 1e-5)
  expect_false(spline1Dpen$ortho)  
  expect_null(spline1Dpen$functions) 
  expect_equal(spline1Dpen$settings, list(bs = "ps", k = 10, m= c(3,3))) 
  
  fpca <- MFPCA:::fpcaBasis(f1, pve = 0.95)
  expect_equal(dim(fpca$scores), c(10,5))
  expect_equal(mean(abs(fpca$scores)), 0.71489, tolerance = 1e-5) 
  expect_null(fpca$B) 
  expect_true(fpca$ortho)  
  expect_false(is.null(fpca$functions))  
  expect_equal(nObs(fpca$functions), 5)
  expect_equal(norm(fpca$functions), rep(1,5))
  
  # given
  expect_error(MFPCA:::givenBasis(funDataObject = f1, 
                                 functions = extractObs(s1$trueFuns, argvals = s1$trueFuns@argvals[[1]][1:10])),
              "Basis functions must be defined on the same domain as the observations.")
  
  expect_error(MFPCA:::givenBasis(funDataObject = f1, functions = s1$trueFuns, scores = as.matrix(1:5)),
               "Scores have wrong dimensions. Must be an N x K matrix with N the number of observations and K the number of basis functions.", fixed = FALSE)
  
  expect_warning(MFPCA:::givenBasis(funDataObject = f1, functions = s1$trueFuns, scores = 10*diag(10)),
               "Scores seem to be not demeaned. Please check.")
  
  given <- MFPCA:::givenBasis(funDataObject = f1, functions = s1$trueFuns)
  expect_equal(dim(given$scores), c(10,10))
  expect_equal(sum(given$scores), 8.029624)
  expect_equal(given$scores[1,1], -0.62431558)
  expect_equal(dim(given$B), c(10,10))
  expect_true(isSymmetric(given$B))
  expect_equal(sum(given$B), 10.6327967)
  expect_equal(given$B[1,1], 1)
  expect_false(given$ortho)
  expect_equal(given$functions, s1$trueFuns)
  
  # currently, inputs are not checked for plausibility!
  given2 <- MFPCA:::givenBasis(funDataObject = f1, functions = s1$trueFuns, scores = diag(10), ortho = TRUE)
  expect_true(given2$ortho)
  expect_null(given2$B)
  expect_equal(given2$scores, diag(10))
  
  
  # fda (only if available)
  if(requireNamespace("fda", quietly = TRUE) & utils::packageVersion("funData") > "1.2")
  {
    library("fda")
    
    # check functionality (all errors are checked in funData2fd)
    Fbasis <- create.fourier.basis(range(f1@argvals), nbasis=15) # Fourier basis with 15 basis functions
    fda <- MFPCA:::fdaBasis(f1, Fbasis)
    expect_equal(dim(fda$scores), c(10,15))
    expect_equal(mean(fda$scores),  0.086700, tolerance = 1e-5) 
    expect_equal(dim(fda$B), c(15,15))
    expect_equal(fda$B, diag(15), tolerance = 1e-5)
    expect_false(fda$ortho)  
    expect_false(is.null(fda$functions))  
    expect_equal(nObs(fda$functions), 15)
    expect_equal(norm(fda$functions), rep(1,15))
  }
  
  # wrapper function
  decompSpline1D <- MFPCA:::univDecomp(type = "splines1D", funDataObject = f1, bs = "ps", m = 3, k = 10)
  expect_equal(decompSpline1D, spline1D)
  
  decompSpline1Dpen <- MFPCA:::univDecomp(type = "splines1Dpen", funDataObject = f1, bs = "ps", m = 3, k = 10)
  expect_equal(decompSpline1Dpen, spline1Dpen)
  
  decompFPCA1D <- MFPCA:::univDecomp(type = "uFPCA", funDataObject = f1, pve = 0.95)
  expect_equal(decompFPCA1D, fpca)
  
  decompGiven <- MFPCA:::univDecomp(type = "given", funDataObject = f1, functions = s1$trueFuns)
  expect_equal(decompGiven, given)
  
  decompGiven2 <- MFPCA:::univDecomp(type = "given", funDataObject = f1, functions = s1$trueFuns, scores = diag(10), ortho = TRUE)
  expect_equal(decompGiven2, given2)
  
  if(requireNamespace("fda", quietly = TRUE) & utils::packageVersion("funData") > "1.2")
  {
    decompFDA <- MFPCA:::univDecomp(type = "fda", funDataObject = f1, Fbasis)
    expect_equal(decompFDA, fda)
  }
})

test_that("PACE function", {
  # check inputs
  expect_error(PACE(funDataObject = 1), "Parameter 'funDataObject' must be a funData or irregFunData object.")
  expect_error(PACE(funDataObject = f2), "PACE: Implemented only for funData objects with one-dimensional support.")
  expect_error(PACE(f1, predData = extractObs(f1, argvals = seq(0,0.5, 0.01))),
               "PACE: funDataObject and predData must be defined on the same domains!")
  expect_error(PACE(funDataObject = f1, nbasis = "10"), "Parameter 'nbasis' must be passed as a number > 0.")
  expect_error(PACE(funDataObject = f1, nbasis = 1:10), "Parameter 'nbasis' must be passed as a number > 0.")
  expect_error(PACE(funDataObject = f1, nbasis = -1), "Parameter 'nbasis' must be passed as a number > 0.")
  expect_error(PACE(funDataObject = f1, pve = "0"), "Parameter 'pve' must be passed as a number between 0 and 1.")
  expect_error(PACE(funDataObject = f1, pve = 0:1), "Parameter 'pve' must be passed as a number between 0 and 1.")
  expect_error(PACE(funDataObject = f1, pve = -1), "Parameter 'pve' must be passed as a number between 0 and 1.")
  expect_error(PACE(funDataObject = f1, pve = 2), "Parameter 'pve' must be passed as a number between 0 and 1.")
  expect_error(PACE(funDataObject = f1, npc = "10"), "Parameter 'npc' must be either NULL or passed as a number > 0.")
  expect_error(PACE(funDataObject = f1, npc = 1:10), "Parameter 'npc' must be either NULL or passed as a number > 0.")
  expect_error(PACE(funDataObject = f1, npc = -1), "Parameter 'npc' must be either NULL or passed as a number > 0.")
  expect_error(PACE(funDataObject = f1, makePD = "yes"), "Parameter 'makePD' must be passed as a logical.")
  expect_error(PACE(funDataObject = f1, cov.weight.type = 0), "Parameter 'cov.weight.type' must be passed as a character.")

  # see also 1D decompositions, fpcaBasis
  pca1D <- PACE(f1, pve = 0.95)
  expect_equal(pca1D$npc, 5)
  expect_equal(nObs(pca1D$fit), 10)
  expect_equal(mean(norm(pca1D$fit)), 4.41694, tolerance = 1e-5)
  expect_equal(dim(pca1D$scores), c(10,5))
  expect_equal(mean(abs(pca1D$scores)), 0.71489, tolerance = 1e-5)
  expect_equal(nObs(pca1D$mu), 1)
  expect_equal(norm(pca1D$mu), 0.54802, tolerance = 1e-5)
  expect_equal(nObs(pca1D$functions), 5)
  expect_equal(norm(pca1D$functions), rep(1,5))
  expect_equal(sum(pca1D$values), 3.77554, tolerance = 1e-5)
  expect_equal(pca1D$values[1], 1.35690, tolerance = 1e-5)
  expect_equal(pca1D$sigma2, 0.01311, tolerance = 1e-5)
  
  pcaPD <- PACE(f1, pve = 0.95, makePD = TRUE)
  expect_equal(pcaPD$npc, 5)
  expect_equal(nObs(pcaPD$fit), 10)
  expect_equal(mean(norm(pcaPD$fit)), 4.41734, tolerance = 1e-5)
  expect_equal(dim(pcaPD$scores), c(10,5))
  expect_equal(mean(abs(pcaPD$scores)), 0.71487, tolerance = 1e-5)
  expect_equal(nObs(pcaPD$mu), 1)
  expect_equal(norm(pcaPD$mu), 0.54802, tolerance = 1e-5)
  expect_equal(nObs(pcaPD$functions), 5)
  expect_equal(norm(pcaPD$functions), rep(1,5))
  expect_equal(sum(pcaPD$values), 3.77567, tolerance = 1e-5)
  expect_equal(pcaPD$values[1], 1.35694, tolerance = 1e-5)
  expect_equal(pcaPD$sigma2, 0.01102, tolerance = 1e-5)
  
  # test also for irregular data
  # suppress warnings in transition to new RNG, as proposed by CRAN maintainers
  suppressWarnings(RNGversion("3.5.0")) 
  set.seed(2)
  f1sparse <- sparsify(f1, minObs=  20, maxObs = 50)
  i1 <- irregFunData(argvals = apply(f1sparse@X,1, function(x){f1sparse@argvals[[1]][which(!is.na(x))]}), X = apply(f1sparse@X, 1, na.omit))
  
  pca1Dirreg <- PACE(i1, pve = 0.95)
  expect_equal(pca1Dirreg$npc, 4)
  expect_equal(nObs(pca1Dirreg$fit), 10)
  expect_equal(mean(norm(pca1Dirreg$fit)), 4.12973, tolerance = 1e-5)
  expect_equal(dim(pca1Dirreg$scores), c(10,4))
  expect_equal(mean(abs(pca1Dirreg$scores)), 0.71889, tolerance = 1e-5)
  expect_equal(nObs(pca1Dirreg$mu), 1)
  expect_equal(norm(pca1Dirreg$mu), 0.99891, tolerance = 1e-5)
  expect_equal(nObs(pca1Dirreg$functions), 4)
  expect_equal(norm(pca1Dirreg$functions), rep(1,4), tolerance = 1e-5)
  expect_equal(sum(pca1Dirreg$values), 3.65081, tolerance = 1e-5)
  expect_equal(pca1Dirreg$values[1], 2.10068, tolerance = 1e-5)
  expect_equal(pca1Dirreg$sigma2, 0.48830, tolerance = 1e-5)
})

test_that("test UMPCA functionality", {
  A <- array(1:24, dim = c(3,4,2))
  
  # check ttv errors
  expect_error(MFPCA:::ttv(NULL, list(rep(1,3)), 1), "Parameter 'A' must be an array.")
  expect_error(MFPCA:::ttv(A, rep(1,3), 1), "Parameter 'v' must be passed as a list.")
  expect_error(MFPCA:::ttv(A, list(rep(1,3)), "Test"), "Parameter 'dim' must be passed as a vector of numerics.")
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
  expect_equal(as.vector(MFPCA:::maxeig(diag(1:5))$x), c(0,0,0,0,1), tolerance = 1e-5)
  
  # check UMPCA inputs
  expect_error(MFPCA:::UMPCA(1, numP = 3), "Parameter 'TX' must be passed as an array.")
  expect_error(MFPCA:::UMPCA(A, numP = "Test"), "Parameter 'numP' must be passed as a number.")
  expect_error(MFPCA:::UMPCA(A, numP = 1:3), "Parameter 'numP' must be passed as a number.")
  
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
  # splines2D
  # check error
  expect_error(MFPCA:::splineBasis2D(funData(x1, rnorm(10)%o%x1), bs = "ps", m = 3, k = 10),
               "splines2D is implemented for 2D functional data (images) only.", fixed = TRUE)
  # check functionality
  spline2D <- MFPCA:::splineBasis2D(f2, bs = "ps", m = c(2,3), k = c(8,10))
  expect_equal(dim(spline2D$scores), c(10,80))
  expect_equal(mean(spline2D$scores),  -0.04045, tolerance = 1e-5) 
  expect_equal(dim(spline2D$B), c(80,80))
  expect_equal(sum(spline2D$B), 5.85241, tolerance = 1e-5)
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
  expect_equal(sum(spline2Dpen$B), 2.00354, tolerance = 1e-5)
  expect_false(spline2Dpen$ortho)  
  expect_null(spline2Dpen$functions)
  expect_equal(spline2Dpen$settings, list(bs = "ps", k = c(8,10), m =list(c(2,2), c(3,3)))) 
  
  expect_warning(umpca2D <- MFPCA:::umpcaBasis(f2, npc = 4))
  expect_error(MFPCA:::umpcaBasis(funData(x1, t(sapply(1:5, function(x){x*x1}))), npc = 4), "UMPCA is implemented for (2D) image data only!", fixed = TRUE)
  expect_equal(dim(umpca2D$scores), c(10, 4))
  expect_equal(mean(umpca2D$scores),  0) 
  expect_equal(dim(umpca2D$B), c(4,4))
  expect_equal(sum(umpca2D$B), 0.00219, tolerance = 1e-5)
  expect_equal(nObs(umpca2D$functions), 4)
  expect_equal(norm(umpca2D$functions)[1], 0.000523, tolerance = 1e-6) # very small value -> check with higher tolerance
   
  set.seed(2)
  fcptpa2D <- MFPCA:::fcptpaBasis(f2, npc = 4, alphaRange = list(v = c(1e-4, 1e4), w = c(1e-4, 1e4)))
  expect_error(MFPCA:::fcptpaBasis(funData(x1, t(sapply(1:5, function(x){x*x1}))), npc = 4), "FCP_TPA is implemented for (2D) image data only!", fixed = TRUE)
  expect_equal(dim(fcptpa2D$scores), c(10, 4))
  expect_equal(mean(fcptpa2D$scores),  -6.34674, tolerance = 1e-5) 
  expect_equal(dim(fcptpa2D$B), c(4,4))
  expect_equal(sum(fcptpa2D$B), 0.00214, tolerance = 1e-5)
  expect_equal(nObs(fcptpa2D$functions), 4)
  expect_equal(norm(fcptpa2D$functions)[1], 0.000521, tolerance = 1e-6) # very small value -> check with higher tolerance
  expect_equal(fcptpa2D$values[1], 643.5194)
  expect_equal(sum(fcptpa2D$values), 658.5576, tolerance = 1e-5)
  
  set.seed(2)
  fcptpa2Dnorm <- MFPCA:::fcptpaBasis(f2, npc = 4, alphaRange = list(v = c(1e-4, 1e4), w = c(1e-4, 1e4)), normalize = TRUE)
  expect_equal(mean(fcptpa2Dnorm$scores), -278.1599, tolerance = 1e-5) 
  expect_true(fcptpa2Dnorm$ortho)
  expect_null(fcptpa2Dnorm$B)
  expect_equal(norm(fcptpa2Dnorm$functions), rep(1, nObs(fcptpa2Dnorm$functions)))
  expect_equal(fcptpa2Dnorm$values[1],  0.33500, tolerance = 1e-5)
  expect_equal(sum(fcptpa2Dnorm$values), 0.34261, tolerance = 1e-5)
  
  # wrapper function
  decompSpline2D <- MFPCA:::univDecomp(type = "splines2D", funDataObject = f2, bs = "ps", m = c(2,3), k = c(8,10))
  expect_equal(decompSpline2D, spline2D, tolerance = 1e-5)
  
  decompSpline2Dpen <- MFPCA:::univDecomp(type = "splines2Dpen", funDataObject = extractObs(f2,1:2), bs = "ps", m = c(2,3), k = c(8,10))
  expect_equal(decompSpline2Dpen, spline2Dpen, tolerance = 1e-5)
  
  expect_warning(decompUMPCA2D <- MFPCA:::univDecomp(type = "UMPCA", funDataObject = f2, npc = 4))
  expect_equal(decompUMPCA2D, umpca2D, tolerance = 1e-5)
  
  set.seed(2)
  decompFCPTPA2D <- MFPCA:::univDecomp(type = "FCP_TPA", funDataObject = f2, npc = 4, alphaRange = list(v = c(1e-4, 1e4), w = c(1e-4, 1e4)))
  expect_equal(decompFCPTPA2D, fcptpa2D, tolerance = 1e-5)
  
  set.seed(2)
  decompFCPTPA2Dnorm <- MFPCA:::univDecomp(type = "FCP_TPA", funDataObject = f2, npc = 4, alphaRange = list(v = c(1e-4, 1e4), w = c(1e-4, 1e4)), normalize = TRUE)
  expect_equal(decompFCPTPA2Dnorm, fcptpa2Dnorm, tolerance = 1e-5)
})