context("Testing functions in MFCPA_calculation.R")

# rather plausibility checks
test_that("test calcBasisIntegrals", {
  basis1D <- univExpansion("splines1D", scores = matrix(1:30, nrow = 3), argvals = list(seq(0, 1, 0.01)), 
                           functions = NULL, params = list(k = 10, m = 2, bs = "ps"))
  calcBasis1D <- MFPCA:::calcBasisIntegrals(basis1D@X, dimSupp = dimSupp(basis1D), argvals = basis1D@argvals)
  
  basis2D <- univExpansion("splines2D", scores = matrix(1:360, nrow = 3), argvals = list(seq(0, 1, 0.01), seq(-2, 2, 0.02)),
                                       functions = NULL, params = list(k = c(10,12), m = c(2,3), bs = "ps"))
  calcBasis2D <- MFPCA:::calcBasisIntegrals(basis2D@X, dimSupp = dimSupp(basis2D), argvals = basis2D@argvals)
  
  # check 1D row sums
  expect_equal(rowSums(calcBasis1D), c(139.7564, 146.9259, 154.0953), tolerance = 1e-5)
  # check symmetry
  expect_equal(calcBasis1D, t(calcBasis1D))
  
  # check 2D row sums
  expect_equal(rowSums(calcBasis2D), c(3552908, 3572102, 3591295), tolerance = 1e-5)
  # check symmetry
  expect_equal(calcBasis2D, t(calcBasis2D))
})

# rather plausibility checks
test_that("test MFPCA main function", {
  # see also MFPCA examples
  set.seed(1)
  sim <-  simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
                          M = 5, eFunType = "Poly", eValType = "linear", N = 100)
  uniExpansions <- list(list(type = "uFPCA"), list(type = "uFPCA"))
  
  # check errors
  expect_error(MFPCA(sim$simData[[1]], M = 5, uniExpansions = uniExpansions),
               "Parameter 'mFData' must be passed as a multiFunData object.")
  expect_error(MFPCA(sim$simData, M = "5", uniExpansions = uniExpansions),
               "Parameter 'M' must be passed as a number > 0.")
  expect_error(MFPCA(sim$simData, M = 1:5, uniExpansions = uniExpansions),
               "Parameter 'M' must be passed as a number > 0.")
  expect_error(MFPCA(sim$simData, M = -5, uniExpansions = uniExpansions),
               "Parameter 'M' must be passed as a number > 0.")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions[[1]]),
               "Parameter 'uniExpansions' must be passed as a list with the same length as 'mFData'.")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"))), 
               "Parameter 'uniExpansions' must be passed as a list with the same length as 'mFData'.")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions, weights = "1"),
               "Parameter 'weights' must be passed as a vector with the same length as 'mFData'.")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions, weights = 1),
               "Parameter 'weights' must be passed as a vector with the same length as 'mFData'.")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions, fit = "Yes"),
               "Parameter 'fit' must be passed as a logical.")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions, approx.eigen = "Yes"),
               "Parameter 'approx.eigen' must be passed as a logical.")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions, bootstrap = "Yes"),
               "Parameter 'bootstrap' must be passed as a logical.")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions, bootstrap = TRUE), 
               "Specify number of bootstrap iterations.")
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions, 
                     bootstrap = TRUE, nBootstrap = 10, bootstrapAlpha = -.1), 
               "Significance level for bootstrap confidence bands must be in (0,1).", fixed = TRUE) # fixed: do not interprete as reg. exp.
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions, 
                     bootstrap = TRUE, nBootstrap = 10, bootstrapAlpha = 1.5), 
               "Significance level for bootstrap confidence bands must be in (0,1).", fixed = TRUE) # fixed: do not interprete as reg. exp.
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions, 
                     bootstrap = TRUE, nBootstrap = 10, bootstrapStrat = 1:nObs(sim$simData)), 
               "bootstrapStrat must be either NULL or a factor.", fixed = TRUE)
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions, 
                     bootstrap = TRUE, nBootstrap = 10, bootstrapStrat = as.factor(1:5)), 
               "bootstrapStrat must have the same length as the number of observations in the mFData object.", fixed = TRUE)
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = uniExpansions, verbose = "Yes"),
               "Parameter 'verbose' must be passed as a logical.")
  
  # check warning
  expect_warning(MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA", npc = 2), list(type = "uFPCA", npc = 2))),
                 "Function MFPCA: total number of univariate basis functions is smaller than given M. M was set to 4")
  
  # check functionality
  expect_warning(uFPCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
                                                          list(type = "uFPCA")), approx.eigen = TRUE, fit = TRUE), 
                 "Calculating a large percentage of principal components, approximation may not be appropriate.
            'approx.eigen' set to FALSE.")
  splines <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "splines1D", k = 10),
                                                            list(type = "splines1D", k = 10)),
                   approx.eigen = TRUE, fit = TRUE)
  mixed <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
                                                            list(type = "splines1D", k = 10)),
                   approx.eigen = TRUE, fit = TRUE)
  
  # values
  expect_equal(length(uFPCA$values), length(splines$values))
  expect_equal(length(uFPCA$values), length(mixed$values))
  expect_equal(sum(uFPCA$values), sum(splines$values), tol = 5e-2)
  expect_equal(sum(uFPCA$values), sum(mixed$values), tol = 5e-2)
  expect_equal(uFPCA$values[1], 1.05174, tolerance = 1e-5)
  expect_equal(splines$values[1], 1.05266, tolerance = 1e-5)
  expect_equal(mixed$values[1], 1.05248, tolerance = 1e-5)
  
  # functions
  expect_equal(nObs(uFPCA$functions), nObs(splines$functions))
  expect_equal(nObs(uFPCA$functions), nObs(mixed$functions))
  expect_equal(norm(uFPCA$functions[[1]])[1], 0.57955, tolerance = 1e-5)
  expect_equal(norm(splines$functions[[1]])[1], 0.57957, tolerance = 1e-5)
  expect_equal(norm(mixed$functions[[1]])[1], 0.57970, tolerance = 1e-5)
  
  # fits
  expect_equal(sum(norm(uFPCA$fit - sim$simData)), 0.86272, tolerance = 1e-5)
  expect_equal(sum(norm(splines$fit - sim$simData)), 3.135592e-07, tolerance = 1e-5) ###
  expect_equal(sum(norm(mixed$fit - sim$simData)), 0.40767, tolerance = 1e-5)
  
  # mean function
  expect_equal(uFPCA$meanFunction[[1]], mixed$meanFunction[[1]])
  expect_equal(splines$meanFunction[[2]], mixed$meanFunction[[2]])
  
  # vectors
  expect_is(uFPCA$vectors, "matrix")
  expect_s4_class(splines$vectors, "dgeMatrix")
  expect_s4_class(mixed$vectors, "dgeMatrix")
  
  expect_equal(dim(uFPCA$vectors), c(6,5))
  expect_equal(dim(splines$vectors), c(20,5))
  expect_equal(dim(mixed$vectors), c(13,5))
  
  expect_equal(sum(abs(uFPCA$scores)), 291.50116, tolerance = 1e-5)
  expect_equal(sum(abs(splines$scores)), 292.12835, tolerance = 1e-5)
  expect_equal(sum(abs(mixed$scores)),  291.83424, tolerance = 1e-5)
  
  expect_equal(abs(uFPCA$scores[1,1]),  0.16779, tolerance = 1e-5)
  expect_equal(abs(splines$scores[1,1]), 0.17022, tolerance = 1e-5)
  expect_equal(abs(mixed$scores[1,1]), 0.17057, tolerance = 1e-5)
  
  # norm factors
  expect_length(uFPCA$normFactors, 5)
  expect_length(splines$normFactors, 5)
  expect_length(mixed$normFactors, 5)
  
  expect_equal(sum(uFPCA$normFactors), 7.18374, tolerance = 1e-5)
  expect_equal(sum(splines$normFactors), 7.16753, tolerance = 1e-5)
  expect_equal(sum(mixed$normFactors), 7.17516, tolerance = 1e-5)
  
  expect_equal(uFPCA$normFactors[1], 0.97509, tolerance = 1e-5)
  expect_equal(splines$normFactors[1], 0.97467, tolerance = 1e-5)
  expect_equal(mixed$normFactors[1], 0.97475, tolerance = 1e-5)
  
  ### Test bootstrap
  # suppress warnings in transition to new RNG, as proposed by CRAN maintainers
  suppressWarnings(RNGversion("3.5.0")) 
  set.seed(2)
  splinesBoot <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "splines1D", k = 10),
                                                                list(type = "splines1D", k = 10)),
                       approx.eigen = FALSE, bootstrap = TRUE, nBootstrap = 100, bootstrapAlpha = 0.05, verbose = FALSE)
  expect_equal(splinesBoot$CIvalues$alpha_0.05$upper[1], 1.39716, tolerance = 1e-5)
  expect_equal(sum(splinesBoot$CIvalues$alpha_0.05$upper), 3.51558, tolerance = 1e-5)
  expect_equal(splinesBoot$CIvalues$alpha_0.05$lower[1], 0.82537, tolerance = 1e-5)
  expect_equal(sum(splinesBoot$CIvalues$alpha_0.05$lower), 2.24919, tolerance = 1e-5)
  expect_equal(nObs(splinesBoot$CI$alpha_0.05$upper), 5)
  expect_equal(nObs(splinesBoot$CI$alpha_0.05$lower), 5)
  expect_equal(norm(splinesBoot$CI$alpha_0.05$upper - splinesBoot$CI$alpha_0.05$lower)[1], 0.980646, tolerance = 1e-5)
  expect_equal(sum(norm(splinesBoot$CI$alpha_0.05$upper - splinesBoot$CI$alpha_0.05$lower)), 23.6923, tolerance = 1e-5)
  
  uFPCABoot <- MFPCA(sim$simData, M = 3, uniExpansions = list(list(type = "uFPCA", npc = 3),
                                                                list(type = "uFPCA", npc = 3)),
                       approx.eigen = FALSE, bootstrap = TRUE, nBootstrap = 10, bootstrapAlpha = 0.1, verbose = FALSE)
  
  # check values
  expect_equal(uFPCABoot$CIvalues$alpha_0.1$upper[1], 1.14779, tolerance = 1e-5)
  expect_equal(sum(uFPCABoot$CIvalues$alpha_0.1$upper), 2.40996, tolerance = 1e-5)
  expect_equal(uFPCABoot$CIvalues$alpha_0.1$lower[1], 0.84857, tolerance = 1e-5)
  expect_equal(sum(uFPCABoot$CIvalues$alpha_0.1$lower), 1.8092398, tolerance = 1e-5)
  
  # check CI
  expect_equal(norm(uFPCABoot$CI$alpha_0.1$upper - uFPCABoot$CI$alpha_0.1$lower)[1], 0.52392, tolerance = 1e-5)
  expect_equal(sum(norm(uFPCABoot$CI$alpha_0.1$upper - uFPCABoot$CI$alpha_0.1$lower)), 8.98431, tolerance = 1e-5)
})

test_that("MFPCA calculation function", {
  set.seed(1)
  sim <-  simFunData(argvals = seq(0,1,0.01), M = 5, eFunType = "Poly", eValType = "linear", N = 100)
  uniB <- univDecomp("uFPCA", sim$simData, npc = 3)
  
  res <- MFPCA:::calcMFPCA(N = 100, p = 1, Bchol = diag(3), M = 3, type = "uFPCA", weights = 1, npc = 3, 
            argvals = list(sim$simData@argvals), uniBasis = list(uniB))
  
  # for p = 1, the univariate and multivariate results should coincide...
  expect_equal(abs(res$scores), abs(uniB$scores), tol = 2e-3, check.attributes = F)
  expect_equal(flipFuns(uniB$functions, res$functions[[1]]), uniB$functions, tol = 2e-3)
  expect_equal(res$values, c(1.30788, 0.79002, 0.42145), tolerance = 1e-5)
})