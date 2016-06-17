context("Testing functions in MFCPA_calculation.R")

# rather plausibility checks
test_that("test calcBasisIntegrals", {
  basis1D <- univExpansion("splines1D", scores = matrix(1:30, nrow = 3), argvals = seq(0, 1, 0.01), 
                           functions = NULL, params = list(k = 10, m = 2, bs = "ps"))
  calcBasis1D <- MFPCA:::calcBasisIntegrals(basis1D@X, dimSupp = dimSupp(basis1D), argvals = basis1D@argvals)
  
  basis2D <- univExpansion("splines2D", scores = matrix(1:360, nrow = 3), argvals = list(seq(0, 1, 0.01), seq(-2, 2, 0.02)),
                                       functions = NULL, params = list(k = c(10,12), m = c(2,3), bs = "ps"))
  calcBasis2D <- MFPCA:::calcBasisIntegrals(basis2D@X, dimSupp = dimSupp(basis2D), argvals = basis2D@argvals)
  
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
               "Significance level for bootstrap confidence bands must be in (0,1).", fixed = TRUE) # fixed: do not interprete as reg. exp.
  expect_error(MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"), list(type = "uFPCA")), 
                     bootstrap = TRUE, nBootstrap = 10, bootstrapAlpha = 1.5), 
               "Significance level for bootstrap confidence bands must be in (0,1).", fixed = TRUE) # fixed: do not interprete as reg. exp.
  
  # check functionality
  expect_warning(uFPCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
                                                          list(type = "uFPCA")), approx.eigen = TRUE), 
                 "Calculating a large percentage of principal components, approximation may not be appropriate.
            'approx.eigen' set to FALSE.")
  splines <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "splines1D", k = 10),
                                                            list(type = "splines1D", k = 10)),
                   approx.eigen = TRUE)
  
  # values
  expect_equal(length(uFPCA$values), length(splines$values))
  expect_equal(sum(uFPCA$values), sum(uFPCA$values))
  expect_equal(uFPCA$values[1], 1.05174404)
  expect_equal(splines$values[1], 1.05266096)
  
  # functions
  expect_equal(nObs(uFPCA$functions), nObs(splines$functions))
  expect_equal(norm(uFPCA$functions[[1]])[1], 0.579550598)
  expect_equal(norm(splines$functions[[1]])[1], 0.57956679)
})

test_that("MFPCA calculation function", {
  set.seed(1)
  sim <-  simFunData(argvals = seq(0,1,0.01), M = 5, eFunType = "Poly", eValType = "linear", N = 100)
  uniB <- univDecomp("uFPCA", sim$simData, npc = 3)
  
  res <- MFPCA:::calcMFPCA(N = 100, p = 1, Bchol = diag(3), M = 3, type = "uFPCA", weights = 1, npc = 3, 
            argvals = sim$simData@argvals, uniBasis = list(uniB))
  
  # for p = 1, the univariate and multivariate results should coincide...
  expect_equal(abs(res$scores), abs(uniB$scores), tol = 2e-3, check.attributes = F)
  expect_equal(flipFuns(uniB$functions, res$functions[[1]]), uniB$functions, tol = 2e-3)
  expect_equal(res$values, c(1.307880364, 0.790022175, 0.421445368))
    
            
          })
