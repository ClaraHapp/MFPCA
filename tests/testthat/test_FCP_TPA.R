context("Testing FCP_TPA functionality")

test_that("normVec", {
  expect_equal(MFPCA:::normVec(c(1,0,0,0)), 1)
  
  x <- runif(4)
  expect_equal(MFPCA:::normVec(x), sqrt(sum(x^2)))
})

test_that("FPC_TPA",{
  u <- 1:30
  v <- sin(seq(-pi, pi, 0.03))
  w <- exp(seq(-0.5, 1, 0.01)) 
  
  X <- u %o% v %o% w
  penMat <- list(v = crossprod(MFPCA:::makeDiffOp(degree = 2, dim = 210)), w = crossprod(MFPCA:::makeDiffOp(degree = 2, dim = 151)))
  alphaRange <- list(v = c(1e-4, 1e4), w = c(1e-4, 1e4))
  
  # check errors
  expect_error(FCP_TPA(X = 1, K = 1, penMat = penMat, alphaRange = alphaRange),
               "Parameter 'X' must be passed as an array.")
  expect_error(FCP_TPA(X = X[1,,], K = 1, penMat = penMat, alphaRange = alphaRange),
               "Parameter 'X' must have three dimensions.")
  expect_error(FCP_TPA(X = X, K = "Test", penMat = penMat, alphaRange = alphaRange),
               "Parameter 'K' must be passed as a number.")
  expect_error(FCP_TPA(X = X, K = 1:5, penMat = penMat, alphaRange = alphaRange),
               "Parameter 'K' must be passed as a number.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = c(u = 1), alphaRange = alphaRange),
               "Parameter 'penMat' must be passed as a list.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = list(u = 1, w = penMat$w), alphaRange = alphaRange),
               "Function FCP_TPA: penMat must be a list of matrices with entries 'v' and 'w'.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = list(v = diag(1), w = penMat$w), alphaRange = alphaRange),
               "Function FCP_TPA: the penalization matrix for dimension v must be 210 x 210.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = list(v = matrix(1:210^2, nrow = 210), w = penMat$w), alphaRange = alphaRange),
               "Function FCP_TPA: the penalization matrix for dimension v must be symmetric.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = list(v = penMat$v, w = diag(1)), alphaRange = alphaRange),
               "Function FCP_TPA: the penalization matrix for dimension w must be 151 x 151.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = list(v = penMat$v, w = matrix(1:151^2, nrow = 151)), alphaRange = alphaRange),
               "Function FCP_TPA: the penalization matrix for dimension w must be symmetric.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = 1),
               "Parameter 'alphaRange' must be passed as a list.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = list(u = 1, v = 2)),
               "Function FCP_TPA: alphaRange must be a list of vectors with entries 'v' and 'w'.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = list(v = 1, w = c(0,1))),
               "Function FCP_TPA: alphaRange$v must be a vector of length 2.", fixed = TRUE) # do not interprete $ as special character
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = list(v = c(-1,1), w = c(0,1))),
               "Function FCP_TPA: Values for alphaV must not be negative.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = list(v = c(0,1), w = 1)),
               "Function FCP_TPA: alphaRange$w must be a vector of length 2.", fixed = TRUE)
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = list(v = c(0,1), w = c(-1,1))),
               "Function FCP_TPA: Values for alphaW must not be negative.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = alphaRange, verbose = "Yes"),
               "Parameter 'verbose' must be passed as a logical.")
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = alphaRange, tol = "No"),
               "Parameter 'tol' must be passed as a number")
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = alphaRange, tol = 1:2),
               "Parameter 'tol' must be passed as a number")
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = alphaRange, maxIter = "100"),
               "Parameter 'maxIter' must be passed as a number")
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = alphaRange, maxIter = 10:12),
               "Parameter 'maxIter' must be passed as a number")
  expect_error(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = alphaRange, adaptTol = "Yes"),
               "Parameter 'adaptTol' must be passed as a logical.")
  
  # check warning
  expect_warning(FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = alphaRange, maxIter = 2, tol = 1e-18),
                 "FCP-TPA algorithm did not converge; iteration 1 stopped.")
  
  # check results 
  set.seed(1)
  res <- FCP_TPA(X = X, K = 1, penMat = penMat, alphaRange = alphaRange)
  
  expect_equal(res$d, MFPCA:::normVec(u) * MFPCA:::normVec(v)  * MFPCA:::normVec(w))
  expect_equal(abs(res$U[,1]), abs(u/MFPCA:::normVec(u)), tolerance = 1e-9)
  expect_equal(abs(res$V[,1]), abs(v/MFPCA:::normVec(v)), tolerance = 1e-6)
  expect_equal(abs(res$W[,1]), abs(w/MFPCA:::normVec(w)), tolerance = 1e-5)
})


test_that("findAlphaOpt", {
  data <- (1:3) %o% 1:2 %o% 1:5
  u <- c(1,0,0)
  w <- c(0,0,0,0,1)
  GammaV <- diag(2)
  alphaW <- 0
  OmegaW <- diag(5)
  lambdaV <- c(1,2)
  
  # no minimum here... (rather plausibility check...)
  expect_equal(MFPCA:::findAlphaVopt(alphaRange = c(1e-4, 1e4), data = data, u = u, w = w, 
                             alphaW = alphaW, OmegaW = OmegaW, GammaV = GammaV, lambdaV = lambdaV),
               1e4, tol = .Machine$double.eps^0.25)
  expect_equal(MFPCA:::findAlphaWopt(alphaRange = c(1e-4, 1e4), data = data, u = u, v = c(1,0), 
                                     alphaV = alphaW, OmegaV = diag(2), GammaW = diag(5), lambdaW = 1:5),
               1e4, tol = .Machine$double.eps^0.25)
})

test_that("gcv", {

  # rather plausibility check...
  expect_equal(MFPCA:::gcv(alpha = 1, n=  5, z = c(1,0,0,0,0), eta = 0.5, lambda = 5:1),
               0.22989, tolerance = 1e-5)
})


test_that("makeDiffOp", {
  expect_equal(MFPCA:::makeDiffOp(degree = 1, dim = 3),
               rbind(c(-1,1,0), c(0,-1,1)))
  expect_equal(MFPCA:::makeDiffOp(degree = 2, dim = 4),
               rbind(c(1,-2,1, 0), c(0, 1, -2, 1)))
})