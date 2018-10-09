#### Test functions for MFPCAfit class ####

context("Testing functions for MFPCAfit class")

pdf(file=NULL) # do not save plots

#### Calculate MFPCA once (see MFPCA help) ####
set.seed(1)
### simulate data (one-dimensional domains)
sim <-  funData::simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
                        M = 5, eFunType = "Poly", eValType = "linear", N = 100)
# MFPCA based on univariate FPCA
pca1D <- MFPCA::MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
                                                        list(type = "uFPCA")))

set.seed(2)
sim2 <-  funData::simMultiFunData(type = "weighted",
                        argvals = list(list(seq(0,1,0.01), seq(-1,1,0.02)), list(seq(-0.5,0.5,0.01))),
                        M = list(c(4,5), 20), eFunType = list(c("Fourier", "Fourier"), "Poly"),
                        eValType = "exponential", N = 150)

# MFPCA based on univariate spline expansions (for images) and univariate FPCA (for functions)
pca2D <- MFPCA::MFPCA(sim2$simData, M = 10,
             uniExpansions = list(list(type = "splines2D", k = c(10,12)),
                                  list(type = "uFPCA")))


#### print function ####
test_that("Test print function", {
  expect_error(MFPCA:::print.MFPCAfit(1:4),
               "Argument is not of class 'MFPCAfit'.")
  
  expect_known_output(print(pca1D), file = "outputs/print_MFPCAfit.txt")
  expect_known_output(print(pca2D), file = "outputs/print_MFPCAfit_2D.txt")
})

#### summary & print function ####
test_that("Test summary and associated print function", {
  expect_error(MFPCA:::summary.MFPCAfit(1:4),
               "Argument is not of class 'MFPCAfit'.")
  
  s1 <- summary(pca1D)
  s2 <- summary(pca2D)
  
  expect_s3_class(s1, "summary.MFPCAfit")
  expect_equal(dim(s1), c(3,5))
  expect_equal(attr(s1, "npc"), 5)
  expect_equal(attr(s1, "nel"), 2)
  
  expect_s3_class(s2, "summary.MFPCAfit")
  expect_equal(dim(s2), c(3,10))
  expect_equal(attr(s2, "npc"), 10)
  expect_equal(attr(s2, "nel"), 2)
  
  
  expect_known_output(print(s1), file = "outputs/summary_MFPCAfit.txt")
  expect_known_output(print(s2), file = "outputs/summary_MFPCAfit_2D.txt")
})

#### predict function ####
test_that("Test predict function", {
  expect_error(MFPCA:::predict.MFPCAfit(1:4),
               "Argument is not of class 'MFPCAfit'.")
  expect_error(MFPCA:::predict.MFPCAfit(pca1D, scores = "Test"),
               "Argument 'scores' must be a matrix with 5 columns.")
  expect_error(MFPCA:::predict.MFPCAfit(pca1D, scores = pca1D$scores[,1:2]),
               "Argument 'scores' must be a matrix with 5 columns.")
  
  pred <- predict(pca1D)
  
  expect_s4_class(pred, "multiFunData")
  expect_equal(dimSupp(pred), c(1,1))
  expect_equal(nObs(pred), 100)
  expect_equal(norm(pred)[[1]], 2.29476, tolerance = 1e-5)
  expect_equal(mean(norm(pred)),  2.88570, tolerance = 1e-5)
})


#### plot function ####
test_that("Test plot function", {
  expect_error(MFPCA:::plot.MFPCAfit(1:4),
               "Argument is not of class 'MFPCAfit'.")
  expect_error(MFPCA:::plot.MFPCAfit(pca1D, plotPCs = "1"),
               "Parameter 'plotPCs' must be a vector with values between 1 and 5.")
  expect_error(MFPCA:::plot.MFPCAfit(pca1D, plotPCs = 1:6),
               "Parameter 'plotPCs' must be a vector with values between 1 and 5.")
  expect_error(MFPCA:::plot.MFPCAfit(pca1D, plotPCs = -1:2),
               "Parameter 'plotPCs' must be a vector with values between 1 and 5.")
  expect_error(MFPCA:::plot.MFPCAfit(pca1D, plotPCs = 4:6),
               "Parameter 'plotPCs' must be a vector with values between 1 and 5.")
  expect_error(MFPCA:::plot.MFPCAfit(pca1D, stretchFactor = "1"),
               "Parameter 'stretchFactor' must be either NULL or a positive number.")
  expect_error(MFPCA:::plot.MFPCAfit(pca1D, stretchFactor = 1:2),
               "Parameter 'stretchFactor' must be either NULL or a positive number.")
  expect_error(MFPCA:::plot.MFPCAfit(pca1D, stretchFactor = -1),
               "Parameter 'stretchFactor' must be either NULL or a positive number.")
  expect_error(MFPCA:::plot.MFPCAfit(pca1D, combined = "Yes"),
               "Parameter 'combined' must be passed as a logical.")
  
  # create eigenfunction of 3-dimensional domain to generate error
  p <- pca1D
  p$functions[[1]] <- funData::tensorProduct(funData::eFun(argvals = c(1:5), M = 2, type = "Fourier"), funData::eFun(argvals = c(2:5), M = 2, type = "Fourier"), funData::eFun(argvals = c(3:5), M = 2, type = "Fourier"))
  expect_error(plot(p, plotPCs = 1:2), 
               "Cannot plot principal components having a 3- or higher dimensional domain.")
  expect_warning(plot(pca2D, combined = TRUE),
               "Cannot combine plots for two-dimensional elements. Will use separate plots (combined = FALSE)",
               fixed = TRUE)
  
  expect_null(plot(pca1D))
  expect_null(plot(pca2D))
})


#### scoreplot function ####
test_that("Test scoreplot function", {
  expect_error(MFPCA:::scoreplot.MFPCAfit(1:4),
               "Argument is not of class 'MFPCAfit'.")
  # create eigenfunction of 3-dimensional domain to generate error
  expect_error(MFPCA:::scoreplot.MFPCAfit(pca1D, choices = 1),
               "Parameter 'choices' must be a vector of length 2 with positive entries.")
  expect_error(MFPCA:::scoreplot.MFPCAfit(pca1D, choices = 1:3),
               "Parameter 'choices' must be a vector of length 2 with positive entries.")
  expect_error(MFPCA:::scoreplot.MFPCAfit(pca1D, choices = c(-1,1)),
               "Parameter 'choices' must be a vector of length 2 with positive entries.")
  expect_error(MFPCA:::scoreplot.MFPCAfit(pca1D, choices = c(1,6)),
                "Argument choices requires 6 scores, MFPCA object contains only 5.")
  expect_error(MFPCA:::scoreplot.MFPCAfit(pca1D, scale = "Yes"),
               "Parameter 'scale' must be passed as a logical.")

  expect_null(scoreplot(pca1D))
  expect_null(scoreplot(pca2D))
  })


#### screeplot function ####
test_that("Test screeplot function", {
  expect_error(MFPCA:::screeplot.MFPCAfit(1:4),
               "Argument is not of class 'MFPCAfit'.")
  expect_error(MFPCA:::screeplot.MFPCAfit(pca1D, npcs = "2"),
               "Parameter 'npcs' must be a number between 1 and 5.")
  expect_error(MFPCA:::screeplot.MFPCAfit(pca1D, npcs = 1:2),
               "Parameter 'npcs' must be a number between 1 and 5.")
  expect_error(MFPCA:::screeplot.MFPCAfit(pca1D, npcs = 0),
               "Parameter 'npcs' must be a number between 1 and 5.")
  expect_error(MFPCA:::screeplot.MFPCAfit(pca1D, npcs = 6),
               "Parameter 'npcs' must be a number between 1 and 5.")
  expect_error(MFPCA:::screeplot.MFPCAfit(pca1D, type = 1),
               "Parameter 'type' must be passed as a character.")
  expect_error(MFPCA:::screeplot.MFPCAfit(pca1D, type = "xxx"),
               "Type xxx not defined in screeplot.")
  expect_error(MFPCA:::screeplot.MFPCAfit(pca1D, main = 1),
               "Parameter 'main' must be either NULL or passed as a character.")
  
  expect_null(screeplot(pca1D))
  expect_null(screeplot(pca2D))
})

dev.off()