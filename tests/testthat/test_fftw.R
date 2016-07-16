context("Testing fftw functionality")

#### DCT ####

test_that("Test fftw: DCT", {
  set.seed(4)
  image2D <- outer(seq(0,1,0.01), sin(seq(-1,1,0.02))) + matrix(rnorm(101^2, sd = 0.2), nrow= 101)
  image3D <- aperm(sapply(seq(0,10,0.5), 
                          function(x){x*outer(seq(0,1,0.01), sin(seq(-1,1,0.02))) + matrix(rnorm(101^2, sd = 0.2), nrow= 101)}, 
                          simplify = "array"), c(3,1,2))
  
  expect_error(MFPCA:::dct2D(image3D, qThresh = 0.9), "dct2D can handle only 2D images")
  expect_error(MFPCA:::dct3D(image2D, qThresh = 0.9), "dct3D can handle only 3D images")
  
  fftw2D <- try(MFPCA:::dct2D(image2D, qThresh = 0.95), silent = TRUE)
  if(class(fftw2D) == "try-error")
    expect_error(stop(fftw2D), "dctBasis2D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.")
  else
  {
    expect_equal(length(fftw2D$ind), 510)
    expect_equal(fftw2D$ind[1:3], c(2,92,102))
    expect_equal(mean(fftw2D$val), -0.00049897171)
    expect_equal(fftw2D$val[1], 0.018282019)
  }
  
  
  fftw3D <- try(MFPCA:::dct3D(image3D, qThresh = 0.95), silent = TRUE)
  if(class(fftw3D) == "try-error")
    expect_error(stop(fftw3D), "dctBasis3D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.")
  else
  {
    expect_equal(length(fftw3D$ind), 10711)
    expect_equal(fftw3D$ind[1:3], c(102, 103, 110))
    expect_equal(mean(fftw3D$val), 0.00027522909)
    expect_equal(fftw3D$val[1], -0.32200299)
  }
})


test_that("test univariate DCT 2D", {
  set.seed(1)
  x1 <- seq(0,1,length.out=50)
  x2 <- seq(-1,1, length.out=75)
  f2 <- funData(argvals = list(x1, x2),
                X = aperm(replicate(10, outer(x1, cos(pi*x2))+matrix(rnorm(50*75, sd = 0.1), nrow = 50)), c(3,1,2)))
  
  dct2D <- try(MFPCA:::dctBasis2D(f2, qThresh = 0.95), silent = TRUE)
  if(class(dct2D) == "try-error")
    expect_error(stop(dct2D), "dctBasis2D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.")
  else
  {
    expect_equal(dim(dct2D$scores), c(10, 3750))
    expect_equal(length(dct2D$scores@i),  1880) 
    expect_equal(dct2D$scores@x[1],  -0.0196469613) 
    expect_equal(dim(dct2D$B), c(3750, 3750))
    expect_equal(dct2D$B@x[1], 0.202642367)
    expect_equal(var(diff(dct2D$B@x)), 0)
    expect_false(dct2D$ortho)  
    expect_null(dct2D$functions)
    
    # wrapper function
    decompDCT2D <- MFPCA:::univDecomp(type = "DCT2D", funDataObject = f2, qThresh = 0.95)
    expect_equal(decompDCT2D, dct2D)
  }
  
  
})


test_that("test univariate DCT 3D", {
  set.seed(3)
  x1 <- seq(0, 1, length.out = 40)
  x2 <- seq(-1, 1, length.out = 30)
  x3 <- seq(0, 0.5, length.out = 20)
  f3 <- funData(argvals = list(x1, x2, x3), X = replicate(20, array(rnorm(10*40*30), dim = c(10, 40, 30))))
  
  dct3D <- try(MFPCA:::dctBasis3D(f3, qThresh = 0.95), silent = TRUE)
  if(class(dct3D) == "try-error")
    expect_error(stop(dct3D), "dctBasis3D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.")
  else
  {
    expect_equal(dim(dct3D$scores), c(10, 23998))
    expect_equal(length(dct3D$scores@i),  12000) 
    expect_equal(dct3D$scores@x[1],  -0.071277532) 
    expect_equal(dim(dct3D$B), c(23998, 23998))
    expect_equal(dct3D$B@x[1], 0.032251534)
    expect_equal(var(diff(dct3D$B@x)), 0)
    expect_false(dct3D$ortho)  
    expect_null(dct3D$functions)
    
    # wrapper function
    decompDCT3D <- MFPCA:::univDecomp(type = "DCT3D", funDataObject = f3, qThresh = 0.95)
    expect_equal(decompDCT3D, dct3D)
  }
})


##### Inverse DCT #####

test_that("Test fftw: IDCT", {
  set.seed(4)
  scores <- rnorm(25, sd = 25:1/25)
  
  expect_error(MFPCA:::idct2D(scores = scores, ind = 1:25, dim = 50), "Function idct2D can handle only 2D images.")
  expect_error(MFPCA:::idct2D(scores = scores, ind = 0:24, dim = c(10, 20)), "Indices must be positive.")
  expect_error(MFPCA:::idct2D(scores = scores, ind = 200+0:24, dim = c(10, 20)), "Index exceeds image dimensions.")
  
  expect_error(MFPCA:::idct3D(scores = scores, ind = 1:25, dim = 50), "Function idct3D can handle only 3D images.")
  expect_error(MFPCA:::idct3D(scores = scores, ind = 0:24, dim = c(10, 20, 30)), "Indices must be positive.")
  expect_error(MFPCA:::idct3D(scores = scores, ind = 2000+0:24, dim = c(10, 20, 10)), "Index exceeds image dimensions.")
  
  idct2D <- try(MFPCA:::idct2D(scores = scores, ind = sample(200, 25), dim = c(10, 20)), silent = TRUE)
  if(class(idct2D) == "try-error")
    expect_error(stop(idct2D), "dctBasis2D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.")
  else
  {
    expect_equal(dim(idct2D), c(10, 20))
    expect_equal(mean(idct2D), 0)
    expect_equal(idct2D[1,1], 1.08970077)
  }
  
  idct3D <- try(MFPCA:::idct3D(scores = scores, ind = sample(2000, 25), dim = c(10, 20, 10)), silent = TRUE)
  if(class(idct3D) == "try-error")
    expect_error(stop(idct3D), "dctBasis3D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.")
  else 
  {
    expect_equal(dim(idct3D), c(10, 20, 10))
    expect_equal(mean(idct3D), 0)
    expect_equal(idct3D[1,1,1], 0.877946026)
  }
})


test_that("test univariate IDCT 2D", {
  set.seed(2)
  scores <- sapply(25:1, function(x){rnorm(20, sd = x/25)})
  argvals <- list(seq(0, 1, 0.01), seq(-1, 1, 0.02))
  
  dct2D <- try(MFPCA:::dctFunction2D(scores = scores, argvals = argvals), silent = TRUE)
  if(class(dct2D) == "try-error")
    expect_error(stop(dct2D), "dctBasis2D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.\"")
  else
  {
    expect_equal(nObs(dct2D), 20)
    expect_equal(nObsPoints(dct2D),  c(101,101))
    expect_equal(mean(norm(dct2D)),  2.06938459) 
    expect_equal(norm(dct2D)[1], 2.45252941)
    
    # wrapper function
    expandDCT2D <- MFPCA:::univExpansion(type = "DCT2D", scores = scores, argvals = argvals, functions = NULL)
    expect_equal(expandDCT2D, dct2D)
  }
})


test_that("test univariate IDCT 3D", {
  set.seed(3)
  scores <- sapply(60:1, function(x){rnorm(20, sd = x/25)})
  argvals <- list(seq(0, 1, 0.01), seq(-1, 1, 0.02), seq(-0.5, 0.5, 0.05))
  
  dct3D <- try(MFPCA:::dctFunction3D(scores = scores, argvals = argvals), silent = TRUE)
  if(class(dct3D) == "try-error")
    expect_error(stop(dct3D), "dctBasis3D requires C-library fftw3 to be installed. Check http://www.fftw.org/ for more information.")
  else
  {
    expect_equal(nObs(dct3D), 20)
    expect_equal(nObsPoints(dct3D),  c(101, 101, 21))
    expect_equal(mean(norm(dct3D)),  7.57861399) 
    expect_equal(norm(extractObs(dct3D, obs = 1)), 7.64625103)
    
    # wrapper function
    expandDCT3D <- MFPCA:::univExpansion(type = "DCT3D", scores = scores, argvals = argvals, functions = NULL)
    expect_equal(expandDCT3D, dct3D)
  }
})
