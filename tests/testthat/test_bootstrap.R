context("Testing bootstrap.R")

test_that("test stratSample", {
  f <- as.factor(c(1,1,1,2,4,4,4,4,6,6)) # contains also a stratum with only one entry -> important  to check sample!
  
  for(i in 1:10) # repeat sampling 10 times
    expect_equal(table(f, dnn = NULL), table(f[stratSample(f)], dnn = NULL)) # dnn: no names for variable that is 'tabled'
})
