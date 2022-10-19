context("options specification for models and tests")


#### try/catch errors about m (maximum breaks) for both model estimation and testing ####
## Model estimation
test_that("Maximum breaks cannot be larger than sample size allowed", {

  T = 50
  y = rnorm(T)
  y[c(25:50)] = y[c(25:50)] + 1
  data = data.frame(y)
  m = 10
  h = round(0.15*T)
  expect_warning(dosequa('y',data=data,m = m),paste('\nNot enough observations for',m+1,'segments with minimum length per segment =',
                                                   h,'.The total required observations for such',m,'breaks would be ',(m+1)*h,'>T=',T,'\n'))

})

test_that("Maximum breaks cannot be negative", {
  T = 50
  y = rnorm(T)
  y[c(25:50)] = y[c(25:50)] + 1
  data = data.frame(y)
  m = -1
  expect_warning(dosequa('y',data=data,m = m),'\\nMaximum number of breaks cannot be less than 1')
})

## Testing
test_that("SupF test maximum breaks cannot be 0", {
  T = 10
  y = rnorm(T)
  #y[c(51:100)] = y[c(51:100)] + 1
  data = data.frame(y)
  m = 0
  expect_warning(dotest('y',data=data,m = m),'The test is undefined for no breaks model')
})

test_that("Sequential test maximum breaks equals to Sup F test of 1 vs 0 for m=1", {
  T = 10
  y = rnorm(T)
  #y[c(51:100)] = y[c(51:100)] + 1
  data = data.frame(y)
  m = 1
  expect_message(doseqtests('y',data=data,m = m),'The test is exactly 0 versus 1 break, hence the sequential test is not repeated')
})

#### try/catch errors about eps1 (trimming level)
test_that("Invalid trimming level eps1 (too small or too large)", {
  T = 50
  y = rnorm(T)
  #y[c(51:100)] = y[c(51:100)] + 1
  data = data.frame(y)
  eps1 = 1
  expect_warning(doseqtests('y',data=data,eps1=eps1),'Invalid trimming level, set trimming level to 15%')
})


