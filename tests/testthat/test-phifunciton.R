
homo1data=readRDS("homoOneStageData.rds")
homo1result=readRDS("homoOneStageResult.rds")
homo2data=readRDS("homoTwoStageData.rds")
homo2result=readRDS("homoTwoStageResult.rds")
heter1data=readRDS("heterOneStageData.rds")
heter1result=readRDS("heterOneStageResult.rds")
heter2data=readRDS("heterTwoStageData.rds")
heter2result=readRDS("heterTwoStageResult.rds")
set.seed(1)
test_that("lglasso", {
  expect_equal(lglasso(data=homo1data[[1]],lambda = 0.01,trace=T), homo1result)
})
set.seed(1)
test_that("lglasso", {
  expect_equal(lglasso(data=homo2data[[1]],lambda = c(0.1,0.1),group = homo2data[[2]],trace=T),
               homo2result)
})

set.seed(1)
test_that("lglasso", {

  expect_equal(lglasso(data=heter1data[[1]],lambda = 0.01,random=T,trace=T,N=100)$wi,
               heter1result$wi)
})


set.seed(1)
test_that("lglasso", {
  expect_equal(lglasso(data=heter2data[[1]],lambda = c(0.1,0.1),random=T,trace=T,
                       group = heter2data[[2]],N=100)$wi,
               heter2result$wi)
})







