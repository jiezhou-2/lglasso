expected_matrix1=readRDS("phidata1.rds")
expected_matrix2=readRDS("phidata2.rds")
test_that("phifunction generate phi matrix correctly", {
  expect_equal(phifunction(t=1:15,tau = c(1.4,1)), expected_matrix1)
})


test_that("phifunction generate phi matrix correctly", {
  expect_equal(phifunction(t=1:15,tau = 1.4), expected_matrix2)
})
