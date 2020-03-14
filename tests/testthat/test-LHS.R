test_that("LHS", {
  a <- simple.LHS(n=10, d=3)
  expect_is(a, 'matrix')
  expect_equal(nrow(a), 10)
  expect_equal(ncol(a), 3)
  expect_true(is.LHS(a))

  a <- simple.LHS(n=10, d=3, centered=TRUE)
  expect_is(a, 'matrix')
  expect_equal(nrow(a), 10)
  expect_equal(ncol(a), 3)
  expect_true(is.LHS(a))



  g <- simple.grid(n=10, d=3)
  expect_true(is.matrix(g))
  expect_equal(dim(g), c(1000, 3))
})

# test_that("phi_p", {
#   phi_p(
#     matrix(runif(2*10), ncol=2)
#   )
# })
