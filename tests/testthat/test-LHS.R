test_that("multiplication works", {
  a <- simple.LHS(n=10, d=3)
  expect_is(a, 'matrix')
  expect_equal(nrow(a), 10)
  expect_equal(ncol(a), d)
  expect_true(is.LHS(a))
})
