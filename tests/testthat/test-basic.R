test_that("adapt works", {
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=TestFunctions::gaussian1,obj="desirability",
                                    des_func=des_func_relmax, n0=20, take_until_maxpvar_below=.9,
                                    package="GauPro", design='sFFLHD', selection_method="max_des_red",
                                    alpha_des=1)
  a$run(2)
  expect_is(a, "adapt.concept2.sFFLHD.seq")
  expect_is(a, "R6")
  # Check other functions that are easy to run
  expect_error(a$plot_des_v_acc(1,1), NA)
  expect_error(a$plot_iwe(), NA)
})

test_that("compare adapt works", {

  ca1 <- compare.adaptR6$new(func=TestFunctions::gaussian1, D=2, L=3,
                             batches=2, reps=2,
                             n0=6, obj="desirability",
                             selection_method=c('max_des', 'SMED'),
                             des_func=c('des_func_relmax', 'des_func_relmax')
  )

  expect_is(ca1, "compare.adaptR6")
  expect_is(ca1, "R6")
})
