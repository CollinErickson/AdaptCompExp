test_that("multiplication works", {
  a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=TestFunctions::gaussian1,obj="desirability",
                                    des_func=des_func_relmax, n0=20, take_until_maxpvar_below=.9,
                                    package="GauPro", design='sFFLHD', selection_method="max_des_red",
                                    alpha_des=1)
  a$run(2)
  expect_is(a, "adapt.concept2.sFFLHD.seq")
  expect_is(a, "R6")
})
