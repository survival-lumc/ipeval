test_that("rhs_is_one() works", {
  expect_equal(rhs_is_one(~ 1), TRUE)
  expect_equal(rhs_is_one(X ~ 1), TRUE)
  expect_equal(rhs_is_one(1 ~ X + 1), FALSE)
  expect_equal(rhs_is_one( ~ X + 1), FALSE)
  expect_equal(rhs_is_one(1 ~ X*Y), FALSE)
  expect_equal(rhs_is_one(1 ~ X*Y + 1), FALSE)
})



test_that("make_named_list works", {
  myfunct <- function(object) {
    # ip_score uses this fct in this way
    make_named_list(object, substitute(object))
  }

  data <- data.frame(x = c(1,0,0), y = c(0,1,0))
  model1 <- glm(y ~ x, family = "binomial", data =data)
  model2 <- glm(x ~ y, family = "binomial", data =data)

  unnamed_list_of_models <- list(model1, model2)
  partial_named_list_of_models <- list("mod1" = model1, model2)
  named_list_of_models <- list("mod1" = model1, "mod2" = model2)

  # passing an object which has previously been defined as a list
  expect_identical(
    myfunct(unnamed_list_of_models),
    list("model.1" = model1, "model.2" = model2)
  )
  expect_identical(
    myfunct(partial_named_list_of_models),
    list("mod1" = model1, "model.2" = model2)
  )
  expect_identical(
    myfunct(named_list_of_models),
    named_list_of_models
  )

  # passing object = list(..)
  expect_identical(
    myfunct(list(model1, model2)),
    list("model1" = model1, "model2" = model2)
  )
  expect_identical(
    myfunct(list("mymod1" = model1, model2)),
    list("mymod1" = model1, "model2" = model2)
  )

  # passing object = model1
  expect_identical(myfunct(model1), list("model1" = model1))

})
