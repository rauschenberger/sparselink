

#--- logit and sigmoid transforms ---

testthat::test_that("sigmoid function returns values in unit interval",{
  eta <- stats::rnorm(n=100)
  prob <- sigmoid(eta)
  cond <- all(prob>=0 & prob<=1)
  testthat::expect_true(cond)
})

testthat::test_that("sigmoid followed by logit returns identity",{
  eta <- stats::rnorm(n=100)
  prob <- sigmoid(eta)
  eta_new <- logit(prob)
  cond <- all(abs(eta-eta_new)<=1e-06)
  testthat::expect_true(cond)
})

testthat::test_that("logit followed by sigmoid returns identity",{
  prob <- stats::runif(n=100)
  eta <- logit(prob)
  prob_new <- sigmoid(eta)
  cond <- all(abs(prob-prob_new)<=1e-06)
  testthat::expect_true(cond)
})

#--- mean and link functions ---

testthat::test_that("Gaussian mean and link functions return identity",{
  eta <- stats::runif(n=100)
  mu <- mean.function(eta=eta,family="gaussian")
  eta_new <- link.function(mu=mu,family="gaussian")
  cond <- all(eta==mu) & all(eta==eta_new)
  testthat::expect_true(cond)
})

testthat::test_that("binomial mean followed by link function returns identity",{
  eta <- stats::rnorm(n=100)
  mu <- mean.function(eta=eta,family="binomial")
  eta_new <- link.function(mu=mu,family="binomial")
  cond <- all(abs(eta-eta_new)<=1e-06)
  testthat::expect_true(cond)
})

testthat::test_that("binomial link followed by mean function returns identity",{
  mu <- stats::runif(n=100)
  eta <- link.function(mu=mu,family="binomial")
  mu_new <- mean.function(eta=eta,family="binomial")
  cond <- all(abs(mu-mu_new)<=1e-06)
  testthat::expect_true(cond)
})
