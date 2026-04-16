test_that("predict link matches fitted values on training data", {
  set.seed(10)
  n <- 100
  u <- rnorm(n)
  x <- pnorm((qnorm(0.03) + 0.1 * u + sqrt(0.02) * rnorm(n)) / sqrt(0.98))
  d <- data.frame(y = x, u = u)
  fit <- vasicekfit(y ~ u, data = d)

  expect_equal(predict(fit, type = "link"), fitted(fit))
})

test_that("predict response is pnorm of link", {
  set.seed(11)
  n <- 100
  u <- rnorm(n)
  x <- pnorm((qnorm(0.03) + 0.1 * u + sqrt(0.02) * rnorm(n)) / sqrt(0.98))
  d <- data.frame(y = x, u = u)
  fit <- vasicekfit(y ~ u, data = d)

  link <- predict(fit, type = "link")
  resp <- predict(fit, type = "response")
  expect_equal(resp, pnorm(link))
})

test_that("predict with newdata works", {
  set.seed(12)
  n <- 100
  u <- rnorm(n)
  x <- pnorm((qnorm(0.03) + 0.1 * u + sqrt(0.02) * rnorm(n)) / sqrt(0.98))
  d <- data.frame(y = x, u = u)
  fit <- vasicekfit(y ~ u, data = d)

  nd <- data.frame(u = c(-1, 0, 1))
  pred <- predict(fit, newdata = nd, type = "response")
  expect_length(pred, 3)
  expect_true(all(pred > 0 & pred < 1))

  # Higher u (positive kappa) should give higher PD
  expect_true(pred[3] > pred[1])
})

test_that("predict with scalar alpha returns vector", {
  set.seed(13)
  n <- 100
  u <- rnorm(n)
  x <- pnorm((qnorm(0.03) + 0.1 * u + sqrt(0.02) * rnorm(n)) / sqrt(0.98))
  d <- data.frame(y = x, u = u)
  fit <- vasicekfit(y ~ u, data = d)

  nd <- data.frame(u = c(0, 1))
  pred <- predict(fit, newdata = nd, alpha = 0.99, type = "response")
  expect_length(pred, 2)
  expect_true(all(pred > 0 & pred < 1))
})

test_that("predict with vector alpha returns matrix", {
  set.seed(14)
  n <- 100
  u <- rnorm(n)
  x <- pnorm((qnorm(0.03) + 0.1 * u + sqrt(0.02) * rnorm(n)) / sqrt(0.98))
  d <- data.frame(y = x, u = u)
  fit <- vasicekfit(y ~ u, data = d)

  nd <- data.frame(u = c(0, 1, 2))
  pred <- predict(fit, newdata = nd, alpha = c(0.95, 0.99), type = "response")
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(3, 2))

  # Higher alpha should give higher loss rate
  expect_true(all(pred[, 2] > pred[, 1]))
})

test_that("predict alpha matches qvasicek for intercept-only model", {
  set.seed(15)
  n <- 5000
  p_true <- 0.04
  rho_true <- 0.08
  z <- rnorm(n)
  x <- pnorm((qnorm(p_true) + sqrt(rho_true) * z) / sqrt(1 - rho_true))
  d <- data.frame(y = x)
  fit <- vasicekfit(y ~ 1, data = d)

  # predict at alpha=0.99 with no covariates should approximately
  # match qvasicek with fitted p and rho
  pred99 <- predict(fit, newdata = data.frame(row.names = 1),
                    alpha = 0.99, type = "response")
  q99 <- qvasicek(0.99, p = fit$p, rho = fit$rho)
  expect_equal(pred99, q99, tolerance = 1e-10)
})
