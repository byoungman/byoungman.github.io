# functions and data for lectures

hilbert <- function(n) {
  # Function to evaluate n by n Hilbert matrix.
  # n is an integer
  # returns n by n matrix.
  ind <- 1:n
  1 / (outer(ind, ind, FUN = '+') - 1)
}

# operating temperatures
temp <- c(35.3, 29.7, 30.8, 58.8, 61.4, 71.3, 74.4, 76.7, 70.7, 57.5,
          46.4, 28.9, 28.1, 39.1, 46.8, 48.5, 59.3, 70, 70, 74.5, 72.1,
          58.1, 44.6, 33.4, 28.6)
# number of operational days
days <- c(20, 20, 23, 20, 21, 22, 11, 23, 21, 20, 20, 21, 21, 19, 23,
          20, 22, 22, 11, 23, 20, 21, 20, 20, 22)
# output from factory
output <- c(10.98, 11.13, 12.51, 8.4, 9.27, 8.73, 6.36, 8.5, 7.82, 9.14,
            8.24, 12.19, 11.88, 9.57, 10.94, 9.58, 10.09, 8.11, 6.83, 8.88,
            7.68, 8.47, 8.86, 10.36, 11.08)
# data frame
prod <- data.frame(temp = temp, days = days, output = output)

y <- c(.7, 1.3, 2.6)
mu <- 1:3
Sigma <- matrix(c(4, 2, 1, 2, 3, 2, 1, 2, 2), 3, 3)

dmvn3 <- function(y, mu, Sigma, log = TRUE) {
  # Function to evaluate MVN pdf
  # y, mu vectors
  # Sigma matrix
  # log is logical indicating whether logarithm
  # returns a scalar
  p <- length(y)
  res <- y - mu
  L <- t(chol(Sigma))
  out <- - sum(log(diag(L))) - 0.5 * p * log(2 * pi) -
    0.5 * sum(forwardsolve(L, res)^2)
  if (!log)
    out <- exp(out)
  out
}

fd <- function(x, f, delta = 1e-6, ...) {
  # Function to evaluate derivative w.r.t. vector by finite-differencing
  # x is a p-vector
  # fn is the function for which the derivative is being calculated
  # delta is the finite-differencing step, which defaults to 10ˆ{-6}
  # returns a vector of length x
  f0 <- f(x, ...)
  p <- length(x)
  f1 <- numeric(p)
  for (i in 1:p) {
    ei <- replace(numeric(p), i, 1)
    f1[i] <- f(x + delta * ei, ...)
  }
}

# Weibull log-likelihood

# Week 10 lecture 1

y0 <- c(3.52, 1.95, 0.62, 0.02, 5.13, 0.02, 0.01, 0.34, 0.43, 15.5, 
        4.99, 6.01, 0.28, 1.83, 0.14, 0.97, 0.22, 0.02, 1.87, 0.13, 0.01, 
        4.81, 0.37, 8.61, 3.48, 1.81, 37.21, 1.85, 0.04, 2.32, 1.06)

weib_d1 <- function(pars, y, mult = 1) {
  # Function to evaluate first derivative of Weibull log-likelihood
  # pars is a vector
  # y can be scalar or vector
  # mult is a scalar defaulting to 1; so -1 returns neg. gradient
  # returns a vector
  n <- length(y)
  z1 <- y / pars[1]
  z2 <- z1^pars[2]
  out <- numeric(2)
  out[1] <- (sum(z2) - n) * pars[2] / pars[1] # derivative w.r.t. lambda
  out[2] <- n * (1 / pars[2] - log(pars[1])) + 
    sum(log(y)) - sum(z2 * log(z1)) # w.r.t k
  mult * out
}

weib_d2 <- function(pars, y, mult = 1) {
  # Function to evaluate second derivative of Weibull log-likelihood
  # pars is a vector
  # y can be scalar or vector
  # mult is a scalar defaulting to 1; so -1 returns neg. Hessian
  # returns a matrix
  n <- length(y)
  z1 <- y / pars[1]
  z2 <- z1^pars[2]
  z3 <- sum(z2)
  z4 <- log(z1)
  out <- matrix(0, 2, 2)
  out[1, 1] <- (pars[2] / pars[1]^2) * (n - (1 + pars[2]) * z3) # w.r.t. (lambda^2)
  out[1, 2] <- out[2, 1] <- (1 / pars[1]) * ((z3 - n) + 
                                               pars[2] * sum(z2 * z4)) # w.r.t. (lambda, k)
  out[2, 2] <- -n/pars[2]^2 - sum(z2 * z4^2) # w.r.t. k^2
  mult * out
}

weib_d0 <- function(pars, y, mult = 1) {
  # Function to evaluate Weibull log-likelihood
  # pars is a vector
  # y can be scalar or vector
  # mult is a scalar defaulting to 1; so -1 returns neg. log likelihood
  # returns a scalar
  n <- length(y)
  if (min(pars) <= 0) {
    out <- -1e8
  } else {
    out <- n * (log(pars[2]) - pars[2] * log(pars[1])) + 
      (pars[2] - 1) * sum(log(y)) - sum((y / pars[1])^pars[2])
  }
  mult * out
}
