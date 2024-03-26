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
        
lk0 <- c(1.6, .6)

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

line_search <- function(theta, p, f, alpha0 = 1, rho = .5, ...) {
  best <- f(theta, ...)
  cond <- TRUE
  while (cond & alpha0 > .Machine$double.eps) {
    prop <- f(theta + alpha0 * p, ...)
    cond <- prop >= best
    if (!cond)
      best <- prop
    alpha0 <- alpha0 * rho
  }
  alpha <- alpha0 / rho
  alpha
}

iH1 <- function(x0, x1, g0, g1, iH0) {
  # Function to update Hessian matrix
  # x0 and x1 are p-vectors of second to last and last estimates, respectively
  # g0 and g1 are p-vectors of second to last and last gradients, respectively
  # iH0 is previous estimate of p x p Hessian matrix
  # returns a p x p matrix
  s0 <- x1 - x0
  y0 <- g1 - g0
  denom <- sum(y0 * s0)
  I <- diag(rep(1, 2))
  pre <- I - tcrossprod(s0, y0) / denom
  post <- I - tcrossprod(y0, s0) / denom
  last <- tcrossprod(s0) / denom
  pre %*% iH0 %*% post + last
}

# simulated annealing

update_T <- function(i, t0 = 10, t1 = 10) {
  # Function to update simulated annealing temperature
  # i is an integer giving the current iteration
  # t0 is a scalar giving the initial temperature
  # t1 is a integer giving how many iterations of each temperature to use
  # returns a scalar
  t0 / log(((i - 1) %/% t1) * t1 + exp(1))
}

q_fn <- function(x) {
  # Function to generate Gaussian proposals with standard deviation 0.1
  # x is the Gaussian mean as either a scalar or vector
  # returns a scalar or vector, as x
  rnorm(length(x), x, .1)
}

sa <- function(p0, h, N, q, T1, ...) {
  # Function to perform simulated annealing
  # p0 p-vector of initial parameters
  # h() function to be minimised
  # N number of iterations
  # q proposal function
  # T1 initial temperature
  # ... arguments to pass to h()
  # returns p x N matrix of parameter estimates at each iteration
  out <- matrix(0, N, length(p0)) # matrix to store estimates at each iteration
  out[1, ] <- p0 # fill first row with initial parameter estimates
  for (i in 2:N) { # N iterations
    T <- update_T(i, T1) # update temperature
    U <- runif(1) # generate U
    out[i, ] <- out[i - 1,] # carry over last parameter estimate, by default
    proposal <- q(out[i - 1,]) # generate proposal
    if (min(proposal) >= 0) { # ensure proposal valid
      h0 <- h(out[i - 1, ], ...) # evaluate h for current theta
      h1 <- h(proposal, ...) # evaluate h for proposed theta
      alpha <- min(exp(- (h1 - h0) / T), 1) # calculate M-H ratio
      if (alpha >= U) # accept if ratio sufficiently high
        out[i, ] <- proposal # swap last with proposal
    }
  }
  out # return all parameter estimates
}

# T_vals <- c(10, 1, .1)
# # loop over values, and plot
# for (j in 1:length(T_vals)) {
#   T1 <- T_vals[j]
#   sa_result <- sa(c(1.6, 0.6), weib_d0, 1e3, q_fn, T1, y = y0, mult = -1)
#   if (j == 1) {
#     plot(sa_result, col = j, pch = 20, xlab = 'lambda', ylab = 'k')
#   } else {
#     points(sa_result, col = j, pch = 20)
#   }
# }
# legend('bottomright', pch = 20, col = 1:length(T_vals),
#        legend = paste("t_0 =", T_vals), bg = 'white')
