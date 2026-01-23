## Week 9 lecture 1

### Challenges I

N <- 1e3
V <- 2
mc <- V * mean(runif(N, -1, 1)^2)
true <- 2/3

### Challenges II

N <- 10^4
d <- 4
x <- matrix(runif(N * d, 1, 3), N)
V <- 2^d
I_hat <- V * mean(apply(1 + (x - 2)^2, 1, prod))
I_hat2 <- V * mean(exp(rowSums(log1p((x - 2)^2))))
true <- (8/3)^d
rel_err <- abs((true - I_hat) / true)
rel_err

## Week 9 lecture 2

### Challenges I

N <- 10^4
d <- 4
x <- matrix(runif(N * d, 1, 3), N)
V <- 2^d
I_hat <- V * mean(apply(1 + (x - 2)^2, 1, prod))
I_hat2 <- V * mean(exp(rowSums(log1p((x - 2)^2))))
true <- (8/3)^d
rel_err <- abs((true - I_hat) / true)
rel_err

### Challenges II

negloglik <- function(theta, y) {
  # Function to evaluate negative log-likelihood
  # theta is a scalar
  # y is a vector
  # returns a scalar
  n <- length(y)
  loglik <- n * (lgamma(2 * theta) - 2 * lgamma(theta))
  loglik <- loglik + (theta - 1) * sum(log(y) + log(1 - y))
  -loglik
}

y <- c(.16, .35, .37, .46, .67, .69, .79, .8, .91, .92)
negloglik(2, y)
optimize(negloglik, c(1, 3), y = y)

## Week 9 lecture 3

### Challenges I

### Q1

nd1 <- function(lambda, y) {
  # Function to evaluate first derivative of exponential
  # negative log-likelihood w.r.t. lambda
  # lambda is a scalar
  # y is a vector
  # returns a scalar
  n <- length(y)
  sum(y) - n / lambda
}

nd2 <- function(lambda, y) {
  # Function to evaluate second derivative of exponential
  # negative log-likelihood w.r.t. lambda
  # lambda is a scalar
  # y is a vector
  # returns a scalar
  n <- length(y)
  n / lambda^2
}

### Q2

y <- c(0.2, 0.5, 0.6, 1.0, 1.8, 2.5)
-nd1(1, y) / nd2(1, y)
