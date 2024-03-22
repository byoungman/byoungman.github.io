## Week 10 lecture 1

### Challenges I

# Q1

nd0 <- function(lambda, y) {
  # Function to evaluate exponential negative log-likelihood
  # lambda is a scalar
  # y is a vector
  # returns a scalar
  n <- length(y)
  -n * log(lambda) + lambda * sum(y)
}

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
  matrix(n / lambda^2, 1, 1)
}

# Q2

y <- c(.2, .5, .6, 1, 1.8, 2.5)
lambda_0 <- 1
nlminb(lambda_0, nd0, nd1, nd2, y = y)

# Q3

optimize(nd0, c(.5, 1.5), y = y)
optim(1, nd0, lower = .5, upper = 1.5, method = 'Brent', y = y)

# Q4

lambda0_seq <- 10^seq(-4, 4)
sapply(lambda0_seq, function(x) nlminb(x, nd0, nd1, nd2, y = y)$par)

## Week 10 lecture 2

### Challenges I

# Q1

nd0_tau <- function(tau, y) {
  # function to evaluate Exp(exp(tau)) neg. log lik.
  # tau is a scalar
  # y is a vector
  # returns a scalar
  n <- length(y)
  - n * tau + exp(tau) * sum(y)
}

y0 <- c(.2, .5, .6, 1, 1.8, 2.5)
tau0 <- 0
nd0_tau(tau0, y0)

nd1_tau <- function(tau, y) {
  # function to 1st deriv. w.r.t. to tau
  # of Exp(exp(tau)) neg. log lik.
  # tau is a scalar
  # y is a vector
  # returns a scalar
  n <- length(y)
  - n + exp(tau) * sum(y)
}

nd2_tau <- function(tau, y) {
  # function to 2nd deriv. w.r.t. to tau
  # of Exp(exp(tau)) neg. log lik.
  # tau is a scalar
  # y is a vector
  # returns a scalar
  as.matrix(exp(tau) * sum(y))
}

nd1_tau(tau0, y0)
nd2_tau(tau0, y0)

fit <- nlminb(tau0, nd0_tau, nd1_tau, nd2_tau, y = y0)
tau_hat <- fit$par
exp(tau_hat)

# Q2

optimize(nd0_tau, c(-1, 1), y = y0)

## Week 10 lecture 3

### Challenges I

# Q1

y <- c(1.8, 2.1, 2.3, 2.4, 2.5, 2.7, 3.1, 4.3, 4.3, 4.4)
sum(dnorm(y, 1, 1, log = TRUE))

# Q2

negloglik <- function(mu, y) {
 # Function to evaluate derivative negative log-likelihood
  # mu is a scalar
  # y is a vector
  # returns a scalar
  n <- length(y)
  n * log(2 * pi) / 2 + sum((y - mu)^2) / 2
}

negloglik(1, y)

# Q3

negloglik_d1 <- function(mu, y) {
  # Function to evaluate derivative of 
  # negative log-likelihood w.r.t. mu
  # mu is a scalar
  # y is a vector
  # returns a scalar
  -sum(y - mu)
}

negloglik_d1(1, y)

# Q4

fit_bfgs <- optim(1, negloglik, negloglik_d1, y = y, method = 'BFGS')
mu_hat <- fit_bfgs$par

# Q5

negloglik_d1(mu_hat, y)
mean(y)
