## Week 11 lecture 1

### Extra example

### Q1

negloglik <- function(pars, y) {
 # Function to evaluate derivative negative log-likelihood
  # pars is a 2-vector, where
  # pars[1] is mu and
  # pars[2] is sigma
  # y is a vector
  # returns a scalar
  mu <- pars[1]
  sigma <- pars[2]
  if (sigma <= 0)
    return(1e6)
  n <- length(y)
  n * log(2 * pi) / 2 + n * log(sigma) + sum((y - mu)^2) / 2 / sigma^2
}

y <- c(1.8, 2.1, 2.3, 2.4, 2.5, 2.7, 3.1, 4.3, 4.3, 4.4)

negloglik(c(1.5, 1.5), y)

-sum(dnorm(y, 1.5, 1.5, log = TRUE))

### Q2

negloglik_d1 <- function(pars, y) {
  # Function to evaluate derivative of 
  # negative log-likelihood w.r.t. mu and sigma
  # pars is a 2-vector, where
  # pars[1] is mu and
  # pars[2] is sigma
  # y is a vector
  # returns a 2-vector
  mu <- pars[1]
  sigma <- pars[2]
  n <- length(y)
  out <- numeric(2)
  res <- y - mu
  out[1] <- -sum(res) / sigma^2
  out[2] <- n / sigma - sum(res^2) / sigma^3
  out
}

negloglik_d1(c(1.5, 1.5), y)

### Q3

fit_bfgs <- optim(c(1.5, 1.5), negloglik, negloglik_d1, y = y, method = 'BFGS')
mles <- fit_bfgs$par

### Q4

negloglik_d1(mles, y)
mean(y)
sqrt(var(y) * (length(y) - 1) / length(y))

### Challenges I

### Q1

weib2_d0 <- function(pars, y, mult = 1) {
  # Function to evaluate re-parameterised Weibull log-likelihood
  # pars is a vector
  # pars[1] is $\tilde \lambda$
  # pars[2] is $\tilde k$
  # y can be scalar or vector
  # mult is a scalar defaulting to 1; so -1 returns neg. log likelihood
  # returns a scalar
  n <- length(y)
  epars <- exp(pars)
  out <- n * (pars[2] - epars[2] * pars[1]) + 
    (epars[2] - 1) * sum(log(y)) - sum((y / epars[1])^epars[2])
  mult * out
}

### Q2

fit_nm2 <- optim(log(c(1.6, 0.6)), weib2_d0, y = y0, mult = -1)
fit_nm2$par

### Q3

fit_nm$par
exp(fit_nm2$par)

## Week 11 lecture 2

### Challenges I

### Q1

weib2_d0 <- function(pars, y, mult = 1) {
  # Function to evaluate re-parameterised Weibull log-likelihood
  # pars is a vector
  # pars[1] is $\tilde \lambda$
  # pars[2] is $\tilde k$
  # y can be scalar or vector
  # mult is a scalar defaulting to 1; so -1 returns neg. log likelihood
  # returns a scalar
  n <- length(y)
  epars <- exp(pars)
  out <- n * (pars[2] - epars[2] * pars[1]) + 
    (epars[2] - 1) * sum(log(y)) - sum((y / epars[1])^epars[2])
  mult * out
}

y0 <- c(3.52, 1.95, 0.62, 0.02, 5.13, 0.02, 0.01, 0.34, 0.43, 15.5, 
        4.99, 6.01, 0.28, 1.83, 0.14, 0.97, 0.22, 0.02, 1.87, 0.13, 0.01, 
        4.81, 0.37, 8.61, 3.48, 1.81, 37.21, 1.85, 0.04, 2.32, 1.06)

optim(c(log(1.6), log(0.6)), weib2_d0, y = y0, mult = -1, method = 'SANN',
      control = list(maxit = 1e3, temp = 10))

### Q2

pars1 <- replicate(20,
                   optim(c(log(1.6), log(0.6)), weib2_d0, y = y0, mult = -1, method = 'SANN',
      control = list(maxit = 1e3, temp = 10))$par)

apply(pars1, 1, sd)

### Q3

pars2 <- replicate(20,
                   optim(c(log(1.6), log(0.6)), weib2_d0, y = y0, mult = -1, method = 'SANN',
                         control = list(maxit = 1e3, temp = 0.1))$par)

apply(pars2, 1, sd)

### Q4

rbind(apply(pars1, 1, sd),
      apply(pars2, 1, sd))
