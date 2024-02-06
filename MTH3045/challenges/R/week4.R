## Week 4 lecture 1

## Challenges I

# Q1

dmvn1 <- function(y, mu, Sigma, log = TRUE) {
  # Function to evaluate multivariate Normal pdf
  # at y given mean mu and variance-covariance
  # matrix Sigma.
  # Returns 1x1 matrix, on log scale, if log == TRUE.
  p <- length(y)
  res <- y - mu
  out <- - 0.5 * determinant(Sigma)$modulus - 0.5 * p * log(2 * pi) -
    0.5 * t(res) %*% solve(Sigma) %*% res
  if (!log) 
    out <- exp(out)
  attributes(out) <- NULL
  out
}

y <- c(.7, 1.3, 2.6)
mu <- 1:3
Sigma <- matrix(c(4, 2, 1, 2, 3, 2, 1, 2, 2), 3, 3)
dmvn1(y, mu, Sigma)

## Week 4 lecture 2

# Question 1
L <- cbind(
  c(2, 2.3, .8, .3),
  c(0, 1.2, 2.1, 1.4),
  c(0, 0, 1.8, 1.1),
  c(0, 0, 0, .5)
)
U <- rbind(
  c(.3, 1.9, 1.9, 1.4),
  c(0, 2.5, .6, .1),
  c(0, 0, 2.5, .5),
  c(0, 0, 0, 2.1)
)
z <- c(-.47, .07, -.48, 1.45)

x1 <- forwardsolve(L, z)
x1
solve(L, z)

x2 <- backsolve(U, z)
x2
solve(U, z)

## Challenges II

# Question 1

A <- rsympdmat(5)
U <- chol(A)
all.equal(crossprod(U), A)
all(U[lower.tri(U)] == 0)
all(diag(U) > 0)

# Inverse of A
iA <- solve(A)
all.equal(chol2inv(U), iA)

# log determinant of A
determinant(A)
2 * sum(log(diag(U)))

# Question 2

