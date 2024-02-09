## Week 4 lecture 1

### Challenges I

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

## Challenges II

# Q2

A <- cbind(
  c(2, 4, -6, 4),
  c(3, 7, -10, 6),
  c(0, 2, 0, 4),
  c(0, 0, 1, 5)
)
b <- c(1, 2, 1, 0)
x <- c(23/2, -22/3, 11/3, -10/3)

all.equal(solve(A, b), x)

## Week 4 lecture 2

### Challenges I

# Q1

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

### Challenges II

# Q1

rsympdmat <- function(n) {
  # function to generate symmetric positive definite matrix
  # of dimension n x n
  # n is an integer
  A <- matrix(rnorm(n * n), n)
  crossprod(A, diag(abs(rnorm(n))) %*% A)
}

# generate matrix
A <- rsympdmat(5)
# upper-triangular Cholesky decomposition
U <- chol(A)
# cross product check
all.equal(crossprod(U), A)
# inverse through Cholesky
all.equal(solve(A), chol2inv(U))
# log determinant through Cholesky
all.equal(as.vector(determinant(A)$modulus), 2 * sum(log(diag(U))))

# as above, but in lower triangular form
L <- t(U)
all.equal(tcrossprod(L), A)
all.equal(solve(A), chol2inv(t(L)))
all.equal(as.vector(determinant(A)$modulus), 2 * sum(log(diag(L))))

# Q2

make_sym <- function(X) {
  # function to generate a symmetric matrix
  # X is a square matrix
  .5 * (X + t(X))
}

A <- matrix(rnorm(16) , 4, 4)
A2 <- make_sym(A)
# will almost always generate an error
chol(A2)
# and so is a quick way of checking whether a matrix
# is positive definite

## Week 4 lecture 3

### Challenges I

# Q1

hilbert <- function(n) {
  # Function to evaluate n by n Hilbert matrix.
  # Returns n by n matrix.
  ind <- 1:n
  1 / (outer(ind, ind, FUN = '+') - 1)
}
H4 <- hilbert(4)
b <- c(.04, 3.19, 3.26, 1.93)
x1 <- solve(H4, b)

U <- chol(H4)
L <- t(U)
x2 <- backsolve(t(L), forwardsolve(L, b))
all.equal(x1, x2)

