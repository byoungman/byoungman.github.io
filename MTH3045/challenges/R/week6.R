# week 6 challenges

## lecture 1

### Q1

m <- 5
n <- 3
A <- matrix(rexp(m * n), m , n)
qrA <- qr(A)
Q <- qr.Q(qrA)
R <- qr.R(qrA)
all.equal(Q %*% R, A)

### Q2

# no code given

### Q3

svdA <- svd(A)
all.equal(prod(svdA$d), abs(prod(diag(R))))

### Q4

crossA <- crossprod(A)
all.equal(crossprod(R), crossA)

### Q5
L <- t(chol(crossA))
all.equal(abs(R), t(L))

## Challenge II

### Q1

hilbert <- function(n) {
  # Function to evaluate n by n Hilbert matrix.
  # n is an integer
  # returns n by n matrix.
  ind <- 1:n
  1 / (outer(ind, ind, FUN = '+') - 1)
}

A <- hilbert(8)

U <- cbind(c(4.7, 3.2, 0.2, 1.1, 5, 8.3, 1.1, 5.6), 
           c(10.6, 22.5, 14.1, 11.5, 3.3, 12.9, 9, 3.9))

wsm <- function(iA, U, C = diag(ncol(U)), V) {
  # function to compute Woodbury-Sherman-Morrison formula
  # iA is a m x m matrix
  # U and V are m x n matrices
  # C is a n x n matrix, and defaults to nrow(U) identity matrix
  # returns m x m matrix
  iC <- solve(C)
  iAU <- iA %*% U
  iA - iAU  %*% solve(iC + crossprod(V, iAU), crossprod(V, iA))
}

wsm(solve(A), U, V = U)

### Q2

B <- A + tcrossprod(U)
iB <- solve(B)
all.equal(iB, wsm(solve(A), U, V = U))

### Q3

C <- cbind(
  c(3.9, -1.7, -1.3, -.9),
  c(-1.7, 2.9, 1, 1.2),
  c(-1.3, 1, 9.4, 1.3),
  c(-.9, 1.2, 1.3, .7)
)

wsm(solve(A), U, C, U)

### Q4

D <- A + U %*% tcrossprod(C, U)
all.equal(solve(D), wsm(solve(A), U, C, U))

# Week 6 lecture 2

## Challenges I
# finite-differencing of exp(-x^2)

x <- seq(-1, 1, by = .1)
deriv0 <- - 2 * x * exp(-x^2)
plot(x, deriv0, type = 'l')
delta1 <- sqrt(.Machine$double.eps)
f <- function(x) exp(-x^2)
deriv1 <- (f(x + delta1) - f(x)) / delta1
lines(x, deriv1, col = 2, lwd = 3, lty = 2)

## Challenges II

### Q1

loglik <- function(lambda, y) {
  n <- length(y)
  n * log(lambda) - lambda * sum(y)
}

y0 <- rexp(20, 3)
loglik(2, y0)

### Q2

# gradient, i.e. first derivative
loglik_d1 <- function(lambda, y) {
  n <- length(y)
  n / lambda - sum(y)
}
loglik_d1(2, y0)

# Hessian, i.e. second derivative
loglik_d2 <- function(lambda, y) {
  n <- length(y)
  matrix(- n / lambda^2, 1, 1)
}
loglik_d2(2, y0)

### Q3
(loglik(2 + 1e-6, y0) - loglik(2, y0)) / 1e-6
loglik_d1(2, y0)

# Week 6 lecture 3

## Challenges I

### Q1

N <- 10
a <- 0
b <- 2
bins <- seq(a, b, length = N + 1)
h <- bins[2] - bins[1]
mids <- bins[-1] - .5 * h
f <- function(x) exp(-x^2/2)
f_mids <- f(mids)
I_hat <- h * sum(f_mids)

### Q2
# true
I0 <- sqrt(2 * pi) * (pnorm(2) - .5)
rel_err <- abs((I0 - I_hat) / I0)
100 * rel_err

## Challenges II

### Q1

f <- function(x) exp(-x^2/2)
N <- 10
a <- 0
b <- 2
h <- (b - a) / N
simpson <- f(a) + f(b)
x1i <- a + h * (2 * (1:N) - 1) / 2
simpson <- simpson + 4 * sum(f(x1i))
x2i <- a + h * (1:(N - 1))
simpson <- simpson + 2 * sum(f(x2i))
simpson <- h * simpson / 6

### Q2

true <- sqrt(2 * pi) * (pnorm(b) - .5)
rel_err <- abs((true - simpson) / true)
rel_err

## Challenges III

### Q1

f <- function(x) exp(-x^2/2)
N <- 7
gq <- pracma::gaussLegendre(N, 0, 2)
I_hat <- sum(gq$w * f(gq$x))

### Q2

true <- sqrt(2 * pi) * (pnorm(2) - .5)
rel_err <- abs((true - I_hat) / true)
