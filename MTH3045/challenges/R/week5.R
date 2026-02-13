## Week 5 lecture 1

### Challenges I

### Q1

Sigma <- matrix(c(4, 2, 1, 2, 3, 2, 1, 2, 2), 3, 3)
eS <- eigen(Sigma, symmetric = TRUE)
eS

### Q2

U <- eS$vectors
all.equal(crossprod(U), diag(3))

### Q3

lambda <- eS$values
detS <- prod(lambda)
all.equal(detS, det(Sigma))
log(detS)
sum(log(lambda))
determinant(Sigma)

### Challenges II

# one loop-based approach
matpow1 <- function(A, k) {
  # function to calculate matrix power of symmetric matrix A
  # based on a for loop
  # k is an integer
  # return matrix of same dimension as A
  Ak <- A
  if (k > 1) {
    for (i in 2:k) {
      Ak <- Ak %*% A
    }
  }
  Ak
}

# an alternative loop-based approach
matpow1 <- function(A, k) {
  # function to calculate matrix power of symmetric matrix A
  # based on a for loop
  # k is an integer
  # return matrix of same dimension as A
  Ak <- diag(nrow(A))
  for (i in 1:k) {
    Ak <- Ak %*% A
  }
  Ak
}

H <- hilbert(4)
Hpow5 <- matpow1(H, 5)
all.equal(matpow1(H, 5), Hpow5)

matpow2 <- function(A, k) {
  # as above, but based on the eigen-decomposition of A
  eA <- eigen(A, symmetric = TRUE)
  U <- eA$vectors
  lambda <- eA$values
  U %*% tcrossprod(diag(lambda^k), U)
}

all.equal(matpow2(H, 5), Hpow5)

### Challenges I

y <- c(.7, 1.3, 2.6)
mu <- 1:3
Sigma <- matrix(c(4, 2, 1, 2, 3, 2, 1, 2, 2), 3, 3)

S.svd <- svd(Sigma)
S.d <- S.svd$d
S.V <- S.svd$v

b2 <- crossprod(S.V, y - mu)
x2 <- b2 / S.d
z <- S.V %*% x2

solve(Sigma, y - mu)

### Challenges II

H <- hilbert(7)
svdH <- svd(H)
r <- 5
U_r <- svdH$u[, 1:r]
V_r <- svdH$v[, 1:r]
D_r <- svdH$d[1:r]
H_r <- U_r %*% (D_r * t(V_r))

# Week 5 lecture 3

## Challenges I

### Q1

hilbert <- function(n) {
  # Function to evaluate n by n Hilbert matrix.
  # Returns n by n matrix.
  ind <- 1:n
  1 / (outer(ind, ind, FUN = '+') - 1)
}
H4 <- hilbert(4)

qrH4 <- qr(H4)
Q <- qr.Q(qrH4)
R <- qr.R(qrH4)
all.equal(Q %*% R, H4)

### Q2

prod(abs(diag(R)))
det(H4)

### Q3

b <- c(.04, 3.19, 3.26, 1.93)
backsolve(R, crossprod(Q, b))
qr.solve(H4, b)
qr.solve(qrH4, b)


## Challenges II

### Q1

m <- 5
n <- 3
A <- matrix(rexp(m * n), m , n)
qrA <- qr(A)
Q <- qr.Q(qrA)
R <- qr.R(qrA)
all.equal(Q %*% R, A)

### Q2

### Q3

svdA <- svd(A)
all.equal(prod(svdA$d), abs(prod(diag(R))))

### Q4

crossA <- crossprod(A)
all.equal(crossprod(R), crossA)

### Q5

L <- t(chol(crossA))
all.equal(abs(R), t(L))

