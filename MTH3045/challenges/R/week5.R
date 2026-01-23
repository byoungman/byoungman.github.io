## Week 5 lecture 1

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
