## Week 3 lecture 1

### Challenges I

# Q1

n <- 35
y <- rgamma(n, 3, 2)
z <- rep(NA, n)
for (i in 4:(n - 3)) z[i] <- mean(y[i + seq(-3, 3)])
matplot(seq_len(n), cbind(y, z))

# Q2

z2 <- filter(y, rep(1/7, 7))
all.equal(z, z2)

z3 <- as.vector(z2)

all.equal(z, z3)

# Q3

all(runif(10) > .5)
any(runif(10) > .5)

### Challenges II

# Q1

browser_test <- function(x) {
  x1 <- 1
  x2 <- 2
  x3 <- 3
  x4 <- 4
  x5 <- 5
  x6 <- 6
  browser()
  x7 <- 7
  x8 <- log(x)
  c(x1, x2, x3, x4, x5, x6, x7, x8)
}


# Q2

# No sketch solutions for Q2-5 given.

## Week 3 lecture 2

### Challenges I

# Q1

n <- 1e3
y <- runif(n)
sort(y)[1]
min(y)

library(microbenchmark)

microbenchmark(
  apply(A, 1, sum),
  rowSums(A)
)

# minima
n <- 1e3
y <- runif(1e3)
microbenchmark(
  sort(y)[1],
  min(y)
)

# Q2

n <- 1e3
a <- rnorm(n)
b <- rnorm(n)
ab <- function(a, b) {
  s <- 0
  for (i in 1:length(a))
    s <- s + a[i] * b[i]
  s
}
ab(a, b)

# Q3

sum(a * b)
microbenchmark(
  ab(a, b),
  sum(a * b)
)

## Week 3 lecture 3

### Challenges I

# Q1

make_sym <- function(X) {
  .5 * (X + t(X))
}

m <- 1e3
n <- 2e2
p <- 1e2

n <- 3
A <- matrix(rnorm(n * n), n, n)
A2 <- make_sym(A)
all.equal(A2, t(A2))

# Q2

make_sym_check <- function(X) {
  # function to generate a symmetric matrix
  # X is a square matrix
  if (nrow(X) != ncol(X))
    stop('Supplied matrix X is not square')
  .5 * (X + t(X))
}

make_sym_check(A)
B <- matrix(rnorm(2 * n), n, 2)
make_sym_check(B)

# Q3

pdmatrix <- function(n) {
  L <- matrix(0, n, n)
  L[!lower.tri(L)] <- abs(rnorm(n * (n + 1) / 2))
  tcrossprod(L)
}
dp <- 2
A <- matrix(c(.91, .32, .62, .32, .31, .15, .62, .15, .47), 3, 3)#round(pdmatrix(3), 2)
x_1 <- as.matrix(c(-.12, .51, .22))#round(rnorm(3), 2))
x_2 <- as.matrix(c(-.32, -1.07, -1.36))#round(rnorm(3), 2))
x_3 <- as.matrix(c(-.97, -.16, -.36))#round(rnorm(3), 2))

A <- matrix(rnorm(1e3 * 2e2), 1e3)
B <- matrix(rnorm(1e3 * 1e2), 1e3)
all.equal(crossprod(A, B), t(A) %*% B)

microbenchmark::microbenchmark(
  crossprod(A, B), 
  t(A) %*% B)

# Q4

<!-- 4. Confirm that $\mathbf{x}_i^{\text{T}} \mathbf{Ax}_i > 0$ for -->
<!-- \[ -->
<!-- \mathbf{A} = `r write_matex(A)` -->
<!-- \] -->
<!-- and -->
<!-- \[ -->
<!-- \mathbf{x}_1 = `r write_matex(x_1)`,~\mathbf{x}_2 = `r write_matex(x_2)`~\text{and}~\mathbf{x}_3 = `r write_matex(x_3)`. -->
<!-- \] -->

<!-- <details><summary>**Solution**</summary> -->
<!-- ```{r, eval = FALSE, echo = week > w0} -->
<!-- ``` -->
<!-- </details> -->

<!-- 5. If $\mathbf{x}_i$ is stored in `R` as `x_i`, how can all $\mathbf{x}_i$ be simultaneously checked by forming `X = cbind(x_1, x_2, x_3)`? -->

<!-- <details><summary>**Solution**</summary> -->
<!-- ```{r, eval = FALSE, echo = week > w0} -->
<!-- ``` -->
<!-- </details> -->

