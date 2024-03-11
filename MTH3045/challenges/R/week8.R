## Week 8 lecture 1

### Challenges I

# Q1

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

# Q2

true <- sqrt(2 * pi) * (pnorm(b) - .5)
rel_err <- abs((true - simpson) / true)
rel_err

### Challenges II

# Q1
f <- function(x) exp(-x^2/2)
N <- 7
gq <- pracma::gaussLegendre(N, 0, 2)
I_hat <- sum(gq$w * f(gq$x))

# Q2
true <- sqrt(2 * pi) * (pnorm(2) - .5)
rel_err <- abs((true - I_hat) / true)

## Week 8 lecture 2

### Challenges I

# Q1
N <- 10
x1 <- x2 <- (1:N - .5) / N
h <- 1/N
w <- h^2
lambda <- 2
midpoint <- 0
for (i in 1:N) for (j in 1:N)
  midpoint <- midpoint + w * lambda^2 * exp(-lambda*(x1[i] + x2[j]))
midpoint  

# Q2
true <- pexp(1, lambda)^2
rel_err <- abs((true - midpoint) / true)

### Challenges II

# Q1
N <- 7
xw1 <- pracma::gaussLegendre(N, 0, 1)
xw2 <- pracma::gaussLegendre(N, 0, 2)
gq <- 0
f <- function(x1, x2, lambda = 2)
  lambda^2 * exp(-lambda * (x1 + x2)) 
for (i in 1:N) for (j in 1:N) 
  gq <- gq + xw1$w[i] * xw2$w[j] * f(xw1$x[i], xw2$x[j])
gq

# Q2
true <- pexp(1, 2) * pexp(2, 2)
rel_err <- abs((true - gq) / true)

## Week 8 lecture 3

### Challenge I

# Q1
N <- 5
d <- 4
b <- 3
a <- 1
h <- (b - a) / N
x <- a + h * (1:N - .5)
xx <- lapply(1:d, function(i) x)
X <- expand.grid(xx)
w <- h^d
f <- apply(1 + (X - 2)^2, 1, prod)
midpoint <- w * sum(f)
true <- (8/3)^d
rel_err <- abs((true - midpoint) / true)
rel_err
