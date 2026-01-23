# Week 1 lecture 1

## Q1

n <- 20
y <- rnorm(n, 1, 3)

mean2 <- function(x) sum(x) / length(x)
mean2(y)
mean(y)

var2 <- function(x) sum((x - mean2(x))^2) / (length(x) - 1)
var2(y)
var(y)

sd2 <- function(x) sqrt(var2(x))
sd2(y)
sd(y)

## Q2

pnorms <- function(z) {
  cbind(
    pnorm(z, lower.tail = FALSE),
    1 - pnorm(z),
    pnorm(-z)
  )
}

vals <- seq(0, 10, by = .5)
pnorms(vals)

## Exploratory and refresher exercises II

## Q3

y1 <- 1:10
y2 <- y1 + 1e9

bvar1 <- function(x) {
  mean(x^2) - mean(x)^2
}

bvar2 <- function(x) {
  mean((x - mean(x))^2)
}

bvar1(y1)
bvar2(y1)

bvar1(y2)
bvar2(y2)

# Week 1 lecture 2

## Challenge

## Q1

# confirm with result below

## Q2

conv2dec <- function(x, B) {
  a <- as.numeric(strsplit(x, '')[[1]])
  pows <- (length(a) - 1):0
  sum(a * B^pows)
}

x <- '101'
conv2dec(x, 2)
conv2dec(x, 3)
conv2dec(x, 4)
conv2dec(x, 5)

## Q3

conv2dec2 <- function(x, B) {
  a0 <- strsplit(x, '[.]')[[1]]
  a_left <- as.numeric(strsplit(a0[1], '')[[1]])
  a_right <- as.numeric(strsplit(a0[2], '')[[1]])
  a <- c(a_left, a_right)
  pows <- (length(a_left) - 1):(-length(a_right))
  sum(a * B^pows)
}

x2 <- '101.101'
conv2dec2(x2, 2)
conv2dec2(x2, 3)
conv2dec2(x2, 4)
conv2dec2(x2, 5)

# Week 1 lecture 3

## Challenge I

## Q1

sum(2^c(0:7))

## Q2

2^(-23)

## Q3

sum(rep(1, 8) * 2^c(0:7)) 

## Challenge II

## Q1

sum(2^c(0:10))

## Q2

2^(-52)

## Q3

sum(rep(1, 11) * 2^c(0:10)) 

## Challenge III

a <- .1
for (i in 1:330) {
  a <- .1 * a
  print(a)
}
