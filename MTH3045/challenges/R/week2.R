# Week 2 lecture 1

## Challenge I

a <- .1
for (i in 1:330) {
  a <- .1 * a
  print(a)
}

## Challenges II

## Q1

x <- 1:12
# %% 3 divides by three and then gives the remainder
x %% 3
# %/% 3 divides by three and then gives the integer part
x %/% 3

## Q2

x <- c(1, 2, NA, Inf)
range(x)                 # gives NA because x includes NA
range(x, na.rm = TRUE)   # recognises that Inf is big and ignores the NA
range(x, finite = TRUE)  # ignores the NA and Inf, i.e. anything non-finite

## Q3

x <- c(1, 1:10)
mean(x)
mean(x, trim = .1) # calculate the mean from the middle 90% of the data

# Week 2 lecture 2

## Challenges I

## Q1

1:3^2
(1:3)^2

1 + 2^3

1 + 2*3
(1 + 2)*3

## Q2

a <- c(2, 4, 6)
B <- matrix(c(2, 4, 4, 3, 3, 2, 1, 1, 3), 3, 3)
c <- matrix(c(2, 4, 6), 1)

B %*% a     # Ba
c %*% B     # cB
c <- t(a)   # creates c as the tranpose of a
t(a) %*% B  # avoids creating c

## Q3

x <- list(
  array(runif(36, -1, 1), c(3, 6, 2)),
  runif(4, -1, 1),
  matrix(runif(12, -1, 1), 4, 3)
)
x

## Challenges II

## Q1

x <- matrix(rnorm(200), 10, 20)
apply(x, 1, median)

## Q2

x <- lapply(1:4, function(i) rnorm(1 + rpois(1, 2), 3, 2))
lapply(1 + rpois(4, 2), function(x) rnorm(x, 3, 2))

## Q3

sapply(x, length)

# Week 2 lecture 3

## Challenge I

## Q1

n <- 15

# as a loop
x <- integer(n)
for (i in 1:n) {
  x[i] <- sample(1:6, 1, prob = c(.1, .05, .05, .2, .2, .4))
}
x

# more tidily with sample()
x <- sample(1:6, n, prob = c(.1, .05, .05, .2, .2, .4), replace = TRUE)

## Challenge 2

## Q1

x <- integer(0)
while(sum(x == 1) < 2) {
  x <- c(x, sample(1:6, 1, prob = c(.1, .05, .05, .2, .2, .4)))
}
x

