# functions and data for lectures

hilbert <- function(n) {
  # Function to evaluate n by n Hilbert matrix.
  # n is an integer
  # returns n by n matrix.
  ind <- 1:n
  1 / (outer(ind, ind, FUN = '+') - 1)
}

# operating temperatures
temp <- c(35.3, 29.7, 30.8, 58.8, 61.4, 71.3, 74.4, 76.7, 70.7, 57.5,
46.4, 28.9, 28.1, 39.1, 46.8, 48.5, 59.3, 70, 70, 74.5, 72.1,
58.1, 44.6, 33.4, 28.6)
# number of operational days
days <- c(20, 20, 23, 20, 21, 22, 11, 23, 21, 20, 20, 21, 21, 19, 23,
20, 22, 22, 11, 23, 20, 21, 20, 20, 22)
# output from factory
output <- c(10.98, 11.13, 12.51, 8.4, 9.27, 8.73, 6.36, 8.5, 7.82, 9.14,
8.24, 12.19, 11.88, 9.57, 10.94, 9.58, 10.09, 8.11, 6.83, 8.88,
7.68, 8.47, 8.86, 10.36, 11.08)
# data frame
prod <- data.frame(temp = temp, days = days, output = output)

y <- c(.7, 1.3, 2.6)
mu <- 1:3
Sigma <- matrix(c(4, 2, 1, 2, 3, 2, 1, 2, 2), 3, 3)
