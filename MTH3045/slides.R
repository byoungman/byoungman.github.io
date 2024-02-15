# functions and data for lectures

hilbert <- function(n) {
  # Function to evaluate n by n Hilbert matrix.
  # n is an integer
  # returns n by n matrix.
  ind <- 1:n
  1 / (outer(ind, ind, FUN = '+') - 1)
}