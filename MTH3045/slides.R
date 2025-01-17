# functions and data for lectures

bit2decimal <- function(x, e, dp = 20) {
  # function to convert bits to decimal form
  # x: the bits as a character string, with appropriate spaces
  # e: the excess
  # dp: the decimal places to report the answer to
  bl <- strsplit(x, ' ')[[1]] # split x into S, E and F components by spaces
  # and then into a list of three character vectors, each element one bit
  bl <- lapply(bl, function(z) as.integer(strsplit(z, '')[[1]]))
  names(bl) <- c('S', 'E', 'F') # give names, to simplify next few lines
  S <- (-1)^bl$S # calculate sign, S
  E <- sum(bl$E * 2^c((length(bl$E) - 1):0)) # ditto for exponent, E
  F <- sum(bl$F * 2^(-c(1:length(bl$F)))) # and ditto to fraction, F
  z <- S * 2^(E - e) * (1 + F) # calculate z
  out <- format(z, nsmall = dp) # use format() for specific dp
  # add (S, E, F) as attributes, for reference
  attr(out, '(S,E,F)') <- c(S = S, E = E, F = F) 
  out
}

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

dmvn3 <- function(y, mu, Sigma, log = TRUE) {
  # Function to evaluate MVN pdf
  # y, mu vectors
  # Sigma matrix
  # log is logical indicating whether logarithm
  # returns a scalar
  p <- length(y)
  res <- y - mu
  L <- t(chol(Sigma))
  out <- - sum(log(diag(L))) - 0.5 * p * log(2 * pi) -
    0.5 * sum(forwardsolve(L, res)^2)
  if (!log)
    out <- exp(out)
  out
}

fd <- function(x, f, delta = 1e-6, ...) {
  # Function to evaluate derivative w.r.t. vector by finite-differencing
  # x is a p-vector
  # fn is the function for which the derivative is being calculated
  # delta is the finite-differencing step, which defaults to 10ˆ{-6}
  # returns a vector of length x
  f0 <- f(x, ...)
  p <- length(x)
  f1 <- numeric(p)
  for (i in 1:p) {
    ei <- replace(numeric(p), i, 1)
    f1[i] <- f(x + delta * ei, ...)
  }
}
