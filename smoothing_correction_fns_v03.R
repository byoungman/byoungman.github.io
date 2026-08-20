mat2arr <- function(x, n_split) {
  # x is a matrix
  array(x, c(nrow(x), n_split, ncol(x) / n_split))
}

blockmax <- function(array, margin = c(1, 3)) {
  apply(array, margin, max)
}

# conv2revexp0 <- function(x, margin = 1) {
#   # x is an array
#   out <- 0 * x
#   d <- nrow(x)
#   ny <- ncol(x)
#   for (i in 1:d) {
#     samp <- sort(runif(length(x[i, ])))
#     temp <- samp[rank(x[i, ])]
#     out[i, ] <- -log(temp)
#   }
#   out
# }
# 
# conv2revexp1 <- function(x, margin = 1) {
#   # x is an array
#   out <- 0 * x
#   d <- nrow(x)
#   ny <- ncol(x)
#   ind2 <- 1:ny
#   ind1 <- ind2 - 1
#   for (i in 1:d) {
#     temp <- runif(ny, ind1 / ny, ind2 / ny)[rank(x[i, ])]
#     out[i, ] <- -log(temp)
#   }
#   out
# }

conv2revexp <- function(x, margin = 1) {
  # x is an array
  out <- 0 * x
  d <- nrow(x)
  ny <- ncol(x)
  ind2 <- 1:ny
  ind1 <- ind2 - 1
  for (i in 1:d) {
    temp <- runif(ny, ind1 / ny, ind2 / ny)[rank(- x[i, ])]
    out[i, ] <- -log(temp)
  }
  out
}

ess_exp <- function(x, margin = 2, cap = FALSE) {
  # x is a matrix
  ny <- ncol(x)
  d <- nrow(x)
  mins <- sapply(1:ny, function(i) min(x[, i]))
  ess <- 1 / mean(mins)
  if (cap)
    ess <- min(ess, d)
  ess
}

weight_exp <- function(x, margin = 2, cap = FALSE, mult = 1) {
  # x is a matrix
  d <- nrow(x)
  ess <- ess_exp(x, margin, cap)
  w <- ess / d
  out <- mult * w
  attr(out, 'ess') <- ess
  out
}

# is.block.max <- function(x, perm = c(2, 1, 3)) {
#   # x is an array
#   x_ind <- apply(x, c(1, 3), function(x) replace(numeric(length(x)), which.max(x), 1))
#   x_ind <- aperm(x_ind, perm)
#   if (!all(dim(x) == dim(x_ind)))
#     stop('Dimension mismatch')
#   x_ind
# }

is.block.max <- function(x, perm = c(2, 1, 3)) {
  # x is an array
  out <- 0 * x
  d <- dim(x)[1]
  ny <- dim(x)[3]
  for (i in 1:d) for (j in 1:ny) {
    out[i, , j] <- replace(numeric(365), which.max(x[i, , j]), 1)
  }
  out
}

ess_daily <- function(x) {
  # x is an array
  d <- dim(x)[1]
  ny <- dim(x)[3]
  nd <- prod(dim(x)) / d / ny
  x <- matrix(x, d)
  R <- tcrossprod(x - 1 / nd) / ny
  #R[R < 2 / nd] <- 0
  d / mean(rowSums(R))
}

weight_daily <- function(x, mult = 1) {
  # x is an array
  d <- dim(x)[1]
  ess <- ess_daily(x)
  w <- ess / d
  out <- mult * w
  attr(out, 'ess') <- ess
  out
}

generate_scaling <- function(n, kappa = 5) {
t <- 1:n
Ct <- exp(-abs(kappa * outer(t / n, t / n, FUN = '-'))^1.9999)
cCt <- chol(Ct)
R <- matrix(c(1, .9, .9, 1), 2, 2) %*% diag(c(1, .1))
cR <- chol(R)
ct <- cR %x% cCt
pars <- matrix(crossprod(ct, rnorm(2 * n)), ncol = 2)
cbind(pars[, 1], exp(pars[, 2]))
}

weigh <- function(D, phi, delta, type = c('exp', 'daily'), cap = FALSE) {
  
  if (as.numeric(as.character(phi)) > 1) {
    C <- diag(nrow(D))
  } else {
    C <- exp(-(as.numeric(as.character(phi)) * D)^as.numeric(as.character(delta)))
  }
  cC <- chol(C)
  
  z0 <- crossprod(cC, matrix(rnorm(n * nrow(D)), nrow(D)))
  z_arr <- mat2arr(z0, 365)
  
  out <- list(z = z_arr)
  
  if ('exp' %in% type) {
    z_max <- blockmax(z_arr, c(1, 3))
    z_rexp <- conv2revexp(z_max)
    wt_exp <- weight_exp(z_rexp, cap = cap)
    n_eff_exp <- attr(wt_exp, 'ess')
    out$exp <- list(ess = n_eff_exp, wt = wt_exp)
  }
  
  if ('daily' %in% type) {
    z_ind <- is.block.max(z_arr)
    wt_daily <- weight_daily(z_ind)
    n_eff_daily <- attr(wt_daily, 'ess')
    out$daily <- list(ess = n_eff_daily, wt = wt_daily)
  }
  
  out
  
}

fit_evgams <- function(data, wts, p0, q99, more.knots = FALSE) {
  if (!more.knots) {
    m_evgam <- lapply(wts, function(wt) evgam(list(z ~ s(t, bs = 'cr'), ~ s(t, bs = 'cr'), ~ 1), data = data, family = 'gev2', gamma = 1/wt))
  } else {
    m_evgam <- lapply(wts, function(wt) evgam(list(z ~ s(t, k = 20, bs = 'cr'), ~ s(t, k = 20, bs = 'cr'), ~ 1), data = data, family = 'gev2', gamma = 1/wt))
  }
s_evgam <- lapply(m_evgam, simulate, nsim = 1e3, prob = p0)
# p_evgam <- sapply(m_evgam, function(x) predict(x, prob = p0)[, 1])
# cov_evgam <- sapply(s_evgam, function(x) apply(x, 1, quantile, .05) < q99 & apply(x, 1, quantile, .95) > q99)
p_evgam <- lapply(m_evgam, function(x) predict(x, prob = p0, se.fit = TRUE))
cov_evgam <- sapply(p_evgam, function(x) x[[1]][, 1] - 2 * x[[2]][, 1] < q99 & x[[1]][, 1] + 2 * x[[2]][, 1] > q99)
p_evgam <- sapply(p_evgam, function(x) x[[1]][, 1])
list(p = p_evgam, cov = cov_evgam)
}

fit_evgmrfs <- function(data, wts, p0, W, q99) {
m_evgmrf <- lapply(wts, function(wt) evgmrf(data, W = W, weights = wt, model = c('icar', 'icar', NA), inits = 'same', outer = 'bfgs'))
s_evgmrf <- lapply(m_evgmrf, simulate, nsim = 1e3, prob = p0)
p_evgmrf <- sapply(m_evgmrf, function(x) predict(x, prob = p0)[[1]])
cov_evgmrf <- sapply(s_evgmrf, function(x) apply(x, 1, quantile, .05) < q99 & apply(x, 1, quantile, .95) > q99)
#p_evgmrf <- lapply(m_evgmrf, function(x) predict(x, prob = p0, se.fit = TRUE, progress = FALSE))
#cov_evgmrf <- sapply(p_evgmrf, function(x) x[[1]][[1]][, 1] - 2 * x[[2]][[1]][, 1] < q99 & x[[1]][[1]][, 1] + 2 * x[[2]][[1]][, 1] > q99)
#p_evgmrf <- sapply(p_evgmrf, function(x) x[[1]][[1]])
list(p = p_evgmrf, cov = cov_evgmrf)
}

make_df_cov <- function(tvals, evgams, evgmrfs, j, phi, q99, info) {
df_cov <- expand.grid(t = tvals, weight = c('standard', 'exponential', 'daily'), model = c('evgam', 'evgmrf'))
df_cov$coverage <- as.vector(cbind(evgams$cov, evgmrfs$cov))
df_cov$simulation <- j
df_cov$phi <- as.character(phi)
df_cov$bias <- abs(as.vector(cbind(evgams$p, evgmrfs$p) - q99))
df_cov$wt_exp <- info$exp$wt
df_cov$wt_daily <- info$daily$wt
df_cov
}

save(blockmax, conv2revexp, ess_daily, ess_exp, fit_evgams, fit_evgmrfs, generate_scaling, is.block.max, make_df_cov, mat2arr, weigh, weight_daily, weight_exp, file = '/home/ben/onedrive/Exeter/Research/Ayu/notes/smoothing_correction_fns.rda')

