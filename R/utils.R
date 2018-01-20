is.wholenumber <-  function(x, tol = .Machine$double.eps^0.5)  {
  return(abs(x - round(x)) < tol)
}

rsimplex <- function(n, m) {
  res <- matrix(NA, nrow = n, ncol = m)
  for(i in 1:n) {
    t <- 1
    for(j in 1:(m-1)) {
      res[i,j] <- runif(1, min = 0, max = t)
      t <- t - res[i,j]
    }
    res[i,m] <- t
  }
  return(res)
}
