CFTP <- function(k, roll.fun, update.fun, monotonic = FALSE, min = NA, max = NA,...) {
  if(monotonic) {
    return(MonotonicCFTP(k = k, roll.fun = rol.fun, update.fun = update.fun, min = min, max = max))
  }
  T <- 2
  X <- 1:k
  B <- roll.fun(n=1,...)
  U <- runif(n=1)
  while(length(unique(X)) > 1) {
    B <- c(roll.fun(n=T/2,...),B) #NOT EFFICIENT
    U <- c(runif(n=T/2),U)
    for(i in 1:k) {
      X[i] <- update.fun(i,B,U)
    }
    T <- 2*T
  }
  return(X[1])
}

MonotonicCFTP <- function(k, roll.fun, update.fun, min, max,...) {
  stop("Not implemented yet.")
}
