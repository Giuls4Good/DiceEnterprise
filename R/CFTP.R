CFTP <- function(k, roll.fun, update.fun, monotonic = FALSE, min = NA, max = NA, verbose = FALSE,...) {
  #if verbose = true, a list is returned containing the time required in terms of iterations
  #min, max are used only in montonoic case, otherwise they are ignored
  if(monotonic) {
    X_span <- c(min,max) #minimum and maximum state
  } else {
    X_span <- 1:k #all states
  }
  T <- 2
  X <- X_span
  B <- roll.fun(n=1,...)
  U <- runif(n=1)
  while(length(unique(X)) > 1) {
    B <- c(roll.fun(n=T/2,...),B) #NOT EFFICIENT
    U <- c(runif(n=T/2),U)
    counter <- 1
    for(i in X_span) {
      X[counter] <- update.fun(i,B,U)
      counter <- counter + 1
    }
    T <- 2*T
  }
  if(verbose) {
    return(list(X=X[1], T=T))
  } else {
    return(X[1])
  }

}
