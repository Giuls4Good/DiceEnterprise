#' Coupling From The Past
#'
#' Sample from the stationary distribution of a discrete Markov chain using
#' the Coupling From the Past algorithm.
#'
#' @details The function is designed to sample from a fine and connected ladder using
#' rolls of a given die and uniform random variables. It requires an update function to be
#' defined as well as a function to roll the original die. It is mostly used internally by
#' the \code{Ladder} class.
#' @return If \code{verbose=FALSE}, returns a state sampled from the stationary distribution.
#' If \code{verbose=TRUE}, returns a list where the first element is a state from the stationary
#' distribution and the second element are the time steps required by the algorithm to finish.
#' @param k number of states of the Markov chain
#' @param roll.fun user defined function that rolls the original die. Must have as first input the size
#' of the returned sample. Additional inputs can be passed.
#' @param update.fun update function of the Markov chain. See the method \code{update.fun} of
#' class \code{Ladder} to check how it should be defined.
#' @param monotonic if \code{TRUE}, monotonic CFTP is implemented
#' @param min minimum state used in the monotonic CFTP. Ignored if \code{monotonic = FALSE}.
#' @param max maximum state used in the monotonic CFTP. Ignored if \code{monotonic = FALSE}.
#' @param verbose if \code{TRUE}, the time steps required by the algorithm is returned as well.
#' @param double_time if \code{TRUE} at each iteration of the algorithm, the time is doubled instead
#' of increased by one. Can lead to slower implementation if rolling the die is computationally costly.
CFTP <- function(k, roll.fun, update.fun, monotonic = FALSE, min = NA, max = NA, verbose = FALSE, double_time = FALSE,...) {
  #if double_time = true, it doubles the time span at each step. Otherwise it is just increased by one.
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
  #Move already forward (as coalesce may occur in 1 time step)
  counter <- 1
  for(i in X_span) {
    X[counter] <- update.fun(i,B,U)
    counter <- counter + 1
  }
  #Loop
  while(length(unique(X)) > 1) {
    if(double_time) {
      B <- c(roll.fun(n=T/2,...),B) #NOT EFFICIENT
      U <- c(runif(n=T/2),U)
    } else {
      B <- c(roll.fun(n=1,...),B) #NOT EFFICIENT
      U <- c(runif(n=1),U)
    }
    counter <- 1
    for(i in X_span) {
      X[counter] <- update.fun(i,B,U)
      counter <- counter + 1
    }
    if(double_time) {
      T <- 2*T
    } else {
      T <- T+1
    }
  }
  if(verbose) {
    if(double_time){
      return(list(X=X[1], T=T/2))
    } else {
      return(list(X=X[1], T=T-1))
    }
  } else {
    return(X[1])
  }

}
