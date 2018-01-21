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

sample.vec <- function(x, ...) {
  #Safer! Works as sample() but also for vector of length 1.
  x[sample.int(length(x), ...)]
}

disaggregationSample <- function(x, origin = c("original", "disaggregation"),
                           disaggDist, A) {
  origin <- match.arg(origin)
  stopifnot("DiscreteDistribution" %in% class(disaggDist),
            is.list(A),
            length(unlist(A)) == disaggDist$get_k(),
            min(x) > 0,
            max(x) <= disaggDist$get_k())
  res <- rep(NA, length(x))
  counter <- 1
  for(s in x) {
    if(origin == "disaggregation") {
      #It is directly mapped into one
      for(i in 1:length(A)) {
        if(s %in% A[[i]]) {
          res[counter] <- i
          counter <- counter +1
          break
        }
      }
    } else {
      #Proportionally assign to one class
      res[counter] <- sample.vec(A[[s]], size = 1, prob = disaggDist$get_R()[A[[s]]])
      counter <- counter + 1
    }
  }
  return(res)
}

discreteSimplex <- function(d,m) {
  #Generates all possible vectors of length m such that the
  #sum of its entries is equal to d.
  #Returns a matrix with m columns
  partitions <- restrictedparts(n=d,m=m,include.zero = TRUE) #all possible partitions. Still need to permutate them.
  l_permutations <- lapply(split(t(partitions), 1:ncol(partitions)), function(vec) {
    data.frame(unique(t(apply(permutations(n = m, r = m), 1, function(x) vec[x]))))
  }) #list of all the permutations (by row)
  return(as.matrix(rbindlist(l_permutations)))
}

multinomialCoeff <- function(d,n) {
  #Compute d!/(n[1]!...n[length(n)]!)
  exp(lfactorial(d)-sum(sapply(n, function(i) {lfactorial(i)})))
}

indicesUniqueMat <- function(M) {
  #Returns the indices of unique elements in the matrix (row-wise)
  #THERE MAY BE A MORE EFFICIENT IMPLEMENTATION
  indices <- rep(0, nrow(M))
  counter <- 1
  for(i in 1:nrow(M)) {
    if(indices[i] == 0) {
      indices[i] <- counter
      if(i < nrow(M)) {
        for(j in (i+1):nrow(M)) {
          if(isTRUE(all.equal(M[i,],M[j,]))) {
            indices[j] <- counter
          }
        }
        counter <- counter + 1
        }
    }
  }
  return(indices)
}
