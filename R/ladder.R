DiscreteDist <- R6::R6Class("DiscreteDistribution",
  public = list(
    get_k = function() {private$k},
    get_R = function() {private$R}
  ),
  private = list(
    k = NULL, #Number of categories
    R = NULL #k long vector of R_i. This are not the associated probabilities, but it is used for disaggregations.
  )
)

Ladder <- R6::R6Class("Ladder",
  inherit = DiscreteDist,
  public = list(
    initialize = function(M,R) {
      #CONSTRUCTOR
      private$m <- ncol(M)
      private$k <- nrow(M)
      private$M <- M
      private$R <- R
      private$degree <- sum(M[1,])
      private$is.ladder() #Check it is valid
      if(private$valid) {
        private$is.fine() #Check fineness
        private$define.neighbourhoods() #Compute neighbourhoods
        private$is.connected() #Check connected
        private$compute.constant() #Compute constant a
      } else {
        stop("The ladder is not valid.")
      }
    },
    print = function() {
      #PRINT
      cat("Ladder Delta^",private$m," -> Delta^",private$k," of degree ",private$degree,"\n", sep="")
      cat("Valid ladder:",private$valid,"\n")
      cat("Fine ladder: ",private$fine,"\n", sep="")
      cat("Connected ladder: ",private$connected,"\n", sep="")
      cat("Constant a = ",private$a,"\n",sep="")
    },
    get_connected = function() {private$connected},
    update.fun = function(i,B,U) {
      #Update function for the ladder
      stopifnot(private$connected, private$fine, length(B)==length(U))
      currentState <- i
      for(c in 1:length(B)) {
        #Find the move
        move <- findInterval(private$a*U[c],cumsum(private$R[private$neigh[[currentState]][[B[c]]]]))+1 #+1 cause they start form 0
        if(move > length(private$neigh[[currentState]][[B[c]]])) {
          currentState <- currentState #Stay still
        } else {
          currentState <- private$neigh[[currentState]][[B[c]]][move]
        }
      }
      return(currentState)
    },
    sample = function(n,roll.fun = NULL, true_p = NULL,...) {
      #Get a sample from the ladder using CFTP
      if(is.null(roll.fun) && is.null(true_p)) {stop("Either declare roll.fun or the trye probabilities.")}
      if(is.null(roll.fun)) {roll.fun <- function(n) {sample(1:private$m, size = n, replace = TRUE, prob = true_p)}}
      res <- rep(NA, n)
      for(i in 1:n) {
        res[i] <- CFTP(k = private$k, roll.fun = roll.fun, update.fun = self$update.fun,
                    monotonic = FALSE,...)
      }
      return(res)
    },
    evalute = function(p) {
      #Return the values of the ladder for a fixed
      #vector of probabilities of rolling each face
      stopifnot(is.numeric(p), length(p) == private$m, isTRUE(all.equal(1,sum(p))))
      aux <- rep(NA, private$k)
      for(i in 1:private$k) {
        aux[i] <- prod(p^private$M[i,])*private$R[i]
      }
      return(aux/sum(aux))
    },
    impose.fineness = function() {
      if(private$fine) {
        return(list(obj = self$clone(), A = lapply(1:private$k, function(i) {i}))) #The ladder is already fine -> returns a copy of it
      }
      new_k <- nrow(unique(private$M))
      A <- vector("list", new_k) #initialize partition
      #WRONG: indices <- match(data.frame(t(private$M)), data.frame(unique(t(private$M)))) #It may have gaps!
      indices <- indicesUniqueMat(private$M) #Find unique rows of M and assign them an index.
      new_R <- numeric(length = new_k)
      new_M <- matrix(NA, ncol = private$m, nrow = new_k)
      for(i in 1:new_k) {
        #Populate A, M, R
        aux <- which(indices == i)
        A[[i]] <- aux
        new_R[i] <- sum(private$R[aux])
        new_M[i,] <- private$M[aux[1],]
      }
      return(list(obj = Ladder$new(M = new_M, R = new_R), A = A))
    },
    impose.connected = function() {
      if(!is.na(private$connected) && private$connected)  {
        return(list(obj = self$clone(), A = lapply(1:private$k, function(i) {i}))) #The ladder is already connected -> returns a copy of it
      }
      #We extend 1 degree at the time until the ladder is connected.
      #NOT EFFICIENT: we could precompute the minimum degree necessary
      #to have a connected ladder
      for(d in 1:private$degree) {
        discrete_simplex <- discreteSimplex(d=d,m=private$m) #all possible combinations
        new_M_list <- vector("list", private$k)
        for(i in 1:private$k) {
          new_M_list[[i]] <- data.frame(t(apply(discrete_simplex, 1, function(vec) {vec+private$M[i,]})))
        }
        #Construct a ladder with the new M and check if it is connected
        new_M <- as.matrix(rbindlist(new_M_list))
        test_ladder <- Ladder$new(M = new_M, R=rep(1,nrow(new_M)))
        if(test_ladder$get_connected()) {
          #The ladder is connected!
          #Compute the new vector R
          discrete_simplex_binomial <- as.numeric(apply(discrete_simplex,1,function(vec) {multinomialCoeff(d=d,n=vec)}))
          new_R <- rep(private$R, each=length(discrete_simplex_binomial))*discrete_simplex_binomial
          A <- split(1:nrow(new_M), rep(1:private$k, each=length(discrete_simplex_binomial)))
          return(list(obj = Ladder$new(M=new_M, R=new_R), A=A))
        }
      }
      stop("Impossible to create a connected ladder.") #Should never happen.
    }
  ),
  private = list(
    #FIELDS
    m = NULL, #Number of faces of the given die
    M = NULL, #k x m matrix describing the power of each p_i. The rows must all sum to the same value.
    degree = NULL, #degree of the ladder
    valid = FALSE, #FALSE if it is not a valid ladder. It is checked after construction.
    fine = NA,
    connected = NA,
    neigh = NULL, #neighbourhood. It is a double indexed list, the first index is the state, the second the roll of the die.
    a = NA, #constant for the Markov chain
    #METHODS
    is.ladder = function() {
      private$valid = TRUE
      #Check if all the fields are properly initialized
      #SPLIT IT TO GIVE INFORMATIVE ERRORS
      if(is.null(private$m) || is.null(private$k) ||
         is.null(private$M) || is.null(private$R) ||
         !is.wholenumber(private$m) || !is.wholenumber(private$k) ||
         !is.matrix(private$M) || !is.numeric(private$R) ||
         nrow(private$M) != private$k || ncol(private$M) != private$m ||
         length(private$R) != private$k ||
         length(unique(rowSums(private$M))) > 1 ||
         !all(private$R > 0)) {
        private$valid = FALSE
      }
      invisible(self)
    },
    is.fine = function() {
      #Check for fineneness condition by checking if the rows of M are unique
      private$fine = TRUE
      if(nrow(unique(private$M)) != nrow(private$M)) {
        private$fine = FALSE
      }
      invisible(self)
    },
    define.neighbourhoods = function() {
      #Define neighbourhoods
      #Initialize
      private$neigh <- vector("list", private$k)
      for(i in 1:private$k) {
        private$neigh[[i]] <- vector("list", private$m)
        for(b in 1:private$m) {
          private$neigh[[i]][[b]] <- numeric()
        }
      }
      #Scroll through matrix M and fill the neighbourhoods dynamically (not efficient)
      for(i in 1:nrow(private$M)) {
        dist <- apply(private$M, 1, function(r) {sum(abs(private$M[i,] - r))}) #Compute distance
        for(j in which(dist <= 2)) {
          if(j != i) { #Different state
            if(!isTRUE(all.equal(private$M[i,], private$M[j,]))) {
              #The two rows are not equal. If they are equal it means that the ladder is not fine.
              aux <- private$M[i,] - private$M[j,]
              private$neigh[[i]][[which(aux == -1)]] <- append(private$neigh[[i]][[which(aux == -1)]],j) #SLOW AND TERRIBLE
            } else {
              #The two rows are equal -> ladder is not fine.
              #Neighbourhoods are computed anyway cause they are used to check for connected condition.
              #We just put the one with the same number in the first neighbourhood by convention.
              #Since for not fine ladders neighbourhoods are used only to check for connected condition
              #this does not affect it.
              private$neigh[[i]][[1]] <- append(private$neigh[[i]][[1]],j)
            }
          }
        }
      }
      invisible(self)
    },
    is.connected = function() {
      #Check for connected condition
      private$connected <- all(private$is.connected.aux(1,rep(FALSE,private$k)))
      invisible(self)
    },
    is.connected.aux = function(x,visited) {
      #Recursively check if it is connected. Starting from point x, remove its neighbourhoods in C and keep on researching.
      if(visited[x]) {
        return(visited)
      }
      visited[x] <- TRUE
      for(y in unlist(private$neigh[[x]])) {
        visited <- private$is.connected.aux(y,visited)
      }
      return(visited)
    },
    compute.constant = function() {
      #Compute the constant a = max(max(sum(R_i)))
      if(private$fine && private$connected) {
        aux <- rep(NA,private$k*private$m)
        counter <- 1
        for(i in 1:private$k) {
          for(b in 1:private$m) {
            aux[counter] <- sum(private$R[private$neigh[[i]][[b]]])
            counter <- counter+1
          }
        }
        private$a <- max(aux)
      }
      invisible(self)
    }
  ))
