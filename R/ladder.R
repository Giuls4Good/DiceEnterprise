Ladder <- R6::R6Class("Ladder",
  public = list(
    initialize = function(m,k,M,R) {
      #CONSTRUCTOR
      private$m <- m
      private$k <- k
      private$M <- M
      private$R <- R
      private$degree <- sum(M[1,])
      private$is.ladder() #Check it is valid
      if(private$valid) {
        private$is.fine() #Check fineness
        if(private$fine) {
          private$define.neighbourhoods() #Compute neighbourhoods
          private$is.connected() #Check connected
        }
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
    }
  ),
  private = list(
    #FIELDS
    m = NULL, #Number of faces of the given die
    k = NULL, #Number of states of the ladder
    M = NULL, #k x m matrix describing the power of each p_i. The rows must all sum to the same value.
    R = NULL, #k long vector of R_i
    degree = NULL, #degree of the ladder
    valid = FALSE, #FALSE if it is not a valid ladder. It is checked after construction.
    fine = NULL,
    connected = NULL,
    neigh = NULL, #neighbourhood. It is a double indexed list, the first index is the state, the second the roll of the die.
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
      if(nrow(unique(M)) != nrow(M)) {
        private$fine = FALSE
      }
      invisible(self)
    },
    define.neighbourhoods = function() {
      #Define neighbourhoods
      if(!private$fine) {
        stop("Neighbourhoods can be computed only for fine ladders.")
      }
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
        dist <- apply(M, 1, function(r) {sum(abs(M[i,] - r))}) #Compute distance
        for(j in which(dist == 2)) {
          aux <- M[i,] - M[j,]
          private$neigh[[i]][[which(aux == -1)]] <- append(private$neigh[[i]][[which(aux == -1)]],j) #SLOW AND TERRIBLE
        }
      }
      invisible(self)
    },
    is.connected = function() {
      #Check for connected condition
      private$connected = FALSE
      invisible(self)
    }
  ))
