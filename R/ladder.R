#' Categorical distributions
#'
#' Class to describe cateogrical distributions such as multi-variate ladders.
#'
#' @docType class
#' @return Object of \code{\link{R6Class}} with no methods.
#' @format \code{\link{R6Class}} object.
#' @field k number of categories
#' @field R vector of length \code{k} with the associated probabilities (may not be normalized)
DiscreteDist <- R6::R6Class("DiscreteDistribution",
  public = list(
    get.k = function() {private$k},
    get.R = function() {private$R}
  ),
  private = list(
    k = NULL, #Number of categories
    R = NULL #k long vector of R_i. These are not the associated probabilities, but it is used for disaggregations.
  )
)

#' Ladders
#'
#' Class to describe both univariate and multivariate ladders.
#'
#' @docType class
#' @return Object of \code{\link{R6Class}} with methods to construct connected and fine ladders as well as sampling from them.
#' @format \code{\link{R6Class}} object.
#' @details A multivariate ladder is a categorical distribution where the probability of each outcome is of the form
#' \deqn{\pi_i(p) = R_i\frac{\prod_{j=1}^m p_j^{n_{i,j}}}{C(p)}}
#' where \eqn{C(p)} is a polynomial with real coefficients not divisible by any \eqn{p_j},
#' \eqn{R_i} is a strictly positive constant, \eqn{n_{i,j}} are possibly null positive integers and \eqn{|n_i| = d}.
#' A ladder is fine if there are no redundant \eqn{n_{i,j}} and it is connected if for any \eqn{\pi_i} there exists
#' a \eqn{\pi_k} where only one of the \eqn{n_{i,j}} is increased by one and another one is decreased by one.
#' @examples
#' Ladder$new(M=matrix(c(0,3,2,1,3,0), ncol = 2, byrow=TRUE),R=c(3,2,sqrt(2)))
#' @section Methods:
#' \describe{
#'   \item{\code{new(M,R)}}{Construct an object of the class \code{Ladder}.
#'   M is a \eqn{m \times k} matrix describing the powers of \eqn{p_1,p_2,...,p_k} and
#'   R is a \eqn{k}-long vector of coefficients.}
#'   \item{\code{print}}{Print information about the ladder, such as if it is a valid/connected/fine ladder.}
#'   \item{\code{get.connected}}{Return \code{TRUE} if the ladder is connected.}
#'   \item{\code{update.fun(i,B,U)}}{Update function for the Markov chain. \code{i} is the current state,
#'   \code{B} is a vector of rolls of the original die, \code{U} is a vector of uniform random variables.}
#'   \item{\code{update.fun.R(i,B,U)}}{Same as \code{update.fun}, but implemented in R rather than in C++. Generally slower.}
#'   \item{\code{update.fun.global(i,B,U)}}{Update function for the Markov chain. This update function is defined using
#'   global properties of the chain. Albeit being valid, it leads to a slow CFTP implementation and it is thus deprecated.
#'   Use \code{update.fun} instead.}
#'   \item{\code{update.fun.slow(i,B,U)}}{Update function for the Markov chain. Albeit being valid, it's a slower
#'   implementation than \code{update.fun.R} and \code{update.fun} and it is thus deprecated. Use \code{update.fun} instead.}
#'   \item{\code{sample(n,roll.fun,true_p=NULL,num_cores=1,verbose=FALSE,global=FALSE,double_time=FALSE,...)}}{Get a sample from a
#'   fine and connected ladder via Coupling From The Past (see \code{CFTP}).
#'   \code{n} is required size of the sample, \code{roll.fun} is a user defined R function that rolls the original die (optional parameters
#'   can be passed). Instead of \code{roll.fun}, the user can defined a fixed value of probabilites of the original die via \code{true_p} for
#'   debug purposes. \code{num_cores} sets the numbers of cores used (supported only on Linux and Mac OS). If \code{verbose = TRUE} a list is
#'   returned where the first element is the sample and the second element are the rolls required to get such sample. When
#'   \code{global = TRUE}, \code{update.fun.global} is used instead of \code{update.fun}, leading to a valid but slower implementation and should not
#'   be used. When \code{double_time = TRUE} at each iteration of the CFTP algorithm, time is doubled instead of increased by 1 leading
#'   to a slower implementation if rolling the die is costly.}
#'   \item{\code{evalute(p)}}{Evaluate the ladder when the true probability of the original die is given. Useful for debug purposes.}
#'   \item{\code{impose.fineness}}{Returns a list where the first object is a new fine ladder and the second object is
#'   a vector that can be used to transform a sample from the new ladder in a sample from the original one. If the original ladder is
#'   connected, so it is the new one. See \code{disaggregation.sample}.}
#'   \item{\code{impose.connected}}{Returns a list where the first object is a new connected ladder and the second object is
#'   a vector that can be used to transform a sample from the new ladder in a sample from the original one. See \code{disaggregation.sample}.}
#'   \item{\code{get.a}}{Returns the global constant \eqn{a} used to define \code{update.fun.global}. Useful for debug purposes.}
#'   \item{\code{get.P}}{Returns the transition matrix of the chain. Contains only the coefficients and not the values of the \eqn{p_i}s.}
#'   \item{\code{get.P.moves}}{Returna a \eqn{k \timex k} matrix, where \eqn{k} is the number of states, describing which roll is necessary
#'   for the chain to move from state \eqn{i} to \eqn{j}.}
#'   \item{\code{get.P.cumsum}}{For internal use. Returns a a doubled indexed list.
#'   The first index is the current state, the second is the current roll. Returns the cumulative sum of the coefficients.}
#'   \item{\code{get.P.moves.list}}{For internal use. Returns a a doubled indexed list.
#'   The first index is the current state, the second is the current roll. Returns the index of the possible moves.}
#'   }
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
        private$compute.transition.matrix() #Compute transition matrix for the update function
        if(private$m == 2){ #Save minimum and maximum states
          private$min_state <- which.min(private$M[,1])
          private$max_state <- which.max(private$M[,1])
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
      cat("Constant a = ",private$a,"\n",sep="")
    },
    get.connected = function() {private$connected},
    get.fine = function() {private$fine},
    update.fun.global = function(i,B,U) {
      #Update function for the ladder using the global constant a (NOT EFFICIENT)
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
    update.fun.R = function(i,B,U) { #This is an R version of the C++ code of update.fun
      #Update function for the ladder using local moves (more efficient than the global bersion)
      stopifnot(private$connected, private$fine, length(B)==length(U))
      currentState <- i
      for(c in 1:length(B)) {
        coeff_cumsum <- private$P_cumsum[[currentState]][[B[c]]]
        if(length(coeff_cumsum) == 0) {
          currentState <- currentState #Stay still -> no possible moves
        } else {
          nextState_index <- findInterval(U[c], coeff_cumsum) + 1 #+1 cause they start from 0
          nextState_possibilities <- private$P_moves_list[[currentState]][[B[c]]]
          if(nextState_index > length(nextState_possibilities)) {
            currentState <- currentState #Stay still -> U > coeff
          } else {
            currentState <- nextState_possibilities[nextState_index]
          }
        }
      }
      return(currentState)
    },
    update.fun = function(i,B,U) {
      return(updateFunCpp(currentState = i,B = B,U = U, connected = private$connected, fine = private$fine,
                          P_cumsum = private$P_cumsum, P_moves_list = private$P_moves_list))
    },
    update.fun.slow = function(i,B,U) {
      #Update function for the ladder using local moves (more efficient than the global bersion)
      #Work as update.fun but the implementation is less efficent
      stopifnot(private$connected, private$fine, length(B)==length(U))
      currentState <- i
      for(c in 1:length(B)) {
        #Find which move
        move <- which(private$P_moves[currentState,] == B[c])
        if(length(move) == 0) {
          currentState <- currentState #There are no ways to go given the current roll
        } else {
          move_coeff <- private$P[currentState,move] #Coefficients of the P_matrix
          nextState_index <- findInterval(U[c], cumsum(move_coeff)) + 1 #+1 cause they start from 0
          if(nextState_index > length(move)) {
            currentState <- currentState #Stay still -> U > coeff
          } else {
            currentState <- move[nextState_index]
          }
        }
      }
      return(currentState)
    },
    debug.update.fun = function(reps,time_span,roll.fun = NULL, true_p = NULL) {
      #This is useful to check if the update function works correctly (just for debug purposes)
      #We just track the chain forward for n steps and then we see if it is close
      #to the stationary distribution
      if(is.null(roll.fun) && is.null(true_p)) {stop("Either declare roll.fun or the true probabilities.")}
      if(is.null(roll.fun)) {
        stopifnot(isTRUE(all.equal(1,sum(true_p))))
        cat("True equilibrium: ",self$evaluate(true_p),"\n")
        roll.fun <- function(n) {sample(1:private$m, size = n, replace = TRUE, prob = true_p)}
      }
      res <- numeric(reps)
      for(i in 1:reps) {
        res[i] <- self$update.fun(1,roll.fun(time_span),runif(time_span))
      }
      cat("Empirical equilibrium: ",table(res)/reps,"\n")
    },
    sample.AR = function(n,roll.fun = NULL, true_p = NULL, num_cores = 1, verbose = FALSE,...) {
      #Sample using accept-reject algorithm
      if(is.na(private$fine) || !private$fine) {
        stop("Sampling is possible only for fine ladders.")
      }
      #Define rolling function
      if(is.null(roll.fun) && is.null(true_p)) {stop("Either declare roll.fun or the true probabilities.")}
      if(is.null(roll.fun)) {
        stopifnot(isTRUE(all.equal(1,sum(true_p))))
        roll.fun <- function(n) {sample(1:private$m, size = n, replace = TRUE, prob = true_p)}
      }

      #Compute constant Q = max R_i/choose(d,n) where n is any member of the discrete simplex
      discrete_simplex <- discrete.simplex(d=private$degree,m=private$m) #all possible combinations
      discrete_simplex_binomial_unique <- unique(as.numeric(apply(discrete_simplex,1,function(vec) {multinomial.coeff(d=private$degree,n=vec)})))
      Q <- max(sapply(private$R, function(r) {max(r/discrete_simplex_binomial_unique)}))
      if(verbose && !is.null(true_p)) { #Compute C_p is true_p is available
        aux <- rep(NA, nrow(private$M))
        for(i in 1:nrow(private$M)) {
          aux[i] <- prod(true_p^private$M[i,])*private$R[i]
        }
        C_p <- sum(aux)
      } else {
        C_p <- NA
      }

      #Sample using A-R
      res_AR <- mclapply(1:n, function(rep) {
        num_rolls <- 0
        num_iter <- 0
        while(TRUE) {
          num_iter <- num_iter + 1
          #Roll the die d times
          toss_dice <- roll.fun(n = private$degree)
          num_rolls <- num_rolls + private$degree
          #Construct vector with the result
          result_dice <- numeric(private$m)
          for(i in 1:private$m) {
            result_dice[i] <- length(which(toss_dice == i))
          }
          #Find which R corresponds to the obtained result
          R_res <- 0 #if there is no row corresponding -> R is 0
          for(i in 1:nrow(private$M)) {
            if(isTRUE(base::all.equal(private$M[i,],result_dice, check.attributes = FALSE))) {
              R_res <- private$R[i]
              categ_output <- i #index of the category of the output
              break
            }
          }
          prob_accept <- R_res/(multinomial.coeff(d=private$degree, n=result_dice)*Q)
          if(prob_accept > 1) { stop("Something's wrong with the algorithm.")} #Debug
          accept <- sample(1:2, size = 1, prob = c(prob_accept, 1-prob_accept))
          if(accept == 1) { #Accept drawn point
            if(verbose) {
              return(list(res = categ_output, rolls = num_rolls, iter = num_iter))
            } else {
              return(categ_output)
            }
          }
        }
      }, mc.cores = num_cores)

      res_AR_sample <- unlist(lapply(res_AR, function(x) {x[[1]]}))
      if(verbose) {
        res_AR_rolls <- unlist(lapply(res_AR, function(x) {x[[2]]}))
        res_AR_iter <- unlist(lapply(res_AR, function(x) {x[[3]]}))
      }

      #Produce output
      if(verbose) {
        return(list(res = res_AR_sample,
                    empirical_tosses = res_AR_rolls,
                    empirical_iter = res_AR_iter,
                    theor_iter = Q/C_p,
                    theor_tosses = private$degree*Q/C_p,
                    C_p = C_p,
                    Q = Q))
      } else {
        return(res_AR_sample)
      }

    },
    sample = function(n,roll.fun = NULL, true_p = NULL, num_cores = 1, verbose = FALSE, global = FALSE, double_time = FALSE,...) {
      #Get a sample from the ladder using CFTP
      #If global = TRUE, uses a different update function that makes use of a global constant -> less efficient!
      if(is.null(roll.fun) && is.null(true_p)) {stop("Either declare roll.fun or the true probabilities.")}
      if(is.null(roll.fun)) {
        stopifnot(isTRUE(all.equal(1,sum(true_p))))
        roll.fun <- function(n) {sample(1:private$m, size = n, replace = TRUE, prob = true_p)}
      }
      if(private$m > 2) {
        monotonic_CFTP <- FALSE
      } else {
        monotonic_CFTP <- TRUE #Univariate case -> monotonic implementation is possible!
      }
      res <- mclapply(1:n, function(i) {
        if(global) {
          stop("global is not supported anymore as it is less efficient. Use global = FALSE")
          #CFTP(k = private$k, roll.fun = roll.fun, update.fun = self$update.fun.global,
          #     monotonic = monotonic_CFTP, min = private$min_state, max = private$max_state,verbose=verbose, double_time = double_time,...) #min, max are used only in monotonic case, otherwise they are ignored
        } else {
          CFTP(k = private$k, roll.fun = roll.fun, update.fun = self$update.fun,
               monotonic = monotonic_CFTP, min = private$min_state, max = private$max_state,verbose=verbose, double_time = double_time,...) #min, max are used only in monotonic case, otherwise they are ignored
        }
      }, mc.cores = num_cores)

      if(verbose) {
        return(list(unlist(lapply(res, function(x) {x[[1]]})), exp_rolls = (unlist(lapply(res, function(x) {x[[2]]})))))
      } else {
        return(unlist(res))
      }

    },
    evaluate = function(p) {
      #Return the values of the ladder for a fixed
      #vector of probabilities of rolling each face
      stopifnot(is.numeric(p), length(p) == private$m, isTRUE(all.equal(1,sum(p))))
      aux <- rep(NA, private$k)
      for(i in 1:private$k) {
        aux[i] <- prod(p^private$M[i,])*private$R[i]
      }
      return(aux/sum(aux))
    },
    increase.degree = function(d = 1) {
      #d = of how many degrees the ladder has to be increased
      if(!is.wholenumber(d) || d < 1) {
        stop("d must be an integer greater than 1.")
      }
      discrete_simplex <- discrete.simplex(d=d,m=private$m) #all possible combinations
      new_M_list <- vector("list", private$k)
      for(i in 1:private$k) {
        new_M_list[[i]] <- data.frame(t(apply(discrete_simplex, 1, function(vec) {vec+private$M[i,]})))
      }
      #Construct a ladder with the new M
      new_M <- as.matrix(rbindlist(new_M_list))
      #Compute the new vector R
      discrete_simplex_binomial <- as.numeric(apply(discrete_simplex,1,function(vec) {multinomial.coeff(d=d,n=vec)}))
      new_R <- rep(private$R, each=length(discrete_simplex_binomial))*discrete_simplex_binomial
      A <- split(1:nrow(new_M), rep(1:private$k, each=length(discrete_simplex_binomial)))
      return(list(obj = Ladder$new(M=new_M, R=new_R), A=A))
    },
    impose.fineness = function() {
      if(private$fine) {
        return(list(obj = self$clone(), A = lapply(1:private$k, function(i) {i}))) #The ladder is already fine -> returns a copy of it
      }
      new_k <- nrow(unique(private$M))
      A <- vector("list", new_k) #initialize partition
      #WRONG: indices <- match(data.frame(t(private$M)), data.frame(unique(t(private$M)))) #It may have gaps!
      indices <- indices.unique.mat(private$M) #Find unique rows of M and assign them an index.
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
    impose.connected = function(increase_degree=NULL) {
      if(!is.na(private$connected) && private$connected && is.null(increase_degree))  {
        return(list(obj = self$clone(), A = lapply(1:private$k, function(i) {i}))) #The ladder is already connected -> returns a copy of it
      }
      #We extend 1 degree at the time until the ladder is connected.
      #NOT EFFICIENT: we could precompute the minimum degree necessary
      #to have a connected ladder
      #If increase_degree is specified, then we just increase of d degree
      #Notice that it may not be connected, in case an error is return
      if(!is.null(increase_degree)) {
        test_ladder_list <- self$increase.degree(d=increase_degree) #Create new ladder that is a disaggregation and has degree increased by increase_degree
        #Check if it is connected
        if(test_ladder_list$obj$get.connected()) {
          return(test_ladder_list)
        }
      } else {
        for(d in 1:private$degree) {
          test_ladder_list <- self$increase.degree(d=d) #Create new ladder that is a disaggregation and has degree increased by 1
          #Check if it is connected
          if(test_ladder_list$obj$get.connected()) {
            return(test_ladder_list)
          }
        }
      }
      stop("Impossible to create a connected ladder.") #Should never happen.
    },
    get.a = function() {private$a},
    get.P = function()  {private$P},
    get.P.moves = function() {private$P_moves},
    get.P.cumsum = function() {private$P_cumsum},
    get.P.moves.list = function() {private$P_moves_list},
    get.M = function() {private$M},
    get.degree = function() {private$degree}
  ),
  private = list(
    #FIELDS
    m = NULL, #Number of faces of the given die
    M = NULL, #k x m matrix describing the power of each p_i. The rows must all sum to the same value.
    degree = NULL, #degree of the ladder
    valid = FALSE, #FALSE if it is not a valid ladder. It is checked after construction.
    fine = NA,
    connected = NA,
    min_state = NA, #Minimum state (used for monotonic CFTP when m=2)
    max_state = NA, #Maximum state (used for monotonic CFTP when m=2)
    neigh = NULL, #neighbourhood. It is a double indexed list, the first index is the state, the second the roll of the die.
    a = NA, #constant for the Markov chain. It is used in global update function (DEPRECATED)
    P = NA, #Transition matrix of the Markov chain. It just contains the cooefficients without the p_i's
    P_moves = NA, #Moves in the transition matrix
    P_cumsum = NA, #List double indexed. The first index is the current state, the second is the current roll. Returns the cumsum of the coefficients
    P_moves_list = NA, #List double index. The first index is the current state, the second is the current roll. Returns the index of the possible moves
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
    },
    compute.transition.matrix = function() {
      #Fill the entries of the transition matrix.
      #P_ij is equal to R_j/max(sum h in neighbour of j given the roll sum R_h,
      #sum h in neighbour of i given the roll and that contains j R_h)
      private$P <- matrix(0, nrow = private$k, ncol = private$k) #Initialize
      private$P_moves <- matrix(0, nrow = private$k, ncol = private$k)
      for(i in 1:private$k) {
        for(b in 1:private$m) {
          neigh <- private$neigh[[i]][[b]]
          private$P_moves[i,neigh] <- b
          for(j in neigh) {
            num <- private$R[j] #numerator
            #Find in which neighbour of j there is i
            for(b_aux in 1:private$m) {
              if(i %in% private$neigh[[j]][[b_aux]]) {
                den_aux <- sum(private$R[private$neigh[[j]][[b_aux]]])
                break
              }
            }
            den <- max(den_aux,sum(private$R[neigh]))
            private$P[i,j] <- num/den
          }
        }
      }
      diag(private$P) <- rep(NA,private$k)
      #Create the list of the cumulative sums to speed up the update function
      private$P_cumsum <- vector("list", length = private$k)
      private$P_moves_list <- vector("list", length = private$k)
      for(i in 1:private$k) {
        private$P_cumsum[[i]] <- vector("list", length = private$m)
        private$P_moves_list[[i]] <- vector("list", length = private$m)
        for(b in 1:private$m) {
          move <- which(private$P_moves[i,] == b) #Find where it can move from i having rolled b
          if(length(move) == 0) {
            private$P_cumsum[[i]][[b]] <- numeric(0)
            private$P_moves_list[[i]][[b]] <- numeric(0)
          } else {
            move_coeff <- private$P[i,move]
            private$P_cumsum[[i]][[b]] <- cumsum(move_coeff)
            private$P_moves_list[[i]][[b]] <- move
          }
        }
      }
    }
  ))
