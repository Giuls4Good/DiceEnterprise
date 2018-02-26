BernoulliFactory <- R6::R6Class("BernoulliFactory",
  inherit = DiceEnterprise,
  public = list(
    initialize = function(f_1,f_2,verbose = FALSE) {
      #f_1 and f_2 are lists, each containing
      #the coefficients of each term (as vector)
      #and a vector of powers of p for f(p) and 1-f(p).
      stopifnot(is.list(f_1), is.list(f_2), length(f_1) == 2, length(f_2) == 2,
                is.numeric(f_1[[1]]), is.numeric(f_1[[2]]), length(f_1[[1]]) == length(f_1[[2]]),
                all(f_1[[1]] != 0), all(is.wholenumber(f_1[[2]])),
                is.numeric(f_2[[1]]), is.numeric(f_2[[2]]), length(f_2[[1]]) == length(f_2[[2]]),
                all(f_2[[1]] != 0), all(is.wholenumber(f_2[[2]]))
                )
      if(verbose) {cat("Converting into a function over simplices...")}
      G <- vector("list", length = 2)
      #f(p)
      G[[1]] <- vector("list", length = 2)
      coeff <- f_1[[1]]
      power <- f_1[[2]]
      G[[1]][[1]] <- coeff
      G[[1]][[2]] <- matrix(c(power,rep(0,length(power))), ncol = 2, byrow = FALSE)
      #1-f(p) = den - num
      G[[2]] <- vector("list", length = 2)
      coeff <- f_2[[1]]
      power <- f_2[[2]]
      G[[2]][[1]] <- coeff
      G[[2]][[2]] <- matrix(c(power,rep(0,length(power))), ncol = 2, byrow = FALSE)
      if(verbose) {cat("DONE!\n")}
      super$initialize(G=G, verbose = verbose)
    },
    print = function() {
      cat("Bernoulli Factory\n")
      cat("Fine and connected ladder: Delta^",private$m," -> Delta^",private$ladder_fine_connected$get.k(),"\n", sep="")
    },
    evaluate = function(p) {
      super$evaluate(c(p,1-p))
    }
  ),
  private = list(

  )
)

DiceEnterprise <- R6::R6Class("DiceEnterprise",
  public = list(
    initialize = function(G, verbose = FALSE) {
      #G is a list of lists (G_1,G_2,...,G_v). Each list must contain two element:
      #1) a vector of coefficients (can not be null)
      #2) a matrix of powers of p_1p_2...p_m OR a vector of strings [THIS CAN WORK ONLY FOR POWER UP TO 9!].
      #for instance p_1^2p_2,p_1p_3^3 becomes either M = [2,1,0;1,0,3] or ("210","103")
      if(verbose) {cat("Converting polynomials in internal representation...")}
      private$generate.G_poly(G) #Generate a list containing the polynomials
      if(verbose) {cat("DONE! \nGenerating initial ladder...")}
      private$generate.ladder_initial() #Generate the initial ladder
      if(verbose) {cat("DONE! \nConstructing a connected ladder...")}
      aux <- private$ladder_initial$impose.connected()
      private$ladder_connected <- aux[["obj"]]
      private$A_initial_connected <- aux[["A"]]
      if(verbose) {cat("DONE! \nConstructing a fine and connected ladder...")}
      aux <- private$ladder_connected$impose.fineness()
      private$ladder_fine_connected <- aux[["obj"]]
      private$A_fine_connected <- aux[["A"]]
      if(verbose) {cat("DONE! \n")}
    },
    print = function() {
      cat("Original function: Delta^",private$m," -> Delta^",private$v,"\n",sep = "")
      cat("Fine and connected ladder: Delta^",private$m," -> Delta^",private$ladder_fine_connected$get.k(),"\n", sep="")
    },
    evaluate = function(p) {
      #This function evaluates the polynomials given fixed probabilites.
      #Useful for testing.
      stopifnot(is.numeric(p), length(p) == private$m, isTRUE(all.equal(1,sum(p))))
      p_out <- sapply(private$G_poly, function(G) {
        G_coeff <- G[[1]] #Coefficients
        G_power <- G[[2]] #Power of p_i
        res <- 0
        for(i in 1:nrow(G_power)) {
          res <- res + prod(p^G_power[i,])*G_coeff[i]
        }
        return(res)
      })
      return(p_out/sum(p_out))
    },
    sample = function(n,roll.fun = NULL, true_p = NULL, num_cores = 1, verbose = FALSE, global = FALSE, double_time = FALSE,...) {
      stopifnot(!is.null(private$ladder_fine_connected),
                !is.null(private$ladder_initial),
                !is.null(private$ladder_connected),
                is.list(private$A_f_initial),
                is.list(private$A_initial_connected),
                is.list(private$A_fine_connected))
      #Get a sample from the fine and connected ladder
      sample_res <- private$ladder_fine_connected$sample(n = n, roll.fun = roll.fun, true_p = true_p,
                                                                    num_cores = num_cores, verbose = verbose,global = global, double_time = double_time,...)
      if(verbose) {
        sample_fine_connected <- sample_res[[1]]
      } else {
        sample_fine_connected <- sample_res
      }
      #Transform sample to the connected ladder
      sample_connected <- disaggregation.sample(sample_fine_connected, origin = "original",
                                               disaggDist = private$ladder_connected,
                                               A = private$A_fine_connected)
      #Transform sample to the initial ladder
      sample_initial <- disaggregation.sample(sample_connected, origin = "disaggregation",
                                             disaggDist = private$ladder_connected,
                                             A = private$A_initial_connected)
      #Transform sample to initial function
      sample_f <- disaggregation.sample(sample_initial, origin = "disaggregation",
                                       disaggDist = private$ladder_initial,
                                       A = private$A_f_initial)
      if(verbose) {
        return(list(sample_f, exp_rolls = sample_res[[2]]))
      } else {
        return(sample_f)
      }

    },
    get.G.poly = function() {private$G_poly},
    get.ladder.initial = function() {private$ladder_initial$clone()},
    get.ladder.fine.connected = function() {private$ladder_fine_connected$clone()}
  ),
  private = list(
    m = NA, #faces of the given die
    v = NA, #faces of the returned die
    d = NA, #degree of the polynomial
    n_coeff = NA, #number of different terms in all the polynomials
    G_poly = NULL, #polynomials
    ladder_initial = NULL, #Initial ladder obtained from the rational function
    ladder_connected = NULL, #Imposing connected condition
    ladder_fine_connected = NULL, #Imposing fine condition
    A_f_initial = NULL, #Disaggregation of f into initial ladder
    A_initial_connected = NULL, #Disaggregation of initial ladder in a connected one
    A_fine_connected = NULL, #Disagreggation of fine&connected ladder in a connected ladder
    generate.G_poly = function(G) {
      stopifnot(is.list(G),
                is.list(G[[1]]))
      private$v <- length(G)
      private$d <- 0
      private$n_coeff <- 0
      G_new <- vector("list", private$v)
      for(i in 1:private$v) {
        stopifnot(is.numeric(G[[i]][[1]]),
                  all(G[[i]][[1]] != 0))
        G_new[[i]] <- vector("list", 2)
        G_new[[i]][[1]] <- G[[i]][[1]]
        if(is.matrix(G[[i]][[2]])) {
          #The user provided a matrix
          if(i == 1) {
            private$m <- ncol(G[[i]][[2]])
          } else if(ncol(G[[i]][[2]]) != private$m) {
            stop("The matrices of power must all have the same number of columns")
          } else if(nrow(G[[i]][[2]]) != length(G[[i]][[1]])) {
            stop("The matrix of power must have the same number of rows as the length of the coefficients.")
          }
          G_new[[i]][[2]] <- G[[i]][[2]]
        } else if(is.character(G[[i]][[2]])) {
          #The user provided a vector of string
          if(i == 1 && length(unique(nchar(G[[i]][[2]]))) == 1) {
            private$m <- unique(nchar(G[[i]][[2]]))
          } else if(length(unique(nchar(G[[i]][[2]]))) != 1 || unique(nchar(G[[i]][[2]])) != private$m) {
            stop("The string of powers must all have the same length.")
          } else if(length(G[[i]][[2]]) != length(G[[i]][[1]])) {
            stop("The string of power must have the same number of elements as the length of the coefficients.")
          }
          G_new[[i]][[2]] <- matrix(as.numeric(unlist(strsplit(G[[i]][[2]],""))), byrow = TRUE, ncol = private$m)
        } else {stop("Provide either a matrix of power or a vector of strings.")}
        #Update number of coefficients
        private$n_coeff <- private$n_coeff + length(G_new[[i]][[1]])
        #Save the degree of the polynomial
        G_new[[i]][[3]] <- rowSums(G_new[[i]][[2]]) #Save the degree of each coefficients
        private$d <- max(private$d,max(rowSums(G_new[[i]][[2]])))
      }
      private$G_poly <- G_new #Save result
    },
    generate.ladder_initial = function() {
      #Using the Theorem, convert the polynomial to the requested form
      #so that a ladder is generated which is a disaggregation of f.
      M_list <- vector("list",private$v) #decompose each polynomial in the required powers
      R_list <- vector("list",private$v)
      private$A_f_initial <- vector("list", private$v)
      M_rows <- 0 #counting the number of rows of the final ladder
      counter <- 0
      for(i in 1:private$v) {
        for(j in 1:length(private$G_poly[[i]][[3]])) { #For each coefficient
          deg <- private$G_poly[[i]][[3]][j] #degree of the a_n
          if(private$d > deg) {
            M_list[[i]] <- rbind(M_list[[i]], t(apply(discrete.simplex(private$d-deg,private$m), 1, function(x) {
              x + private$G_poly[[i]][[2]][j,]
            }))) #NOT EFFICIENT
            R_list[[i]] <- c(R_list[[i]], apply(discrete.simplex(private$d-deg,private$m), 1, function(x) {
              private$G_poly[[i]][[1]][j]*multinomial.coeff(private$d-deg,x)
            })) #NOT EFFICIENT
          } else {
            M_list[[i]] <- rbind(M_list[[i]], private$G_poly[[i]][[2]][j,]) #NOT EFFICIENT
            R_list[[i]] <- c(R_list[[i]], private$G_poly[[i]][[1]][j]) #NOT EFFICIENT
          }
        }
        #Sum all the coefficients in R that correspond to the same powers
        indices <- indices.unique.mat(M_list[[i]])
        M_list[[i]] <- unique(M_list[[i]])
        R_aux <- rep(NA, length = nrow(M_list[[i]]))
        for(l in 1:nrow(M_list[[i]])) {
          R_aux[l] <- round(sum(R_list[[i]][which(indices == l)]), digits = 8) #round to avoid machine error
        }
        R_list[[i]] <- R_aux
        if(any(R_list[[i]] < 0)) {stop("Generated negative coefficient. Are you sure it is a valid function?")}
        #Remove zero coefficients
        if(any(R_list[[i]] == 0)) { #Need to check or M_list[[i]][numeric(0),] is empty!
          M_list[[i]] <- M_list[[i]][-which(R_list[[i]] == 0),] #NOT GREAT, but it should work cause we rounded up the values of R
          R_list[[i]] <- R_list[[i]][-which(R_list[[i]] == 0)]
        }
        M_rows <- M_rows + nrow(M_list[[i]])
        private$A_f_initial[[i]] <- (counter+1):(counter+nrow(M_list[[i]]))
        counter <- counter+nrow(M_list[[i]])
      }
      #Construct ladder
      M_new <- do.call(rbind, M_list)
      R_new <- do.call(c, R_list)
      # M_new <- do.call(rbind, sapply(M_list, unlist))
      # R_new <- do.call(c, sapply(R_list, unlist))
      #Remove zero coefficients

      private$ladder_initial <- Ladder$new(M = M_new, R=R_new)
    }
  )
)

CoinsEnterprise <- R6::R6Class("CoinsEnterprise",
                               inherit = DiceEnterprise,
                               public = list(
                                 initialize = function(G,toss.coins,die_type = c("toss_all","first_heads"), verbose = FALSE) {
                                   private$die_type <- match.arg(die_type)
                                   private$toss.coins.fun <- toss.coins
                                   super$initialize(G=G,verbose=verbose)
                                 },
                                 roll.die = function(n,...) {
                                   #If die_type = "toss_all"
                                   #This function constructs an m+2 sided die given m independent coins.
                                   #It returns 1 if all the tosses are heads
                                   #It returns i+1 if all the tosses are heads, except for the ith coin
                                   #It returns m+2 otherwise
                                   #A function f(p) of the probabilities of the independent coins can be converted in this
                                   #representation by substituing p_i = q_1/(q_1+q_{i+1})
                                   #If die_type = "first_heads"
                                   #This function constructs an m+1 sided die given m independent coins.
                                   #It returns the position of the first heads and returns m+1 if all tails
                                   #A function f(p) of the probabilities of the independent coins can be converted in this
                                   #representation by substituting p_i = q_i/(1-sum_{k=1}^{i-1} q_k)
                                   res <- numeric(n)
                                   for(i in 1:n) {
                                     toss_res <- private$toss.coins.fun(...) #Toss the m coins
                                     m <- length(toss_res)
                                     if(private$die_type == "toss_all") {
                                       if(isTRUE(all.equal(toss_res,rep(1,m)))) {
                                         res[i] <- 1
                                       } else if(length(which(toss_res == 2)) > 1) {
                                         res[i] <- m+2
                                       } else { #Only one tails
                                         res[i] <- which(toss_res == 2)+1
                                       }
                                     } else if(private$die_type == "first_heads") {
                                       if(any(toss_res == 1, na.rm = TRUE)) {
                                         res[i] <- which(toss_res == 1)[1]
                                       } else {
                                         res[i] <- m+1 #All tails
                                       }
                                     } else {
                                       stop("Unknown specified type of die.")
                                     }
                                   }
                                   return(res)
                                 },
                                 sample = function(n, num_cores = 1, verbose = FALSE, global = FALSE, double_time = FALSE,...) {
                                   super$sample(n=n,roll.fun = self$roll.die, num_cores = num_cores, verbose = verbose, global = global,double_time = double_time,...)
                                 }
                               ),
                               private = list(
                                die_type = NULL,
                                toss.coins.fun = NULL
                               ))
