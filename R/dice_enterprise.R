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
    initialize = function(G, d = NULL, verbose = FALSE) {
      #G is a list of lists (G_1,G_2,...,G_v). Each list must contain two element:
      #1) a vector of coefficients (can not be null)
      #2) a matrix of powers of p_1p_2...p_m OR a vector of strings [THIS CAN WORK ONLY FOR POWER UP TO 9!].
      #for instance p_1^2p_2,p_1p_3^3 becomes either M = [2,1,0;1,0,3] or ("210","103")
      private$G_orig <- G
      if(verbose) {cat("Converting polynomials in internal representation...")}
      private$generate.G_poly(G) #Generate a list containing the polynomials
      if(verbose) {cat("DONE! \nGenerating initial ladder...")}
      private$generate.ladder_initial() #Generate the initial ladder
      if(verbose) {cat("DONE! \nConstructing a connected ladder...")}
      aux <- private$ladder_initial$impose.connected(increase_degree = d)
      private$ladder_connected <- aux[["obj"]]
      private$A_initial_connected <- aux[["A"]]
      if(verbose) {cat("DONE! \nConstructing a fine and connected ladder...")}
      aux <- private$ladder_connected$impose.fineness()
      private$ladder_fine_connected <- aux[["obj"]]
      private$A_fine_connected <- aux[["A"]]
      private$check.efficiency.condition() #Check if efficiency condition is satisfied
      if(verbose) {cat("DONE! \n")}
    },
    print = function() {
      cat("Original function: Delta^",private$m," -> Delta^",private$v,"\n",sep = "")
      cat("Fine and connected ladder: Delta^",private$m," -> Delta^",private$ladder_fine_connected$get.k(),"\n", sep="")
      cat("Efficiency condition is satisfied: ",private$efficiency_condition,"\n", sep ="")
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
    get.ladder.fine.connected = function() {private$ladder_fine_connected$clone()},
    get.efficiency.condition = function() {private$efficiency_condition},
    expected.tosses.bound = function(true_p,threshold = 1e6) {
      #This method computed the expected tosses required by CFTP for a given p.
      #Works only when m = 2.
      #true_p = (p,1-p)
      if(private$m > 2) {
        stop("The expected number of tosses is computable only when the original die has 2 faces.")
      }
      stopifnot(is.numeric(true_p), length(true_p) == private$m, isTRUE(all.equal(1,sum(true_p))))
      #The expected number of tosses is computed by using a telescopic sum
      #E[D_t^(1,k)] = sum E[D_t^(n,n+1)]
      #where E[D_t^(i,j)] is the expected distance at time between two coupled
      #particles started at i and j.
      #We then use Markov property to compute the distance at time t+1 given the
      #expected distance at time t.
      p <- true_p[1]
      k <- private$ladder_fine_connected$get.k()
      R <- private$ladder_fine_connected$get.R()
      exp_dist_old <- rep(1,k-1) #Distances between particles at time 0
      exp_dist_new <- exp_dist_old #Initialization
      t <- 1
      while(TRUE) {
        if(t > threshold) { #If t gets greater than 10^6, we just interrupt
          t <- Inf
          warning(paste0("Time required is over threshold (",threshold,")"))
          break
        } else if(sum(exp_dist_new) < exp(-1)) { #The expected distance is smaller than 1,
          #so coalescence has occured
          break
        }
        #Compute the expected distance at time t given the distances at previous time
        for(i in 1:(k-1)) {
          #Notice that some of the following may be NA/NULL/0/weird
          #as we exceed the entries of R
          a <- R[i+1]/max(R[i],R[i+1])
          b <- R[i-1]/max(R[i-1],R[i])
          c <- R[i+2]/max(R[i+1],R[i+2])
          d <- R[i]/max(R[i],R[i+1])
          if(i == 1 && i+1 == k) { #There are only two states
            b <- NA; c <- NA;
          } else if(i == 1) {
            b <- NA;
          } else if(i == k-1) {
            c <- NA;
          }
          #Compute the expected distance at time t between particles started at i and i+1
          exp_dist_new[i] <- 0
          if(!is.na(c) && !is.na(a) && i+1 < k) {
            exp_dist_new[i] <- exp_dist_new[i] + max((c-a)*p*(exp_dist_old[i]+exp_dist_old[i+1]),0)
          }
          if(!is.na(b) && !is.na(d) && i-1 >= 1) {
            exp_dist_new[i] <- exp_dist_new[i] + max((b-d)*(1-p)*(exp_dist_old[i-1]+exp_dist_old[i]),0)
          }
          if(i+1 < k) {
            exp_dist_new[i] <- exp_dist_new[i] + min(a,c,na.rm = TRUE)*p*exp_dist_old[i+1]
          }
          if(i-1 >= 1) {
            exp_dist_new[i] <- exp_dist_new[i] + min(b,d,na.rm = TRUE)*(1-p)*exp_dist_old[i-1]
          }
          exp_dist_new[i] <- exp_dist_new[i] + min(1-a,1-c,na.rm = TRUE)*p*exp_dist_old[i] +
            min(1-b,1-d,na.rm = TRUE)*(1-p)*exp_dist_old[i]
        }
        exp_dist_old <- exp_dist_new #Update
        t <- t+1
      }
      return(19.2*t)
    },
    increase.degree = function(d = 1,verbose=FALSE) {
      #We increase the degree of the dice enterprise by just applying the
      #multinomial theorem and create a new Dice Enterprise
      increase_degree <- private$ladder_fine_connected$get.degree() - private$ladder_initial$get.degree() + d
      incresed_de <- DiceEnterprise$new(private$G_orig,d=increase_degree,verbose = verbose)
      return(incresed_de)
    }
  ),
  private = list(
    m = NA, #faces of the given die
    v = NA, #faces of the returned die
    d = NA, #degree of the polynomial
    n_coeff = NA, #number of different terms in all the polynomials,
    G_orig = NULL, #f(p) as defined by the user. Still need to be converted in internal representation (see G_poly)
    G_poly = NULL, #polynomials
    ladder_initial = NULL, #Initial ladder obtained from the rational function
    ladder_connected = NULL, #Imposing connected condition
    ladder_fine_connected = NULL, #Imposing fine condition
    A_f_initial = NULL, #Disaggregation of f into initial ladder
    A_initial_connected = NULL, #Disaggregation of initial ladder in a connected one
    A_fine_connected = NULL, #Disagreggation of fine&connected ladder in a connected ladder
    efficiency_condition = NA, #Efficiency of the dice enterprise. If TRUE the expected time for coalescence is bounded exponentially in time for every p. Available only for fine and connected.
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
        if(any(R_list[[i]] < 0)) {stop("Generated negative coefficient. Are you sure it is a valid function? Check that numerator and
                                       denominator do not share any common roots. Also, check that they are both positive polynomials
                                       (in case, switch signs of the coefficients accordingly).")}
        #Remove zero coefficients
        if(any(R_list[[i]] == 0)) { #Need to check or M_list[[i]][numeric(0),] is empty!
          M_list[[i]] <- matrix(M_list[[i]][-which(R_list[[i]] == 0),], ncol = ncol(M_list[[i]]), byrow = FALSE) #NOT GREAT, but it should work cause we rounded up the values of R
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
    },
    check.efficiency.condition = function() {
      #This method checks if an efficiency condition is satisfied,
      #that is if the expected number of tosses is bounded exponentially
      #in time
      if(private$m == 2) {
        private$efficiency_condition <- TRUE
        R <- private$ladder_fine_connected$get.R()
        M <- private$ladder_fine_connected$get.M() #matrix of coefficients
        R <- R[order(private$ladder_fine_connected$get.M()[,1]+1)]  #Reorder the coefficients in the right order starting from power 0 to power k-1
        k <- private$ladder_fine_connected$get.k()
        for(i in 1:k) {
          if(i+2 <= k) {
            if(R[i+1]/max(R[i],R[i+1]) < R[i+2]/max(R[i+1],R[i+2])) {
              private$efficiency_condition <- FALSE
              break
            }
          }
          if(i-1 >= 1 && i+1 <= k) {
            if(R[i]/max(R[i],R[i+1]) < R[i-1]/max(R[i-1],R[i])) {
              private$efficiency_condition <- FALSE
              break
            }
          }
        }
      } else { #Not supported for now in the not monotonic case
        #TO BE STUDIED AND FILLED
      }
    }
  )
)

CoinsEnterprise <- R6::R6Class("CoinsEnterprise",
                               inherit = DiceEnterprise,
                               public = list(
                                 initialize = function(G,toss.coins,num_coins,die_type = c("uniform","toss_all","first_heads"), d = NULL, verbose = FALSE) {
                                   private$num_coins <- num_coins
                                   private$die_type <- match.arg(die_type)
                                   private$toss.coins.fun <- toss.coins
                                   super$initialize(G=G,d=d,verbose=verbose)
                                 },
                                 roll.die = function(n,...) {
                                   #If die_type = "uniform"
                                   #This function constructs an m+1 sided die given m independent coins.
                                   #It selects uniformly which coin to toss, say the ith, and returns
                                   #i if the result is heads, and m+1 otherwise.
                                   #A function f(p) of the probabilities of the independent coins can be converted in this
                                   #representation by substituing p_i = mq_i
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
                                     if(private$die_type == "uniform") {
                                       which_coin <- sample(1:private$num_coins, size = 1) #Select which coin to flip
                                       if(private$toss.coins.fun(which_coin,...) == 1) {
                                         res[i] <- which_coin
                                       } else {
                                         res[i] <- private$num_coins+1
                                       }
                                     } else if(private$die_type == "toss_all") {
                                       toss_res <- private$toss.coins.fun(...) #Toss the coins
                                       m <- length(toss_res)
                                       if(isTRUE(all.equal(toss_res,rep(1,m)))) {
                                         res[i] <- 1
                                       } else if(length(which(toss_res == 2)) > 1) {
                                         res[i] <- m+2
                                       } else { #Only one tails
                                         res[i] <- which(toss_res == 2)+1
                                       }
                                     } else if(private$die_type == "first_heads") {
                                       toss_res <- private$toss.coins.fun(...) #Toss the coins
                                       m <- length(toss_res)
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
                                 },
                                 increase.degree = function(d = 1,verbose=FALSE) {
                                   increase_degree <- private$ladder_fine_connected$get.degree() - private$ladder_initial$get.degree() + d
                                   incresed_de <- CoinsEnterprise$new(private$G_orig,private$toss.coins.fun,
                                                                      private$num_coins, private$die_type,
                                                                      d=increase_degree,verbose = verbose)
                                   return(incresed_de)
                                 }
                               ),
                               private = list(
                                die_type = NULL,
                                toss.coins.fun = NULL,
                                num_coins = NULL
                               ))
