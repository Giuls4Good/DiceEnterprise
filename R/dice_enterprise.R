BernoulliFactory <- R6::R6Class("BernoulliFactory",
  inherit = DiceEnterprise,
  public = list(
    initialize = function(f_1,f_2,verbose = FALSE, threshold = 1e2) {
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
      super$initialize(G=G, verbose = verbose, threshold = threshold)
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
    initialize = function(G, d = NULL, verbose = FALSE, threshold = 1e2) {
      #G is a list of lists (G_1,G_2,...,G_v). Each list must contain two element:
      #1) a vector of coefficients (can not be null)
      #2) a matrix of powers of p_1p_2...p_m OR a vector of strings [THIS CAN WORK ONLY FOR POWER UP TO 9!].
      #for instance p_1^2p_2,p_1p_3^3 becomes either M = [2,1,0;1,0,3] or ("210","103")
      private$G_orig <- G
      if(verbose) {cat("Converting polynomials in internal representation...")}
      private$generate.G_poly(G) #Generate a list containing the polynomials
      if(verbose) {cat("DONE! \nGenerating initial ladder...")}
      private$generate.ladder_initial(threshold = threshold) #Generate the initial ladder
      if(verbose) {cat("DONE! \nConstructing a connected ladder...")}
      aux <- private$ladder_initial$impose.connected(increase_degree = d)
      private$ladder_connected <- aux[["obj"]]
      private$A_initial_connected <- aux[["A"]]
      if(verbose) {cat("DONE! \nConstructing a fine and connected ladder...")}
      aux <- private$ladder_connected$impose.fineness()
      private$ladder_fine_connected <- aux[["obj"]]
      private$A_fine_connected <- aux[["A"]]
      private$check.efficiency.condition() #Check if efficiency condition is satisfied
      private$d <- private$ladder_fine_connected$get.degree()
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
    sample.AR = function(n,roll.fun = NULL, true_p = NULL, num_cores = 1, verbose = FALSE,...) {
      #This function produces samples from the original ladder (not necessarily connected)
      #using Acceptance-Rejection
      stopifnot(!is.null(private$ladder_initial),
                is.list(private$A_f_initial))

      #Construct fine ladder from the initial one
      aux <- private$ladder_initial$impose.fineness()
      ladder_initial_fine <- aux[["obj"]]
      A_initial_fine <- aux[["A"]]

      #Get a sample from the fine ladder
      sample_res <- ladder_initial_fine$sample.AR(n = n, roll.fun = roll.fun, true_p = true_p,
                                                         num_cores = num_cores, verbose = verbose,...)

      if(verbose) {
        sample_initial_fine <- sample_res[[1]]
      } else {
        sample_initial_fine <- sample_res
      }

      #Transform sample to initial ladder
      sample_initial <- disaggregation.sample(sample_initial_fine, origin = "original",
                                                disaggDist = private$ladder_initial,
                                                A = A_initial_fine)
      #Transform sample to initial function
      sample_f <- disaggregation.sample(sample_initial, origin = "disaggregation",
                                        disaggDist = private$ladder_initial,
                                        A = private$A_f_initial)

      if(verbose) {
        return(list(sample_f, empirical_tosses = sample_res[["empirical_tosses"]],
                    empirical_iter = sample_res[["empirical_iter"]],
                    theor_iter = sample_res[["theor_iter"]],
                    theor_tosses = sample_res[["theor_tosses"]],
                    C_p = sample_res[["C_p"]],
                    Q = sample_res[["Q"]]))
      } else {
        return(sample_f)
      }

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
    # bound.expected.tosses = function(p) {
    #   #This function returns a bound on the expected # tosses required
    #   #by CFTP when the ladder is fine and connected and monotonic
    #   #CFTP can be used (i.e. m=2 original faces of the die)
    #   #and the efficiency condition is satisfied
    #   if(private$m != 2 || !private$efficiency_condition) {
    #     stop("A bound on the expected number of tosses is available only
    #          when the efficiency condition is satisfied and the original
    #          number of faces is 2. Try to increase the degree of the ladder.")
    #   }
    #   R <- private$ladder_fine_connected$get.R()
    #   M <- private$ladder_fine_connected$get.M() #matrix of coefficients
    #   R <- R[order(private$ladder_fine_connected$get.M()[,1]+1)]
    #   k <-  private$ladder_fine_connected$get.k()
    #
    #   #Generate matrix of coefficients (with right order)
    #   P <- matrix(0, nrow = k, ncol = k)
    #   for(i in 1:k) {
    #     if(i > 1) {P[i,i-1] <- R[i-1]/max(R[i-1],R[i])}
    #     if(i < k) {P[i,i+1] <- R[i+1]/max(R[i],R[i+1])}
    #   }
    #   rho <- 0
    #   for(i in 1:(k-1)) {
    #     if(i == 1) {
    #       new_rho <- 1 - (P[i,i+1]-P[i+1,i+2])*p - (P[i+1,i]-0)*(1-p)
    #     } else if(i == k-1) {
    #       new_rho <- 1 - (P[i,i+1]-0)*p - (P[i+1,i]-P[i,i-1])*(1-p)
    #     } else {
    #       new_rho <- 1 - (P[i,i+1]-P[i+1,i+2])*p - (P[i+1,i]-P[i,i-1])*(1-p)
    #     }
    #     if(new_rho > rho) {
    #       rho <- new_rho
    #     }
    #   }
    #   return((k-1)/(1-rho))
    # },
    get.G.poly = function() {private$G_poly},
    get.ladder.initial = function() {private$ladder_initial$clone()},
    get.ladder.fine.connected = function() {private$ladder_fine_connected$clone()},
    get.efficiency.condition = function() {private$efficiency_condition},
    expected.tosses.bound.generic = function() {
      #This bound is independent on p
      warning("This bound is not theoretically valid.")
      if(private$m > 2) {
        stop("The expected number of tosses is computable only when the original die has 2 faces.")
      }
      P <- private$ladder_fine_connected$get.P()
      a <- min(P[which(P > 0)], na.rm = TRUE)
      k <-  private$ladder_fine_connected$get.k()
      return(private$m^(private$d-1)/(a^(private$d)))
    },
    expected.tosses.bound = function(true_p) {
      #This method computed the expected tosses required by CFTP for a given p.
      #Works only when m = 2 and the efficiency condition is satisfied.
      #true_p = (p,1-p)
      if(private$m > 2) {
        stop("The expected number of tosses is computable only when the original die has 2 faces.")
      }
      if(!private$efficiency_condition) {
        stop("The bound can be obtained only when the efficiency condition is satisfied. Try to increase the degree of the ladder.")
      }
      stopifnot(is.numeric(true_p), length(true_p) == private$m, isTRUE(all.equal(1,sum(true_p))))
      #The expected number of tosses is computed by using a telescopic sum
      #E[D_t^(1,k)] = sum E[D_t^(n,n+1)]
      #where E[D_t^(i,j)] is the expected distance at time between two coupled
      #particles started at i and j.
      #We then use Markov property to compute the distance at time t+1 given the
      #expected distance at time t.
      p <- true_p[1]
      R <- private$ladder_fine_connected$get.R()
      M <- private$ladder_fine_connected$get.M() #matrix of coefficients
      R <- R[order(private$ladder_fine_connected$get.M()[,1]+1)] #Reorder coefficients
      k <-  private$ladder_fine_connected$get.k()

      #Generate matrix of coefficients (with right order)
      P <- matrix(0, nrow = k, ncol = k)
      P_coeff <- matrix(0, nrow = k, ncol = k) #just coefficients without p
      for(i in 1:k) {
        if(i > 1) {
          P[i,i-1] <- (1-p)*R[i-1]/max(R[i-1],R[i])
          P_coeff[i,i-1] <- R[i-1]/max(R[i-1],R[i])
        }
        if(i < k) {
          P[i,i+1] <- p*R[i+1]/max(R[i],R[i+1])
          P_coeff[i,i+1] <- R[i+1]/max(R[i],R[i+1])
        }
      }

      #Construct tridiagonal matrix
      if(k == 2) {
        #there are only two states
        A <- matrix(1-P[1,2]-P[2,1], nrow=1, ncol = 1)
      } else {
        A <- diag(sapply(1:(k-1), function(i) {(1-P[i,i+1]-P[i+1,i])}))
        #Fill upper off-diagonal
        A[row(A) == (col(A) - 1)] <- sapply(1:(k-2), function(i) {P[i+1,i+2]})
        #Fill lower off-diagonal
        A[row(A) == (col(A) + 1)] <- sapply(2:(k-1), function(i) {P[i,i-1]})
      }

      #Spectral decomposition of the matrix (to compute exponentiation
      # more quickly)
      # eigen_decomp <- eigen(A)
      # Q <- Re(eigen_decomp$vectors)
      # Qinv <- solve(Q)
      # D <- Re(eigen_decomp$values)
      # exp_tosses <- 0
      # new_exp_tosses <- 1
      # t <- 1
      # while(new_exp_tosses - exp_tosses > 1e-4) {
      #   if(t > 1) {exp_tosses <- new_exp_tosses}
      #   new_exp_tosses <- exp_tosses + sum(Q%*%diag(D^t, nrow = length(D))%*%Qinv)
      #   t <- t+1
      # }
      new_exp_tosses <- expected_tosses_bound_cpp(A) #C++ implementation

      return(new_exp_tosses)

      # OLD:
      # exp_distance <- numeric(k-1) #Vector of expected distances between coupled particles
      # new_exp_distance <- numeric(k-1)
      # new_exp_distance2 <- numeric(k-1)
      # exp_tosses <- 0
      # new_exp_tosses <- Inf
      # t <- 1
      #
      # while(new_exp_tosses - exp_tosses > 1e-6) {
      #   if(t > 1) {exp_tosses <- new_exp_tosses}
      #   for(i in 1:(k-1)) {
      #     if(t == 1) {
      #       #Expected distance after first time step
      #
      #       # THIS IS VALID ONLY IF THE EFFICIENCY CONDITION IS SATISFIED:
      #       # if(i == 1 && k >= 3) {new_exp_distance[i] <- 1 - (P[i,i+1] - P[i+1,i+2]) - (P[i+1,i] - 0) }
      #       # else if(i == 1) {new_exp_distance[i] <- 1 - (P[i,i+1] - 0) - (P[i+1,i] - 0) }
      #       # else if(i == k-1) {new_exp_distance[i] <- 1 - (P[i,i+1] - 0) - (P[i+1,i] - P[i,i-1])}
      #       # else {new_exp_distance[i] <- 1 - (P[i,i+1] - P[i+1,i+2]) - (P[i+1,i] - P[i,i-1])}
      #
      #       # THIS IS ALWAYS VALID:
      #       if(i == 1 && k >= 3) { new_exp_distance[i] <- 2*(max(0,P[i+1,i+2] - P[i,i+1]) + max(0,0 - P[i+1,i])) +
      #         1*(1-(max(0,P[i+1,i+2] - P[i,i+1]) + max(0,0 - P[i+1,i]) +
      #                 max(0, P[i,i+1] - P[i+1,i+2]) + max(0,P[i+1,i] - 0))) }
      #       else if(i == 1) {new_exp_distance[i] <- 2*(max(0,0 - P[i,i+1]) + max(0,0 - P[i+1,i])) +
      #         1*(1-(max(0,0 - P[i,i+1]) + max(0,0 - P[i+1,i]) +
      #                 max(0, P[i,i+1] - 0) + max(0,P[i+1,i] - 0)))}
      #       else if(i == k-1) {new_exp_distance[i] <- 2*(max(0,0 - P[i,i+1]) + max(0,P[i,i-1] - P[i+1,i])) +
      #         1*(1-(max(0,0 - P[i,i+1]) + max(0,P[i,i-1] - P[i+1,i]) +
      #                 max(0, P[i,i+1] - 0) + max(0,P[i+1,i] - P[i,i-1])))}
      #       else {new_exp_distance[i] <- 2*(max(0,P[i+1,i+2] - P[i,i+1]) + max(0,P[i,i-1] - P[i+1,i])) +
      #         1*(1-(max(0,P[i+1,i+2] - P[i,i+1]) + max(0,P[i,i-1] - P[i+1,i]) +
      #                 max(0, P[i,i+1] - P[i+1,i+2]) + max(0,P[i+1,i] - P[i,i-1])))
      #       }
      #
      #     } else {
      #       #Expected distance after more than one time step
      #
      #       # THIS IS VALID ONLY IF THE EFFICIENCY CONDITION IS SATISFIED:
      #       # if(i == 1 && k>= 3) {new_exp_distance[i] <- P[i+1,i+2]*exp_distance[i+1] +
      #       #   (1-P[i,i+1]-P[i+1,i])*exp_distance[i] }
      #       # else if(i == 1) {new_exp_distance[i] <-(1-P[i,i+1]+P[i+1,i])*exp_distance[i]}
      #       # else if(i == k-1) {new_exp_distance[i] <- P[i,i-1]*exp_distance[i-1] +
      #       #   (1-P[i,i+1]-P[i+1,i])*exp_distance[i]}
      #       # else {new_exp_distance[i] <- P[i+1,i+2]*exp_distance[i+1] +
      #       #   P[i,i-1]*exp_distance[i-1] +
      #       #   (1-P[i,i+1]-P[i+1,i])*exp_distance[i]}
      #
      #       # THIS IS ALWAYS VALID:
      #       if(i == 1 && k >= 3) {new_exp_distance[i] <- max(P[i+1,i+2]-P[i,i+1],0)*(exp_distance[i]+exp_distance[i+1]) +
      #         max(0-P[i+1,i],0)*(0+exp_distance[i]) +
      #         min(P[i,i+1],P[i+1,i+2])*exp_distance[i+1] +
      #         min(0,P[i+1,i])*0 +
      #         (1-max(P[i+1,i+2]-P[i,i+1],0)-max(0-P[i+1,i],0)-min(P[i,i+1],P[i+1,i+2])-min(0,P[i+1,i])-max(P[i,i+1]-P[i+1,i+2],0)-max(P[i+1,i]-0,0))*exp_distance[i]}
      #       else if(i == 1) {new_exp_distance[i] <- max(0-P[i,i+1],0)*(exp_distance[i]+exp_distance[i+1]) +
      #         max(0-P[i+1,i],0)*(0+exp_distance[i]) +
      #         min(P[i,i+1],0)*exp_distance[i+1] +
      #         min(0,P[i+1,i])*0 +
      #         (1-max(0-P[i,i+1],0)-max(0-P[i+1,i],0)-min(P[i,i+1],0)-min(0,P[i+1,i])-max(P[i,i+1]-0,0)-max(P[i+1,i]-0,0))*exp_distance[i]}
      #       else if(i == k-1) {new_exp_distance[i] <- max(0-P[i,i+1],0)*(exp_distance[i]+0) +
      #         max(P[i,i-1]-P[i+1,i],0)*(exp_distance[i-1]+exp_distance[i]) +
      #         min(P[i,i+1],0)*0 +
      #         min(P[i,i-1],P[i+1,i])*exp_distance[i-1] +
      #         (1-max(0-P[i,i+1],0)-max(P[i,i-1]-P[i+1,i],0)-min(P[i,i+1],0)-min(P[i,i-1],P[i+1,i])-max(P[i,i+1]-0,0)-max(P[i+1,i]-P[i,i-1],0))*exp_distance[i]}
      #       else { new_exp_distance[i] <- max(P[i+1,i+2]-P[i,i+1],0)*(exp_distance[i]+exp_distance[i+1]) +
      #         max(P[i,i-1]-P[i+1,i],0)*(exp_distance[i-1]+exp_distance[i]) +
      #         min(P[i,i+1],P[i+1,i+2])*exp_distance[i+1] +
      #         min(P[i,i-1],P[i+1,i])*exp_distance[i-1] +
      #         (1-max(P[i+1,i+2]-P[i,i+1],0)-max(P[i,i-1]-P[i+1,i],0)-min(P[i,i+1],P[i+1,i+2])-min(P[i,i-1],P[i+1,i])-max(P[i,i+1]-P[i+1,i+2],0)-max(P[i+1,i]-P[i,i-1],0))*exp_distance[i]}
      #
      #     }
      #
      #   }
      #   #Update exp_distance
      #   exp_distance <- new_exp_distance
      #   #Compute new expected # tosses
      #   new_exp_tosses <- exp_tosses + sum(new_exp_distance)
      #   t <- t+1
      # }
      # if(new_exp_tosses == 0) {new_exp_tosses <- 1}
      # return(new_exp_tosses)
    },
    expected.tosses.bound.loose = function() {
      #This method computes a loose bound for the expected number
      #of tosses when m=2 and the efficiency condition is satisfied
      #independent on p

      if(private$m > 2) {
        stop("The expected number of tosses is computable only when the original die has 2 faces.")
      }
      if(!private$efficiency_condition) {
        stop("The bound can be obtained only when the efficiency condition is satisfied. Try to increase the degree of the ladder.")
      }

      R <- private$ladder_fine_connected$get.R()
      M <- private$ladder_fine_connected$get.M() #matrix of coefficients
      R <- R[order(private$ladder_fine_connected$get.M()[,1]+1)] #Reorder coefficients
      k <-  private$ladder_fine_connected$get.k()

      #Generate matrix of coefficients (with right order)
      P_coeff <- matrix(0, nrow = k, ncol = k) #just coefficients without p
      for(i in 1:k) {
        if(i > 1) {
          P_coeff[i,i-1] <- R[i-1]/max(R[i-1],R[i])
        }
        if(i < k) {
          P_coeff[i,i+1] <- R[i+1]/max(R[i],R[i+1])
        }
      }
      #Find the i that minimises min(P_{i,i+1}-P_{i+1,i+2}, P_{i+1,i} - P_{i,i-1})
      min_i <- 1
      rho <- 1
      for(i in 1:(k-1)) {
        if(i == 1 && k > 2) {
          val <- 1-min(P_coeff[i,i+1]-P_coeff[i+1,i+2], P_coeff[i+1,i])
        } else if((i == 1 && k==2) || (i == k-1 && k == 2)) {
          val <- 1-min(P_coeff[i,i+1], P_coeff[i+1,i])
        } else if(i == k-1 && k > 2) {
          val <- 1-min(P_coeff[i,i+1], P_coeff[i+1,i]-P_coeff[i,i-1])
        } else {
          val <- 1-min(P_coeff[i,i+1]-P_coeff[i+1,i+2], P_coeff[i+1,i]-P_coeff[i,i-1])
        }

        if(val < rho) {
          min_i <- i
          rho <- val
        }
      }

      #The bound is (k-1)/(1-rho) where rho = 1-min(P_{i,i+1}-P_{i+1,i+2} + P_{i+1,i} - P_{i,i-1})
      return((k-1)/(1-rho))

    },
    expected.tosses.extreme = function(p) {
      browser()
      if(private$m > 2) {
        stop("The expected number of tosses is computable only when the original die has 2 faces.")
      }
      if(p != 0 && p != 1) {
        stop("This method provides the exact expected value of tosses only for p=0 or p=1")
      }

      R <- private$ladder_fine_connected$get.R()
      M <- private$ladder_fine_connected$get.M() #matrix of coefficients
      R <- R[order(private$ladder_fine_connected$get.M()[,1]+1)] #Reorder coefficients
      k <-  private$ladder_fine_connected$get.k()

      #Generate matrix of coefficients (with right order)
      P <- matrix(0, nrow = k, ncol = k)
      P_coeff <- matrix(0, nrow = k, ncol = k) #just coefficients without p
      for(i in 1:k) {
        if(i > 1) {
          P[i,i-1] <- (1-p)*R[i-1]/max(R[i-1],R[i])
          P_coeff[i,i-1] <- R[i-1]/max(R[i-1],R[i])
        }
        if(i < k) {
          P[i,i+1] <- p*R[i+1]/max(R[i],R[i+1])
          P_coeff[i,i+1] <- R[i+1]/max(R[i],R[i+1])
        }
      }

      if(p == 1) {
        #Get the upper offdiagonal of P
        P_offdiagonal <- P[row(P) == (col(P) - 1)]
      } else {
        #Get the lower offdiagonal of P
        P_offdiagonal <- P[row(P) == (col(P) + 1)]
      }
      return(sum(1/P_offdiagonal))
    },
    increase.degree = function(d = 1,verbose=FALSE) {
      #We increase the degree of the dice enterprise by just applying the
      #multinomial theorem and create a new Dice Enterprise
      if(d == 0) {
        return(self$clone())
      }
      increase_degree <- private$ladder_fine_connected$get.degree() - private$ladder_initial$get.degree() + d
      incresed_de <- DiceEnterprise$new(private$G_orig,d=increase_degree,verbose = verbose)
      return(incresed_de)
    },
    impose.efficiency = function(method = c("efficiency_condition","bound"), threshold = 1e2, true_p = NULL) {
      #This function keeps on increasing the degree of the DE
      #until the efficiency condition is satisfied (if method is efficiency_condition)
      #or until the bound on the expected number of tosses is minimised.
      method <- match.arg(method)
      if(method == "bound") {original_bound <- self$expected.tosses.bound(true_p = true_p)}
      for(d in 1:threshold) {
        new_de <- self$increase.degree(d=d)
        if(method == "efficiency_condition" && isTRUE(new_de$get.efficiency.condition())) {
          return(new_de)
        } else if(method == "bound") {
          new_bound <- new_de$expected.tosses.bound(true_p = true_p)
          if(new_bound >= original_bound) {
            return(self$increase.degree(d=d-1))
          } else {
            original_bound <- new_bound
          }
        }
      }
      stop("Impossible to construct efficienct dice enterprise. Try to increase the threshold.")
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
    efficiency_condition = NA, #Efficiency of the dice enterprise. If TRUE the expected time for coalescence is bounded exponentially in time for every p. Available only for fine and connected.,
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
    generate.ladder_initial2.increase_degree = function(R_coeff, M_row, incr_degree, m) {
      #Auxiliary function of generate.ladder_initial2
      #Given R and M, increases the degree of one coefficient of the ladder
      #of the amount specified in incr_degree

      M_row <- as.numeric(M_row)
      R_coeff <- as.numeric(R_coeff)
      #Increase the degree
      if(incr_degree == 0) {
        new_R <- R_coeff
        new_M <- matrix(M_row, nrow = 1, byrow = TRUE)
      } else {
        discrete_simplex <- discrete.simplex(incr_degree,m) #discrete simplex
        new_R <- t(apply(discrete_simplex, 1, function(vec_simplex) {
          multinomial.coeff(incr_degree,vec_simplex)*R_coeff
        }))
        new_M <- t(apply(discrete_simplex, 1, function(vec_simplex) {
          vec_simplex + M_row
        }))
      }
      return(list(R=new_R,M=new_M))
    },
    generate.ladder_initial2.reduce = function(R,M) {
      #Auxiliary function of generate.ladder_initial2
      #This function sum all the coeffcieints in R with the same power and remove zero coeffiicients
      if(!is.numeric(R) || !is.matrix(M)) {
        stop("R must be a vector and M a matrix.")
      }
      #Sum all the coefficients in R that correspond to the same powers
      indices <- indices.unique.mat(M)
      M <- unique(M) #Remove repeated rows
      R_aux <- rep(NA, length = nrow(M)) #will remove duplicates
      for(l in 1:nrow(M)) {
        R_aux[l] <- round(sum(R[which(indices == l)]), digits = 8) #round to avoid machine error
      }
      #Remove all zeros coefficients
      if(any(R_aux == 0)) { #Need to check for the next two lines to work properly
        new_M <- matrix(M[-which(R_aux == 0),], ncol = ncol(M), byrow = FALSE) #NOT GREAT, but it should work cause we rounded up the values of R
        new_R <- R_aux[-which(R_aux == 0)]
      } else {
        new_M <- M
        new_R <- R_aux
      }
      return(list(R=new_R, M=new_M))
    },
    generate.ladder_initial = function(threshold) {
      #This method makes sure that:
      #1) The polynomials are homogeneous (same degree)
      #2) The coefficients of the polynomials are all positive
      #3) All the polynomials have the same degree

      #G_poly is a list of v elements. Each element is itself a list of 3 elements
      #where the first elements are the coefficients R, the second element
      #is a matrix of powers of the p_i, the 3rd element is the degree associated
      #to each coefficients (it's just the sum of the rows of the second element)
      #Prepare output
      M_list <- vector("list",private$v)
      R_list <- vector("list",private$v)
      private$A_f_initial <- vector("list", private$v) #To store disaggregation

      #Store degree of the generate members
      new_degree <- numeric(private$v)
      threshold_counter <- 0 #To control how much we keep on increasing the degree

      for(i in 1:private$v) {
        R_list[[i]] <- private$G_poly[[i]][[1]] #Original Coefficients
        M_list[[i]] <- private$G_poly[[i]][[2]] #Original Powers
        degree_R <- private$G_poly[[i]][[3]] #Degree of each coefficient
        #STEP 1-2: Convert to homogeneous polynomial
        #by increasing the degree and until all the coefficients are positive
        while(new_degree[i] < private$d || any(R_list[[i]] < 0)) {
          temp_R_list <- vector("list", length(R_list[[i]]))
          temp_M_list <- vector("list", length(R_list[[i]]))
          for(j in 1:length(R_list[[i]])) {
            if(new_degree[i] == 0) {
              #First step -> create homogeneous polynomial
              aux <- private$generate.ladder_initial2.increase_degree(R_coeff = R_list[[i]][j],
                    M_row = M_list[[i]][j,], incr_degree = private$d-degree_R[j], m = private$m)
            } else {
              #Polynomal is already homogeneous, but some coefficients are negative
              #we increase the degree by 1
              aux <- private$generate.ladder_initial2.increase_degree(R_coeff = R_list[[i]][j],
                     M_row = M_list[[i]][j,], incr_degree = 1, m = private$m)
            }
            temp_R_list[[j]] <- aux$R
            temp_M_list[[j]] <- as.data.frame(aux$M)
          }
          #Update with new degree
          if(new_degree[i] == 0) {
            new_degree[i] <- private$d
          } else {
            new_degree[i] <- new_degree[i] + 1
          }
          #Put them all together
          R_list[[i]] <- do.call(c, temp_R_list) #rbindlist is extremely efficient
          M_list[[i]] <- as.matrix(rbindlist(temp_M_list))
          #Sum all the coefficients in R that corresponds to same power and remove zero coefficients
          aux2 <- private$generate.ladder_initial2.reduce(R=R_list[[i]], M=M_list[[i]])
          R_list[[i]] <- aux2$R
          M_list[[i]] <- as.data.frame(aux2$M)
          threshold_counter <- threshold_counter + 1
          if(threshold_counter > threshold) {
            stop("I am having a hard time constructing a valid ladder. Are you sure
                 it is a valid function? Check that numerator and
                 denominator do not share any common roots. Also, check that they are both positive polynomials
                 (in case, switch signs of the coefficients accordingly). In case, increase the threshold.")
          }
        }

      }

      #STEP 3: make sure that all the members of the ladder have the same degree
      if(length(unique(new_degree)) > 1) {
        #Polynomials have different degrees
        private$d <- max(new_degree) #update new degree of the ladder
        for(i in 1:private$v) {
          if(new_degree[i] != private$d) {
            #Need to increase the degree
            temp_R_list <- vector("list", length(R_list[[i]]))
            temp_M_list <- vector("list", length(R_list[[i]]))
            for(j in 1:length(R_list[[i]])) {
              aux <- private$generate.ladder_initial2.increase_degree(R_coeff = R_list[[i]][j],
                          M_row = M_list[[i]][j,], incr_degree = private$d-new_degree[i], m = private$m)
              temp_R_list[[j]] <- aux$R
              temp_M_list[[j]] <- as.data.frame(aux$M)
            }
            #Put them all together
            R_list[[i]] <- do.call(c, temp_R_list) #rbindlist is extremely efficient
            M_list[[i]] <- as.matrix(rbindlist(temp_M_list))
            #Sum all the coefficients in R that corresponds to same power and remove zero coefficients
            aux2 <- private$generate.ladder_initial2.reduce(R=R_list[[i]], M=M_list[[i]])
            R_list[[i]] <- aux2$R
            M_list[[i]] <- as.data.frame(aux2$M)
          }
        }
      }

      #Put them all together to construct the new ladder
      R_final <- do.call(c, R_list)
      M_final <- as.matrix(rbindlist(M_list))
      private$ladder_initial <- Ladder$new(M = M_final, R=R_final)
      #Prepare for disaggregation
      counter <- 0
      for(i in 1:private$v) {
        private$A_f_initial[[i]] <- (counter+1):(counter+nrow(M_list[[i]]))
        counter <- counter+nrow(M_list[[i]])
      }
    },
    generate.ladder_initialOLD = function() {
      #Using the Theorem, convert the polynomial to the requested form
      #so that a ladder is generated which is a disaggregation of f.

      #CONVERT TO BERNSTEIN REPRESENTATION (NOT EFFICIENT)
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
      browser()
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
                                 initialize = function(G,toss.coins,num_coins,die_type = c("uniform","toss_all","first_heads"), d = NULL, verbose = FALSE, threshold = 1e2) {
                                   private$num_coins <- num_coins
                                   private$die_type <- match.arg(die_type)
                                   private$toss.coins.fun <- toss.coins
                                   super$initialize(G=G,d=d,verbose=verbose,threshold = threshold)
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
