BernoulliFactory <- R6::R6Class("BernoulliFactory",
  inherit = DiceEnterprise,
  public = list(
    # initialize = function(num_string, den_string, variable = "p") {
    #   #f(p) = num_string/den_string. The variable name can be changed.
    #   #num_string and den_string are strings defining the function, in LaTeX style.
    #   #for instance: p^5+5p+4p^1-3p^0+14p^{15}+7
    #   #NOTICE: the only allowed symbols are: numbers, ^, {, }, +, -
    #   #The coefficients must be in front of the variable. I.e. 5p, not p5.
    #   #This means that after a p there can only be a ^, +, -, end of string.
    #   if(!(variable %in% letters) && !(variable %in% LETTERS)) {
    #     stop("The variable must be a single lower case letter.")
    #   }
    #   if(any(grepl("[^0-9p+\\-\\^{}]",c(num_string,den_string)))) {
    #     stop("The only allowed symbols are: 0-9, ^, {, }, +, -")
    #   }
    #   if(any(grepl(paste0(variable,"[^\\^+\\-]"), c(num_string,den_string)))) {
    #     stop(paste0(variable," can only be followed by ^,+ or -"))
    #   }
    #
    #   #Check if all the R_i are positive. Discard null ones. Error if there are negatives (function not in a simplex).
    #}
    sample = function(n,roll.fun,...) {
      sample_f <- super$sample(n,roll.fun,...) #Sample from the dice enterprise
      stop("Do some other stuff")
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
      if(verbose) {cat("DONE! \n")}
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
    sample = function(n,roll.fun = NULL, true_p = NULL,...) {
      stopifnot(!is.null(private$ladder_fine_connected),
                !is.null(private$ladder_initial),
                !is.null(private$ladder_connected),
                is.list(private$A_f_initial),
                is.list(private$A_initial_connected),
                is.list(private$A_fine_connected))
      if(is.null(roll.fun) && is.null(true_p)) {stop("Either declare roll.fun or the trye probabilities.")}
      if(is.null(roll.fun)) {roll.fun <- function(n) {sample(1:private$m, size = n, replace = TRUE, prob = true_p)}}
      #Get a sample from the fine and connected ladder
      sample_fine_connected <- private$ladder_fine_connected$sample(n = n, roll.fun = roll.fun,...)
      #Transform sample to the connected ladder
      sample_connected <- disaggregationSample(sample_fine_connected, origin = "original",
                                               disaggDist = private$ladder_connected,
                                               A = private$A_fine_connected)
      #Transform sample to the initial ladder
      sample_initial <- disaggregationSample(sample_connected, origin = "disaggregation",
                                             disaggDist = private$ladder_connected,
                                             A = private$A_initial_connected)
      #Transform sample to initial function
      sample_f <- disaggregationSample(sample_initial, origin = "disaggregation",
                                       disaggDist = private$ladder_initial,
                                       A = private$A_f_initial)
      return(sample_f)
    },
    get_G_poly = function() {private$G_poly},
    get_ladder_initial = function() {private$ladder_initial$clone()}
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
      private$A_initial_connected <- vector("list", private$v)
      M_rows <- 0 #counting the number of rows of the final ladder
      counter <- 0
      for(i in 1:private$v) {
        for(j in 1:length(private$G_poly[[i]][[3]])) { #For each coefficient
          deg <- private$G_poly[[i]][[3]][j] #degree of the a_n
          if(private$d > deg) {
            M_list[[i]] <- rbind(M_list[[i]], t(apply(discreteSimplex(private$d-deg,private$m), 1, function(x) {
              x + private$G_poly[[i]][[2]][j,]
            }))) #NOT EFFICIENT
            R_list[[i]] <- c(R_list[[i]], apply(discreteSimplex(private$d-deg,private$m), 1, function(x) {
              private$G_poly[[i]][[1]][j]*multinomialCoeff(private$d-deg,x)
            })) #NOT EFFICIENT
          } else {
            M_list[[i]] <- rbind(M_list[[i]], private$G_poly[[i]][[2]][j,]) #NOT EFFICIENT
            R_list[[i]] <- c(R_list[[i]], private$G_poly[[i]][[1]][j]) #NOT EFFICIENT
          }
        }
        M_rows <- M_rows + nrow(M_list[[i]])
        private$A_initial_connected[[i]] <- (counter+1):(counter+nrow(M_list[[i]]))
        counter <- counter+nrow(M_list[[i]])
      }
      #Construct ladder
      M_new <- do.call(rbind, sapply(M_list, unlist))
      R_new <- do.call(c, sapply(R_list, unlist))
      private$ladder_initial <- Ladder$new(M = M_new, R=R_new)
    }
  )
)
