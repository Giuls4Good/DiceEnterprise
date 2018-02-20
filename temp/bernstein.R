require(matrixcalc)

A <- matrix(c(0,0,0,0,
              0,0,0,0,
              0,0,0,0,
              sqrt(2),0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
P <- pascal.matrix(nrow(A))
res <- t(P%*%t(P%*%A))
coeff <- rev(diag(res[nrow(res):1,])*P[nrow(P),])

B <- matrix(c(3,0,0,0,
              -9,0,0,0,
              11,0,0,0,
              -5,0,0,0), nrow = 4, ncol = 4, byrow = TRUE)
C <- matrix(NA, nrow= nrow(B), ncol = ncol(B))
k <- nrow(B)-1 #degree
for(i in 1:nrow(B)) {
  for(j in 1:ncol(B)) {
    if(i+j <= k+2) {
      C[i,j] <- B[i,j]*factorial(k-i-j+2)*factorial(i-1)*factorial(j-1)/factorial(k)
    } else {
      C[i,j] <- B[i,j]
    }
  }
}
P <- pascal.matrix(nrow(C))
res <- t(P%*%t(P%*%C))
coeff <- rev(diag(res[nrow(res):1,])*P[nrow(P),])

###
# MULTIVARIATE - WRONG, THIS IS A DIFFERENT REPRESENTATION
# f(p) = sqrt(2)p_1^3+p_1^2p_3+1/4p_1p_2^2+2p_1p_2p_3+1/2p_1p_3^2+3/4p_^2p_3
rm(list=ls())
deg <- 3 #degree
m <- 3 #number of variables
coeff_power <- array(data = 0, dim = rep(deg+1,deg)) #We store the coefficients in a bigger matrix
# coeff_power[3+1,0+1,0+1] <- sqrt(2) #sqrt(2)p_1^3
# coeff_power[2+1,0+1,1+1] <- 1 #p_1^2p_3
# coeff_power[1+1,2+1,0+1] <- 1/4 #1/4p_1p_2^2
# coeff_power[1+1,1+1,1+1] <- 2 #2p_1p_2p_3
# coeff_power[1+1,0+1,2+1] <- 1/2 #1/2p_1p_3^2
# coeff_power[0+1,2+1,1+1] <- 3/4 #3/4p_2^2p_3
coeff_power[3+1,0+1,0+1] <- 6
coeff_power[1+1,1+1,1+1] <- 36

A <- matrix(0, nrow = deg+1, ncol = (deg+1)^(m-1))
C <- matrix(0, nrow = deg+1, ncol = (deg+1)^(m-1))
indices <- vector("list", deg+1) #Store how to convert from array to matrix
for(aux in 1:(deg+1)) {
  indices[[aux]] <- vector("list", (deg+1)^(m-1))
}
for(i_1 in 0:(dim(coeff_power)[1]-1)) { #To generalize for generic array
  for(i_2 in 0:(dim(coeff_power)[2]-1)) {
    for(i_3 in 0:(dim(coeff_power)[3]-1)) {
      if(i_1+i_2+i_3 <= deg) {
        l <- i_1 + 1
        kappa <- i_2 + 1 + i_3*(deg+1)^(3-2) #To generalize
        A[l,kappa] <- coeff_power[i_1+1,i_2+1,i_3+1]
        C[l,kappa] <- A[l,kappa]*factorial(deg-i_1-i_2-i_3)*factorial(i_1)*factorial(i_2)*factorial(i_3)/factorial(deg) #To generalize
        indices[[l]][[kappa]] <- c(i_1+1,i_2+1,i_3+1)
      }
    }
  }
}

mat2arr <- function(C,indices,deg) { #useless
  coeff <- array(data = 0, dim = rep(deg+1,deg))
  for(i in 1:(deg+1)) {
    for(j in 1:((deg+1)^(m-1))) {
      aux <- indices[[i]][[j]]
      coeff[aux[1],aux[2],aux[3]] <- C[i,j] #to generalize
    }
  }
  return(coeff)
}

cyclic_ordering <- function(M,deg) {
  res <- matrix(NA, nrow=nrow(M), ncol = ncol(M))
  for(i_1 in 0:deg) { #To generalize for generic matrix
    for(i_2 in 0:deg) {
      for(i_3 in 0:deg) {
        l <- i_1 + 1
        kappa <- i_2 + 1 + i_3*(deg+1)^(3-2) #To generalize
        #Cyclic ordering - to generalize
        new_i_1 <- i_3
        new_i_2 <- i_1
        new_i_3 <- i_2
        new_l <- new_i_1 + 1
        new_kappa <- new_i_2 + 1 + new_i_3*(deg+1)^(3-2) #To generalize
        res[new_l,new_kappa] <- M[l,kappa]
      }
    }
  }
  return(res)
}

all.equal(cyclic_ordering(cyclic_ordering(cyclic_ordering(A, deg), deg), deg), A)
C_t <- vector("list",m+1)
C_t[[1]] <- C
P <- pascal.matrix(deg+1)
for(t in 2:(m+1)) {
  C_t[[t]] <- cyclic_ordering(P%*%C_t[[t-1]], deg)
}
bernstein_patch <- C_t[[m+1]]

prova <- vector("list")
counter <- 1
for(i_1 in 0:deg) { #To generalize for generic matrix
  for(i_2 in 0:deg) {
    for(i_3 in 0:deg) {
      if(i_1+i_2+i_3 <= deg) {
        l <- i_1 + 1
        kappa <- i_2 + 1 + i_3*(deg+1)^(3-2) #To generalize
        if(round(bernstein_patch[l,kappa],4) > 1e-6) {
          bern_coeff <- bernstein_patch[l,kappa]*factorial(deg)/(factorial(i_1)*factorial(i_2)*factorial(i_3))
          prova[[counter]] <- paste0(bern_coeff,"*x^",i_1,"*(1-x)^",deg-i_1,"*y^",i_2,"*(1-y)^",deg-i_2,"*z^",i_3,"*(1-z)^",deg-i_3)
          print(prova[[counter]])
          if(counter == 1) {
            res_print <- prova[[counter]]
          }
          else {
            res_print <- paste0(res_print,"+",prova[[counter]])
          }
          counter <- counter + 1
        }
      }
    }
  }
}

print(res_print)

res_print <- str_replace_all(res_print,"x","0.2")
res_print <- str_replace_all(res_print,"y","0.5")
res_print <- str_replace_all(res_print,"z","0.3")
eval(parse(text=res_print))

# MULTIVARIATE
# This is correct and it is the classical method described in
# "A matrix method for efficient computation of Bern. coefficients"
# Its complexity is O(n^(2m)) where n is the maximum degree of a single variable (not of the polynomial)
# and m is the number of variables.
# More efficient methods are available, such as Garloff's method (complexity O(n^(m+1))) or their
# proposed method (complexity O(n^(m+1))). However I didn't quite get them and this package
# doesn't "care" about efficieny for now
rm(list=ls())
N <- c(2,3,2) #Maximum degree of each variables, i.e. p_1 has max. degree equal to 2, p_2 equal to 3...
coeff_power <- array(data = 0, dim = N+1) #+1 cause the degree may be equal to 0. coeff_power stores the usual coefficients
coeff_power[3,4,3] <- -366 #-366x^2y^3z^2. Notice that the degree and the indices have a difference of 1.
coeff_power[3,4,2] <- 396
coeff_power[3,4,1] <- -30
coeff_power[3,3,3] <- 672
coeff_power[3,3,2] <- -672
coeff_power[3,2,3] <- -336
coeff_power[3,2,2] <- 336
coeff_power[2,4,3] <- 30
coeff_power[2,4,2] <- -60
coeff_power[2,4,1] <- 30
#This should correspond to bernstein coefficients: b_130 = 15 and b_211 = 56. All the others are 0.

compute_coeff_bern <- function(N,I,a) {
  res <- 0
  for(j_1 in 0:I[1]) { #To generalize
    for(j_2 in 0:I[2]) {
      for(j_3 in 0:I[3]) {
        if(!is.na(a[j_1+1,j_2+1,j_3+1]) && a[j_1+1,j_2+1,j_3+1] != 0) {
          res <- res + (choose(I[1],j_1)*choose(I[2],j_2)*choose(I[3],j_3))/
            (choose(N[1],j_1)*choose(N[2],j_2)*choose(N[3],j_3))*a[j_1+1,j_2+1,j_3+1]
        }
      }
    }
  }
  return(res)
}

coeff_bern <- array(data = 0, dim = N+1)
my_coeff_bern <- array(data = 0, dim = N+1) #The coefficients I need are not exactly the Bernstein ones
#as we want to incorporate the binomial coefficients
for(i_1 in 0:N[1]) { #To generalize
  for(i_2 in 0:N[2]) {
    for(i_3 in 0:N[3]) {
      coeff_bern[i_1+1,i_2+1,i_3+1] <- compute_coeff_bern(N,c(i_1,i_2,i_3),coeff_power)
      my_coeff_bern[i_1+1,i_2+1,i_3+1] <- coeff_bern[i_1+1,i_2+1,i_3+1]*(choose(N[1],i_1)*choose(N[2],i_2)*choose(N[3],i_3))
    }
  }
}


