Ladder$debug("is.ladder")
m <- 3
k <- 2
M <- matrix(c(0,0,2,1,1,0), ncol = 3, byrow = TRUE)
R <- c(0.5, sqrt(2))
ciao <- Ladder$new(m = m, k = k, M = M, R = R)

p <- rsimplex(5,m)

for(i in 1:nrow(p)) {
  R[i]*prod(p[i,]^M[i,])
}

#Valid fine and connected ladder
m <- 3
k <- 6
M <- matrix(c(3,0,0,
              2,0,1,
              1,2,0,
              1,1,1,
              1,0,2,
              0,2,1), byrow = TRUE, ncol = 3)
R <- c(sqrt(2),1,1/4,2,1/2,3/4)
l1 <- Ladder$new(m = m, k = k, M = M, R = R)

#Valid connected ladder (not fine)
m <- 3
k <- 7
M <- matrix(c(3,0,0,
              2,0,1,
              1,2,0,
              1,1,1,
              1,0,2,
              0,2,1,
              1,1,1), byrow = TRUE, ncol = 3)
R <- c(sqrt(2),1,1/4,2,1/2,3/4,0.7)
l2 <- Ladder$new(m = m, k = k, M = M, R = R)

#Valid fine ladder (not connected)
m <- 3
k <- 5
M <- matrix(c(3,0,0,
              1,2,0,
              1,1,1,
              1,0,2,
              0,2,1), byrow = TRUE, ncol = 3)
R <- c(sqrt(2),1,1/4,2,1/2)
l3 <- Ladder$new(m = m, k = k, M = M, R = R)

#Valid ladder (not fine not connected)
m <- 3
k <- 6
M <- matrix(c(3,0,0,
              1,1,1,
              1,2,0,
              1,1,1,
              1,0,2,
              0,2,1), byrow = TRUE, ncol = 3)
R <- c(sqrt(2),1,1/4,2,1/2,0.44)
l4 <- Ladder$new(m = m, k = k, M = M, R = R)
