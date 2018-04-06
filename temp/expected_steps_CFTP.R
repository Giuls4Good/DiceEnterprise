# This is to test to find upper bounds on the number of
# steps required by CFTP.
# Lemma 3.4 in Huber book states that in Monotonic CFTP if
# P(X_t \neq Y_t) <= e^-1 where X_t is the chain starting at the
# minimum  state and Y_t is the state starting at the maximum state,
# then the expected number of steps required by Monotonic CFTP
# is less than 19.2t.

#Consider the example of the Bernoulli Factory given by
# f(p) = \frac{\sqrt{2}p^3}{(\sqrt{2}-5)p^3+11p^2-9p+3}
bf <- BernoulliFactory$new(f_1 = list(coeff = c(sqrt(2)), power = c(3)),
                           f_2 = list(coeff = c(-5,11,-9,3), power = c(3,2,1,0)))
coeff_P <- bf$get.ladder.fine.connected()$get.P() #Matrix P without the powers of p

#Let's consider the (extreme) case where p=1
p <- 1
#Compute the transition matrix
P <- coeff_P
P[1,] <- P[1,]*c(0,p,0,0,0)
P[2,] <- P[2,]*c(1-p,0,p,0,0)
P[3,] <- P[3,]*c(0,1-p,0,p,0)
P[4,] <- P[4,]*c(0,0,1-p,0,p)
P[5,] <- P[5,]*c(0,0,0,1-p,0)
for(i in 1:nrow(P)) {
  P[i,i] <- 1-sum(P[i,],na.rm = TRUE)
}

#Since p=1, the random walk can only go right and Y_t = 5 for all t
#Then P(X_t \neq Y_t) = 1-P(X_t = Y_t) = 1-P(X_t = 5)
#Let's find the smallest t such that this probability is smaller than e^-1
t <- 1
while(TRUE) {
  prob <- matrix(c(1,rep(0,nrow(P)-1)),nrow = 1)%*%(P%^%t)
  if(1-prob[1,5] <= exp(-1)) {
    break
  }
  t <- t+1
}
