Ladder$debug("is.ladder")
Ladder$debug("define.neighbourhoods")
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
M_ciao <- matrix(c(3,0,0,
              2,0,1,
              1,2,0,
              1,1,1,
              1,0,2,
              0,2,1), byrow = TRUE, ncol = 3)
R_ciao <- c(sqrt(2),1,1/4,2,1/2,3/4)
l1 <- Ladder$new(M = M_ciao, R = R_ciao)
print(l1)

#Valid connected ladder (not fine)
M_ciao <- matrix(c(3,0,0,
              2,0,1,
              1,2,0,
              1,1,1,
              1,0,2,
              0,2,1,
              1,1,1), byrow = TRUE, ncol = 3)
R_ciao <- c(sqrt(2),1,1/4,2,1/2,3/4,0.7)
l2 <- Ladder$new(M = M_ciao, R = R_ciao)
print(l2)

#Valid fine ladder (not connected)
M_ciao <- matrix(c(3,0,0,
              1,2,0,
              1,1,1,
              1,0,2,
              0,2,1), byrow = TRUE, ncol = 3)
R_ciao <- c(sqrt(2),1,1/4,2,1/2)
Ladder$debug("impose.connected")
l3 <- Ladder$new(M = M_ciao, R = R_ciao)
l3$impose.connected()
print(l3)

#Valid ladder (not fine not connected)
M_ciao <- matrix(c(3,0,0,
              1,1,1,
              1,2,0,
              1,1,1,
              1,0,2,
              0,2,1), byrow = TRUE, ncol = 3)
R_ciao <- c(sqrt(2),1,1/4,2,1/2,0.44)
l4 <- Ladder$new(M = M_ciao, R = R_ciao)
print(l4)

###
# BF
###

rm(list=ls())
pi_prime <- Ladder$new(M=matrix(c(0,3,2,1,3,0), ncol = 2, byrow=TRUE),
                       R=c(3,2,sqrt(2)))
print(pi_prime) #not connected
#1) Impose connected condition
aux <- pi_prime$impose.connected()
pi_tilde <- aux[[1]]
A_prime_tilde <- aux[[2]]
print(pi_tilde)
#2) Impose fineness condition
aux <- pi_tilde$impose.fineness()
pi <- aux[[1]]
A_pi_tilde <- aux[[2]]
print(pi)
#3) Sample from pi
true_p <- c(0.6,0.4)
sample_size <- 10000
sample_pi <- pi$sample(n = sample_size, roll.fun = function(n) {sample(1:length(true_p), size = n, replace = TRUE, prob = c(0.6,0.4))})
# check correctness
pi$evalute(p = true_p)
table(sample_pi)/sample_size
#4) Transform sample to pi_tilde
sample_pi_tilde <- disaggregationSample(sample_pi, origin = "original",
                                        disaggDist = pi_tilde, A = A_pi_tilde)
# check correctness
pi_tilde$evalute(p = true_p)
table(sample_pi_tilde)/sample_size
#5) Transform sample to pi_prime
sample_pi_prime <- disaggregationSample(sample_pi_tilde, origin = "disaggregation",
                                        disaggDist = pi_tilde, A = A_prime_tilde)
pi_prime$evalute(p = true_p)
table(sample_pi_prime)/sample_size

##
# Dice enterprise

rm(list=ls())
DiceEnterprise$debug("generate.ladder_initial")
de <- DiceEnterprise$new(G=list(
  list(c(1,sqrt(2)),c("120","022")),
  list(c(4,1/2,3),c("470","000","123")),
  list(c(7,2),matrix(c(1,3,4,0,0,2),byrow=TRUE,ncol=3))
), verbose = TRUE)

##
# Bernoulli Factory

rm(list=ls())
BernoulliFactory$debug("initialize")
BernoulliFactory$debug("sample")
#f(p) = sqrt(2)p^3/((sqrt(2)-5)p^3+11p^2-9p+3)
#1-f(p) = (-5p^3+11p^2-9p+3)/((sqrt(2)-5)p^3+11p^2-9p+3)
bf <- BernoulliFactory$new(f_1 = list(coeff = c(sqrt(2)), power = c(3)),
                       f_2 = list(coeff = c(-5,11,-9,3), power = c(3,2,1,0)),
                       verbose = TRUE)
true_p <- c(0.7,0.3)
sample_f <- bf$sample(n = 100, true_p = true_p, num_cores = 4)
print(bf$evaluate(true_p))
print(table(sample_f)/length(sample_f))
plotConfidenceInterval(sample_f,bf$evaluate(true_p))
