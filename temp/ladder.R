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
Ladder$debug("initialize")
de <- DiceEnterprise$new(G=list(
  list(c(1,sqrt(2)),c("120","022")),
  list(c(4,1/2,3),c("470","000","123")),
  list(c(7,2),matrix(c(1,3,4,0,0,2),byrow=TRUE,ncol=3))
), verbose = TRUE)
true_p <- c(0.6,0.2,0.2)
sample_f <- de$sample(n = 1000, true_p = true_p, num_cores = 4, verbose = TRUE, global = FALSE)
sample_f_inefficient <- de$sample(n = 10, true_p = true_p, num_cores = 1, verbose = TRUE, global = TRUE) #Uses global consatnt -> inefficient
print(de$evaluate(true_p))
print(table(sample_f)/length(sample_f))
plotConfidenceInterval(sample_f,de$evaluate(true_p))
require(profvis)
profvis({de$sample(n = 1000, true_p = true_p, num_cores = 1, verbose = TRUE)})
require(rbenchmark)
benchmark(de$sample(n = 1000, true_p = true_p, num_cores = 1, verbose = TRUE, global = FALSE),
          de$sample(n = 1000, true_p = true_p, num_cores = 2, verbose = TRUE, global = FALSE),
          de$sample(n = 1000, true_p = true_p, num_cores = 4, verbose = TRUE, global = FALSE),
          de$sample(n = 1000, true_p = true_p, num_cores = 8, verbose = TRUE, global = FALSE),
          order = "relative", replications = 1)
##
# Dice enterprise example
de_ladder <- DiceEnterprise$new(G=list(
  list(sqrt(2),"300"), #G_1
  list(1,"201"), #G_1
  list(1/4,"120"), #G_1
  list(2,"111"), #G_1
  list(1/2,"102"), #G_1
  list(3/4,"021") #G_6
), verbose = TRUE)
true_p <- c(0.6,0.2,0.2)
sample_f <- de_ladder$sample(n = 10000, true_p = true_p, num_cores = 1, verbose = TRUE, global = FALSE)
sample_f_global <- de_ladder$sample(n = 10000, true_p = true_p, num_cores = 1, verbose = TRUE, global = TRUE) #Less efficient (slightly in this case, but can be a lot less efficient!)
print(de_ladder$evaluate(true_p))
print(table(sample_f)/length(sample_f))
plotConfidenceInterval(sample_f,de_ladder$evaluate(true_p))

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
sample_f <- bf$sample(n = 1000, true_p = true_p, num_cores = 1, verbose = TRUE)
print(bf$evaluate(true_p))
print(table(sample_f)/length(sample_f))
plotConfidenceInterval(sample_f,bf$evaluate(true_p))

##
# Dice enterprise example 2
rm(list=ls())
de2 <- DiceEnterprise$new(G=list(
  list(c(7,sqrt(3)),c("0700","2040")),
  list(c(4,5),matrix(c(1,1,1,1,10,0,0,0),byrow=TRUE,ncol=4)),
  list(c(pi/8),"0000")
), verbose = TRUE)
print(de2)
true_p <- c(0.2,0.7,0.06,0.04)
sample_f <- de2$sample(n = 1000, true_p = true_p, num_cores = 4, verbose = TRUE)
print(de2$evaluate(true_p))
print(table(sample_f)/length(sample_f))
plotConfidenceInterval(sample_f,de2$evaluate(true_p))

##
# 2 coin algorithm
c1 <- 1
c2 <- 1
p1 <- 0.00001
p2 <- 0.00005
p_dice <- c(p1*p2,(1-p1)*(1-p2),p1*(1-p2),p2*(1-p1)) #convert to Bernstein
print((c1+c2)/(c1*p1+c2*p2)) #Theoretical of 2 coin
f_2coin <- list(
  list(c(c1,c1),c("1000","0010")),
  list(c(c2,c2),c("1000","0001"))
)
de_2coin <- DiceEnterprise$new(f_2coin, verbose = TRUE)
sample_2coin <- de_2coin$sample(n = 100, true_p = p_dice, num_cores = 2, verbose = TRUE) #Slower :(
print(de_2coin$evaluate(p_dice))
print(table(sample_2coin)/length(sample_2coin))
plotConfidenceInterval(sample_2coin,de_2coin$evaluate(p_dice))

##
# INDEPENDENT COINS
#Given 3 independent coins, this simulates with probs
#f(p) \propto (p_1,p_2,p_3)

rm(list=ls())
p_coins <- c(0.16,0.18,0.06)
de_indep <- DiceEnterprise$new(G=list(
  list(rep(1,4),c("20000","10100","10010","00110")),
  list(rep(1,4),c("20000","11000","10010","01010")),
  list(rep(1,4),c("20000","11000","10100","01100"))
), verbose = TRUE)

toss.coins <- function(true_p) { #tosses the three coins
  return(sapply(true_p, function(p) {sample(1:2, size = 1, prob = c(p,1-p))})) #1 or 2 (not 0))
}
roll.die <- function(n,toss.fun) { #roll the die
  res <- numeric(n)
  for(i in 1:n) {
    while(TRUE) {
     toss_res <- toss.fun(p_coins)
     if(isTRUE(all.equal(toss_res,c(1,1,1)))) {
       res[i] <- 1 #q0
       break
     } else if(isTRUE(all.equal(toss_res,c(2,1,1)))) {
       res[i] <- 2 #q1
       break
     } else if(isTRUE(all.equal(toss_res,c(1,2,1)))) {
       res[i] <- 3 #q2
       break
     } else if(isTRUE(all.equal(toss_res,c(1,1,2)))) {
       res[i] <- 4 #q3
       break
     } else {
       res[i] <- 5 #none of the above
       break
     }
    }
  }
  return(res)
}

sample_size <- 1000
print(paste0("True prob: ",prod(p_coins)," - ",
             (1-p_coins[1])*p_coins[2]*p_coins[3]," - ",
             p_coins[1]*(1-p_coins[2])*p_coins[3]," - ",
             p_coins[1]*p_coins[2]*(1-p_coins[3])," - "))
print(table(roll.die(sample_size, toss.fun = toss.coins))/sample_size)

res <- de_indep$sample(n = sample_size, roll.fun = roll.die, verbose = TRUE, toss.fun = toss.coins, num_cores = 4)
print(paste0("True prob:",p_coins[1]/sum(p_coins)," - ",p_coins[2]/sum(p_coins)," - ",p_coins[3]/sum(p_coins)))
print(table(res)/sample_size)
plotConfidenceInterval(res,p_coins/sum(p_coins))
