## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7, fig.height=4,
  cache=TRUE
)
require("DiceEnterprise")
set.seed(10, "L'Ecuyer-CMRG") 

## ----toss_coin-----------------------------------------------------------
toss_coin <- function(n) {
  sample(1:2, size = n, replace = TRUE, prob = c(3/4,1/4))
}

## ----bf_poly-------------------------------------------------------------
f_1 <- list(coeff = c(sqrt(2)), power = c(3)) #f(p)
f_2 <- list(coeff = c(-5,11,-9,3), power = c(3,2,1,0)) #1-f(p)

## ----bf_construction-----------------------------------------------------
bf <- BernoulliFactory$new(f_1 = f_1, f_2 = f_2) #f_1 = f(p), f_2 = 1-f(p)

## ----bf_toss-------------------------------------------------------------
fp_tosses <- bf$sample(n = 10, roll.fun = toss_coin) #Produces 10 tosses of the f(p)-coin
print(fp_tosses)

## ----bf_diagnosis--------------------------------------------------------
bf <- BernoulliFactory$new(f_1 = f_1, f_2 = f_2, verbose = TRUE)
print(bf)

## ----bf_evaluate---------------------------------------------------------
print(bf$evaluate(3/4))

## ----bf_toss_2-----------------------------------------------------------
fp_tosses <- bf$sample(n = 1000, roll.fun = toss_coin, num_cores = 2, verbose = TRUE, double_time= TRUE) #Produces 1000 tosses of the f(p)-coin, using 2 cores and doubling the time step at each iteration of CFTP.
print(table(fp_tosses[[1]])/1000) #Empirical probabilities. Notice that the theoretical ones are given by print(bf$evaluate(3/4))
print(paste0("Average number of tosses required: ",mean(fp_tosses[[2]])))

## ----bf_plot-------------------------------------------------------------
plot.confidence.interval(fp_tosses[[1]],print(bf$evaluate(3/4)))

## ----bf_error, error=TRUE------------------------------------------------
bf_amp <- BernoulliFactory$new(f_1 = list(2,1), f_2 = list(c(1,-2),c(0,1)))

## ----bf_threshold, error=TRUE--------------------------------------------
bf_error <- BernoulliFactory$new(f_1 = list(c(-1,1,0.749),c(2,1,0)), 
                                 f_2 = list(c(1,-1,0.251),c(2,1,0)))
bf_no_error <- BernoulliFactory$new(f_1 = list(c(-1,1,0.749),c(2,1,0)), 
                                 f_2 = list(c(1,-1,0.251),c(2,1,0)),
                                 threshold = 500)

## ----bf_AR---------------------------------------------------------------
fp_tosses_AR <- bf$sample.AR(n = 1000, roll.fun = toss_coin)
plot.confidence.interval(fp_tosses_AR,bf$evaluate(3/4))

## ----bf_AR_verbose1------------------------------------------------------
fp_tosses_AR2 <- bf$sample.AR(n = 1000, roll.fun = toss_coin, verbose = TRUE)
print(paste0("Empirical number of tosses required: ",mean(fp_tosses_AR2$empirical_tosses)))

## ----bf_AR_verbose2------------------------------------------------------
fp_tosses_AR3 <- bf$sample.AR(n = 1000, true_p = c(3/4,1/4), verbose = TRUE)
print(paste0("Theoretical number of tosses required: ",mean(fp_tosses_AR3$theor_tosses)))

## ----de_polynomials------------------------------------------------------
f_dice <- list(
  list(c(1,sqrt(2)),c("120","022")),
  list(c(4,1/2,3),c("470","000","123")),
  list(c(7,2),matrix(c(1,3,4,0,0,2),byrow=TRUE,ncol=3))
)

## ----de_construction-----------------------------------------------------
de <- DiceEnterprise$new(f_dice)

## ----de_roll-------------------------------------------------------------
roll_die <- function(n) {
  sample(1:3, size = n, replace = TRUE, prob = c(1/5,1/4,1-1/5-1/4))
} #The original die has probability 1/5, 1/4, 11/20
sample_die <- de$sample(n = 10, roll.fun = roll_die)
print(sample_die)

## ----de_diagnosis--------------------------------------------------------
sample_die <- de$sample(n = 1000, roll.fun = roll_die, num_cores = 2, verbose = TRUE, double_time = TRUE) #Produces 1000 rolls of the f(p)-die, using 2 cores
print(table(sample_die[[1]])/1000) #Empirical probabilities. Notice that the theoretical ones are given by print(de$evaluate(c(1/5,1/4,1-1/5-1/4)))
print(de$evaluate(c(1/5,1/4,1-1/5-1/4)))
print(paste0("Average number of rolls required: ", mean(sample_die[[2]])))
plot.confidence.interval(sample_die[[1]],de$evaluate(c(1/5,1/4,1-1/5-1/4)))


## ----example_2_paper-----------------------------------------------------
de_ex2 <- DiceEnterprise$new(G=list(
  list(sqrt(2),"300"),
  list(1,"201"),
  list(1/4,"120"),
  list(2,"111"),
  list(1/2,"102"),
  list(3/4,"021")
))
#Define the original die
true_prob_original_die <- c(1/5,1/4,11/20)
roll_die <- function(n) {
  sample(1:3, size = n, replace = TRUE, prob = true_prob_original_die)
} 
#Get a sample of size 1000 from the multivariate ladder
#and plot the estimates with confidence intervals
set.seed(17)
sample_ex2 <- de_ex2$sample(n=1000, roll.fun = roll_die, verbose = TRUE)
print(paste0("Average number of rolls required: ", mean(sample_ex2[[2]])))
plot.confidence.interval(sample_ex2[[1]],de_ex2$evaluate(true_prob_original_die))


## ----ce_toss_all---------------------------------------------------------
toss.all.coins <- function(probs) { 
  return(sapply(probs, function(p) {sample(1:2, size = 1, prob = c(p,1-p))})) #1 or 2 (not 0))
}

## ----ce_toss_all_poly----------------------------------------------------
f_indep_coins1 <- list(
  list(rep(1,4),c("20000","10100","10010","00110")),
  list(rep(1,4),c("20000","11000","10010","01010")),
  list(rep(1,4),c("20000","11000","10100","01100"))
)

## ----ce_toss_all_initialize----------------------------------------------
ce1 <- CoinsEnterprise$new(f_indep_coins1, toss.coins = toss.all.coins, num_coins = 3, die_type = "toss_all")

## ----ce_toss_all_sample--------------------------------------------------
indep_coins_probs <- c(0.4,0.7,0.55)
result <- ce1$sample(n = 1000, num_cores = 2, verbose = TRUE,  double_time = FALSE, probs = indep_coins_probs) #the argument probs is passed to toss.coins
print(table(result[[1]])/1000) #Empirical probabilities. Notice that the theoretical ones are given by
print(indep_coins_probs/sum(indep_coins_probs))
print(paste0("Average number of rolls required: ", mean(result[[2]])))
print(paste0("Average number of tosses required: ", length(indep_coins_probs)*mean(result[[2]])))
plot.confidence.interval(result[[1]],indep_coins_probs/sum(indep_coins_probs))


## ----ce_toss_until_heads-------------------------------------------------
toss.until.heads <- function(probs) {
  m <- length(probs)
  res <- rep(NA, m)
  for(i in 1:m) {
    res[i] <- sample(1:2, size = 1, prob = c(probs[i],1-probs[i]))
    if(res[i] == 1) { break } #Stop loop if heads is obtained
  }
  return(res)
}

## ----ce_toss_until_heads_poly--------------------------------------------
f_indep_coins2 <- list(
  list(c(1,-2,1,-1,1),c("1000","2000","3000","1100","2100")),
  list(c(1,-1,-1),c("0100","1100","0200")),
  list(c(1,-1),c("0010","1010"))
)

## ----ce_toss_until_heads_initialize--------------------------------------
ce2 <- CoinsEnterprise$new(f_indep_coins2, toss.coins = toss.until.heads, num_coins = 3, die_type = "first_heads")

## ----ce_toss_until_heads_sample------------------------------------------
result <- ce2$sample(n = 1000, num_cores = 2, verbose = TRUE, double_time = FALSE, probs = indep_coins_probs) #the argument probs is passed to toss.coins
print(table(result[[1]])/1000) #Empirical probabilities. Notice that the theoretical ones are given by
print(indep_coins_probs/sum(indep_coins_probs))
print(paste0("Average number of rolls required: ", mean(result[[2]])))
print(paste0("Average number of tosses required: ", (indep_coins_probs[1]+
                                                       2*(1-indep_coins_probs[1])*indep_coins_probs[2]+
                                                       3*(1-indep_coins_probs[1])*(1-indep_coins_probs[2]))*
                                                       mean(result[[2]])))
plot.confidence.interval(result[[1]],indep_coins_probs/sum(indep_coins_probs))


## ----ce_unif_toss--------------------------------------------------------
toss.coins.single <- function(which_coin, probs) {
  return(sample(c(1,2), size = 1, prob = c(probs[which_coin], 1-probs[which_coin])))
}

## ----ce_unif_poly--------------------------------------------------------
f_indep_coins3 <- list(
  list(1, "1000"),
  list(1, "0100"),
  list(1, "0010")
)

## ----ce_unif_def---------------------------------------------------------
ce3 <- CoinsEnterprise$new(f_indep_coins3, toss.coins = toss.coins.single, num_coins = 3, die_type = "uniform")

## ----ce_unif_sample------------------------------------------------------
result <- ce3$sample(n = 1000, num_cores = 2, verbose = TRUE, double_time = FALSE, probs = indep_coins_probs) #the argument probs is passed to toss.coins.single
print(table(result[[1]])/1000) #Empirical probabilities. Notice that the theoretical ones are given by
print(indep_coins_probs/sum(indep_coins_probs))
print(paste0("Average number of tosses required: ", mean(result[[2]])))
plot.confidence.interval(result[[1]],indep_coins_probs/sum(indep_coins_probs))

## ----eliminiami, eval=FALSE,echo=FALSE-----------------------------------
#  roll.uniform <- function(n,probs_coins,probs_unif = rep(1,length(probs_coins))) {
#    m <- length(probs_coins)
#    res <- numeric(n)
#    for(k in 1:n) {
#    i <- sample(1:m, size = 1, prob = probs_unif)
#    res[k] <- (sample(c(i,m+1), size = 1, prob = c(probs_coins[i],1-probs_coins[i])))
#    }
#    return(res)
#  }
#  
#  indep_coins_probs <- seq(0.01,0.09,by=0.005)
#  f_unif <- vector("list", length = length(indep_coins_probs))
#  for(i in 1:length(indep_coins_probs)) {
#    zeros <- rep(0, length(indep_coins_probs)+1)
#    zeros[i] <- 1
#    f_unif[[i]] <- list(c(1), paste0(zeros,collapse=""))
#  }
#  
#  de_uniform <- DiceEnterprise$new(f_unif)
#  
#  size_sample <- 10000
#  sample_coins_uniform <- de_uniform$sample(n = size_sample, roll.fun = roll.uniform, probs_coins = indep_coins_probs, verbose = TRUE)
#  table(sample_coins_uniform[[1]])/size_sample
#  mean(sample_coins_uniform[[2]])
#  length(indep_coins_probs)/sum(indep_coins_probs)

## ----bug_to_correct, error=TRUE, echo = FALSE----------------------------
#(0.2501 - 1 x + 1 x^2)/((0.3 - 1 x + 1 x^2)+(0.2501 - 1 x + 1 x^2))
# bf_bug <- BernoulliFactory$new(
#   f_1 = list(c(0.2501,-1,1), c(0,1,2)),
#   f_2 = list(c(0.3,-1,1), c(0,1,2)),
# )
#(2-5x+5x^2)/(10x+10)
bf_bug <- DiceEnterprise$new(G = list(
  list(c(2,-5,5),c("00","10","20")),
  list(c(10,10),c("10","00"))
))

