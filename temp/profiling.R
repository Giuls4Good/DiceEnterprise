rm(list=ls())
require(DiceEnterprise)
require(profvis)
require(rbenchmark)
toss.until.heads <- function(probs) {
  m <- length(probs)
  res <- rep(NA, m)
  for(i in 1:m) {
    res[i] <- sample(1:2, size = 1, prob = c(probs[i],1-probs[i]))
    if(res[i] == 1) { break } #Stop loop if heads is obtained
  }
  return(res)
}

set.seed(17)
indep_coin_probs <- matrix(c(0.00001,0.00005), ncol = 2, byrow= TRUE)
const <- matrix(c(1,1), ncol=2, byrow = TRUE) #c1, c2
sample_size <- 100 #Number of tosses
num_cores <- 1 #Number of cores used (parallel computing is not supported on Windows)

i <- 1
ce_2coins_until_heads <- CoinsEnterprise$new(list(
  list(c(const[i,1],-const[i,1]), c("100","200")),
  list(c(const[i,2]), c("010"))
), toss.coins = toss.until.heads, die_type = "first_heads")

profvis({
  sample_2coins_untilheads <- ce_2coins_until_heads$sample(n = sample_size, num_cores = num_cores, verbose = TRUE,
                                                         double_time = FALSE, probs = indep_coin_probs[i,])
}) #Most of the time is spent in the update function

#Test findInterval
for(j in 1:20) {
  print(findInterval(j,c(5,10,11,12)) == findIntervalSingle(j,c(5,10,11,12)))
}

ce_2coins_until_heads <- CoinsEnterprise$new(list(
  list(c(const[i,1],-const[i,1]), c("100","200")),
  list(c(const[i,2]), c("010"))
), toss.coins = toss.until.heads, die_type = "first_heads")

set.seed(17)
size <- 100
for(k in 1:1000) {
  B <- sample(1:3, size = size, replace = TRUE)
  U <- runif(size)
  res_cpp <- ce_2coins_until_heads$get.ladder.fine.connected()$update.fun.Cpp(1,B,U)
  res_R <- ce_2coins_until_heads$get.ladder.fine.connected()$update.fun(1,B,U)
  if(res_cpp != res_R) {
    print("Error")
    break
  }
}

size <- 10000
B <- sample(1:3, size = size, replace = TRUE)
U <- runif(size)
benchmark(
  ce_2coins_until_heads$get.ladder.fine.connected()$update.fun.Cpp(1,B,U),
  ce_2coins_until_heads$get.ladder.fine.connected()$update.fun(1,B,U), replications = 100
)


#################
rm(list=ls())
require(DiceEnterprise)
require(profvis)
require(rbenchmark)
set.seed(1007)

de <- DiceEnterprise$new(list(
  list(c(1,-1), c("100","200")),
  list(c(1), c("010"))
))
ladder <- de$get.ladder.fine.connected()

CFTPCpp(4, c(0.1,0.4,0.5), TRUE, TRUE,
        ladder$get.P.cumsum(), ladder$get.P.moves.list(), FALSE, 0, 0, TRUE)

######
rm(list=ls())
require(DiceEnterprise)
require(profvis)
require(rbenchmark)
set.seed(1007)

true_p <- c(3/4,1/4)
sample_size <- 10000
num_cores <- 2

toss_coin <- function(n, prob) {
  sample(1:2, size = n, replace = TRUE, prob = prob)
}
f_1 <- list(coeff = c(sqrt(2)), power = c(3)) #f(p)
f_2 <- list(coeff = c(-5,11,-9,3), power = c(3,2,1,0)) #1-f(p)
bf <- BernoulliFactory$new(f_1 = f_1, f_2 = f_2) #f_1 = f(p), f_2 = 1-f(p)
fp_tosses <- bf$sample(n = sample_size, roll.fun = toss_coin, verbose = TRUE, num_cores = num_cores,
                       prob = true_p)
aux <- bf$get.ladder.fine.connected()
CFTPCpp(aux$get.k(), true_p, TRUE, TRUE,
        aux$get.P.cumsum(), aux$get.P.moves.list(), monotonic = TRUE, min = 1, max = aux$get.k(), verbose = TRUE)

print(bf$evaluate(true_p[1]))
print(table(fp_tosses[[1]])/sample_size)
print(table(fp_tosses_cpp[[1]])/sample_size)
print(mean(fp_tosses[[2]]))
print(mean(fp_tosses_cpp[[2]]))

benchmark(bf$sample(n = sample_size, roll.fun = toss_coin, verbose = TRUE, num_cores = num_cores,
                    prob = true_p),
         bf$sample(n = sample_size, true_p = true_p, verbose = TRUE, num_cores = num_cores) )
