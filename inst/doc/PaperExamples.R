## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7, fig.height=4,
  cache=TRUE
)
require("DiceEnterprise")
require("pander")
require("ggplot2")
require("xtable")
set.seed(10,"L'Ecuyer-CMRG")

## ----bf_settings---------------------------------------------------------
sample_size <- 1000 #Number of tosses
probs_vec <- c(0.01,0.1,0.25,0.5,0.75,0.9,0.99) #Vector of true probabilities of the p-coin
num_cores <- 2 #Number of cores used (parallel computing is not supported on Windows)
conf_level <- 0.95 #Confidence interval level

## ----bf_toy, echo=FALSE, results='asis'----------------------------------
toss.coin <- function(n, p) {
  sample(1:2, size = n, replace = TRUE, prob = c(p,1-p))
}
bf <- BernoulliFactory$new(f_1 = list(coeff = c(sqrt(2)), power = c(3)), 
                           f_2 = list(coeff = c(-5,11,-9,3), power = c(3,2,1,0))) 

bf.res <- function(toss.coin, bf) {
  bf_res <- matrix(NA, nrow = 6, ncol = length(probs_vec))
  rownames(bf_res) <- c("f(p)", "LowerCI", "Emp.f(p)", "UppCI", "Exp.TossesNoDoubling", "Exp.TossesDoubling")
  colnames(bf_res) <- as.character(probs_vec)
  for(i in 1:length(probs_vec)) {
    sample_res <- bf$sample(n = sample_size, roll.fun = toss.coin, p = probs_vec[i], num_cores = num_cores, verbose = TRUE, double_time = FALSE)
    sample_res_double_time <- bf$sample(n = sample_size, roll.fun = toss.coin, p = probs_vec[i], num_cores = num_cores, verbose = TRUE, double_time = TRUE)
    #Get lower and upper CI
    if(isTRUE(all.equal(sample_res[[1]],rep(1,sample_size)))) {
      #All equal to 1 (Heads)
      lower_ci <- 1
      upper_ci <- 1
    } else if(isTRUE(all.equal(sample_res[[1]],rep(2,sample_size)))) {
      #All equal to 2 (Tails)
      lower_ci <- 0
      upper_ci <- 0
    } else {
      ci <- MultinomCI(table(sample_res[[1]]),conf.level = conf_level)
      lower_ci <- ci[1,2]
      upper_ci <- ci[1,3]
    }
    bf_res[1,i] <- round(bf$evaluate(p = probs_vec[i])[1],2)
    bf_res[2,i] <- round(lower_ci,2)
    bf_res[3,i] <- round(1-mean(sample_res[[1]]-1),2)
    bf_res[4,i] <- round(upper_ci,2)
    bf_res[5,i] <- round(mean(sample_res[[2]]),2)
    bf_res[6,i] <- round(mean(sample_res_double_time[[2]]),2)
  }
  return(bf_res)
}
bf_res <- bf.res(toss.coin,bf)
pandoc.table(bf_res)

## ----bf_latex, eval = FALSE, echo = FALSE--------------------------------
#    xtable(bf_res, align = c("l|",rep("c",ncol(bf_res))))

## ----bf_efficiency, echo=FALSE, results='asis'---------------------------
R <- c(1,1000,1,500,1)
max_increase_degree <- 2 #How much the degree of the original ladder is increased by
original_ladder <- DiceEnterprise$new(list(
  list(R[1], "04"),
  list(R[2], "13"),
  list(R[3], "22"),
  list(R[4], "31"),
  list(R[5], "40")
))
#Sample and get the empirical number of tosses required
bf.res.eff <- function(toss.coin, original_ladder,max_increase_degree) {
  #Prepare increased degree ladders
  increased_ladders <- vector("list", length= max_increase_degree)
  for(i in 1:max_increase_degree) {
    increased_ladders[[i]] <- original_ladder$increase.degree(d=i)
  }
  #Prepare output
  bf_res <- matrix(NA, nrow = max_increase_degree+1, ncol = length(probs_vec)+1)
  rownames(bf_res) <- paste0("+",0:max_increase_degree)
  colnames(bf_res) <- c(as.character(probs_vec),"Eff.Condition") #The last column is true if the efficiency condition is met
  tot_rows <- length(probs_vec)*(max_increase_degree+1)
  df_res <- data.frame(p=rep(NA,tot_rows),
                     degree = rep(NA, tot_rows),
                     empirical_tosses = rep(NA, tot_rows))
  counter <- 1
  #Compute
  for(i in 1:length(probs_vec)) {
    #Check if the result already exists
    bf_res[1,i] <- mean(original_ladder$sample(n=sample_size, roll.fun = toss.coin, p = probs_vec[i], num_cores = num_cores, verbose = TRUE, double_time = FALSE)$exp_rolls)
    bf_res[1,ncol(bf_res)] <- as.numeric(original_ladder$get.efficiency.condition())
    df_res[counter,] <- c(probs_vec[i],0,bf_res[1,i])
    counter <- counter + 1
    for(deg in 1:max_increase_degree) {
      bf_res[deg+1,i] <- mean(increased_ladders[[deg]]$sample(n=sample_size, roll.fun = toss.coin, p = probs_vec[i], num_cores = num_cores, verbose = TRUE, double_time = FALSE)$exp_rolls)
      bf_res[deg+1,ncol(bf_res)] <- as.numeric(increased_ladders[[deg]]$get.efficiency.condition())
      df_res[counter,] <- c(probs_vec[i],deg,bf_res[deg+1,i])
      counter <- counter + 1
    }
  }
  return(list(bf_res,df_res))
}
if(!file.exists(paste0("saved_data/efficiency_mono_cftp_example_n_",sample_size,".rds"))) {
  aux <- bf.res.eff(toss.coin,original_ladder,max_increase_degree)
  #Save the result
  saveRDS(list(aux=aux),file=paste0("saved_data/efficiency_mono_cftp_example_n_",sample_size,".rds"))
} else {
  aux <- readRDS(paste0("saved_data/efficiency_mono_cftp_example_n_",sample_size,".rds"))$aux
}
bf_res <- aux[[1]]
print(bf_res)
pandoc.table(bf_res)

## ----de_setup------------------------------------------------------------
sample_size <- 1000 #Number of rolls
probs_mat <- matrix(c(1/5,1/4,11/20,
                      1/10,4/10,5/10), ncol = 3, byrow = TRUE) #Matrix of true probabilities of the die
num_cores <- 4 #Number of cores used (parallel computing is not supported on Windows)
conf_level <- 0.95 #Confidence interval level

## ----de, eval=TRUE, echo = FALSE, result='asis'--------------------------

roll_die <- function(n, probs_die) {
  sample(1:3, size = n, replace = TRUE, prob = probs_die)
} 
de <-  DiceEnterprise$new(G=list(
  list(sqrt(2),"300"),
  list(1,"201"),
  list(1/4,"120"),
  list(2,"111"),
  list(1/2,"102"),
  list(3/4,"021")
))

de_res <- matrix(NA, nrow = 6, ncol = nrow(probs_mat))
rownames(de_res) <- c("f(p)", "LowerCI", "Emp.f(p)", "UppCI", "Exp.RollsNoDoubling", "Exp.RollsDoubling")
colnames(de_res) <- apply(probs_mat, 1, function(x) { paste0("(",paste0(x, collapse = ", "),")")} )

for(i in 1:nrow(probs_mat)) {
  de_sample_double_time <- de$sample(n=sample_size, roll.fun = roll_die, double_time = TRUE, num_cores = num_cores, verbose = TRUE, probs_die = probs_mat[i,])
  de_sample <- de$sample(n=sample_size, roll.fun = roll_die, double_time = FALSE, num_cores = num_cores, verbose = TRUE, probs_die = probs_mat[i,])
  if(length(unique(de_sample[[1]])) > 1) {
    #There are more outcomes
    ci <- MultinomCI(table(de_sample[[1]]),conf.level = conf_level)
    ci_aux <- matrix(NA, nrow = 2, ncol = 6) #6 possible outcomes of the die
    for(j in 1:nrow(ci)) {
      ci_aux[1,as.integer(rownames(ci)[j])] <- ci[j,2]
      ci_aux[2,as.integer(rownames(ci)[j])] <- ci[j,3]
    }
    de_res[2,i] <- paste0("(",paste0(round(ci_aux[1,],2), collapse = ", "),")") #Lower confidence interval
    de_res[4,i] <- paste0("(",paste0(round(ci_aux[2,],2), collapse = ", "),")")  #Uppter confidence interval
  }
  de_res[1,i] <- paste0("(",paste0(round(de$evaluate(p = probs_mat[i,]),2), collapse = ", "),")")
  de_res[3,i] <- paste0("(",paste0(round(table(de_sample[[1]])/sample_size,2), collapse = ", "),")")
  de_res[5,i] <- round(mean(de_sample[[2]]),2)
  de_res[6,i] <- round(mean(de_sample_double_time[[2]]),2)
}
pandoc.table(de_res,split.tables=Inf)

## ----high_power_univariate_ladder, results='hide', cache=FALSE-----------
rm(list=ls())
require("DiceEnterprise")

p <- 1/2 #True value of p
a <- 10 #Exponent
n <- 1000 #Size sample CFTP
num_cores <- 4 #Number of cores used
initial_seed <- 17 #Seed for RNG
increase_degree <- 260 #How much the degree is increased
de_list <- vector("list", length = increase_degree+1) #Store all the de objects

set.seed(initial_seed, "L'Ecuyer-CMRG") 
#Construct the ladder
cat("Constructing the dice enterprise \n")
  de <- DiceEnterprise$new(G=list(
    list(1,matrix(c(2*a,0),nrow=1,ncol=2)),
    list(1,matrix(c(a,a),nrow=1,ncol=2)),
    list(1,matrix(c(0,2*a),nrow=1,ncol=2))
  ))
  de_list[[1]] <- de
  cat("Constructing the augmented dice enterprises \n")
  for(i in 1:increase_degree) {
    #Check if it already exists
    if(file.exists(paste0("efficiency_uni_data/ladder_a_",a,"_degree_",i,".rds"))) {
      cat("Restoring previous work for degree +",i," \n")
      prev_work <- readRDS(paste0("efficiency_uni_data/ladder_a_",a,"_degree_",i,".rds"))
      de_list[[i+1]]  <- prev_work$de_aug
    } else {
      stop("WHAT?!")
      cat("Constructing the augmented dice enterprises for +",i," \n")
      if(i == 1) {
        de_list[[i+1]] <- de$increase.degree(1)
      } else {
        de_list[[i+1]] <- de_list[[i-1]]$increase.degree(1)
      }
      #Save object
      saveRDS(list(de_aug=de_list[[i+1]],i=i), file = paste0("efficiency_uni_data/ladder_a_",a,"_degree_",i,".rds"))
    }
  
}

#For each ladder get empirical number of tosses and bound
cat("Computing empirical number of tosses and bound \n")
emp_tosses <- matrix(NA, nrow = n, ncol = increase_degree+1)
emp_res <- matrix(NA, nrow = n, ncol = increase_degree+1)
colnames(emp_tosses) <- paste0("+",0:increase_degree)
colnames(emp_res) <- paste0("+",0:increase_degree)
for(i in 0:increase_degree) {
  #Check if it already exists
  if(file.exists(paste0("efficiency_uni_data/n_",n,"_a_",a,"_degree_",i,".rds"))) {
    cat("Restoring previous work for degree +",i," \n")
    prev_work <- readRDS(paste0("efficiency_uni_data/n_",n,"_a_",a,"_degree_",i,".rds"))
    emp_tosses[,i+1] <- prev_work$emp_tosses
    emp_res[,i+1] <- prev_work$emp_res
  } else {
    cat("Getting empirical tosses for +",i," \n")
    aux <- de_list[[i+1]]$sample(n=n, true_p = c(p,1-p), num_cores = num_cores, verbose = TRUE)
    emp_tosses[,i+1] <- aux[[2]]
    emp_res[,i+1] <- aux[[1]]
    #Save object
    saveRDS(list(emp_tosses=aux[[2]],emp_res=aux[[1]]), file = paste0("efficiency_uni_data/n_",n,"_a_",a,"_degree_",i,".rds"))
  }
}

## ----high_power_univariate_ladder_print----------------------------------
#Prepare output
res_output <- matrix(NA, nrow = 2, ncol = increase_degree+1)
colnames(res_output) <- paste0("+",0:increase_degree)
rownames(res_output) <- c("Exp.Tosses","Efficiency.Condition")
res_output[1,] <- apply(emp_tosses,2,mean)
res_output[2,] <- sapply(de_list, function(x) {x$get.efficiency.condition()})
print(res_output)

## ----high_power_multivariate_ladder, results='hide', cache=FALSE---------
rm(list=ls())
require("DiceEnterprise")
setwd("~/Dropbox/OxWaSP/Dice Enterprise/R/DiceEnterprise/vignettes") #Change accordingly

true_p <- c(1/3,1/3,1/3) #True value of p
a <- 5 #Exponent
n <- 1000 #Size sample CFTP
num_cores <- 4 #Number of cores used
initial_seed <- 7117 #Seed for RNG
increase_degree <- 78 #How much the degree is increased
de_list <- vector("list", length = increase_degree+1) #Store all the de objects

set.seed(initial_seed, "L'Ecuyer-CMRG") 
#Construct the ladder
cat("Constructing the dice enterprise \n")
  de <- DiceEnterprise$new(G=list(
    list(1,matrix(c(3*a,0,0),nrow=1,ncol=3)),
    list(1,matrix(c(a,a,a),nrow=1,ncol=3)),
    list(1,matrix(c(0,0,3*a),nrow=1,ncol=3)),
    list(1,matrix(c(0,3*a,0),nrow=1,ncol=3))
  ))
  de_list[[1]] <- de
  cat("Constructing the augmented dice enterprises \n")
  for(i in 1:increase_degree) {
    #Check if it already exists
    if(file.exists(paste0("efficiency_multi_data/ladder_n_",n,"_a_",a,"_degree_",i,".rds"))) {
      cat("Restoring previous work for degree +",i," \n")
      prev_work <- readRDS(paste0("efficiency_multi_data/ladder_n_",n,"_a_",a,"_degree_",i,".rds"))
      de_list[[i+1]]  <- prev_work$de_aug
    } else {
      cat("Constructing the augmented dice enterprises for +",i," \n")
      if(i == 1) {
        de_list[[i+1]] <- de$increase.degree(1)
      } else {
        de_list[[i+1]] <- de_list[[i-1]]$increase.degree(1)
      }
      #Save object
      saveRDS(list(de_aug=de_list[[i+1]],i=i), file = paste0("efficiency_multi_data/ladder_n_",n,"_a_",a,"_degree_",i,".rds"))
    }
  }
  

#For each ladder get empirical number of tosses
cat("Computing empirical number of tosses \n")
emp_tosses <- matrix(NA, nrow = n, ncol = increase_degree+1)
emp_res <- matrix(NA, nrow = n, ncol = increase_degree+1)
colnames(emp_tosses) <- paste0("+",0:increase_degree)
colnames(emp_res) <- paste0("+",0:increase_degree)
for(i in 0:increase_degree) {
  #Check if it already exists
  if(file.exists(paste0("efficiency_multi_data/n_",n,"_a_",a,"_degree_",i,".rds"))) {
    cat("Restoring previous work for degree +",i," \n")
    prev_work <- readRDS(paste0("efficiency_multi_data/n_",n,"_a_",a,"_degree_",i,".rds"))
    emp_tosses[,i+1] <- prev_work$emp_tosses
    emp_res[,i+1] <- prev_work$emp_res
  } else {
    cat("Getting empirical tosses for +",i," \n")
    aux <- de_list[[i+1]]$sample(n=n, true_p = true_p, num_cores = num_cores, verbose = TRUE)
    emp_tosses[,i+1] <- aux[[2]]
    emp_res[,i+1] <- aux[[1]]
    #Save object
    saveRDS(list(emp_tosses=aux[[2]],emp_res=aux[[1]]), file = paste0("efficiency_multi_data/n_",n,"_a_",a,"_degree_",i,".rds"))
  }
}

apply(emp_tosses,2,mean)

## ----high_power_multivariate_ladder_print--------------------------------
#Prepare output
res_output <- matrix(NA, nrow = 1, ncol = increase_degree+1)
colnames(res_output) <- paste0("+",0:increase_degree)
res_output[1,] <- apply(emp_tosses,2,mean)
print(res_output)

