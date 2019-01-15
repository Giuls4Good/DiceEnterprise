## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.width=7, fig.height=4,
cache=TRUE
)
rm(list=ls())
require("DiceEnterprise")
library(plyr)

## ----high_power_univariate_ladder----------------------------------------
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
bound_tosses <- matrix(NA, nrow = 2, ncol = increase_degree+1)
colnames(emp_tosses) <- paste0("+",0:increase_degree)
colnames(emp_res) <- paste0("+",0:increase_degree)
colnames(bound_tosses) <- paste0("+",0:increase_degree)
rownames(bound_tosses) <- c("p=1","p=0")
for(i in 0:increase_degree) {
  #Check if it already exists
  # if(file.exists(paste0("efficiency_uni_data/bound_a_",a,"_degree_",i,".rds"))) {
  #   cat("Restoring previous bound for degree +",i," \n")
  #   prev_work <- readRDS(paste0("efficiency_uni_data/bound_a_",a,"_degree_",i,".rds"))
  #   bound_tosses[,i+1] <- prev_work$bound_tosses
  # } else {
  #   cat("Getting bound tosses for +",i," \n") 
  #   bound_tosses[1,i+1] <- de_list[[i+1]]$expected.tosses.extreme(1)
  #   bound_tosses[2,i+1] <- de_list[[i+1]]$expected.tosses.extreme(0)
  #   #Save object
  #   saveRDS(list(bound_tosses=bound_tosses[,i+1]), file = paste0("efficiency_uni_data/bound_a_",a,"_degree_",i,".rds"))
  # }
  
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

cat(bound_tosses,"\n ****** \n")
eff_cond <- sapply(de_list, function(x) {x$get.efficiency.condition()})
names(eff_cond) <- paste0("+",0:increase_degree)
print(eff_cond)
apply(emp_tosses,2,mean)

## ----high_power_multivariate_ladder--------------------------------------
rm(list=ls())
require("DiceEnterprise")

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

