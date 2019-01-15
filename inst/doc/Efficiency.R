## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7, fig.height=4,
  cache=TRUE
)
require("DiceEnterprise")
library(plyr)

## ----setup_efficiency, echo = FALSE--------------------------------------
p_seq <- c(0.01,0.1,0.25,0.5) #Sequence of true values of p
a_seq <- 1:10 #Sequence of a
size_sample_CFTP <- 1000 #Size sample CFTP
size_sample_est_p <- 1000 #Size tosses to get an estimate of p
num_cores <- 2 #Number of cores used
initial_seed <- 17 #Seed for RNG

generate.data <- function() {
  set.seed(initial_seed, "L'Ecuyer-CMRG") 
  data_out <- vector("list", length = length(a_seq))
  counter <- 1
  for(a in a_seq) {
    cat("Started a = ",a,"\n")
    if(file.exists(paste0("efficiency_a_",a,".rds"))) {
      #Recover previous work
      cat("Recover previous work for a = ",a,"\n")
      work_prev <- readRDS(paste0("efficiency_a_",a,".rds"))
      data_out[[counter]] <- work_prev$data_out
      counter <- counter+ 1
      seed <- work_prev$seed
      set.seed(seed, "L'Ecuyer-CMRG")  #Restore random seed
    } else {
      #Construct the ladder
      cat("Constructing the ladder\n")
      de <- DiceEnterprise$new(G=list(
        list(1,matrix(c(a,0),nrow=1,ncol=2)),
        list(1,matrix(c(0,a),nrow=1,ncol=2))
      ))
      #Get expected value for A-R
      cat("Get expected tosses of AR\n")
      exp_tosses_AR <- sapply(p_seq, function(p) {
        C_p = p^a+(1-p)^a
        Q <- 1
        return(Q*a/C_p)
      })
      cat("Get empirical expected tosses of naive CFTP\n")
      #Get value for naive CFTP
      tosses_CFTP <- lapply(p_seq, function(p) {
        return(1)
        return(de$sample(n=size_sample_CFTP, true_p = c(p,1-p), num_cores = num_cores, verbose = TRUE)[[2]])
      })
      cat("Get empirical expected tosses of augmented CFTP\n")
      #Get an estimate of p
      p_est_seq <- sapply(p_seq, function(p) {
        mean(sample(c(1,0), prob = c(p,1-p), size = size_sample_est_p, replace = TRUE))
      })
      #Get a value for augmented CFTP
      tosses_augmented_CFTP <- lapply(1:length(p_est_seq), function(i){
        p_est <- p_est_seq[i]
        p <- p_seq[i]
        #Construct augmented ladder
        new_de <- de$impose.efficiency(method="bound",true_p = c(p_est,1-p_est))
        saveRDS(list(new_de=new_de,a=a,p_est=p_est,p=p), file = paste0("augmented_de_a_",a,"_pest_",p_est,".rds"))
        return(new_de$sample(n=size_sample_CFTP, true_p = c(p,1-p), num_cores = num_cores, verbose = TRUE)[[2]])
      })

      #Prepare output
      data_out[[counter]] <- list(a = a,
                                p_seq = p_seq,
                                tosses_CFTP = tosses_CFTP,
                                exp_tosses_AR = exp_tosses_AR,
                                tosses_augmented_CFTP = tosses_augmented_CFTP) 
      #Generate new random seed
      seed <- ceiling(runif(1)*1e8)
      set.seed(seed, "L'Ecuyer-CMRG") 
      
      #Save to outer file
      cat("Saving result for a = ",a,"\n")
      saveRDS(list(data_out=data_out[[counter]],seed=seed), file = paste0("efficiency_a_",a,".rds"))
      counter <- counter + 1
    }
  }
  return(data_out)
}
generate.data.frame <- function(data_out) {
  #Prepare output
  tot_rows <- length(a_seq)*length(p_seq)*3
  res <- data.frame(p=rep(NA,tot_rows),
                    a=rep(NA,tot_rows), 
                    meanTosses=rep(NA,tot_rows),
                    maxMeanTosses=rep(NA,tot_rows),
                    minMeanTosses=rep(NA,tot_rows),
                    method=rep(NA,tot_rows))
  counter <- 1
  for(i in 1:length(data_out)) {
    a <- data_out[[i]]$a
    p_seq_saved <- data_out[[i]]$p_seq
    exp_tosses_AR <- data_out[[i]]$exp_tosses_AR
    exp_tosses_CFTP <- sapply( data_out[[i]]$tosses_CFTP, mean)
    max_exp_tosses_CFTP <- sapply( data_out[[i]]$tosses_CFTP, max)
    min_exp_tosses_CFTP <- sapply( data_out[[i]]$tosses_CFTP, min)
    exp_tosses_augmented_CFTP <- sapply( data_out[[i]]$tosses_augmented_CFTP, mean)
    max_exp_tosses_augmented_CFTP <- sapply( data_out[[i]]$tosses_augmented_CFTP, max)
    min_exp_tosses_augmented_CFTP <- sapply( data_out[[i]]$tosses_augmented_CFTP, min)

    #Sanity check
    if(!(a %in% a_seq)) {
      stop("There's a conflict between the setup and saved results (a is not in a_seq).")
    }
    if(!isTRUE(all.equal(p_seq_saved,p_seq))) {
      stop("There's a conflict between the setup and saved results (saved p_seq is not equal to p_seq).")
    }
    #Save to data frame
    for(k in 1:length(exp_tosses_CFTP)) {
      res[counter,] <- c(p_seq_saved[k],a,exp_tosses_CFTP[k],max_exp_tosses_CFTP[k],min_exp_tosses_CFTP[k],1)
      res[counter+1,] <- c(p_seq_saved[k],a,exp_tosses_augmented_CFTP[k],max_exp_tosses_augmented_CFTP[k],min_exp_tosses_augmented_CFTP[k],2)
      res[counter+2,] <- c(p_seq_saved[k],a,exp_tosses_AR[k],exp_tosses_AR[k],exp_tosses_AR[k],3)
      counter <- counter + 3
    }
  }
  res$method <- as.factor(res$method)
  res$method <- revalue(res$method, c("1" = "NaiveCFTP", "2" = "OptCFTP" , "3" = "AR"))
  return(res)
}

data_out <- generate.data()
res <- generate.data.frame(data_out)

ggplot(data = res[which(res$p == 0.01),], aes(x = factor(a), y = meanTosses, colour = factor(method))) + ggtitle("True p: 0.01") +
  geom_line(aes(group = factor(method))) + geom_point()
ggplot(data = res[which(res$p == 0.1),], aes(x = factor(a), y = meanTosses, colour = factor(method))) + ggtitle("True p: 0.1") +
  geom_line(aes(group = factor(method))) + geom_point()
ggplot(data = res[which(res$p == 0.25),], aes(x = factor(a), y = meanTosses, colour = factor(method))) + ggtitle("True p: 0.25") +
  geom_line(aes(group = factor(method))) + geom_point()
ggplot(data = res[which(res$p == 0.5),], aes(x = factor(a), y = meanTosses, colour = factor(method))) + ggtitle("True p: 0.5") +
  geom_line(aes(group = factor(method))) +  geom_point()
  #geom_ribbon(aes(ymin = minMeanTosses, ymax = maxMeanTosses, fill = factor(method)), alpha = 0.2, color = NA)


## ----old, echo=FALSE, eval=FALSE-----------------------------------------
#  # for(a in a_seq) {
#  #   cat("Started a = ",a,"\n")
#  #   if(file.exists(paste0("efficiency_a_",a,".rds"))) {
#  #     #Recover previous work
#  #     cat("Recover previous work for a = ",a,"\n")
#  #     work_prev <- readRDS(paste0("efficiency_a_",a,".rds"))
#  #     res_prev <- work_prev$res
#  #     seed <- work_prev$seed
#  #     res[counter:(counter+nrow(res_prev)-1),] <- res_prev
#  #     counter <- counter + nrow(res_prev)
#  #     .Random.seed <- seed #Restore random seed
#  #   } else {
#  #     #Construct the ladder
#  #     cat("Constructing the ladder\n")
#  #     de <- DiceEnterprise$new(G=list(
#  #       list(1,matrix(c(a,0),nrow=1,ncol=2)),
#  #       list(1,matrix(c(0,a),nrow=1,ncol=2))
#  #     ))
#  #     #Get expected value for A-R
#  #     cat("Get expected tosses of AR\n")
#  #     exp_tosses_AR <- sapply(p_seq, function(p) {
#  #       C_p = p^a+(1-p)^a
#  #       Q <- 1
#  #       return(Q*a/C_p)
#  #     })
#  #     cat("Get empirical expected tosses of naive CFTP\n")
#  #     #Get value for naive CFTP
#  #     tosses_CFTP <- lapply(p_seq, function(p) {
#  #       return(de$sample(n=size_sample_CFTP, true_p = c(p,1-p), num_cores = num_cores, verbose = TRUE)[[2]])
#  #     })
#  #     exp_tosses_CFTP <- sapply(tosses_CFTP, mean)
#  #     cat("Get empirical expected tosses of augmented CFTP\n")
#  #     #Get an estimate of p
#  #     p_est_seq <- sapply(p_seq, function(p) {
#  #       mean(sample(c(1,0), prob = c(p,1-p), size = size_sample_est_p, replace = TRUE))
#  #     })
#  #     #Get a value for augmented CFTP
#  #     tosses_augmented_CFTP <- sapply(1:length(p_est_seq), function(i){
#  #       p_est <- p_est_seq[i]
#  #       p <- p_seq[i]
#  #       #Construct augmented ladder
#  #       new_de <- de$impose.efficiency(method="bound",true_p = c(p_est,1-p_est))
#  #       return(new_de$sample(n=size_sample_CFTP, true_p = c(p,1-p), num_cores = num_cores, verbose = TRUE)[[2]])
#  #     })
#  #     exp_tosses_augmented_CFTP <- sapply(tosses_augmented_CFTP, mean)
#  #     #Save in dataframe
#  #     initial_counter <- counter
#  #     for(i in 1:length(exp_tosses_CFTP)) {
#  #       res[counter,] <- c(p_seq[i],a,exp_tosses_CFTP[i],1)
#  #       res[counter+1,] <- c(p_seq[i],a,exp_tosses_augmented_CFTP[i],2)
#  #       res[counter+2,] <- c(p_seq[i],a,exp_tosses_AR[i],3)
#  #       counter <- counter + 3
#  #     }
#  #
#  #     #Save to outer file
#  #     cat("Saving result for a = ",a,"\n")
#  #     res_temp <- res[initial_counter:(counter-1),]
#  #     saveRDS(list(res=res_temp,seed=.Random.seed), file = paste0("efficiency_a_",a,".rds"))
#  #   }
#  # }
#  #
#  # res$method <- as.factor(res$method)
#  # res$method <- revalue(res$method, c("1" = "NaiveCFTP", "2" = "OptCFTP" , "3" = "AR"))

