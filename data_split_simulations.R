library('permute')
library('gbm')
library('randomForestSRC')
library('party')
library('mboost')
library("survival")
library('MachineShop')
# to install required packages, use:
# install.packages(c("MachineShop","gbm","randomForestSRC","party","mboost"))

usageMsg <- "Valid usage:\n\tSetup a full run:\tRscript data_split_simulations.R setup\n\tRun simulation number 1:\tRscript data_split_simulations.R 1\n\tAggregate final results of a completed run:\tRscript data_split_simulations.R final"

args <- commandArgs(trailingOnly = TRUE)
if(length(args)==0 | length(args)>1) {
  stop(usageMsg,call.=FALSE)
}

## ML settings:
# these can be modified
model_list <- c("GLMBoostModel","RFSRCModel","CForestModel","BlackBoostModel")
metric_list <- c("cindex","mse","mae")
# these should not be modified
method_list <- c("cv","dmcv","ho","dmho")
finals_list <- c("rse","rb")
fo <- Surv(times,status) ~ .

# Simulation settings:
num_reps <- 30 # repetitions per simulated outcome
num_vars <- 20 # number of variables contributing to simulated outcomes
num_points <- 500 # number of data points to sample in each repetition
num_sims <- 1000 # total number of simulations
cv_folds <- 10 # number of cross validation folds
out_dir <- "simulation_output" # directory to store results and simulated data
sim_data_file <- paste0(out_dir,"/","simulated_features.rds") # file to store simulated data
#mRNA_file = "TCGA-HNSC/gdac.broadinstitute.org_HNSC.mRNAseq_Preprocess.Level_3.2016012800.0.0/HNSC.uncv2.mRNAseq_RSEM_all.txt" # file path to mRNA expr data
mRNA_file = "HNSC.uncv2.mRNAseq_RSEM_all.txt" # file path to mRNA expr data
TCGA-like = TRUE # should data generation use TCGA-derived weibull parameters?

final_filename <- paste0(out_dir,"/final_output.rds")
result_filename <- function(i) {
  return(paste0(out_dir,"/simulation_results_",i,".rds"))
}



# Distribution Matching: needs a better name
#   advanced distribution-mimicking/matching? (adm?) - too close to automated decision making :(
#   "limited data" non-tractable
#   small data sampling (sds) ?
#   test-train sampling for extensive feature engineering in contexts with limited data and high-dimensionality


## function to split data
# n: number of data points
# k: number of folds
# returns indices of k folds
get_fold_indices <- function(n,k){
  ind <- c()
  while (length(ind) < n) { # add indices k at a time
    tempset <- shuffle(k)
    diff <- n-length(ind)
    if (diff < k) { # if we don't need the whole set, trim it down
      tempset <- tempset[1:diff]
    }
    ind <- append(ind,tempset)
  }
  return(ind)
}

# to use indices with your data:
# sort the data on the regression/time-to-event outcome value
# get indices, use them to select from sorted data

# e.g.:
# k <- 5 # num folds
# n <- 29 # num data points
# a <- get_fold_indices(n,k)
# for (fold in 1:k) { # loop for cross validation
#   train_ind <- a == fold
#   test_ind <- a != fold
# }


# Function to calculate standard error of a vector
std_err <- function(x) sd(x)/sqrt(length(x))

# Function to plot density histograms for data generation testing
denshist <- function(data, title) {
  # Get the density estimate
  dens=density(data)
  # Plot y-values scaled by number of observations against x values
  plot(dens$x,length(data)*dens$y,type="l",xlab="Survival Time",ylab="Count estimate",xlim = c(-500,3000),main=title)
}


# Function to normalize a vector
normalize_col <- function(column) {
  return(column/max(column))
}

## SETUP SECTION
if (args[1]=="setup") {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  set.seed(9001)
  if (!file.exists(sim_data_file)) { # check if data has already been generated
    # load in TCGA-HNSC mRNA tumor data
    full_data_load <- read.csv(mRNA_file, sep = "\t")
    full_data <- data.frame(t(full_data_load[-1]))
    colnames(full_data) <- full_data_load[, 1]
    
    # select only rows with primary tumor samples
    tumor_data <- full_data[grep("^TCGA.........01",rownames(full_data)),]
    rownames(tumor_data) <- gsub("^(TCGA........).*","\\1",rownames(tumor_data))
    # matching IDs to data split format
    rownames(tumor_data) <- gsub("\\.","-",rownames(tumor_data))
    
    # normalize tumor data
    norm_data <- apply(tumor_data,2,normalize_col)
    
    # Function to exclude columns that have NaNs or low variance 
    var_thresh <- 0.005
    sim_filter <- function(i) {
      if (sum(is.nan(norm_data[,i])) == 0) {
        return(var(norm_data[,i]) >= var_thresh)
      } else {
        return(FALSE)
      }
    }
    # apply filter
    filt_data <- norm_data[,unlist(lapply(1:dim(norm_data)[2],sim_filter))]
    
    # Function to generate dummy data points from tumor data
    dummy_col <- function(i) {
      return(sample(filt_data[,i],3000,replace=TRUE))
    }
    sim_data <- sapply(1:dim(filt_data)[2],dummy_col) # generation takes <15s to run
    # write data to file
    saveRDS(sim_data,file=sim_data_file)
  }
} else if (args[1] %in% 1:num_sims) { ## INDIVIDUAL RUN SECTION
  
  # Structure of object for storing results for each of the 1k sims
  #   - beta values (1 set)
  #   - var set (1 set)
  #   - point sets (30 sets)
  #   - true extra-sample error:
  #       - 4d array: 30 reps x 4 models x 3 metrics
  #   - est extra-sample error:
  #       - 4d array: 30 reps x 4 models x 3 metrics x 4 methods
  #   - final results for each model/method/metric combo:
  #       - relative standard error (1 value)
  #       - relative bias (1 value)
  
  # Constructors for objects to store results
  con_sim_est_error <- function() {
    x <- array(dim=c(num_reps,length(model_list),length(metric_list),length(method_list)))
    dimnames(x)[2:4] <- list(model_list,metric_list,method_list)
    class(x) <- "SimEstError"
    return(x)
  }
  con_sim_true_error <- function() {
    x <- array(dim=c(num_reps,length(model_list),length(metric_list)))
    dimnames(x)[2:3] <- list(model_list,metric_list)
    class(x) <- "SimTrueError"
    return(x)
  }
  con_sim_finals <- function () {
    x <- array(dim=c(length(finals_list),length(model_list),length(metric_list),length(method_list)))
    dimnames(x) <- list(finals_list,model_list,metric_list,method_list)
    class(x) <- "SimFinals"
    return(x)
  }
  
  # Simulation Results object constructor
  construct_results_obj <- function(i) {
    x <- list(index=i, betas=-1, vars=-1, points=list(),
              true_error=con_sim_true_error(),est_error=con_sim_est_error(),finals=con_sim_finals())
    class(x) <- "SimResults"
    return(x)
  }
  
  sim_data <- readRDS(file=sim_data_file)
  
  i <- args[1]
  results <- construct_results_obj(i)
  results_file <- result_filename(i)
  set.seed(i)

  if (!file.exists(results_file)) {
    # select 20 random vars
    var_set <- sample(1:dim(sim_data)[2],num_vars)
    results$vars <- var_set
    full_dataset <- data.frame(sim_data[,var_set])
    
    if (TCGA-like) {
      lam <- 2062
      alph <- 0.92
    } else {
      # TODO: random lam/alph generation
      lam <- 2062
      alph <- 0.92
    }
    
    betas <- rnorm(num_vars,mean=0,sd=1)
    results$betas <- betas
    #u <- runif(1,min=0,max=1) # <--- examine this?
    u <- 9000
    gen_weib <- function(vector) {
      return((u/(lam*exp(sum(betas*vector))))^(1/alph))
    }
    # generate outcomes and censoring variables, scaled to approx. match TCGA data
    int_outcomes <- apply(full_dataset,1,gen_weib)
    true_outcomes <- ceiling(normalize_col(apply(full_dataset,1,gen_weib))*3000)
    # survival dist is approximately the same regardless of censoring status
    censoring <- sample(c(0,1),size=dim(full_dataset)[1],replace=TRUE,prob = c(357/527,170/527))
    
    # generate distribution comparison plots between dataset and TCGA data (for testing)
    # par(mfrow=c(2,1))
    # denshist(true_outcomes,"Weib Simulated Data")
    # denshist(survival$times,"Actual Data")
    
    
    # complete data frame for generated survival data
    full_dataset["times"] <- true_outcomes
    full_dataset["status"] <- censoring
    
    for (j in 1:num_reps) {
      cat(paste0("Start of repetition ",j,"\n"))
      # select 500 random points to be our "available" data for this copy of the sim
      point_set <- sample(1:dim(full_dataset)[1],num_points)
      results$points[[length(results$points)+1]] <- point_set
      # pull out sampled data from the full set of simulated data
      dataset <- data.frame(full_dataset[point_set,])
      extrasample <- data.frame(full_dataset[-point_set,])
      
      # split data for DM experiments
      
      # split on status
      data0 <- dataset[which(dataset$status==0),]
      data1 <- dataset[which(dataset$status==1),]
      # sort the individual sets on survival time
      data0 <- data0[order(data0$times),]
      data1 <- data1[order(data1$times),]
      # recombine sorted data
      sorted_dataset <- rbind(data0,data1)
      
      # cross validation
      folds_DMCV_0 <- get_fold_indices(nrow(data0),10)
      folds_DMCV_1 <- get_fold_indices(nrow(data1),10)
      folds_DMCV <- c(folds_DMCV_0,folds_DMCV_1)
      
      # hold out
      folds_DMHO_0 <- get_fold_indices(nrow(data0),3)
      folds_DMHO_1 <- get_fold_indices(nrow(data1),3)
      folds_DMHO <- c(folds_DMHO_0,folds_DMHO_1)
      
      for(model in model_list) { # for each model type we're interested in
        # train model on the entire sample of data
        full_fit <- fit(fo, data = dataset, model = model)
        
        # get true extra-sample error on the 2.5k data points (stored/indexed on i,j)
        extra_pred <- predict(full_fit,newdata=extrasample,type="response",dist="weibull")
        
        # par(mfrow=c(2,1))
        # denshist(extra_pred,"Predicted Survival (sim)")
        # denshist(extrasample$times,"True Survival (sim)")
        
        for (metric in metric_list) {
          results$true_error[j,model,metric] <- eval(parse(text=paste0(metric,"(Surv(extrasample$times,extrasample$status), predicted = extra_pred)")))
        }
        
        
        # estimate extra-sample error with the 4 CV/HO methods on the 500 data points (stored/indexed on i,j)
        
        # ho
        ho_fit <- resample(fo,data=dataset,model=model,control=SplitControl(prop=2/3),ndpost=250)
        perf_ho <- performance(ho_fit,metrics=metric_list)
        
        for (metric in metric_list) {
          results$est_error[j,model,metric,"ho"] <- perf_ho@.Data[,metric]
        }
        
        # dmho
        testset <- sorted_dataset[folds_DMHO == 3,]
        trainset <- sorted_dataset[folds_DMHO != 3,]
        dmho_fit <- fit(fo, data = trainset, model = model) 
        dmho_pred <- predict(dmho_fit,newdata=testset,type="response",dist="weibull")
        
        # par(mfrow=c(2,1))
        # denshist(dmho_pred,"Predicted Survival (sim)")
        # denshist(testset$times,"True Survival (sim)")
        
        
        for (metric in metric_list) {
          results$est_error[j,model,metric,"dmho"] <- eval(parse(text=paste0(metric,"(Surv(testset$times,testset$status), predicted = dmho_pred)")))
        }
        
        # cv
        cv_fit <- resample(fo,data=dataset,model=model,control=CVControl(folds=cv_folds,repeats=1))
        perf_cv <- performance(cv_fit,metrics=metric_list)
        
        for (metric in metric_list) {
          results$est_error[j,model,metric,"cv"] <- mean(perf_cv@.Data[,metric])
        }
        
        
        # dmcv
        dmcv_results <- array(dim=c(length(metric_list),cv_folds))
        dimnames(dmcv_results)[1] <- list(metric_list)
        for (k in 1:cv_folds) {
          testset <- sorted_dataset[folds_DMCV == k,]
          trainset <- sorted_dataset[folds_DMCV != k,]
          dmcv_fit <- fit(fo, data = trainset, model = model) 
          dmcv_pred <- predict(dmcv_fit,newdata=testset,type="response",dist="weibull")
          
          for (metric in metric_list) {
            dmcv_results[metric,k] <- eval(parse(text=paste0(metric,"(Surv(testset$times,testset$status), predicted = dmcv_pred)")))
          }
        }
        
        for (metric in metric_list) {
          results$est_error[j,model,metric,"dmcv"] <- mean(dmcv_results[metric,])
        }
      }
    }
    # relative standard error calculation
    rse_calc <- function(true,est) {
      v <- sqrt(sum((est-true)^2)/num_reps)
      p <- sqrt(sum((true-mean(true))^2)/num_reps)
      return(v/p)
    }
    # relative bias calculation
    rb_calc <- function(true,est) {
      return(sum((est-true)/(est+true))/num_reps)
    }
    
    for (model in model_list) { # calculate and store rse, rb values for each model
      for (method in method_list) {
        for (metric in metric_list) {
          results$finals["rse",model,metric,method] <- rse_calc(results$true_error[,model,metric],
                                                                results$est_error[,model,metric,method])
          results$finals["rb",model,metric,method] <- rb_calc(results$true_error[,model,metric],
                                                              results$est_error[,model,metric,method])
        }
      }
    }
    # write full_results to file
    saveRDS(results,file=results_file)
  } else {
    cat(paste0("Tried to run simulation number ",i,", but it is already completed.\n"))
  }
} else if (args[1]=="final") { ## FINAL AGGREGATION SECTION
  if (file.exists(final_filename)) {
    stop(paste0("Final file ",final_filename," already exists."),call.=FALSE)
  }
  
  # check if all expected files exist
  check_exists <- function(i) {
    return(file.exists(result_filename(i)))
  }
  file_check <- sapply(1:num_sims,check_exists)
  
  if (sum(!file_check)>0) {
    str <- paste0("Missing ", sum(!file_check)," files: ",paste(which(file_check %in% FALSE),collapse=", "),"\n")
    stop(str,call.=FALSE)
  }
  
  full_results <- array(dim=c(length(finals_list),length(model_list),length(metric_list),length(method_list),num_sims))
  dimnames(full_results) <- list(finals_list,model_list,metric_list,method_list,1:num_sims)
  read_results <- function(i) {
    file <- result_filename(i)
    if (!file.exists(file)) {
      stop(paste0("Expected file ",file," to exist.\n"),call.=FALSE)
    }
    result <- readRDS(file)
    for (model in model_list) {
      for (metric in metric_list) {
        for (method in method_list) {
          for (final in finals_list) {
            full_results[final,model,metric,method,i] <- result[final,model,metric,method]
          }
        }
      }
    }
    return(TRUE)
  }
  
  lap <- lapply(1:num_sims,read_results)
  
  mean_results <- array(dim(c(length(finals_list),length(model_list),length(metric_list),length(method_list))))
  dimnames(mean_results) <- list(finals_list,model_list,metric_list,method_list)
  for (model in model_list) {
    for (metric in metric_list) {
      for (method in method_list) {
        mean_results["rse",model,metric,method] <- mean(full_results["rse",model,metric,method,])
        mean_results["rb",model,metric,method] <- mean(abs(full_results["rb",model,metric,method,]))
      }
    }
  }
  saveRDS(mean_results,file=final_filename)
  
} else {
  stop(usageMsg,call.=FALSE)
}



