##############################################################
# Metabolite predictability - machine learning pipeline

# The script iterates over a set of metabolites in a specified dataset,
#  and given specific machine learning modelling settings (e.g. model type).
#  Then, for each metabolite, a regression model is trained, using genera 
#  abundances as features. Cross validation is used to evaluate model's
#  performance on out-of-fold samples. 
# The script is meant to be used with a job management system such as HTCondor.

# The script requires the following arguments:
#  (1) Rdata object with all needed data / functions / variables (as created by the main R notebook, "Analysis_Notebook.Rmd")
#  (2) dataset name
#  (3) model type (out of: RF, ENet)
#  (4) runs per task
#  (5) tuning (TRUE/FALSE)
#  (6) tuning.cv.k
#  (7) test.cv.k (0 for loocv)
#  (8) split of metabolites - in order to further narrow each job (FIRST_THIRD, SECOND_THIRD, THIRD_THIRD, ALL)

# Updated: December 2020
# Author: Efrat Muller
# Email: efratmuller@mail.tau.ac.il
##############################################################

##############################################################
# Load libraries
##############################################################

# Search for libraries in below path
.libPaths( c("/usr/local/lib/R-local/R-3.6.0", .libPaths()) )

# Machine learning packages
library(ranger) # Fast implementation of Random Forests
library(compositions)
library(tidymodels)

# Files I/O packages
library(readr)

# Data visualization packages
library(RColorBrewer)
library(ggplot2)     

# Misc. packages
library(tictoc) # Track runtimes

# For parallel processing
library(doParallel)
library(foreach)
all_cores <- parallel::detectCores(logical = FALSE)
print(paste("DEBUG: detected ",all_cores," cores"))
registerDoParallel(cores = 8) #round(all_cores/4))

##############################################################
# Read arguments
##############################################################

# Arguments example: Image_For_Condor.RData FRANZOSA_IBD_HEALTHY RF 10 TRUE 10 0 ALL
args = commandArgs(trailingOnly = TRUE)
n.args <- 8
if(length(args) != n.args) {stop(paste("Expected",n.args,"command line arguments, found", length(args)))}

load(args[1])
print("DEBUG: Loaded R objects")

dataset.arg <- args[2] #"KOSTIC_INFANTS"
print(paste("DEBUG: dataset =",dataset.arg))
model.choice <- args[3] #"RF"
print(paste("DEBUG: model =",model.choice))
num.runs <- as.numeric(args[4]) 
print(paste("DEBUG: runs per task =",num.runs))
param.tuning <- as.logical(args[5]) #FALSE/TRUE
print(paste("DEBUG: parameter tuning =",param.tuning))
tuning.cv.k <- as.numeric(args[6]) 
print(paste("DEBUG: (inner) k =",tuning.cv.k))
test.cv.k <- as.numeric(args[7]) #0 means LOOCV
print(paste("DEBUG: (outer) k =",test.cv.k))
metabs.half <- args[8] # FIRST_THIRD/SECOND_THIRD/THIRD_THIRD/ALL
print(paste("DEBUG: metabs split =",metabs.half))

# Additional debugging info (includes packages versions
print("Session info:")
print(sessionInfo())
print("--------------------------------------------------------------")

##############################################################
# Preparations
##############################################################

# Place holder for performance results per metabolite
main.id <- 0
metab.pred.results <- data.frame(stringsAsFactors = F)
metab.contributors <- data.frame(stringsAsFactors = F)

# Override
genera.trans.options <- c("RELATIVE")

##############################################################
# Pipeline
##############################################################

print(paste("Started pipeline at:", format(Sys.time(), format="%H:%M:%S")))

# Iterate over datasets
for (dataset in c(dataset.arg)) { 
  cat("Working on dataset:", dataset, "\n")
  
  # Get list of metabolites to work on (metabolites with HMDB ID's)
  metabs.for.training <- get.metabs.for.training(dataset, 
                                                 metab.mappings, 
                                                 mtb.trans,
                                                 common.HMDBs)
  print(paste(metabs.for.training$n, "metabolites identified in dataset"))
  print(paste("About to train models for:",metabs.half))
  third.n <- floor(metabs.for.training$n / 3)
  
  # For longitudinal datasets, get a named vector that will serve as a map between sample names to subject names
  # If the dataset is not longitudinal, the map is from sample name to sample name. (so it can be effectively treated the same)
  sample.to.subject <- get.sample.to.subject(dataset, metadatas, mtb.trans)
  tmp.unique.subjects <- unique(unname(sample.to.subject))
  # Iterate over genera transformation
  for (genera.trans.choice in genera.trans.options) { 
    cat("Using genera transformation:", genera.trans.choice, "\n")
    
    # Prepare feature table (genera abundances)
    genera <- get.feature.table(dataset, 
                                genera.trans, 
                                genera.trans.choice)
    
    # Model chosen
    cat("Using model:",model.choice,"\n")
    
    # Shuffle sample id's?
    for (shuffle.flag in c(FALSE, TRUE)) { 
      cat("Shuffling set to",shuffle.flag,"\n")
      
      # Iterate over metabolites (from i.start to i.end)
      i.start <- ifelse(metabs.half %in% c("FIRST_THIRD","ALL"), 1, ifelse(metabs.half == "SECOND_THIRD", third.n + 1, 2 * third.n + 1))
      i.end <- ifelse(metabs.half == "FIRST_THIRD", third.n, ifelse(metabs.half == "SECOND_THIRD", 2 * third.n, metabs.for.training$n))
      
      print(paste("Running over metabolites",i.start,"to",i.end))														
      for (i in i.start:i.end) { 
        cat(".")
        if((i.end - i) %% 50 == 0) {cat(" ", (i.end - i), " metabolites left\n")}
        
        # Get metabolite name
        metabolite <- metabs.for.training$names[i]
        metabolite.hmdb <- unname(metabs.for.training$hmdbs)[i]
        
        # Extract target variable (metabolite values) 
        metab.vector <- 
          mtb.trans$LOG[[dataset]][metabolite, rownames(genera$dt)]
        # Scale (zero-mean unit-variance)
        metab.vector.scaled <- c(scale(metab.vector))
        # Shuffle if needed
        if (shuffle.flag) {metab.vector.scaled <- sample(metab.vector.scaled)}
        # Add to data frame
        dt <- cbind(genera$dt, Metabolite = metab.vector.scaled)
        
        # Run pipeline 'num.runs' times
        Run.Results <- foreach (run = 1:num.runs, 
                                .combine = 'my.combine', 
                                .multicombine = TRUE,
                                .packages = c('tidymodels', 'glmnet', 'tictoc', 'foreach')) %dopar% { 
                                  # Reset timer log
                                  tic.clearlog()
                                  # Init timer log
                                  tic("Timer")
                                  # Increment our global id (TODO: not working with dopar, ignored for now)
                                  main.id <- main.id + 1
                                  
                                  # Place holder for logging warnings
                                  tmp.warnings <- c("")
                                  
                                  if (param.tuning) {
                                    # TODO: support longitudinal datasets or remove non loo/loso validations completely
                                    test.cv.predictions <- data.frame(stringsAsFactors = F)
                                    final.args.string <- ""
                                    
                                    # Iterate over outer folds (LOO/LOSO) and collect predictions for out-of-fold data
                                    # TODO: can be shortened using purrr
                                    for (j in 1:length(tmp.unique.subjects)) {
                                      tmp.samples.to.test.on <- names(sample.to.subject)[sample.to.subject == tmp.unique.subjects[j]]
                                      dt.train <- dt[! rownames(dt) %in% tmp.samples.to.test.on,]
                                      dt.test <- dt[sample(tmp.samples.to.test.on, 1),]
                                      
                                      # Specify model (with place holders for hyperparameters. 
                                      #  They'll be tuned later)
                                      model <- my.get.model(model.choice)
                                      
                                      # Create workflow object (tidymodels convention)
                                      workflow.obj <- my.get.workflow(model)
                                      
                                      # Specify the hyper-parameters grids
                                      #  (sets of hyper-parameters to test in next step)
                                      hyparam.grid <- my.get.hyperparameter.grid(model.choice, 
                                                                                 n.features = genera$n)
                                      
                                      # Use cross-validation to find best set of hyper-parameters
                                      # TODO: catch all warnings
                                      tuning.results <- my.tune.params(dt.train, 
                                                                       workflow.obj, 
                                                                       hyparam.grid,
                                                                       tuning.cv.k)
                                      
                                      # Train final model using entire train dataset
                                      final.model <- 
                                        workflow.obj %>% 
                                        finalize_workflow(parameters = tuning.results$final.hyparams) %>% 
                                        fit(data = dt.train)
                                      
                                      # Get final hyper-parameters
                                      final.args.string <- paste0(final.args.string, ";",
                                                                  get.final.args.pretty(final.model))
                                      
                                      # Evaluate model on train data 
                                      eval.metrics <- my.evaluate(final.model, dt.train)
                                      
                                      # Get predictions on test data
                                      test.preds <- my.predict(final.model, dt.test)
                                      
                                      # Add to previous folds
                                      test.cv.predictions <- bind_rows(test.cv.predictions,
                                                                       test.preds)
                                    
                                    }
                                    
                                    # Calculate test CV performance (outer cv)
                                    test.cv.performance <- my.get.metrics(test.cv.predictions)
                                    
                                    # Record runtime
                                    toc(log = TRUE, quiet = TRUE)
                                    
                                    # Update main results table
                                    return(list(metab.pred.results = 
                                                  data.frame(ID = main.id,
                                                             PipeSettings.Genera.Trans = genera.trans.choice,
                                                             PipeSettings.ML.Model = model.choice,
                                                             PipeSettings.Tuning = TRUE,
                                                             PipeSettings.Shuffled = shuffle.flag,
                                                             Task.Dataset = dataset,
                                                             Task.Metabolite.Raw = metabolite,
                                                             Task.HMDB = metabolite.hmdb,
                                                             Final.Model.Params = final.args.string,
                                                             N.Genera.Features = genera$n,
                                                             N.Train.Samples = nrow(dt),
                                                             N.Unique.Subjects = length(tmp.unique.subjects),
                                                             Run = run,
                                                             Test.CV.RMSE = test.cv.performance$rmse,
                                                             Test.CV.RSq = test.cv.performance$rsq,
                                                             Test.CV.Spearman.rho = test.cv.performance$spearman.rho,
                                                             Test.CV.Spearman.p = test.cv.performance$spearman.p,
                                                             Train.RMSE = eval.metrics$train.rmse,
                                                             Train.RSq = eval.metrics$train.rsq,
                                                             Train.Spearman.rho = eval.metrics$train.spearman.rho,
                                                             Train.Spearman.p = eval.metrics$train.spearman.p,
                                                             Runtime = gsub(" elapsed", "", 
                                                                            gsub("Timer: ", "", 
                                                                                 tic.log(format = TRUE)[[1]])),
                                                             Warnings = paste(tmp.warnings[-1], collapse = "; "),
                                                             stringsAsFactors = F),
                                                # TODO: place holder
                                                metab.contributors = data.frame(Dummy = "")))
                                  } else { # No tuning
                                    # Specify model 
                                    model <- my.get.model(model.choice)
                                    
                                    # Create workflow object (tidymodels convention)
                                    workflow.obj <- my.get.workflow(model)
                                    
                                    # Specify a "dummy" hyper-parameters grid 
                                    #  (only includes a single default options)
                                    hyparam.grid <- my.get.hyperparameter.grid(model.choice, n.features = genera$n, dummy = TRUE)
                                    
                                    # Cross validation
                                    if (test.cv.k == 0) {
                                      loocv.preds <- foreach(j = 1:length(tmp.unique.subjects), .combine = rbind, .packages = c('tidymodels')) %do% {
                                        tmp.samples.to.test.on <- names(sample.to.subject)[sample.to.subject == tmp.unique.subjects[j]]
                                        tmp.trained.model <- workflow.obj %>% 
                                          finalize_workflow(parameters = hyparam.grid) %>% 
                                          fit(data = dt[! rownames(dt) %in% tmp.samples.to.test.on,])
                                        # To reduce biases in performance metrics we choose only 1 random sample per subject
                                        # (if the data is not longitudinal, the sampling does nothing)
                                        # Note: the returned table is collected by the foreach loop
                                        my.predict(tmp.trained.model, dt[sample(tmp.samples.to.test.on, 1),])
                                      }
                                      tuning.results <- my.get.metrics(loocv.preds)
                                      names(tuning.results) <- paste0("cv.", names(tuning.results))
                                    } else {
                                      # Here we're not actually doing any tuning, 
                                      #  just collecting performance metrics on out-of-fold samples
                                      # (for loocv this trick can't be used so we use a straigth forward loop)
                                      # TODO: support longitudinal datasets or remove non loo/loso validations completely
                                      tuning.results <- my.tune.params(dt, 
                                                                       workflow.obj, 
                                                                       hyparam.grid, 
                                                                       ifelse(test.cv.k == 0, 20, test.cv.k))
                                    }
                                    
                                    # Train final model using entire train dataset
                                    final.model <- 
                                      workflow.obj %>% 
                                      finalize_workflow(parameters = hyparam.grid) %>% 
                                      fit(data = dt)
                                    
                                    # Get final hyper-parameters
                                    final.args.string <- get.final.args.pretty(final.model)
                                    
                                    # Evaluate model on entire data (more of a sanity)
                                    eval.metrics <- my.evaluate(final.model, dt.train = dt)
                                    
                                    # Get built-in feature importance (RF only)
                                    feat.importance <- my.get.feat.importance(final.model)
                                    feat.importance$PipeSettings.Genera.Trans <- genera.trans.choice
                                    feat.importance$PipeSettings.ML.Model <- model.choice
                                    feat.importance$PipeSettings.Tuning <- FALSE
                                    feat.importance$PipeSettings.Shuffled <- shuffle.flag
                                    feat.importance$Task.Dataset <- dataset
                                    feat.importance$Task.Metabolite.Raw <- metabolite
                                    feat.importance$Task.HMDB <- metabolite.hmdb
                                    feat.importance$Run <- run
                                    
                                    # Record runtime
                                    toc(log = TRUE, quiet = TRUE)
                                    
                                    # Update main results table
                                    return(list(metab.pred.results = 
                                                  data.frame(ID = main.id,
                                                             PipeSettings.Genera.Trans = genera.trans.choice,
                                                             PipeSettings.ML.Model = model.choice,
                                                             PipeSettings.Tuning = FALSE,
                                                             PipeSettings.Shuffled = shuffle.flag,
                                                             Task.Dataset = dataset,
                                                             Task.Metabolite.Raw = metabolite,
                                                             Task.HMDB = metabolite.hmdb,
                                                             Final.Model.Params = final.args.string,
                                                             N.Genera.Features = genera$n,
                                                             N.Train.Samples = nrow(dt),
                                                             N.Unique.Subjects = length(tmp.unique.subjects),
                                                             Run = run,
                                                             Test.CV.RMSE = tuning.results$cv.rmse,
                                                             Test.CV.RSq = tuning.results$cv.rsq,
                                                             Test.CV.Spearman.rho = tuning.results$cv.spearman.rho,
                                                             Test.CV.Spearman.p = tuning.results$cv.spearman.p,
                                                             Train.RMSE = eval.metrics$train.rmse,
                                                             Train.RSq = eval.metrics$train.rsq,
                                                             Train.Spearman.rho = eval.metrics$train.spearman.rho,
                                                             Train.Spearman.p = eval.metrics$train.spearman.p,
                                                             Runtime = gsub(" elapsed", "", 
                                                                            gsub("Timer: ", "", 
                                                                                 tic.log(format = TRUE)[[1]])),
                                                             Warnings = paste(tmp.warnings[-1], collapse = "; "),
                                                             stringsAsFactors = F),
                                                metab.contributors = feat.importance))
                                    
                                  } # End of tuning / no tuning ifelse
                                } # Done iterating over runs (foreach+dopar loop)
        
        metab.pred.results <- bind_rows(metab.pred.results, Run.Results$metab.pred.results)
        metab.contributors <- bind_rows(metab.contributors, Run.Results$metab.contributors)
        
      } # Done iterating over metabs
      #} # Done iterating over models
    } # Done iterating over shuffling Y/N
  } # Done iterating over genera-transformation choices
} # Done iterating over datasets

print(paste("Completed pipeline at:", format(Sys.time(), format="%H:%M:%S")))

# Wrap up
print(paste("Predictability results table has", nrow(metab.pred.results), "rows"))

metab.pred.results$file.ID = paste(args[2:n.args], collapse = '_')

fname <- paste0('Pipeline_',
                Sys.Date(),'_',
                paste(args[2:n.args], collapse = '_'))
print(paste("Results saved with prefix:", fname))

# Write result tables
write.table(metab.pred.results, 
            file = paste0(fname,".tsv"),
            row.names = FALSE, sep = "\t", 
            quote = FALSE)
write.table(metab.contributors, 
            file = paste0(fname,"_contribs.tsv"),
            row.names = FALSE, sep = "\t", 
            quote = FALSE)
print("Wrote tables")





