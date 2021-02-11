##############################################################
# Utility functions

# Updated: January 2021
# Author: Efrat Muller
# Email: efratmuller@mail.tau.ac.il
##############################################################

##############################################################
# 1. Misc. functions
##############################################################

# Get significant marks given p values, for plotting
get.signif.marks <- function(p.vals) {
  signif.marks <- sapply(p.vals,
                         function(x) {
                           ifelse(x <= 0.001, "***", 
                                  ifelse(x <= 0.01, "**", 
                                         ifelse(x <= 0.05, "*", 
                                                ifelse(x <= 0.1, "^", ""))))
                         })
  signif.marks[is.na(signif.marks)] <- ""
  return(signif.marks)
}

# Get confidence interval and mean of values, given a numeric vector
CI <- function(x, 
               ci=.95) {
  a<-mean(x)
  s<-sd(x)
  n<-length(x)
  error<-qt(ci+(1-ci)/2,df=n-1)*s/sqrt(n)
  return(c(upper=a+error,mean=a,lower=a-error))
}

# Function for downsampling multiple samples for the same subject
# Returns a (reproducable) list of samples to include
downsample.by.subject <- function(dataset, 
                                  metadatas, 
                                  genera.trans, 
                                  max.samples.per.subject = 3) {
  set.seed(27)
  tmp <- metadatas[[dataset]] %>% 
    # We first assure we consider only samples actually in data (metadata table may include additional samples)
    filter(Sample %in% colnames(genera.trans$RELATIVE[[dataset]])) %>%
    group_by(Subject.ID) %>% 
    # We sample 3 per subject (unless there are less than 3 to begin with, in which case we take all samples of that subject)
    sample_n(if(n() < max.samples.per.subject) n() else max.samples.per.subject)
  return(as.character(tmp$Sample))
}

# Returns a map from sample names to subject names, or sample names to sample names if the data is not longitudinal (does not have a "Subject.ID" column)
get.sample.to.subject <- function(dataset, 
                                  metadatas, 
                                  mtb.trans) {
  tmp <- metadatas[[dataset]] %>%
    # Take only relevant samples
    filter(Sample %in% colnames(mtb.trans$LOG[[dataset]]))
  if ("Subject.ID" %in% names(tmp)) {
    sample.to.subject <- tmp$Subject.ID
  } else {
    sample.to.subject <- tmp$Sample
  }
  names(sample.to.subject) <- tmp$Sample
  return(sample.to.subject)
}

# Takes an arbitrary number of lists x all of which much have the same structure
my.combine <- function(x, ...) {
  mapply(rbind, x, ..., SIMPLIFY=FALSE)
}

##############################################################
# 2. Machine learning pipeline functions
##############################################################

# Get list of metabolites to work on (metabolites with HMDB ID's)
get.metabs.for.training <- function(dataset, 
                                    metab.mappings, 
                                    mtb.trans,
                                    inclusion.list.hmdbs) {
  # Create metabolite raw name to HMDB mapping
  mtb2hmdb <- metab.mappings[[dataset]]$HMDB
  names(mtb2hmdb) <- metab.mappings[[dataset]]$Compound
  
  # Get list of metabolites to work on (metabolites with HMDB ID's)
  hmdbs <- sapply(rownames(mtb.trans$LOG[[dataset]]),
                  function(m) {unname(mtb2hmdb[m])})
  # Take only those metabolites that are shared by 3 or more datasets
  hmdbs <- hmdbs[hmdbs %in% inclusion.list.hmdbs]
  metabs.names <- names(hmdbs)
  n.metabs <- length(metabs.names)
  
  return(list("hmdbs" = hmdbs, 
              "names" = metabs.names,
              "n" = n.metabs))
}

# Prepare/organize feature table (genera abundances)
get.feature.table <- function(dataset, 
                              genera.trans, 
                              genera.trans.choice) {
  # Fetch genera abundances
  genera.dt <- genera.trans[[genera.trans.choice]][[dataset]]
  genera.dt <- t(genera.dt) # Columns as microbes
  genera.dt <- as.data.frame(genera.dt)
  names(genera.dt) <- gsub(".*\\|g__", "g__", names(genera.dt))
  
  # Fix column names (save mapping) to avoid various errors in tidymodels
  new.names <- make.names(names(genera.dt))
  tmp.names.mapping <- 
    data.frame(old.names = names(genera.dt), new.names = new.names)
  names(genera.dt) <- new.names
  n.genera <- ncol(genera.dt)
  
  return(list("dt" = genera.dt,
              "n" = n.genera,
              "names.mapping" = tmp.names.mapping))
}

# Return a model object with pre-selected parameters to tune
# Currently supports only random forest or ENet
my.get.model <- function(model.type, 
                         rf.importance = "impurity_corrected") {
  if (model.type == "RF") {
    model <- rand_forest() %>%
      set_engine("ranger", importance = rf.importance) %>% 
      set_mode("regression") %>%
      set_args(mtry = tune(), trees = tune(), min_n = tune()) 
  } else if (model.type == "ENet") {
    model <- linear_reg() %>%
      set_engine("glmnet") %>% 
      set_args(penalty = tune(), mixture = tune())
  } else {
    print("Invalid model type")
    model <- NULL
  }
  return(model)
}

# Return a pre-defined grid with several hyper-parameter options.
# If 'dummy' flag is set to TRUE, this means we return a "dummy" grid that actually contains only a single option (for the sake of using the same pipeline even when we skip tuning)
my.get.hyperparameter.grid <- function(model.type, 
                                       n.features = NA, 
                                       dummy = FALSE) {
  if (model.type == "RF" & ! dummy) {
    hyparam.grid <- grid_regular(        
      mtry(range = c(5, floor(n.features/2))), 
      trees(range = c(500, 1000)),
      min_n(range = c(4, 6)),
      levels = c(4,1,2)
    ) 
  } else if (model.type == "RF" & dummy) {
    # Dummy grid - all default params
    def.mtry <- round(sqrt(n.features))
    def.min.n <- 5
    hyparam.grid <- grid_regular(        
      mtry(range = c(def.mtry, def.mtry)), 
      trees(range = c(500, 500)),
      min_n(range = c(def.min.n, def.min.n)),
      levels = c(1,1,1)
    ) 
  } else if (model.type == "ENet" & ! dummy) {
    hyparam.grid <- grid_regular(        
      mixture(range = c(0.2, 1)), 
      penalty(range = c(-8, -1)),
      levels = c(3,3)
    ) 
  } else if (model.type == "ENet" & dummy) {
    hyparam.grid <- grid_regular(        
      mixture(range = c(1, 1)), 
      penalty(range = c(-3, -3)),
      levels = c(1,1)
    ) 
  } else {
    print("Invalid model type")
    hyparam.grid <- NULL
  }
  return(hyparam.grid)
}

# Returns a simple "workflow" object (tidymodels standard)
my.get.workflow <- function(model) {
  return(workflow() %>% 
           add_formula(Metabolite ~ .) %>% 
           # add_recipe(....) %>%
           add_model(model))
}

# Tune hyper-parameters using a model specification, hyper-parameter grid, and a training dataset.
# Uses a k-fold cross-validation (if k=0, LOOCV).
# Returns the best performing set of hyper-parameters, and performance metrics for this set.
my.tune.params <- function(dt.train, 
                           workflow.obj, 
                           hyparam.grid, 
                           k) {
  # 10-fold stratified CV (splits are stratified to have 
  #   similar metabolite-level distributions)
  # Low sample numbers will generate warnings regarding to 
  #   stratification. We log them instead of printing
  # TODO: log warnings instead of only supressing them
  if (k == 0) {
    dt.train.cv <- loo_cv(dt.train)
  } else {
    suppressWarnings({
      dt.train.cv <- vfold_cv(dt.train, v = k, strata = Metabolite)
    })
  }
  
  # Define metrics of interest
  my.metrics <- metric_set(rmse) #, mae)
  
  # Actual param tuning
  hyparams.tuning <- tune_grid(
    workflow.obj,
    resamples = dt.train.cv,
    grid = hyparam.grid,
    metrics = my.metrics,
    control = control_resamples(save_pred = TRUE)
  )
  
  # Inspect tuning results:
  # hyparams.tuning %>% collect_metrics() %>% arrange(mean)
  # autoplot(hyparams.tuning, metric = "rmse")
  
  # Get best params
  final.hyparams <- hyparams.tuning %>%
    select_best(metric = "rmse")
  
  # Extract actual cross-validation predictions, for the best model only
  cv.preds.best.model <- collect_predictions(hyparams.tuning) %>% 
    inner_join(final.hyparams, by = names(final.hyparams))
  
  # Get CV performance metrics. I aggregate predictions on all out-of-fold samples (each sample appears once), and calculate the metrics on the entire dataset.
  # Another approach would be: calculate per fold, then get average performance over folds. 
  # --> Spearman
  cv.spearman <- cor.test(cv.preds.best.model$.pred, 
                          cv.preds.best.model$Metabolite, 
                          method = "spearman")
  # --> RMSE
  cv.rmse <- cv.preds.best.model %>% 
    rmse(Metabolite, .pred) %>% 
    pull(.estimate)
  
  # --> R^2
  cv.rsq <- cv.preds.best.model %>% 
    rsq(Metabolite, .pred) %>% 
    pull(.estimate)
  
  return(list(final.hyparams = final.hyparams,
              cv.spearman.rho = cv.spearman$estimate,
              cv.spearman.p = cv.spearman$p.value,
              cv.rmse = cv.rmse,
              cv.rsq = cv.rsq))
}

# Returns predictions for a given dataset
my.predict <- function(final.model, 
                       dt) {
  preds <- 
    final.model %>%
    predict(new_data = dt) %>% 
    bind_cols(dt %>% select(Metabolite)) %>%
    rename(actual = Metabolite, predicted = .pred)
  return(preds)
}

# Returns RMSE, Spearman, R^2
# Expects a table with columns: "actual" and "predicted"
my.get.metrics <- function(preds) {
  # --> Spearman
  tmp.spearman <- 
    cor.test(preds$predicted, 
             preds$actual, 
             method = "spearman") 
  # --> RMSE
  tmp.rmse <- preds %>% 
    rmse(actual, predicted) %>% 
    pull(.estimate)
  # --> RSQ
  tmp.rsq <- preds %>% 
    rsq(actual, predicted) %>% 
    pull(.estimate)
  
  return(list(spearman.rho = unname(tmp.spearman$estimate),
              spearman.p = tmp.spearman$p.value,
              rmse = tmp.rmse,
              rsq = tmp.rsq))
}

# Returns a list of performance metrics for both a train and test set
my.evaluate <- function(final.model, 
                        dt.train = NULL, 
                        dt.test = NULL) {
  
  # Get predictions for train dataset
  if (! is.null (dt.train)) {
    train.preds <- my.predict(final.model, dt.train)
    # Calculate performance metrics
    train.metrics <- my.get.metrics(train.preds)
  }
  
  # Get predictions for test dataset
  if (! is.null (dt.test)) {
    test.preds <- my.predict(final.model, dt.test)
    
    # Calculate performance metrics
    test.metrics <- my.get.metrics(test.preds)
  }
  
  return(list(train.spearman.rho = ifelse(is.null(dt.train),NA,
                                          train.metrics$spearman.rho),
              train.spearman.p = ifelse(is.null(dt.train),NA,
                                        train.metrics$spearman.p),
              train.rmse = ifelse(is.null(dt.train),NA,
                                  train.metrics$rmse),
              train.rsq = ifelse(is.null(dt.train),NA,
                                 train.metrics$rsq),
              test.spearman.rho = ifelse(is.null(dt.test),NA,
                                         test.metrics$spearman.rho),
              test.spearman.p = ifelse(is.null(dt.test),NA,
                                       test.metrics$spearman.p),
              test.rmse = ifelse(is.null(dt.test),NA,
                                 test.metrics$rmse),
              test.rsq = ifelse(is.null(dt.test),NA,
                                test.metrics$rsq)))
}

# Extract feature importance from final model
# Currently only RF feature importance is supported
my.get.feat.importance <- function(final.model) {
  # Extract model object from the final workflow
  model.obj <- pull_workflow_fit(final.model)$fit
  model.spec <- pull_workflow_spec(final.model)
  
  # Organize feature importance in a table
  # For RF models:
  if (model.spec$engine == "ranger") {
    feat.importance <- data.frame(model.obj$variable.importance,
                                  stringsAsFactors = F)
    names(feat.importance) <- "Score"
    feat.importance$Feature <- rownames(feat.importance)
    rownames(feat.importance) <- NULL
    
    # Sort
    feat.importance <- feat.importance %>% arrange(-Score)
  } else if (model.spec$engine == "glmnet") {
    # TODO 
    feat.importance = NULL
  }
  
  return(feat.importance)
}

# Get a simple string with all final model hyper-paremeters
get.final.args.pretty <- function(final.model) {
  final.args <- pull_workflow_spec(final.model)$args
  final.args.string <- paste(names(final.args), 
                             lapply(final.args, 
                                    function(a) {rlang::get_expr(a)}), 
                             sep=" = ",collapse="; " )
  return(final.args.string)
}





