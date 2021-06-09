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


##############################################################
# 3. Model comparison functions
##############################################################

# Given a full contributions table, extract only contributions relevant to a specific hmdb id 
#  and add cross-study statistics such as mean importance score.
get.hmdb.contribs <- function(contribs, hmdb, healthy.only = T) {
  tmp <- contribs %>% 
    filter(Task.HMDB == hmdb) 
  if (healthy.only) {
    tmp <- tmp %>% 
      filter(Healthy) 
  }
  tmp <- tmp %>%
    # Add mean score per feature 
    group_by(Feature) %>%
    mutate(Mean.Feat.P = mean(Mean.P), N = n(), N.Signif = sum(Mean.P < 0.1)) %>%
    ungroup() %>%
    arrange(-N.Signif, Mean.Feat.P)
  
  # Order features by number of times that feature is significant / mean P
  tmp$Feature <- factor(tmp$Feature, levels = unique(tmp$Feature))
  
  # Get additional versions of the table
  tmp2 <- tmp %>%
    tidyr::pivot_wider(id_cols = Feature, names_from = Task.Dataset, values_from = Mean.P) 
  
  tmp3 <- as.matrix(tmp2[,-1])
  rownames(tmp3) <- tmp2$Feature
  
  return(list(contrib.hmdb = tmp,
              contrib.hmdb.wide = tmp2,
              contrib.hmdb.m = tmp3))
}

# Returns a table with comparison of feature importances between each pair of datasets
get.feat.importance.comparisons <- function(contrib.hmdb.wide, datasets.ordering) {
  fi.comparisons <- data.frame(stringsAsFactors = F)
  
  # Iterate over pairs of datasets
  n.datasets <- length(datasets.ordering)
  for (i in 1:(n.datasets-1)) {
    for (j in (i+1):n.datasets) {
      dataset1 <- datasets.ordering[i]
      dataset2 <- datasets.ordering[j]
      tmp <- contrib.hmdb.wide[,c("Feature", dataset1, dataset2)]
      
      # Add significance flag for Fisher's test
      tmp <- tmp %>%
        rename(P.D1 = 2, P.D2 = 3) %>%
        filter((! is.na(P.D1)) & (! is.na(P.D2))) %>%
        mutate(Signif.0.05.D1 = (P.D1 <= 0.05)) %>%
        mutate(Signif.0.05.D2 = (P.D2 <= 0.05)) %>%
        mutate(Signif.0.1.D1 = (P.D1 <= 0.1)) %>%
        mutate(Signif.0.1.D2 = (P.D2 <= 0.1))
      
      # Calculate several comparison metrics
      fisher.0.1 <- fisher.test(factor(tmp$Signif.0.1.D1, levels=c("TRUE","FALSE")), factor(tmp$Signif.0.1.D2, levels=c("TRUE","FALSE")))
      fisher.0.05 <- fisher.test(factor(tmp$Signif.0.05.D1, levels=c("TRUE","FALSE")), factor(tmp$Signif.0.05.D2, levels=c("TRUE","FALSE")))
      cor.sp <- cor.test(x = tmp$P.D1, y = tmp$P.D2, method = 'spearman')
      cor.pe <- cor.test(x = tmp$P.D1, y = tmp$P.D2, method = 'pearson')
      num.feat.sh.0.1 <- sum(tmp$Signif.0.1.D1 & tmp$Signif.0.1.D2)
      num.feat.sh.0.05 <- sum(tmp$Signif.0.05.D1 & tmp$Signif.0.05.D2)
      
      # Weighted correlation - here we give higher weights to more important features 
      wcor <- wtd.cor(x = tmp$P.D1, y = tmp$P.D2, 
                      weight = pmax(-log(tmp$P.D1), -log(tmp$P.D2))) 
      
      fi.comparisons <- bind_rows(fi.comparisons,
                                  data.frame(Dataset1 = dataset1,
                                             Dataset2 = dataset2,
                                             N.Feat.In.Common = nrow(tmp),
                                             Cor.Spearman = unname(cor.sp$estimate),
                                             Cor.Spearman.P = cor.sp$p.value,
                                             Cor.Pearson = unname(cor.pe$estimate),
                                             Cor.Pearson.P = cor.pe$p.value,
                                             W.Cor.Pearson = wcor[1,'correlation'],
                                             W.Cor.Pearson.P = wcor[1,'p.value'],
                                             Fisher.0.1.OR = unname(fisher.0.1$estimate),
                                             Fisher.0.1.P = fisher.0.1$p.value,
                                             Fisher.0.05.OR = unname(fisher.0.05$estimate),
                                             Fisher.0.05.P = fisher.0.05$p.value,
                                             Shared.Feats.0.1 = num.feat.sh.0.1,
                                             Shared.Feats.0.05 = num.feat.sh.0.05,
                                             stringsAsFactors = F))
    }
  }
  
  # Convert datasets to factors according to given order
  fi.comparisons <- fi.comparisons %>%
    mutate(Dataset1 = factor(Dataset1, levels = datasets.ordering),
           Dataset2 = factor(Dataset2, levels = datasets.ordering))
  
  return(fi.comparisons)
}


# Returns a table with comparison of feature importances between each pair of datasets,
#  using adjusted models, where models are re-trained with only features shared between the pair of datasets compared.
# TODO: merge with function above
get.adj.feat.importance.comparisons <- function(pairwise.comp.contribs, datasets.ordering, hmdb) {
  fi.comparisons <- data.frame(stringsAsFactors = F)
  
  # Iterate over dataset pairs
  n.datasets <- length(datasets.ordering)
  for (i in 1:(n.datasets-1)) {
    for (j in (i+1):n.datasets) {
      dataset1 <- datasets.ordering[i]
      dataset2 <- datasets.ordering[j]
      
      # First we get the list of feature scores for each pair of datasets, 
      #  using the "adjusted" models (where only shared features were considered)
      tmp1 <- pairwise.comp.contribs %>%
        ungroup() %>%
        filter(HMDB == hmdb) %>%
        filter(Dataset1 == dataset1) %>%
        filter(Dataset.Aligned.To == dataset2) %>%
        rename(P.D1 = Mean.P) %>%
        select(Feature, P.D1) 
      tmp2 <- pairwise.comp.contribs %>%
        ungroup() %>%
        filter(HMDB == hmdb) %>%
        filter(Dataset1 == dataset2) %>%
        filter(Dataset.Aligned.To == dataset1) %>%
        rename(P.D2 = Mean.P) %>%
        select(Feature, P.D2) 
      # Merge together
      tmp <- tmp1 %>%
        left_join(tmp2, by = "Feature") %>%
        mutate(Dataset1 = dataset1, Dataset2 = dataset2) %>%
        # For Fisher's test
        mutate(Signif.0.05.D1 = (P.D1 <= 0.05)) %>%
        mutate(Signif.0.05.D2 = (P.D2 <= 0.05)) %>%
        mutate(Signif.0.1.D1 = (P.D1 <= 0.1)) %>%
        mutate(Signif.0.1.D2 = (P.D2 <= 0.1))
      
      # Calculate several comparison metrics
      fisher.0.1 <- fisher.test(factor(tmp$Signif.0.1.D1, levels=c("TRUE","FALSE")), factor(tmp$Signif.0.1.D2, levels=c("TRUE","FALSE")))
      fisher.0.05 <- fisher.test(factor(tmp$Signif.0.05.D1, levels=c("TRUE","FALSE")), factor(tmp$Signif.0.05.D2, levels=c("TRUE","FALSE")))
      cor.sp <- cor.test(x = tmp$P.D1, y = tmp$P.D2, method = 'spearman')
      cor.pe <- cor.test(x = tmp$P.D1, y = tmp$P.D2, method = 'pearson')
      num.feat.sh.0.1 <- sum(tmp$Signif.0.1.D1 & tmp$Signif.0.1.D2)
      num.feat.sh.0.05 <- sum(tmp$Signif.0.05.D1 & tmp$Signif.0.05.D2)
      
      # Weighted correlation - here we give higher weights to more important features 
      # Weighting option 1: max raw importance (the function normalizes to mean 1)
      wcor <- wtd.cor(x = tmp$P.D1, y = tmp$P.D2, 
                      weight = pmax(-log(tmp$P.D1), -log(tmp$P.D2))) 
      
      # Organize all stats in one table
      fi.comparisons <- bind_rows(fi.comparisons,
                                  data.frame(Dataset1 = dataset1,
                                             Dataset2 = dataset2,
                                             N.Feat.In.Common = nrow(tmp),
                                             Cor.Spearman = unname(cor.sp$estimate),
                                             Cor.Spearman.P = cor.sp$p.value,
                                             Cor.Pearson = unname(cor.pe$estimate),
                                             Cor.Pearson.P = cor.pe$p.value,
                                             W.Cor.Pearson = wcor[1,'correlation'],
                                             W.Cor.Pearson.P = wcor[1,'p.value'],
                                             Fisher.0.1.OR = unname(fisher.0.1$estimate),
                                             Fisher.0.1.P = fisher.0.1$p.value,
                                             Fisher.0.05.OR = unname(fisher.0.05$estimate),
                                             Fisher.0.05.P = fisher.0.05$p.value,
                                             Shared.Feats.0.1 = num.feat.sh.0.1,
                                             Shared.Feats.0.05 = num.feat.sh.0.05,
                                             stringsAsFactors = F))
    }
  }
  
  # Convert datasets to factors according to given order
  fi.comparisons <- fi.comparisons %>%
    mutate(Dataset1 = factor(Dataset1, levels = datasets.ordering),
           Dataset2 = factor(Dataset2, levels = datasets.ordering))
  
  return(fi.comparisons)
}

# Receives an hmdb id and a contributions table (scores + P's for all features in all analyzed models)
# Returns a list of features
get.features.to.include.in.plot <- function(contribs, hmdb, healthy.only = T) {
  features.to.include <- contribs %>%
    filter(Task.HMDB == hmdb) %>%
    filter(Mean.P < 0.1) 
  if (healthy.only) {
    features.to.include <- features.to.include %>% 
      filter(Healthy) 
  }
  return(unique(features.to.include$Feature))
}

plot.ggplot.heatmap.feat.imp2 <- function(contrib.hmdb, features.to.include) {
  fixed.feat.levels <- levels(contrib.hmdb$Feature)
  names(fixed.feat.levels) <- levels(contrib.hmdb$Feature)
  fixed.feat.levels["g__Ruminococcaceae.UCG.013"] <- "g__Ruminococcaceae\nUCG.013"
  fixed.feat.levels["g__Lachnospiraceae.ND3007.group"] <- "g__Lachnospiraceae\nND3007.group"
  fixed.feat.levels["g__Erysipelotrichaceae.UCG.003"] <- "g__Erysipelotrichaceae\nUCG.003"
  fixed.feat.levels["g__Ruminococcaceae.UCG.003"] <- "g__Ruminococcaceae\nUCG.003"
  fixed.feat.levels["g__Phascolarctobacterium"] <- "g__Phascolarcto-\nbacterium"
  fixed.feat.levels["g__Lachnospiraceae.NC2004.group"] <- "g__Lachnospiraceae\nNC2004.group"
  fixed.feat.levels["g__Erysipelatoclostridium"] <- "g__Erysipelato-\nclostridium"
  fixed.feat.levels["g__Ruminococcaceae.NK4A214.group"] <- "g__Ruminococcaceae\nNK4A214.group"
  fixed.feat.levels["g__Ruminococcaceae.UCG.010"] <- "g__Ruminococcaceae\nUCG.010"
  fixed.feat.levels["g__Lachnospiraceae.NK4A136.group"] <- "g__Lachnospiraceae\nNK4A136.group"
  fixed.feat.levels["g__Lachnospiraceae.UCG.001"] <- "g__Lachnospiraceae\nUCG.001"
  fixed.feat.levels["g__Family.XIII.AD3011.group"] <- "g__Family.XIII\nAD3011.group"
  fixed.feat.levels["g__Ruminococcaceae.UCG.002"] <- "g__Ruminococcaceae\nUCG.002"
  fixed.feat.levels["g__Christensenellaceae.R.7.group"] <- "g__Christensenellaceae\nR.7.group"
  fixed.feat.levels["g__Ruminococcaceae.UCG.014"] <- "g__Ruminococcaceae\nUCG.014"
  fixed.feat.levels["g__Ruminococcaceae.UCG.005"] <- "g__Ruminococcaceae\nUCG.005"
  fixed.feat.levels["g__Escherichia-Shigella"] <- "g__Escherichia-\nShigella"
  fixed.feat.levels["g__Lachnoclostridium"] <- "g__Lachno-\nclostridium"
  fixed.feat.levels <- gsub("g__", "", fixed.feat.levels)
  
  contrib.hmdb2 <- contrib.hmdb %>%
    filter(Feature %in% features.to.include) %>%
    mutate(Feature = factor(Feature, levels = names(fixed.feat.levels),
                            labels = unname(fixed.feat.levels))) %>%
    mutate(Mean.P = ifelse(Mean.P > 0.1, NA, Mean.P)) %>%
    mutate(Mean.P.Signif.Mark = get.signif.marks(Mean.P)) %>%
    mutate(Mean.P.Signif.Mark = ifelse(Mean.P.Signif.Mark=="^","",Mean.P.Signif.Mark))
  
  p <- ggplot(contrib.hmdb2, aes(x = Feature,y = Task.Dataset, color = "")) + 
    geom_tile(aes(fill = Mean.P)) +
    scale_fill_gradientn(colors = c("#5f0f40","#C87EAC"),
                         values = c(0,1),
                         breaks = c(0.08, 0.05, 0.01, 0.005),
                         na.value = "gray88") +
    guides(fill = guide_colourbar(frame.colour = "black", 
                                  title = "Feature P Value", 
                                  #barheight = 4,
                                  reverse = TRUE)) +
    ylab(NULL) + 
    xlab(NULL) +
    theme_bw() +
    # Dummy scale just to add to legend
    scale_color_manual(values = NA) +              
    guides(color = guide_legend("P > 0.1", 
                                title.vjust = 0.5,
                                override.aes=list(fill="gray88", color="black"))) +
    geom_text(aes(label = Mean.P.Signif.Mark), size = 4, color = "black") +
    theme(axis.text.x = element_text(angle = 90, size = 7, vjust = 0.2, hjust = 1),
          axis.text.y = element_text(size = 8),
          legend.position = "bottom",
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 9, vjust = 0.8),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(5.5, 0, 5.5, 5.5), "points"))
  
  return(p)
}

