rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

main_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/scripts/")
data_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/data/preprocessed_data/")
out_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/output/")

setwd(main_dir)

require(mixOmics)
require(dplyr)

source("fun_utils.R")

seed <- 1
set.seed(seed)

# For design matrix
bool_borrow_info <- TRUE
if (bool_borrow_info) {
  bool_empirical <- TRUE
} else {
  bool_empirical <- FALSE 
}

tune_keep <- TRUE # CV for number of variables per block 
n_cpus <- 4

bool_training_test <- TRUE

data_names <- c("MS",
                "cell_types",
                "glyco-lipo-proteins",
                "log_ratios",
                "covariates")

vec_col_diablo <- c("mediumaquamarine", 
                    "lightskyblue",
                    "grey75", 
                    "wheat4",
                    "#CC79A7")

names(vec_col_diablo) <- data_names

if ("covariates" %in% data_names) {
  
  df_clin <- df_info
  df_clin$gender <- as.numeric(df_clin$gender) - 1 # 1 = M, 0 = F
  
  mat_indic <- unmap(df_clin$ethnicity) 
  df_clin$asian <- mat_indic[,1]
  df_clin$other <- mat_indic[,2]
  df_clin$white <- mat_indic[,3]
  
  covariates <- c("age", "gender", "bmi", "asian", "other", "white")[c(1:2)]
  df_clin <- df_clin[, c("sample_id_v2", "subject_id", covariates), drop = FALSE]
  covariates <- first_letter_upper_case(covariates)
  colnames(df_clin) <- c("sample_id_v2", "subject_id", covariates)
  nb_fixed_covariates <- length(covariates)
  
}


list_df <- list(df_ms,
                df_ct,
                df_glprot,
                df_all_ratios,
                df_clin)


bool_save <- FALSE
if (bool_save) {
  res_dir <- paste0(out_dir, "pred_gCCA_days_thres_", days_thres, "_", dep_var, 
                    "_", paste0(data_names, collapse = "-"), "_",
                    paste0(covariates, collapse = "-"),
                    ifelse(bool_training_test, "_training_test", "_training"),
                    ifelse(bool_borrow_info, 
                           paste0("_borrow_info", 
                                  ifelse(bool_empirical, "_emp", "_fixed"))
                           , "_no_borrow_info"),
                    ifelse(tune_keep, "_CV", "_not_sparse"),
                    "_seed_", seed, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F,
       split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"),
       type = "message")
} else {
  res_dir <- NULL
}


list_X <- list_df[data_ids]
names(list_X) <- data_names

# Restrict to samples available for each data type
#
inter_samples <- Reduce(intersect, lapply(list_X, rownames)) 
list_X <- lapply(list_X, function(x) { x[inter_samples,] })
stopifnot(length(unique(lapply(list_X, rownames)))==1) # check that all datasets have the same rownames

print(paste0("Number of samples before taking early timepoints: ", length(inter_samples)))
print(paste0("Datasets: ", paste0(data_names, collapse = ", ")))
print(paste0("Number of analysed variables in each dataset: ", paste0(sapply(list_X, ncol), collapse = ", ")))

sub_df_info <- df_info[match(inter_samples, df_info$sample_id_v2),, drop = FALSE]
stopifnot(isTRUE(all.equal(rownames(sub_df_info), rownames(list_X[[1]]))))


# Choose severity groups to be analysed
#
selected_severity_groups <- c("B", "C", "D", "E") 
sub_df_info <- choose_severity_groups(sub_df_info, selected_severity_groups)


# Prepare response (factor if discrete)
#
trunc_days <- 50
load(file.path(out_dir, "all_scores.RData")) 
sub_df_info <- left_join(sub_df_info, df_scores) 
sub_df_info <- sub_df_info[!is.na(sub_df_info[, dep_var]),, drop = FALSE]

if (discrete_resp) {
  sub_df_info[, dep_var] <- factor(sub_df_info[, dep_var])
}


# Restrict to the first sample of each subject within a trunc_days window post symptom onset
#
enforce_subj_w_samples_avail_after_trunc_days <- FALSE # we only use the subjects who have at least one sample taken after trunc_days days # for persisting vs recovering CRP, ensures that the last CRP samples is taken after trunc_days days, i.e., not too early
list_fs <- get_early_or_late_samples(sub_df_info, 
                                     vec_var = NULL, 
                                     bool_early = TRUE,
                                     trunc_days = days_thres, 
                                     single_sample_per_subject = TRUE, 
                                     enforce_subj_w_samples_avail_after_trunc_days = enforce_subj_w_samples_avail_after_trunc_days)

sub_df_info <- list_fs$df_comb
sub_df_info <- sub_df_info[order(sub_df_info[, dep_var]), ]

stopifnot(!all(duplicated(sub_df_info$subject_id)))

# Restrict the variables to the first sample of each subject within the desired time window
#
list_X <- lapply(list_X, function(x) { x[match(sub_df_info$sample_id_v2, rownames(x)),] })


if (bool_training_test) {
  set.seed(seed)
  train_ids <- caret::createDataPartition(sub_df_info[, dep_var], p = 0.7, list = FALSE)
  
  list_X_to_fit <- lapply(list_X, function(ll) ll[train_ids,, drop = FALSE])
  sub_df_info_to_fit <- sub_df_info[train_ids,, drop = FALSE]
  
  list_X_test <- lapply(list_X, function(ll) ll[-train_ids,, drop = FALSE])
  sub_df_info_test <- sub_df_info[-train_ids,, drop = FALSE]
  
  print(paste0("Number of subjects in training set: ", nrow(sub_df_info_to_fit)))
  print(paste0("Number of subjects in test set: ", nrow(sub_df_info_test)))
  
} else {
  list_X_to_fit <- list_X
  sub_df_info_to_fit <- sub_df_info
  
  list_X_test <- sub_df_info <- NULL
}


# Sets the design matrix
#
if (bool_borrow_info) {
  
  if (bool_empirical) {
      design <- sapply(list_X_to_fit, function(X1) {
      sapply(list_X_to_fit, function(X2) {
        print(class(X1))
        print(class(X2))
        res_pls <- pls(X1, X2, ncomp = 1)
        cor(res_pls$variates$X, res_pls$variates$Y)
        # ifelse(cor(res_pls$variates$X, res_pls$variates$Y)>thres, 1, 0)
      })
    })
    
    
  } else {

    design <- matrix(0.1, ncol = length(list_X_to_fit), nrow = length(list_X_to_fit),
                     dimnames = list(names(list_X_to_fit), names(list_X_to_fit)))
  }
  
} else {
  design <- matrix(0, ncol = length(list_X_to_fit), nrow = length(list_X_to_fit),
                   dimnames = list(names(list_X_to_fit), names(list_X_to_fit)))
}
diag(design) <-  0

print("Design matrix")
print(design)

allSame <- function(x) length(unique(x)) == 1
stopifnot(allSame(lapply(list_X_to_fit, function(ll) rownames(ll))))
stopifnot(all.equal(rownames(list_X_to_fit[[1]]), sub_df_info_to_fit$sample_id_v2))


if (discrete_resp & tune_keep) {# it seems that tune.block.spls doesn't exist for continuous responses.... so cv to select keepX not available for continuous responses
  
  set.seed(seed)
  
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_cpus-1)
  
  list_keep_all <- lapply(list_X_to_fit, function(X) {vec <- seq(min(ncol(X), 4), min(ncol(X), 24), by = 2)})
  
  
  tune.splsda.srbct <- tune.block.splsda(list_X_to_fit, 
                                         Y= sub_df_info_to_fit[, dep_var], 
                                         ncomp = 2, # max number of component we suggest to push ncomp a bit more, e.g. 4
                                         validation = 'Mfold',
                                         design = design,
                                         folds = 3, 
                                         dist = 'max.dist', 
                                         progressBar = FALSE,
                                         measure = "BER", 
                                         test.keepX = list_keep_all,
                                         nrepeat = 100,  # increase to 100 ?
                                         BPPARAM = BPPARAM) 
  
  # plot(tune.splsda.srbct)
  
  error <- tune.splsda.srbct$error.rate 
  # print(paste0("CV error rate: "))
  # print(error)
  
  ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
  print(paste0("CV optimal number of components: ", ncomp))
  ncomp <- ifelse(ncomp > 1, ncomp, 2) # if optimal ncomp is 1, increase to 2 (and then we can always restrict to the first - doesn't change the results)
  print(paste0("ncomp used for fitting: ", ncomp))
  
  select.keepX <- lapply(tune.splsda.srbct$choice.keepX, function(keep) keep[1:ncomp])  # optimal number of variables to select
  print("CV optimal number of variables per component and per data type: ")
  print(unlist(select.keepX))
  
  # Fit final model to training set, using the CV-optimised parameters:
  #
  set.seed(seed)
  top_res <- block.splsda(list_X_to_fit, 
                          Y = sub_df_info_to_fit[, dep_var],
                          design = design,
                          keepX = select.keepX,
                          ncomp = ncomp)
  
  
} else {
  
  stopifnot(!tune_keep)
  
  ncomp <- 3
  keepX <- lapply(list_X_to_fit, function(ll) rep(ncol(ll), ncomp))
  names(keepX) <- names(list_X_to_fit)
  
  set.seed(seed)
  if (discrete_resp) {
    top_res <- block.splsda(list_X_to_fit, 
                            Y = sub_df_info_to_fit[, dep_var],
                            design = design,
                            ncomp = ncomp,
                            keepX = keepX)
  } else {
    resp_mat <- as.matrix(sub_df_info_to_fit[, dep_var, drop = FALSE])
    rownames(resp_mat) <- sub_df_info_to_fit$sample_id_v2
    top_res <- block.spls(list_X_to_fit, 
                          Y =  resp_mat,
                          design = design,
                          ncomp = ncomp,
                          keepX = keepX)
  }
  
  
  
}

if (discrete_resp) {
  
  if (ncomp > 1) {
    if (bool_save) {
      pdf(paste0(res_dir, "plot_indiv.pdf"), width = 12, height = 10, paper='special')
    }
    plotIndiv(top_res, legend = TRUE, ind.names = FALSE, ellipse = TRUE, star = F)
    if (bool_save) {
      dev.off()
    }
  }
  
  for (comp in 1:ncomp) {
    if (bool_save) {
      pdf(paste0(res_dir, "loadings_comp_", comp, ".pdf"), width = 15, height = 12, paper='special')
    }
    plotLoadings(top_res,  
                 comp = comp, 
                 legend.col = vec_col_gr_2,
                 contrib = 'max', 
                 size.name = 1)
    if (bool_save) {
      dev.off()
    }
    if (bool_save) {
      pdf(paste0(res_dir, "loadings_comp_", comp, "_bw.pdf"), width = 15, 
          height = 12, paper='special')
    }
    plotLoadings(top_res,  
                 comp = comp, 
                 size.name = 1)
    if (bool_save) {
      dev.off()
    }
    if (bool_save) {
      pdf(paste0(res_dir, "loadings_comp_", comp, "_horiz.pdf"), width = 15, 
          height = 3, paper='special')
    }
    par(mfrow = c(1, length(data_ids)))
    for (bl in seq_along(data_ids)) {
      plotLoadings(top_res,
                   block = bl,
                   comp = comp,
                   legend.col = vec_col_gr_2,
                   contrib = 'max',
                   size.name = 0.5,
                   legend = FALSE)
    }
    if (bool_save) {
      dev.off()
    }
    par(mfrow = c(1, 1), mar = c(5, 5, 5, 5))
    for (bl in seq_along(data_ids)) {
      if (bool_save) {
        pdf(paste0(res_dir, "loadings_comp_", comp, "_bl_", bl, ".pdf"), 
            width = 6, height = 4, paper='special')
      }
      plotLoadings(top_res,
                   block = bl,
                   comp = comp,
                   legend.col = vec_col_gr_2,
                   contrib = 'max',
                   size.name = 1,
                   legend = FALSE)
      if (bool_save) {
        dev.off()
      }
    }
    if (bool_save) {
      pdf(paste0(res_dir, "loadings_comp_", comp, "_horiz_bw.pdf"), 
          width = 22, height = 6, paper='special')
    }
    par(mfrow = c(1, length(data_ids)))
    for (bl in seq_along(data_ids)) {
      
      plotLoadings(top_res,  
                   block = bl,
                   comp = comp, 
                   size.name = 1)
      
    }
    if (bool_save) {
      dev.off()
    }
    
  }
  
  
  vec_col <- seq_along(data_names)
  
  methods <- c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid")
  
  for (comp in 1:2) {
    pal <- colorRampPalette(c("blue", "white", "red"))(1000)
    for (mm in methods) {
      if (bool_save) {
        pdf(paste0(res_dir, "heatmap_", mm, "_comp_", comp, ".pdf"), 
            width = 18, height = 15, paper='special')
      }
      cimDiablo(top_res, 
                color.blocks = vec_col_diablo,
                color.Y = vec_col_gr_2, color = pal, 
                comp = comp, margin=c(15,30), legend.position = "right",
                clust.method = c(mm, mm))
      if (bool_save) {
        dev.off()
      }
      if (bool_save) {
        pdf(paste0(res_dir, "heatmap_", mm, "_no_row_order_comp_", comp, ".pdf"), 
            width = 18, height = 15, paper='special')
      }
      cimDiablo(top_res, 
                color.blocks = vec_col_diablo,
                color.Y = vec_col_gr_2, color = pal, 
                cluster = "column",
                comp = comp, margin=c(15,30), legend.position = "right",
                clust.method = c(mm, mm))
      if (bool_save) {
        dev.off()
      }
    }
  }
  
  par(mfrow = c(1, 1))
 
  for (comp in 1:ncomp) {
    
    if (bool_save) {  
      pdf(paste0(res_dir, "ROC_training_all_data_types_comp_", comp, ".pdf"), 
          width = 6.45, height = 4.25, paper='special')
    }
    my_auroc(top_res, roc.comp = comp, bool_multiple = TRUE, 
             vec_col = vec_col_diablo)  
    
    if (bool_save) {
      dev.off()
      pdf(paste0(res_dir, "ROC_training_all_data_types_comp_", comp, "_plus_avg.pdf"), 
          width = 6.45, height = 4.25, paper='special')
    }
    my_auroc(top_res, roc.comp = comp, bool_multiple = TRUE, 
             bool_add_avg = TRUE, vec_col = vec_col_diablo) 
    
    if (bool_save) {
      dev.off()
      pdf(paste0(res_dir, "ROC_training_all_data_types_comp_", comp, "_plus_weight.pdf"), 
          width = 6.45, height = 4.25, paper='special')
    }
    my_auroc(top_res, roc.comp = comp, bool_multiple = TRUE, 
             bool_add_weight = TRUE, vec_col = vec_col_diablo)
    
    if (bool_save) {
      dev.off()
      pdf(paste0(res_dir, "ROC_training_weighted_y_comp_", comp, ".pdf"), 
          width = 6.45, height = 4.25, paper='special')
    }
    my_auroc(top_res, 
             roc.block = dn, roc.comp = comp, my_predict = "Weighted")
    
    if (bool_save) {
      dev.off()
      pdf(paste0(res_dir, "ROC_training_avg_y_comp_", comp, ".pdf"), 
          width = 6.45, height = 4.25, paper='special')
    }
    my_auroc(top_res, 
             roc.block = dn, roc.comp = comp, my_predict = "Averaged")
    if (bool_save) {
      dev.off()
    }
  }
  print(top_res$prop_expl_var)
  
} else {
  
  print(top_res$prop_expl_var)
  rowSums(do.call(rbind, top_res$prop_expl_var))
  
  for (comp in 1:ncomp) {
    if (bool_save) {
      pdf(paste0(res_dir, "loadings_comp_", comp, ".pdf"), width = 15, height = 12, paper='special')
    }
    plotLoadings(top_res,  comp = comp)
    if (bool_save) {
      dev.off()
    }
  }
  
}


top_res$prop_expl_var

if(bool_save) {
  save.image(file = file.path(res_dir, "output.RData"))
}


dim(sub_df_info_to_fit)

if (bool_training_test) {
  
  # For non discriminant analysis, the predicted values (predict) are returned 
  # for each block and these values are combined by average (AveragedPredict) 
  # or weighted average (WeightedPredict), using the weights of the blocks that 
  # are calculated as the correlation between a block's components and the outcome's components.
  # For discriminant analysis, the predicted class is returned for each block (class) 
  # and each distance (dist) and these predictions are combined by majority vote 
  # (MajorityVote) or weighted majority vote (WeightedVote), using the weights of 
  # the blocks that are calculated as the correlation between a block's components 
  # and the outcome's components. NA means that there is no consensus among the 
  # block. For PLS-DA and sPLS-DA objects, the prediction area can be visualised 
  # in plotIndiv via the background.predict function.
  #
  my_pred <- predict(top_res, newdata = list_X_test)
  
  # pls by default uses 2 components, so the first prediction is only based on 
  # the first component, the second prediction is based on the first two 
  # components (this is how I understand it).
  #
  par(mfrow = c(1,1), pty="s")
  for (ii in seq_along(list_X)) {
    if (bool_save) {
      pdf(paste0(res_dir, "pred_block_", ii, ".pdf"), width = 5.7, height = 6, paper='special')
    }
    train_subj <- rownames(top_res$ind.mat)
    group <- as.factor(df_scores$CRP_hclust_2gr_EF1_EF2[match(train_subj, df_scores$sample_id_v2)])
    test_subj <- rownames(my_pred$variates[[ii]])
    vec_col_test <- vec_col_gr_2[df_scores$CRP_hclust_2gr_EF1_EF2[match(test_subj, df_scores$sample_id_v2)]]
    plotIndiv(top_res, comp = 1:2, blocks = ii, 
              rep.space = "X-variate",
              style="graphics",
              ind.names=FALSE, 
              group = group,
              col.per.group = vec_col_gr_2,
              pch = 19)
    points(my_pred$variates[[ii]][, 1], 
           my_pred$variates[[ii]][, 2], 
           pch = 13, cex = 1.2, col=vec_col_test, )
    if (bool_save) {
      dev.off()
    }
  }
  
  if (bool_save) {
    pdf(paste0(res_dir, "pred_average.pdf"), width = 5.7, height = 6, paper='special')
  }
  train_subj <- rownames(top_res$ind.mat)
  group <- as.factor(df_scores$CRP_hclust_2gr_EF1_EF2[match(train_subj, df_scores$sample_id_v2)])
  test_subj <- rownames(my_pred$variates[[ii]])
  vec_col_test <- vec_col_gr_2[df_scores$CRP_hclust_2gr_EF1_EF2[match(test_subj, df_scores$sample_id_v2)]]
  aa <- plotIndiv(top_res, comp = 1:2, blocks = "average",
                  style="graphics",
                  ind.names=FALSE, 
                  group = group,
                  col.per.group = vec_col_gr_2,
                  pch = 19)
  points(Reduce("+", lapply(my_pred$variates, function(ll) ll[, 1])) / length(my_pred$variates), 
         Reduce("+", lapply(my_pred$variates, function(ll) ll[, 2])) / length(my_pred$variates),
         pch = 13, cex = 1.2, col=vec_col_test)
  if (bool_save) {
    dev.off()
  }
  
  if (discrete_resp) {
  
    
    # print(vec_col)
    for (comp in 1:ncomp) {
      if (bool_save) {  
        pdf(paste0(res_dir, "ROC_test_all_data_types_comp_", comp, ".pdf"), 
            width = 6.9, height = 4, paper='special')
      }
      my_auroc(top_res, newdata = list_X_test, outcome.test = sub_df_info_test[, dep_var],
               roc.comp = comp, bool_multiple = TRUE, vec_col = vec_col_diablo) 
      if (bool_save) {
        dev.off()
        pdf(paste0(res_dir, "ROC_test_all_data_types_comp_", comp, "_plus_avg.pdf"), 
            width = 6.9, height = 4, paper='special')
      }
      my_auroc(top_res, newdata = list_X_test, 
               outcome.test = sub_df_info_test[, dep_var],
               roc.comp = comp, bool_multiple = TRUE, bool_add_avg = TRUE, 
               vec_col = vec_col_diablo)
      if (bool_save) {
        dev.off()
        pdf(paste0(res_dir, "ROC_test_all_data_types_comp_", comp, "_plus_weight.pdf"), 
            width = 6.9, height = 4,  paper='special')
      }
      my_auroc(top_res, newdata = list_X_test, 
               outcome.test = sub_df_info_test[, dep_var],
               roc.comp = comp, bool_multiple = TRUE, bool_add_weight = TRUE, 
               vec_col = vec_col_diablo)
      if (bool_save) {
        dev.off()
        pdf(paste0(res_dir, "ROC_test_weighted_y_comp_", comp, ".pdf"), 
            width = 6.9, height = 4,  paper='special')
      }
      my_auroc(top_res, newdata = list_X_test, 
               outcome.test = sub_df_info_test[, dep_var],
               roc.block = dn, roc.comp = comp, my_predict = "Weighted")
      if (bool_save) {
        dev.off()
        pdf(paste0(res_dir, "ROC_test_avg_y_comp_", comp, ".pdf"), 
            width = 6.9, height = 4,  paper='special')
      }
      my_auroc(top_res, newdata = list_X_test, 
               outcome.test = sub_df_info_test[, dep_var],
               roc.block = dn, roc.comp = comp, my_predict = "Averaged")
      if (bool_save) {
        dev.off()
      }
    }
    
  } else {
    plot(sub_df_info_test[, dep_var], my_pred$predict[[ii]][,,1], pch = 20)
    plot(sub_df_info_test[, dep_var], my_pred$predict[[ii]][,,2], pch = 20)
    print(mean((sub_df_info_test[, dep_var] - my_pred$predict[[ii]][,,1])^2))
    print(mean((sub_df_info_test[, dep_var] - my_pred$predict[[ii]][,,2])^2))
  }
  
} else if (discrete_resp) {
  # not available for continous variables
  
  n_cpus <- 4
  set.seed(seed) # for reproducibility in this vignette
  res_perf <- perf(top_res, validation = 'Mfold', folds = 3, 
                   nrepeat = 50, auc = TRUE, n_cpus = n_cpus) 
  
  
  res_perf$MajorityVote.error.rate
  res_perf$auc
  
  if (bool_save) {
    pdf(paste0(res_dir, "performance_plot.pdf"), 
        width = 6, height = 6, paper='special')
  }
  # performance plot
  #
  plot(res_perf)
  if (bool_save) {
    dev.off()
  }
  
}

if(bool_save) {
  save.image(file = file.path(res_dir, "output.RData"))
}


if (bool_save) {
  pdf(paste0(res_dir, "cor_circle.pdf"), width = 11, height = 10, paper='special')
}
plotVar(top_res, style = 'graphics', legend = TRUE, col = vec_col_diablo)
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "circos.pdf"), width = 6.5, height = 6.5, paper='special')
}
par(mar = c(4,4,4,4))
circosPlot(top_res, cutoff=0.75, color.blocks = vec_col_diablo, 
           line = T, size.labels = 0.000001, size.legend = 0, #showIntraLinks = T,
           size.variables = 0.8, color.Y = vec_col_gr_2, color.cor = c("red", "blue"), 
           mar=c(10,10,2,2)) # showIntraLinks = T)
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "cor_panels.pdf"), width = 10, height = 10, paper='special')
}
Y <- top_res$Y
plotDiablo(top_res, ncomp = 1)
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "network.pdf"), width = 10, height = 10, paper='special')
}
network(top_res, blocks = 1:4,  cutoff = 0.75,  cex.node.name = 0.5)
if (bool_save) {
  dev.off()
}


names(top_res$X) <- c("AA - Tryptophan", "Cell types", 
                          "Glyco- & Lipo-proteins", "Ratios", "Covariates")
top_res$names$blocks <- names(top_res$X) 
names(vec_col_diablo) 
names(list_X_test) <- names(top_res$X) 
names(top_res$ncomp) <- names(top_res$X) 
names(top_res$loadings) <- c(names(top_res$X), "Y")


require(ggpubr)

for (comp in 1:2) {
  all_out <- NULL
  for (bl in c(2, 1, 3, 4, 5)) {
    
    out <- plotLoadings(top_res,
                        block = bl,
                        comp = comp,
                        legend.col = vec_col_gr_2,
                        contrib = 'max',
                        size.name = 0.5,
                        legend = FALSE, 
                        plot = FALSE)
    
    
    out$abs_importance <- abs(out$importance)
    out$names <- rownames(out)
    out$type <- vec_col_diablo[[bl]]
    out <- out[order(out$abs_importance, decreasing = T), ]
    
    if (bool_save) {
      pdf(paste0(res_dir, "plot_dot_loadings_bl_", bl, "_comp_", comp, ".pdf"), 
          width = 3.6, height = 10, paper='special')
    }
    p <- ggdotchart(out, x = "names", y = "abs_importance",
                    color = out$color,
                    sorting = "none",
                    add = "segments",
                    rotate = FALSE,
                    group = "color",
                    dot.size = 6,
                    xlab = "",
                    ylab = "",
                    add.params = list(color = out$color, size = 1.25), 
                    ggtheme = theme_pubr(),
                    ylim = c(0, 0.85))
    print(p)
    
    if (bool_save) {
      dev.off()
    }
    all_out <- rbind(all_out, out)
    
  }
  all_out <- purrr::map_df(all_out, rev)
  if (bool_save) {
    pdf(paste0(res_dir, "plot_dot_loadings_all_blocks_comp_", comp, ".pdf"), 
        width = ifelse(comp == 1, 6.5, 7.3), height = 7.5, paper='special')
  }
  p <- ggdotchart(all_out, x = "names", y = "abs_importance",
                  color = all_out$color,
                  sorting = "none",
                  add = "segments",
                  rotate = TRUE,
                  dot.size = 6,
                  xlab = "",
                  ylab = "",
                  add.params = list(color = "grey20", size = 1), 
                  ggtheme = theme_pubr())
  print(p)
  
  if (bool_save) {
    dev.off()
  }
  
}
