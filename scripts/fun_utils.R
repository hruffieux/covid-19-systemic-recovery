create_named_list <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}

filter_invalid_test_covariates <- function(vec_test_covariates, selected_severity_groups) {
  
  if (!("HC" %in% selected_severity_groups)) {
    vec_test_covariates <- setdiff(vec_test_covariates, "covid_status") 
  }
  
  if (!("A" %in% selected_severity_groups) | !any(c("B", "C", "D", "E") %in% selected_severity_groups)) {
    vec_test_covariates <- setdiff(vec_test_covariates, "asymptomatics_vs_symptomatics") 
  }
  
  if (!("E" %in% selected_severity_groups) | !any(c("A", "B", "C", "D") %in% selected_severity_groups)) {
    vec_test_covariates <- setdiff(vec_test_covariates, 
                                   c("non_severe_vs_severe", "alive_vs_deceased")) 
  }
  
  if (length(intersect(c("A", "B", "C", "D", "E"), selected_severity_groups)) < 2) {
    vec_test_covariates <- setdiff(vec_test_covariates, "sorted_severity") 
  }
  
  
  if (!any(c("C", "D", "E") %in% selected_severity_groups) | 
      !any(c("A", "B", "E") %in% selected_severity_groups) |
      !any(c("A", "B", "C", "D") %in% selected_severity_groups)) {
    vec_test_covariates <- setdiff(vec_test_covariates, 
                                   c("screening_vs_hospital", "alive_vs_deceased")) # cannot be tested as "E" excluded
  }
  
  vec_test_covariates
}


choose_severity_groups <- function(df_comb, severity_groups) {
  
  stopifnot(all(severity_groups %in% c("HC", "A", "B", "C", "D", "E")))
  
  df_comb <- df_comb[df_comb$severity %in% severity_groups,, drop = FALSE]
  df_comb$severity <- factor(df_comb$severity, levels = severity_groups)
  
  print(paste0("Retained levels: ", paste0(levels(df_comb$severity), collapse = ", ")))
  print(table(df_comb$severity))
  
  df_comb
  
}

get_early_or_late_samples <- function(df_comb, vec_var, bool_early, 
                                      trunc_days = NULL, 
                                      drop_HC = FALSE, nzf_thres = NULL,
                                      single_sample_per_subject = TRUE,
                                      enforce_subj_w_samples_avail_after_trunc_days = FALSE) {
  if (!bool_early) {
    stopifnot(!enforce_subj_w_samples_avail_after_trunc_days) # not used for late samples.
  }
  
  # get first sample from each subjects, if it is < trunc_days
  # if two samples from a same subject < trunc_days, only the first is kept.
  if (!is.null(vec_var)) {
    stopifnot(all(vec_var %in% colnames(df_comb)))
  } else {
    stopifnot(is.null(nzf_thres))
  }
  
  stopifnot(all(!is.na(df_comb$days_from_sx_or_swab_imputed[df_comb$severity != "HC"])))
  
  
  if (bool_early & enforce_subj_w_samples_avail_after_trunc_days) { 
    
    df_comb <- filter_out_subj_with_no_late_samples(df_comb, trunc_days)
    
  }
  
  if (bool_early) {
    df_comb <- df_comb[order(df_comb$days_from_sx_or_swab_imputed),]
  } else {
    df_comb <- df_comb[order(df_comb$days_from_sx_or_swab_imputed, decreasing = TRUE),]
  }
  
  if (single_sample_per_subject) {
    df_comb <- df_comb[!duplicated(df_comb$subject_id),, drop = FALSE]
  }
  
  if (!bool_early) { # reorder increasingly
    df_comb <- df_comb[order(df_comb$days_from_sx_or_swab_imputed),]
  }
  
  if (!is.null(trunc_days)) {
    
    if (bool_early) {
      bool_keep <- df_comb$days_from_sx_or_swab_imputed <= trunc_days
    } else {
      bool_keep <- df_comb$days_from_sx_or_swab_imputed > trunc_days
    }
    
  } else {
    bool_keep <- rep(TRUE, length(df_comb$days_from_sx_or_swab_imputed))
  }
  
  if (drop_HC) {
    bool_keep <- bool_keep & df_comb$severity != "HC"
  } else {
    bool_keep <- bool_keep | df_comb$severity == "HC"
  }
  
  df_sub_comb <- df_comb[bool_keep, , drop = FALSE]
  
  if (!is.null(vec_var)) {
    df_sub <- df_sub_comb[, vec_var, drop = FALSE]
    
    na_df_sub <- apply(df_sub, 2, function(vv) all(is.na(vv)))
    if (any(na_df_sub)) {
      na_var <- names(na_df_sub)[na_df_sub]
      warning(paste0("Removed all NA variable(s): ", na_var, "\n"))
      vec_var <- vec_var[!(vec_var %in% na_var)]
      df_sub_comb <- df_sub_comb[, !(colnames(df_sub_comb) %in% na_var), drop = FALSE] 
      df_sub <- df_sub[, vec_var, drop = FALSE] 
    }
    
    
    var_df_sub <- apply(df_sub, 2, function(vv) var(vv, na.rm = TRUE))
    if (any(is.na(var_df_sub)) | any(var_df_sub == 0)) {
      cst_var <- names(var_df_sub)[is.na(var_df_sub) | var_df_sub == 0]
      warning(paste0("Removed constant variable(s): ", cst_var, "\n"))
      vec_var <- vec_var[!(vec_var %in% cst_var)]
      df_sub_comb <- df_sub_comb[, !(colnames(df_sub_comb) %in% cst_var), drop = FALSE] 
      df_sub <- df_sub[, vec_var, drop = FALSE] 
    }
    
    if (!is.null(nzf_thres)) {
      list_exclude <- exclude_low_nzf(df_sub_comb, vec_var, nzf_thres)
      vec_var <- list_exclude$vec_var
      df_sub_comb <- list_exclude$df_comb
      df_sub <- df_sub[, vec_var, drop = FALSE] 
    }
  } else {
    df_sub <- vec_var <- NULL
  }
  
  list("df_comb" = df_sub_comb, "df" = df_sub, "vec_var" = vec_var)
}


run_mixed_models <- function(df_comb, test_covariate, nuis_covariates, var_names, 
                             fdr_thres = 0.05, nzf_thres = NULL, 
                             excluded_severity_levels = NULL) {
  
  
  stopifnot(nuis_covariates %in% 
              c("gender", "age", "ethnicity", "log_CRP", "IL6", "bmi", "log_bmi"))
  
  stopifnot(all(var_names %in% names(df_comb)))
  
  
  if (test_covariate %in% c("log_CRP", "bmi", "log_bmi", "IL6")) {
    nuis_covariates <- setdiff(nuis_covariates, test_covariate)
  } else if (test_covariate == "white_ethnicity_vs_other") {
    nuis_covariates <- setdiff(nuis_covariates, "ethnicity")
  }
  covariates <- c(test_covariate, nuis_covariates)
  
  df_comb <- add_factor_conversion(df_comb, "gender") 
  df_comb <- add_factor_conversion(df_comb, "ethnicity") 
  
  
  if (test_covariate %in% c("log_CRP", "bmi", "log_bmi")) {
    if (!is.null(excluded_severity_levels)) {
      df_comb <- df_comb[!(df_comb$severity %in% excluded_severity_levels),,drop = FALSE] 
    }
    
  } else {
    df_comb <- add_grouped_factor(df_comb, output_var = test_covariate, excluded_severity_levels) 
  }
  
  vv <- unique(df_comb[, test_covariate])
  vv <- vv[!is.na(vv)]
  if (length(vv)>1) {
    
    if (!is.null(nzf_thres)) {
      list_exclude <- exclude_low_nzf(df_comb, var_names, nzf_thres)
      df_comb <- list_exclude$df_comb
      var_names <- list_exclude$vec_var
    }
    
    if (test_covariate %in%  var_names) {
      var_names <- setdiff( var_names, test_covariate)
    }
    
    formula_colnames <- paste0("X", gsub("\\+", ".", 
                                         gsub(":", ".",  
                                              gsub("\\(", ".", 
                                                   gsub("\\)", ".", 
                                                        gsub(" ", ".", 
                                                             gsub("-", ".", var_names))))))) # colnames starting with a digit cause an error
    stopifnot(all.equal(names(df_comb)[match(var_names, names(df_comb))], var_names))
    names(df_comb)[match(var_names, names(df_comb))] <- formula_colnames
    
    # require(nlme) # directly provides p-values, unlike lme4
    require(lmerTest)
    require(dplyr)
    
    
    vec_pval <- vec_est <- vec_se <- NULL
    for (ii in seq_along(formula_colnames)) {
      
      dep_var_name <- formula_colnames[ii]
      plain_dep_var_name <- var_names[ii]
      
      print(paste0("Dependent variable: ", plain_dep_var_name))
      
      formula_string <- paste0(dep_var_name, " ~ ", 
                               paste(covariates, sep = "", collapse = " + "),
                               " + (1 | subject_id)")
      print(substring(formula_string, 2))
      
      
      df_comb_dep <- df_comb[!is.na(df_comb[, dep_var_name]),, drop = FALSE]
      
      out_lmm <- lmer(formula_string, # dep_var ~ covid_status + gender + age + (1|subject_id),
                      data=df_comb_dep,
                      REML=F)
      
      sum_dep_var <- summary(out_lmm)
      
      if ("bmi" %in% covariates) {
        if (sum_dep_var$coefficients["bmi",5]<0.05) {
          print(paste0("SIGNIF FOR", dep_var_name))
          print(sum_dep_var)
        }
      }
      
      an_dep_var <- anova(out_lmm)
      
      pval <- an_dep_var[test_covariate, "Pr(>F)"]
      
      print(test_covariate)
      if (test_covariate == "covid_status") {
        fixed_effect <- "covid_statusPositive"
      } else if (test_covariate %in% c("log_CRP", "bmi", "log_bmi")) {
        fixed_effect <- test_covariate
      } else {
        stop("test_covariate not implemented.")
      }
    
      est <- sum_dep_var$coefficients[fixed_effect, "Estimate"]
      se <- sum_dep_var$coefficients[fixed_effect, "Std. Error"]
      
      
      vec_pval <- c(vec_pval, pval)
      vec_est <- c(vec_est, est)
      vec_se <- c(vec_se, se)
      
    }
    stopifnot(all.equal(names(df_comb)[match(formula_colnames, names(df_comb))], formula_colnames))
    names(df_comb)[match(formula_colnames, names(df_comb))] <- var_names
    names(vec_pval) <- names(vec_est) <- var_names 
    
    vec_adj_pval <- p.adjust(vec_pval, method = "fdr")
    
    sign <- ifelse(vec_est > 0, "Upregulated", "Downregulated")
    sign[vec_adj_pval > fdr_thres] <- "NS"
    
    df_res <- data.frame(var_names, 
                         vec_est, 
                         vec_se,
                         vec_pval,
                         vec_adj_pval,
                         sign, stringsAsFactors = FALSE)
    
    names(df_res) <- c("var_name", "estimate", "se", "pval", "adj_pval", "sign")
    
  } else {
    df_res <- NULL
  }
  
  
  create_named_list(df_comb, df_res, var_names)
}

set_analysis_title <- function(test_covariate, selected_severity_groups) {
  
  if (test_covariate == "covid_status") {
    if (isTRUE(all.equal(sort(selected_severity_groups), 
                         sort(c("HC", "A", "B", "C", "D", "E"))))) {
      title <- "COVID+ vs HC"
    } else if (isTRUE(all.equal(sort(selected_severity_groups), 
                                sort(c("HC", "B", "C", "D", "E"))))) {
      title <- "COVID+ symptomatics vs HC"
    } else if (isTRUE(all.equal(sort(selected_severity_groups), 
                                sort(c("HC", "A", "B", "C", "D"))))) {
      title <- "Non severe vs HC"
    } else if (isTRUE(all.equal(sort(selected_severity_groups), 
                                sort(c("HC", "A"))))) {
      title <- "Asymptomatics vs HC" 
    } else if (isTRUE(all.equal(sort(selected_severity_groups), 
                                sort(c("HC", "E"))))) {
      title <- "Severe vs HC" 
    } else if (isTRUE(all.equal(sort(selected_severity_groups), 
                                sort(c("HC", "A", "B"))))) {
      title <- "Mild vs HC"
    } else if (isTRUE(all.equal(sort(selected_severity_groups), 
                                sort(c("HC", "C", "D"))))) {
      title <- "Moderate vs HC"
    } else {
      stop("Invalid comparison.")
    }
  } else if (test_covariate == "log_CRP") {
    title <- paste0("CRP (", 
                    paste0(intersect(c("A", "B", "C", "D", "E"), 
                                     selected_severity_groups), 
                           collapse = ", "), ")")
  } else if (test_covariate %in% 
             c("bmi", "log_bmi", vec_cytokines, vec_cplt, vec_ig, vec_glyc)) {
    title <- paste0(test_covariate, " (", 
                    paste0(intersect(c("A", "B", "C", "D", "E"), 
                                     selected_severity_groups), 
                           collapse = ", "), ")")
  } else {
    stop("test_covariate not implemented.")
  }
  
  title
}

add_factor_conversion <- function(df, var, 
                                  levels = NULL, ordered = FALSE,
                                  list_groups = NULL, grouped_var_name = NULL,
                                  drop_unlisted_levels = FALSE) {
  
  stopifnot(var %in% colnames(df))
  
  if (!is.null(list_groups)) {
    
    stopifnot(length(list_groups)>1)
    stopifnot(!is.null(names(list_groups))) # this gives the new level names
    stopifnot(is.null(levels)) # the levels should be provided as names of the list_groups
    stopifnot(!is.null(grouped_var_name))
    
    df[, grouped_var_name] <- as.character(df[, var])
    for (ii in seq_along(list_groups)) {
      
      ll <- list_groups[[ii]]
      ll_name <- names(list_groups)[ii]
      df[, grouped_var_name][df[, var] %in% ll] <- ll_name
      
    }
    var <- grouped_var_name
    
    if (drop_unlisted_levels) {
      df <- df[df[, grouped_var_name] %in% names(list_groups),, drop = FALSE]
      levels <- names(list_groups)
    } else {
      levels <- c(names(list_groups), sort(setdiff(unique(df[, grouped_var_name]), names(list_groups))))
    }
    
    
  } else {
    
    stopifnot(is.null(grouped_var_name))
    
  }
  
  if (is.null(levels)) {
    levels <- sort(unique(as.character(df[, var])))
  }
  
  df[, var] <- factor(df[, var], levels = levels, ordered = ordered)
  
  print(table(df[, var]))
  
  df
}

myEnhancedVolcano <- function (toptable, lab, x, y, selectLab = NULL, 
                               xlim = c(min(toptable[, x], na.rm = TRUE), 
                                        max(toptable[, x], na.rm = TRUE)), 
                               ylim = c(0, max(-log10(toptable[, y]), na.rm = TRUE) + 5), 
                               xlab = bquote(~Log[2] ~ "fold change"), 
                               ylab = bquote(~-Log[10] ~ italic(P)), axisLabSize = 18, 
                               title = "Volcano plot", 
                               subtitle = "Bioconductor package EnhancedVolcano", 
                               caption = paste0("Total = ", nrow(toptable), " variables"), 
                               titleLabSize = 18, subtitleLabSize = 14, captionLabSize = 14, 
                               pCutoff = 1e-05, pLabellingCutoff = pCutoff, FCcutoff = 1, 
                               FCLabellingCutoff = FCcutoff,
                               cutoffLineType = "longdash", cutoffLineCol = "black", cutoffLineWidth = 0.4, 
                               transcriptPointSize = 0.8, transcriptLabSize = 3, transcriptLabCol = "black", 
                               transcriptLabFace = "plain", transcriptLabhjust = 0, transcriptLabvjust = 1.5, 
                               boxedlabels = FALSE, shape = 19, shapeCustom = NULL, 
                               col = c("grey30", "grey30", "grey30", "forestgreen", "royalblue", "red2"), 
                               colCustom = NULL, 
                               colAlpha = 1/2, 
                               legend = c("NS", "Log2 FC -", "Log2 FC +", "P", "P & Log2 FC -",  "P & Log2 FC +"), 
                               legendPosition = "top", legendLabSize = 13, legendIconSize = 3.75, 
                               legendVisible = TRUE, shade = NULL, shadeLabel = NULL, shadeAlpha = 1/2, 
                               shadeFill = "grey", shadeSize = 0.01, shadeBins = 2, drawConnectors = FALSE, 
                               widthConnectors = 0.5, typeConnectors = "closed", endsConnectors = "first", 
                               lengthConnectors = unit(0.01, "npc"), colConnectors = "grey10", 
                               hline = NULL, hlineType = "longdash", hlineCol = "black", 
                               hlineWidth = 0.4, vline = NULL, vlineType = "longdash", vlineCol = "black", 
                               vlineWidth = 0.4, gridlines.major = TRUE, gridlines.minor = TRUE, 
                               border = "partial", borderWidth = 0.8, borderColour = "black") 
{
  require(ggplot2)
  require(ggrepel)
  
  if (!requireNamespace("ggplot2")) {
    stop("Please install ggplot2 first.", call. = FALSE)
  }
  if (!requireNamespace("ggrepel")) {
    stop("Please install ggrepel first.", call. = FALSE)
  }
  if (!is.numeric(toptable[, x])) {
    stop(paste(x, " is not numeric!", sep = ""))
  }
  if (!is.numeric(toptable[, y])) {
    stop(paste(y, " is not numeric!", sep = ""))
  }
  i <- xvals <- yvals <- Sig <- NULL
  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[(toptable[, x] > FCcutoff)] <- "FC+"
  toptable$Sig[(toptable[, x] < -FCcutoff)] <- "FC-"
  toptable$Sig[(toptable[, y] < pCutoff)] <- "P"
  toptable$Sig[(toptable[, y] < pCutoff) & (toptable[,x] > FCcutoff)] <- "FC+_P"
  toptable$Sig[(toptable[, y] < pCutoff) & (toptable[,x] < -FCcutoff)] <- "FC-_P"
  toptable$Sig <- factor(toptable$Sig, levels = c("NS", "FC-", "FC+", "P", "FC-_P", "FC+_P"))
  if (min(toptable[, y], na.rm = TRUE) == 0) {
    warning(paste("One or more P values is 0.", "Converting to minimum possible value..."), 
            call. = FALSE)
    toptable[which(toptable[, y] == 0), y] <- .Machine$double.xmin
  }
  toptable$lab <- lab
  toptable$xvals <- toptable[, x]
  toptable$yvals <- toptable[, y]
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }
  th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(), 
                                         plot.title = element_text(angle = 0, size = titleLabSize, 
                                                                   face = "bold", vjust = 1), plot.subtitle = element_text(angle = 0, 
                                                                                                                           size = subtitleLabSize, face = "plain", vjust = 1), 
                                         plot.caption = element_text(angle = 0, size = captionLabSize, 
                                                                     face = "plain", vjust = 1), axis.text.x = element_text(angle = 0, 
                                                                                                                            size = axisLabSize, vjust = 1), axis.text.y = element_text(angle = 0, 
                                                                                                                                                                                       size = axisLabSize, vjust = 1), axis.title = element_text(size = axisLabSize), 
                                         legend.position = legendPosition, legend.key = element_blank(), 
                                         legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = legendLabSize), 
                                         title = element_text(size = legendLabSize), legend.title = element_blank())
  if (!is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = factor(names(shapeCustom))), alpha = colAlpha, 
                 size = transcriptPointSize, na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = shapeCustom)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           1) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom))), 
                 alpha = colAlpha, shape = shape, size = transcriptPointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(guide = TRUE)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 6) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = factor(Sig)), alpha = colAlpha, size = transcriptPointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = c(NS = shape[1], 
                                    "FC-" = shape[2], 
                                    "FC+" = shape[3], 
                                    P = shape[4], 
                                    "FC-_P" = shape[5], 
                                    "FC+_P" = shape[6]), 
                         labels = c(NS = legend[1], 
                                    "FC-" = paste(legend[2], sep = ""), 
                                    "FC+" = paste(legend[3], sep = ""),
                                    P = paste(legend[4], sep = ""), 
                                    "FC-_P" = paste(legend[5], sep = ""),
                                    "FC+_P" = paste(legend[6], sep = "")
                         ), 
                         guide = TRUE)
  }
  else if (is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(Sig), shape = factor(names(shapeCustom))), 
                 alpha = colAlpha, size = transcriptPointSize, 
                 na.rm = TRUE) + scale_color_manual(values = c(NS = col[1], 
                                                               "FC-" = col[2], 
                                                               "FC+" = col[3], 
                                                               P = col[4], 
                                                               "FC-_P" = col[5],
                                                               "FC+_P" = col[6]),
                                                    labels = c(NS = legend[1], 
                                                               "FC-" = paste(legend[2], sep = ""), 
                                                               "FC+" = paste(legend[3], sep = ""), 
                                                               P = paste(legend[4], sep = ""), 
                                                               "FC-_P" = paste(legend[5], sep = ""),
                                                               "FC+_P" = paste(legend[6], sep = ""))) + 
      scale_shape_manual(values = shapeCustom)
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) ==  1) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals)))  + 
      th + 
      guides(colour = guide_legend(order = 1, override.aes = list(shape = shape, 
                                                                  size = legendIconSize))) + geom_point(aes(color = factor(Sig)), 
                                                                                                        alpha = colAlpha, shape = shape,
                                                                                                        size = transcriptPointSize, na.rm = TRUE, 
                                                                                                        show.legend = legendVisible) +
      scale_color_manual(values = c("NS" = col[1],
                                    "FC-" = col[2],
                                    "FC+" = col[3],
                                    "P" = col[4],
                                    "FC-_P" = col[5],
                                    "FC+_P" = col[6]),
                         labels = c("NS" = legend[1],
                                    "FC-" = paste(legend[2], sep = ""),
                                    "FC+" = paste(legend[3], sep = ""),
                                    "P" = paste(legend[4], sep = ""),
                                    "FC-_P" = paste(legend[5], sep = ""),
                                    "FC+_P" = paste(legend[6], sep = ""))
      )
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 6) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(shape = c(NS = shape[1], 
                                                                                 "FC-" = shape[2], 
                                                                                 "FC+" = shape[3], 
                                                                                 P = shape[4], 
                                                                                 "FC-_P" = shape[5], 
                                                                                 "FC+_P" = shape[6]), 
                                                                       size = legendIconSize))) + 
      geom_point(aes(color = factor(Sig), shape = factor(Sig)), 
                 alpha = colAlpha, size = transcriptPointSize, 
                 na.rm = TRUE, show.legend = legendVisible) + 
      scale_color_manual(values = c(NS = col[1], 
                                    "FC-" = col[2], 
                                    "FC+" = col[3], 
                                    P = col[4], 
                                    "FC-_P" = col[5],
                                    "FC+_P" = col[6]),
                         labels = c(NS = legend[1], 
                                    "FC-" = paste(legend[2], sep = ""), 
                                    "FC+" = paste(legend[3], sep = ""),
                                    P = paste(legend[4], sep = ""), 
                                    "FC-_P" = paste(legend[5], sep = ""),
                                    "FC+_P" = paste(legend[6], sep = "")
                         )) + 
      scale_shape_manual(values = c(NS = shape[1], 
                                    "FC-" = shape[2], 
                                    "FC+" = shape[3], 
                                    P = shape[4], 
                                    "FC-_P" = shape[5], 
                                    "FC+_P" = shape[6]), guide = FALSE)
  }
  plot <- plot + xlab(xlab) + ylab(ylab) + xlim(xlim[1], xlim[2]) + 
    ylim(ylim[1], ylim[2]) +
    geom_hline(yintercept = -log10(pCutoff), linetype = cutoffLineType, colour = cutoffLineCol, size = cutoffLineWidth) # +
  # geom_vline(xintercept = c(-FCcutoff, 
  #                                                    FCcutoff), linetype = cutoffLineType, colour = cutoffLineCol, 
  #                                     size = cutoffLineWidth) + 
  
  plot <- plot + labs(title = title, subtitle = subtitle, caption = caption)
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline, linetype = vlineType, 
                              colour = vlineCol, size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline), 
                              linetype = hlineType, colour = hlineCol, size = hlineWidth)
  }
  if (border == "full") {
    plot <- plot + theme(panel.border = element_rect(colour = borderColour, 
                                                     fill = NA, size = borderWidth))
  }
  else if (border == "partial") {
    plot <- plot + theme(axis.line = element_line(size = borderWidth, 
                                                  colour = borderColour), panel.border = element_blank(), 
                         panel.background = element_blank())
  }
  else {
    stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")
  }
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line(size = 0.5))
  }
  else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }
  if (boxedlabels == FALSE) {
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_text_repel(data = subset(toptable, 
                                                   toptable[, y] < pLabellingCutoff & abs(toptable[, 
                                                                                                   x]) > FCLabellingCutoff), aes(label = subset(toptable, 
                                                                                                                                                toptable[, y] < pLabellingCutoff & abs(toptable[, 
                                                                                                                                                                                                x]) > FCLabellingCutoff)[, "lab"]), 
                                     size = transcriptLabSize, 
                                     max.overlaps = Inf,
                                     segment.color = colConnectors, segment.size = widthConnectors, 
                                     arrow = arrow(length = lengthConnectors, type = typeConnectors, 
                                                   ends = endsConnectors), hjust = transcriptLabhjust, 
                                     vjust = transcriptLabvjust, colour = transcriptLabCol, 
                                     fontface = transcriptLabFace, na.rm = TRUE)
    }
    else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_text_repel(data = subset(toptable, 
                                                   !is.na(toptable[, "lab"])), aes(label = subset(toptable, 
                                                                                                  !is.na(toptable[, "lab"]))[, "lab"]), size = transcriptLabSize, 
                                     segment.color = colConnectors, segment.size = widthConnectors, 
                                     max.overlaps = Inf,
                                     arrow = arrow(length = lengthConnectors, type = typeConnectors, 
                                                   ends = endsConnectors), hjust = transcriptLabhjust, 
                                     vjust = transcriptLabvjust, colour = transcriptLabCol, 
                                     fontface = transcriptLabFace, na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_text(data = subset(toptable, 
                                             !is.na(toptable[, "lab"])), aes(label = subset(toptable, 
                                                                                            !is.na(toptable[, "lab"]))[, "lab"]), size = transcriptLabSize, 
                               check_overlap = TRUE, hjust = transcriptLabhjust, 
                               vjust = transcriptLabvjust, colour = transcriptLabCol, 
                               fontface = transcriptLabFace, na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_text(data = subset(toptable, 
                                             toptable[, y] < pLabellingCutoff & abs(toptable[, 
                                                                                             x]) > FCLabellingCutoff), aes(label = subset(toptable, 
                                                                                                                                          toptable[, y] < pLabellingCutoff & abs(toptable[, 
                                                                                                                                                                                          x]) > FCLabellingCutoff)[, "lab"]), size = transcriptLabSize, 
                               check_overlap = TRUE, hjust = transcriptLabhjust, 
                               vjust = transcriptLabvjust, colour = transcriptLabCol, 
                               fontface = transcriptLabFace, na.rm = TRUE)
    }
  }
  else {
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                                    toptable[, y] < pLabellingCutoff & abs(toptable[, 
                                                                                                    x]) > FCLabellingCutoff), aes(label = subset(toptable, 
                                                                                                                                                 toptable[, y] < pLabellingCutoff & abs(toptable[, 
                                                                                                                                                                                                 x]) > FCLabellingCutoff)[, "lab"]), size = transcriptLabSize, 
                                      segment.color = colConnectors, segment.size = widthConnectors, 
                                      max.overlaps = Inf,
                                      arrow = arrow(length = lengthConnectors, type = typeConnectors, 
                                                    ends = endsConnectors), hjust = transcriptLabhjust, 
                                      vjust = transcriptLabvjust, colour = transcriptLabCol, 
                                      fontface = transcriptLabFace, na.rm = TRUE)
    }
    else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                                    !is.na(toptable[, "lab"])), aes(label = subset(toptable, 
                                                                                                   !is.na(toptable[, "lab"]))[, "lab"]), size = transcriptLabSize, 
                                      segment.color = colConnectors, segment.size = widthConnectors, 
                                      max.overlaps = Inf,
                                      arrow = arrow(length = lengthConnectors, type = typeConnectors, 
                                                    ends = endsConnectors), hjust = transcriptLabhjust, 
                                      vjust = transcriptLabvjust, colour = transcriptLabCol, 
                                      fontface = transcriptLabFace, na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                              !is.na(toptable[, "lab"])), aes(label = subset(toptable, 
                                                                                             !is.na(toptable[, "lab"]))[, "lab"]), size = transcriptLabSize, 
                                hjust = transcriptLabhjust, vjust = transcriptLabvjust, 
                                colour = transcriptLabCol, fontface = transcriptLabFace, 
                                na.rm = TRUE)
    }
    else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                              toptable[, y] < pLabellingCutoff & abs(toptable[, 
                                                                                              x]) > FCLabellingCutoff), aes(label = subset(toptable, 
                                                                                                                                           toptable[, y] < pLabellingCutoff & abs(toptable[, 
                                                                                                                                                                                           x]) > FCLabellingCutoff)[, "lab"]), size = transcriptLabSize, 
                                hjust = transcriptLabhjust, vjust = transcriptLabvjust, 
                                colour = transcriptLabCol, fontface = transcriptLabFace, 
                                na.rm = TRUE)
    }
  }
  if (!is.null(shade)) {
    plot <- plot + stat_density2d(data = subset(toptable, 
                                                rownames(toptable) %in% shade), fill = shadeFill, 
                                  alpha = shadeAlpha, geom = "polygon", contour = TRUE, 
                                  size = shadeSize, bins = shadeBins, show.legend = FALSE, 
                                  na.rm = TRUE) + scale_fill_identity(name = shadeLabel, 
                                                                      labels = shadeLabel, guide = "legend")
  }
  
  if (is.null(legend)) {
    plot <- plot + ggpubr::rremove("legend")
  }
  return(plot)
}

add_grouped_factor <- function(df_comb, output_var, excluded_levels = NULL) {
  
  drop_unlisted_levels <- TRUE
  
    
  if (output_var == "covid_status") {
    
    input_var <- "severity"
    grouped_var_name <- output_var
    levels <- NULL
    ordered <- FALSE
    list_groups <- list("Negative" = setdiff("HC", excluded_levels),
                        "Positive" = setdiff(c("A", "B", "C", "D", "E"), excluded_levels))
    
  } else {
    
    stop("Grouping for output_var not implemented.")
    
  }
  
  
  print(list_groups)
  list_groups <- list_groups[sapply(list_groups, length)>0]
  
  
  stopifnot(length(list_groups) > 1)
  
  df_comb <- add_factor_conversion(df_comb, input_var,
                                   levels = levels, 
                                   ordered = ordered,
                                   list_groups = list_groups,
                                   grouped_var_name = grouped_var_name,
                                   drop_unlisted_levels = drop_unlisted_levels)
  
  df_comb
}

get_mface_data <- function(df_comb, vec_var, min_nb_timepoints = 2) {
  
  list_subjects <- lapply(vec_var, function(vv) {
    
    df_sub <- df_comb[!is.na(df_comb[,vv]),, drop = FALSE]
    tb <- table(df_sub$subject_id)
    names(tb)[tb >= min_nb_timepoints] # subjects with enough timepoints
    
  })
  
  common_subjects <- Reduce(intersect, list_subjects)
  
  if (length(common_subjects)==0) {
    stop("Stop. No subjects common to all variables.")
  }
  
  data <- lapply(vec_var, function(vv) {
    
    df_var <- df_comb[!is.na(df_comb[,vv]),, drop = FALSE]
    df_var <- df_var[df_var$subject_id %in% common_subjects,, drop = FALSE]
    
    vec_sub <- c("subject_id", "days_from_sx_or_swab_imputed", vv)
    df_var <- subset(df_var, select = vec_sub)
    rownames(df_var) <- NULL
    names(df_var) <- c("subj", "argvals", "y")
    df_var
  })
  names(data) <- paste0("y", 1:length(data))
  
  # check that all variables have measurements for each subject 
  # (otherwise potential subject coding mismatch!)
  stopifnot(length(unique(lapply(data, function(ll) sort(unique(ll$subj)))))==1)
  
  if (length(vec_var) == 1) {
    data[[1]]
  } else {
    data
  }
  
}

get_mface_scores <- function(fit, list_data, tnew, df_comb, vec_col) {
  
  if (is.data.frame(list_data)) {
    data <- list("y1" = list_data) # in case of univariate fpca
  } else {
    data <- list_data # in case of multivariate fpca
  }
  
  # the same subjects for all metabolites, thanks to the stopifnot check above
  vec_subj <- unique(data[[1]]$subj)
  
  scores <- NULL
  for(i in vec_subj){
    
    sel <- lapply(data, function(x){which(x$subj==i)})
    dat_i <- mapply(function(data, sel){data[sel,]}, 
                    data = data, sel = sel, SIMPLIFY = FALSE)
    dat_i_pred <- lapply(dat_i, function(x){
      data.frame(subj=rep(x$subj[1],nrow(x) + length(tnew)),
                 argvals = c(rep(NA,nrow(x)),tnew),
                 y = rep(NA,nrow(x) + length(tnew)))
    })
    for(j in 1:length(dat_i)){
      dat_i_pred[[j]][1:nrow(dat_i[[j]]), ] <- dat_i[[j]]
    }
    
    if (is.data.frame(list_data)) {
      dat_i_pred <- dat_i_pred[[1]]
    }
    pred <- predict(fit, dat_i_pred)
    scores <- rbind(scores, pred$scores$scores)
  }
  colnames(scores) <- paste0("EF", 1:ncol(scores))
  scores <- as.data.frame(scores)
  scores$subject_id <- vec_subj
  
  df_sub_info <- subset(df_comb, select = c("subject_id", "severity", "hospital_outcome"))
  rownames(df_sub_info) <- NULL
  df_sub_info <- unique(df_sub_info)
  df_sub_info$subject_id <- as.character(df_sub_info$subject_id)
  scores <- dplyr::left_join(scores, df_sub_info)
  
  df_col <- data.frame("severity" = names(vec_col), "color_severity" = vec_col)
  scores <- dplyr::left_join(scores, df_col)
  rownames(scores) <- vec_subj
  
  scores
}

scatterplot_scores <- function(fit, scores, color, pcs = c(1, 2)) {
  
  par(mfrow=c(1,1),mar=c(4.5,5,3,2))
  plot(scores[, paste0("EF", pcs[1])], 
       scores[, paste0("EF", pcs[2])], 
       col = color,
       pch = 20, 
       xlim = c(-max(abs(scores$EF1)), max(abs(scores[, paste0("EF", pcs[1])]))),
       ylim = c(-max(abs(scores$EF2)), max(abs(scores[, paste0("EF", pcs[2])]))),
       xlab = paste0("Scores FPC ", pcs[1], " (", 
                     format(get_pve(fit, pcs[1]), digits = 3), "%)"),
       ylab = paste0(" \n Scores FPC ", pcs[2], " (", 
                     format(get_pve(fit, pcs[2]), digits = 3), "%)"),
       main = "Scores"
  )
  points(scores[!is.na(scores$hospital_outcome) & scores$hospital_outcome == "dead", 
                paste0("EF", pcs[1])], 
         scores[!is.na(scores$hospital_outcome) & scores$hospital_outcome == "dead", 
                paste0("EF", pcs[2])], 
         pch = 13, col = adjustcolor("black", alpha.f = 0.4), cex = 1)
  
}

get_subjects_to_show <- function(list_data, 
                                 scores = NULL, 
                                 nb_to_show = 4, 
                                 show_extreme_subjects = FALSE, 
                                 min_timepoints = 3) { 
  
  
  if (is.data.frame(list_data)) {
    data <- list("y1" = list_data) # in case of univariate fpca
  } else {
    data <- list_data # in case of multivariate fpca
  }
  
  if (show_extreme_subjects) {
    
    stopifnot(!is.null(scores))
    stopifnot((nb_to_show %% 2) == 0)
    
    subjects_to_show <- scores$subject_id[c(order(scores$EF2, decreasing = T)[1:(nb_to_show/2)],
                                            order(scores$EF2)[1:(nb_to_show/2)])]
    print("Extreme scores:")
    print(scores$EF2[match(subjects_to_show,scores$subject_id)])
  } else {
    
    subjects_to_show <- names(which(table(data[[1]]$subj)>=min_timepoints))
    
    if (!is.null(nb_to_show)) {
      subjects_to_show <- subjects_to_show[1:nb_to_show]
    }
    
    if (any(is.na(subjects_to_show))) {
      stop(paste0("Not enough subjects with at least ", min_timepoints, 
                  " samples. Reduce min_timepoints."))
    }
  }
  
  subjects_to_show
}


plot_subject_trajectories <- function(fit, scores, tnew, list_data, 
                                      subjects_to_show, 
                                      vec_var, df_comb, ylim_offset = 1,
                                      var_disp_names = NULL,
                                      bool_trunc = FALSE,
                                      bool_par = TRUE,
                                      main_plus = NULL,
                                      sub = NULL,
                                      ylab = NULL,
                                      cex_main = 1.5,
                                      cex_lab = 1,
                                      col_curve = NULL,
                                      col_ci = NULL,
                                      bool_sub = TRUE
) {
  
  if (is.data.frame(list_data)) {
    data <- list("y1" = list_data) # in case of univariate fpca
  } else {
    data <- list_data # in case of multivariate fpca
  }
  
  if (bool_trunc) {
    df_comb$time <- df_comb$hospital_outcome_days_from_start_sx_or_first_pos_swab
    df_comb$time[!(df_comb$hospital_outcome %in% "dead")] <- df_comb$last_observation_days_from_start_sx_or_first_pos_swab[!(df_comb$hospital_outcome %in% "dead")]
    
    df_comb$hospital_outcome[df_comb$severity %in% c("HC", "A", "B")] <- "alive"
    
    df_comb$status <- ifelse(df_comb$hospital_outcome == "dead", 2, 1) # 1 = cencored, 2 = dead
  }
  
  if (bool_par) {
    par(mar=c(5,4.5,5,7), mfrow = c(length(subjects_to_show), length(vec_var)))
  }
  
  # for y-axis limits
  sel_subj_to_show <- lapply(data, function(x){which(x$subj %in% subjects_to_show)})
  # print(sel_subj_to_show)
  dat_subj_to_show <- mapply(function(data, sel_subj_to_show){data[sel_subj_to_show,]}, 
                             data = data, sel_subj_to_show = sel_subj_to_show, SIMPLIFY = FALSE)
  
  if (is.null(col_curve)) {
    col_curve <- rep("red", length(subjects_to_show))
  } else if (length(col_curve)==1) {
    col_curve <- rep(col_curve, length(subjects_to_show))
  }
  names(col_curve) <- subjects_to_show
  
  if (is.null(col_ci)) {
    col_ci <- "blue"
  }
  
  for(i in subjects_to_show){
    
    sev <- unique(df_comb$severity[df_comb$subject_id %in% i])
    sel <- lapply(data, function(x){which(x$subj==i)})
    dat_i <- mapply(function(data, sel){data[sel,]}, 
                    data = data, sel = sel, SIMPLIFY = FALSE)
    dat_i_pred <- lapply(dat_i, function(x){
      data.frame(subj=rep(x$subj[1],nrow(x) + length(tnew)),
                 argvals = c(rep(NA,nrow(x)),tnew),
                 y = rep(NA,nrow(x) + length(tnew)))
    })
    for(j in 1:length(dat_i)){
      dat_i_pred[[j]][1:nrow(dat_i[[j]]), ] <- dat_i[[j]]
    }
    
    if (is.data.frame(list_data)) {
      dat_i_pred <- dat_i_pred[[1]]
    }
    
    pred <- predict(fit, dat_i_pred)
    y_pred <- mapply(function(pred_y.pred, dat_i){
      pred_y.pred[nrow(dat_i)+1:length(tnew)]}, pred_y.pred = pred$y.pred, 
      dat_i = dat_i, SIMPLIFY = TRUE)
    
    pre <- pred
    
    for (k in seq_along(data)) {
      
      vv <- vec_var[k]
      
      if (is.null(var_disp_names)) {
        var_disp_name <- vv
      } else {
        var_disp_name <- var_disp_names[k]
      }
      
      if (is.data.frame(list_data)) {
        y_pred_k <- pre$y.pred
        se_pred_k <- pre$se.pred
      } else {
        y_pred_k <- pre$y.pred[[paste0("y", k)]] 
        se_pred_k <- pre$se.pred[[paste0("y", k)]]
      }
      
      if (any(df_comb$severity == "HC")) {
        lo <- quantile(df_comb[df_comb$severity == "HC", vv], probs = 0.25, na.rm = TRUE)
        up <- quantile(df_comb[df_comb$severity == "HC", vv], probs = 0.75, na.rm = TRUE)
        
        mmax <- max(max(dat_subj_to_show[[k]]$y, na.rm = T), up, na.rm = T)
        mmin <- min(min(dat_subj_to_show[[k]]$y, na.rm = T), lo, na.rm = T)
      } else {
        mmax <- max(dat_subj_to_show[[k]]$y, na.rm = T)
        mmin <- min(dat_subj_to_show[[k]]$y, na.rm = T)
      }
      
      Ylim = c(max(0.5, mmin) - ylim_offset, mmax + ylim_offset)
      Xlim = c(0,max(data[[k]]$argvals)+1)
      Ylab = ifelse(is.null(ylab), paste0(var_disp_name, " (log-scale)"), ylab)#bquote(y^(1))
      Xlab = "Days from symptom onset"
      
      main <- var_disp_name
      if (!is.null(main_plus)) {
        main <- paste0(main, main_plus)
      } else {
        main <- paste0(main, ", subj. ", i, "\n (sev.: ", sev, 
                       ", gender: ", unique(df_comb$gender[df_comb$subject_id %in% i]), 
                       ", age: ", unique(df_comb$age[df_comb$subject_id %in% i]), 
                       ")")
      }
      
      if (is.null(sub)) {
        sub_i <- paste0("Scores 1st EF: ", format(scores$EF1[scores$subject_id %in% i], digits =2),
                        ", 2nd EF: ", format(scores$EF2[scores$subject_id %in% i], digits =2))
      } else {
        sub_i <- sub
      }
      
      idx = (nrow(dat_i[[k]])+1):(nrow(dat_i[[k]])+length(tnew))
      
      if (bool_sub) {
        plot(dat_i[[k]][,"argvals"],dat_i[[k]][,"y"],ylim=Ylim,xlim=Xlim,ylab=Ylab,xlab=Xlab,
             main=main,cex.lab=cex_lab,cex.axis = 1.0,cex.main = cex_main,pch=1, 
             sub = sub_i)
      } else {
        plot(dat_i[[k]][,"argvals"],dat_i[[k]][,"y"],ylim=Ylim,xlim=Xlim,
             # ylab=paste0(" \n \n", Ylab, "\n", sub_i), 
             ylab = Ylab, xlab=Xlab,
             main=main,cex.lab=cex_lab,cex.axis = 1.0,cex.main = cex_main,pch=1)
        legend("top", sub_i, bty = "n")
      }
      
      
      if (any(df_comb$severity == "HC")) {
        rect(-5, lo, days_thres+5, up, col = adjustcolor("gray", alpha.f=0.5), border = "gray")
      }
      
      status <- unique(df_comb$status[df_comb$subject_id %in% i])
      time <- unique(df_comb$time[df_comb$subject_id %in% i])
      time <- time[!is.na(time)]
      
      stopifnot(length(status) == 1 & length(time) == 1)
      
      
      if (status == 2 & time < 50) { # cut prediction if patient died within the analysis window
        lines(tnew[tnew<=time], y_pred_k[idx][tnew<=time],col=col_curve[i],lwd=2)
        points(time,  y_pred_k[idx][time], pch = 4, cex = 1.25, lwd = 1.5)
        lines(tnew[tnew<=time], y_pred_k[idx][tnew<=time]-1.96*se_pred_k[idx][tnew<=time],
              col=col_ci,lwd=2,lty=2)
        lines(tnew[tnew<=time], y_pred_k[idx][tnew<=time]+1.96*se_pred_k[idx][tnew<=time],
              col=col_ci,lwd=2,lty=2)
      } else {
        lines(tnew, y_pred_k[idx],col=col_curve[i],lwd=2)
        lines(tnew,y_pred_k[idx]-1.96*se_pred_k[idx],col=col_ci,lwd=2,lty=2)
        lines(tnew,y_pred_k[idx]+1.96*se_pred_k[idx],col=col_ci,lwd=2,lty=2)

      }
      
    }
    
  }
  
}


plot_eigenfunctions <- function(eigenfunctions, eigenvalues, list_data, vec_var, days_thres) {
  
  if (is.data.frame(list_data)) {
    data <- list("y1" = list_data) # in case of univariate fpca
  } else {
    data <- list_data # in case of multivariate fpca
  }
  
  par(mfrow = c(1,length(vec_var)))
  for (k in seq_along(data)) {
    
    vv <- vec_var[k]
    lim_eigen <- c(min(c(eigenfunctions[((k-1)*days_thres + 1):(k*days_thres),1:2], 0)),
                   max(c(eigenfunctions[((k-1)*days_thres + 1):(k*days_thres),1:2], 0)))
    plot(eigenfunctions[((k-1)*days_thres + 1):(k*days_thres),1], 
         xlab = "Days from swab or symptom onset",
         ylab = "Eigenfunctions",
         main = vv,
         ylim = lim_eigen, type = "l", col = "black", lwd = 2)
    abline(h = 0, col = "grey", lwd = 2, lty = 3)
    lines(eigenfunctions[((k-1)*days_thres + 1):(k*days_thres),2], 
          type = "l", lwd = 2, col = "black", lty = 2)
    legend("bottomleft", 
           col = c("black", "black"),
           lwd = 2,
           lty = c(1,2),
           bty = "n",
           legend = paste0(c("1st: ", "2nd: "), #, "3rd"), 
                           format(eigenvalues[1:2] / sum(eigenvalues)*100, digits = 3), "%")) 
  }
  
}


plot_variance <- function(fit, list_data, tnew, vec_var, disp_param = TRUE, 
                          var_disp_names = NULL) {
  
  if (is.data.frame(list_data)) {
    data <- list("y1" = list_data) # in case of univariate fpca
  } else {
    data <- list_data # in case of multivariate fpca
  }
  
  Cov <- as.matrix(fit$Chat.new)
  Cov_diag <- diag(Cov)
  
  if (disp_param) {
    par(mfrow=c(1, length(vec_var)),mar=c(4.5,4.1,3,4.5))
  }
  
  Xlab = "Days from symptom onset"
  
  if (!is.null(var_disp_names)) {
    vec_var <- var_disp_names
  }
  
  for (k in seq_along(data)) {
    vv <- vec_var[k]
    plot(tnew,Cov_diag[seq_along(tnew)+length(tnew)*(k-1)],type="l",
         xlab = Xlab, ylab="variance",main=vv,
         cex.axis=1.25,cex.lab=1.25,cex.main=1.25,lwd=2)
  }
  
}


plot_correlation <- function(fit, list_data, tnew, vec_var, nb_spaces = 3, nb_ticks = 4,
                             var_disp_names = NULL, sub_var = NULL) {
  
  if (is.data.frame(list_data)) {
    data <- list("y1" = list_data) # in case of univariate fpca
  } else {
    data <- list_data # in case of multivariate fpca
  }
  
  if (!is.null(var_disp_names)) {
    vec_var <- var_disp_names
  }
  
  if (!is.null(sub_var)) {
    data <- data[sub_var]
    ind_cor <- Reduce(c, sapply(sub_var, function(id) (id-1)*length(tnew) + (tnew+1)))
    fit$Cor.new <- fit$Cor.new[ind_cor, ind_cor, drop = FALSE]
    vec_var <- vec_var[sub_var]
  }
  
  Cor <- as.matrix(fit$Cor.new)
  
  Xlab = "Days from symptom onset"
  
  require(fields)
  par(mar=c(5,4.5,4,7))
  par(mfrow = c(1,1))
  mycols <- colorRampPalette(colors = c("blue","white", "red"))(200)
  
  nb_var <- length(data)
  vec_trunc <- stringr::str_trunc(vec_var, 25)
  vec_spaces <- ceiling(nb_spaces / length(vec_var) - nchar(vec_trunc))
  
  main <- paste0(sapply(1:length(vec_var), 
                        function(i) paste0(paste0(rep(" ", max(2, ceiling(vec_spaces[i]/2))), collapse = ""),
                                           vec_trunc[i], 
                                           paste0(rep(" ", max(2, ceiling(vec_spaces[i]/2))), collapse = ""))),
                 collapse = "")
  
  image(Cor,axes=F, col=mycols, xlab=Xlab, ylab = Xlab,
        main = main, cex.main = 0.95)
  
  ax <- as.vector(sapply(seq_along(data), function(k) {
    seq(0, 1/nb_var, l = nb_ticks) + (k-1)/nb_var
  }))
  axis(1,at=ax,labels=format(rep(seq(0, max(tnew), l = nb_ticks), times = nb_var), digits = 2))
  axis(2,at=ax,labels=format(rep(seq(0, max(tnew), l = nb_ticks), times = nb_var), digits = 2))
  
  image.plot(1:(nb_var*length(tnew)), 
             1:(nb_var*length(tnew)), 
             as.matrix(Cor),
             col=mycols,
             cex.axis=1.25,cex.lab=1,cex.main=1,
             axis.args = list(at = c(-1.0, -0.5, 0, 0.5,1.0)),
             legend.shrink=0.75,legend.line=-1.5, legend.only = T)
  
  
}



my_auroc = function(
    object,
    newdata = object$X,
    outcome.test = as.factor(object$Y),
    multilevel = NULL,
    plot = TRUE,
    roc.block = 1,
    roc.comp = 1,
    my_predict = NULL,
    bool_multiple = FALSE,
    bool_add_avg = FALSE,
    bool_add_weight = FALSE,
    vec_col = NULL,
    ...)
{
  
  data=list()
  auc.mean = graph=list()
  data$outcome=factor(outcome.test)
  
  # note here: the dist does not matter as we used the predicted scores only
  
  list.predict  =  predict(object, newdata = newdata, dist = "max.dist", multilevel = multilevel)
  
  
  if (is.null(my_predict)) { # default function
    
    res.predict <- list.predict$predict
    
    if (is.null(vec_col)) {
      vec_col <- c("#999999", # "#E69F00", 
                   "#56B4E9", "#009E73", # "#F0E442", 
                   # "#0072B2", 
                   "#D55E00", "#CC79A7")[seq_along(res.predict)]
      # print(vec_col)
    }
    
    if (bool_multiple) {
      
      # controls the ordering of colors
      names(res.predict) <- paste0(seq_along(res.predict), ": ", names(res.predict))
      names(vec_col) <- names(res.predict) 
      
      if (bool_add_avg) {
        res.predict <- append(res.predict, list(list.predict$AveragedPredict))
        names(res.predict)[length(res.predict)] <- "t: average contributions"
        vec_col <- c(vec_col, "black")
        names(vec_col)[length(res.predict)] <- "t: average contributions"
      } 
      
      if (bool_add_weight) {
        res.predict <- append(res.predict, list(list.predict$WeightedPredict))
        names(res.predict)[length(res.predict)] <- "t: weighted contributions"
        vec_col <- c(vec_col, "black")
        names(vec_col)[length(res.predict)] <- "t: weighted contributions"
      }
      
      # for (i in 1:object$ncomp[1]) # ii = component id (can differ between data-types, that's why depends on j)
      # {
      # data$data=lapply(res.predict, function(pp) pp[,,i])
      # title=paste("ROC Curve per data type\n", "Comp: ",i, sep="")
      data$data=lapply(res.predict, function(pp) pp[,,roc.comp])
      title=paste("ROC Curve per data type\n", "Comp: ",roc.comp, sep="")
      
      temp = my_statauc_multiple(data, plot = TRUE, title = title, vec_col = vec_col)
      
      # }
      
      out = NULL
      
    } else {
      
      block.all = names(res.predict)
      block.temp = names(res.predict[roc.block])
      
      for(j in 1:length(res.predict)) # j = data-type index
      {
        vec_col <- vec_col[j]
        for (i in 1:object$ncomp[j]) # ii = component id (can differ between data-types, that's why depends on j)
        {
          data$data=res.predict[[j]][,,i]
          title=paste("ROC Curve\nBlock: ", names(res.predict)[j], ", comp: ",i, sep="")
          
          plot.temp = ifelse(i%in%roc.comp && names(res.predict)[j]%in%block.temp, plot, FALSE)
          temp = my_statauc(data, plot = plot.temp, title = title, vec_col = vec_col)
          auc.mean[[names(res.predict)[j]]][[paste("comp",i,sep = "")]] = temp[[1]]
          graph[[names(res.predict)[j]]][[paste("comp",i,sep = "")]] = temp$graph
          
        }
      }
      out = c(auc.mean,graph=graph)
      
    }
    
    
  } else {
    
    vec_col <- "black"
    
    if (my_predict == "Averaged") { # uses the average predicted values over the blocks
      res.predict <- list.predict$AveragedPredict
    } else if (my_predict == "Weighted") { # uses the weighted average of the predicted values over the blocks (using the predict and weights outputs)
      res.predict <- list.predict$WeightedPredict
    } else {
      stop("Predict type not implemented.")
    }
    
    
    # for (i in 1:object$ncomp[1])
    # {
    # data$data=res.predict[,,i]
    # title=paste("ROC Curve\n", my_predict, " block contributions, comp: ",i, sep="")
    data$data=res.predict[,,roc.comp]
    title=paste("ROC Curve\n", my_predict, " block contributions, comp: ",roc.comp, sep="")
    
    plot.temp = plot
    temp = my_statauc(data, plot = plot.temp, title = title, vec_col = vec_col)
    auc.mean[[paste("comp",roc.comp,sep = "")]] = temp[[1]]
    graph[[paste("comp",roc.comp,sep = "")]] = temp$graph
    # }
    out = c(auc.mean,graph=graph) 
  }
  
  print(auc.mean)
  return(invisible(out))
  
}

my_statauc <- function(data = NULL, plot = FALSE, title = NULL, vec_col = NULL){
  res.predict = data$data; outcome = data$outcome
  
  
  ann_text = matrix(ncol=2,nrow=nlevels(outcome))
  colnames(ann_text) = c("AUC", "p-value")
  
  df = NULL; seauc = NULL; zauc = NULL
  for (i in 1 : nlevels(outcome)){
    tempout = outcome
    levels(tempout)[-i] = "Other(s)"
    tempout = factor(tempout, c(levels(outcome)[i], "Other(s)"))
    temp = my_roc.default(response = tempout,
                          predictor = as.matrix(res.predict[, i],ncol=1))
    
    temp_plot <- pROC::roc(
      response = tempout, #outcome, # two_class_example$truth,
      predictor = res.predict[, i],#res.predict[,1], # two_class_example$Class1,
      #levels =rev(levels(outcome)),# rev(levels(two_class_example$truth)),
      #direction = "<",
      smooth = TRUE
    )
    
    # print(coords(temp,x="local maximas" ,ret="threshold"))
    # print(coords(temp,x="best" ,ret="threshold",best.method = "youden"))
    # print(coords(temp,x="best" ,ret="threshold",best.method = "closest.topleft"))
    
    seauc = sqrt((0.25 + (sum(table(tempout)) -2) * 0.083333)/(prod(table(tempout))))
    zauc = (temp$auc-0.5)/seauc
    zauc[zauc < 0] = - zauc[zauc < 0]
    for (k in unique(temp_plot$specificities)){
      temp_plot$sensitivities[which(temp_plot$specificities == k)] = temp_plot$sensitivities[rev(which(temp_plot$specificities == k))]
    }
    if (nlevels(outcome) == 2){
      ann_text = matrix(ncol=2,nrow=1)
      colnames(ann_text) = c("AUC", "p-value")
      df = rbind(df, cbind(temp_plot$specificities, temp_plot$sensitivities, 
                           paste(paste(levels(outcome)[1], levels(outcome)[2], sep = " vs "),"\n",
                                 "AUC: ", signif(temp$auc*100, 3), "%")))
      #ann_text[i , 1] =
      ann_text[i , 1] = signif(temp$auc, 3)
      ann_text[i , 2] = signif((1 - pnorm(zauc, 0, 1))*2, 3)
      break
    } else {
      df = rbind(df, cbind(temp_plot$specificities, temp_plot$sensitivities, 
                           paste(levels(outcome)[i], "vs Other(s)\n",
                                 "AUC: ", signif(temp$auc*100, 3), "%")))
      #ann_text[i , 1] =
      ann_text[i , 1] = as.numeric(signif(temp$auc, 3))
      ann_text[i , 2] = as.numeric(signif((1 - pnorm(zauc, 0, 1))*2, 3))
    }
  }
  #define rownames for ann_text
  if (nlevels(outcome) == 2){
    rownames(ann_text) = paste(levels(outcome)[1], levels(outcome)[2], sep = " vs ")
  }else{
    rownames(ann_text) = paste(levels(outcome), "vs Other(s)")
  }
  
  df = data.frame(df, stringsAsFactors = FALSE)
  names(df) = c("Specificity", "Sensitivity", "Outcome")
  df$Specificity = 100 - 100*as.numeric(df$Specificity)
  df$Sensitivity = 100*as.numeric(df$Sensitivity)
  
  # print(head(outcome))
  # print(head(res.predict))
  
  # print(names(smooth_curv))
  # print(smooth_curv$response)
  # print(head(smooth_curv$specificities))
  # print(head(df$Specificity))
  # 
  # df$Specificity <- smooth_curv$specificities
  # df$Sensitivity <- smooth_curv$sensitivities
  # df$Outcome <- as.factor(smooth_curv$response)
  
  # print(smooth_curv)
  if(plot)
  {
    
    # p=NULL
    # Sensitivity = Specificity = Outcome = NULL #R check
    if(is.null(title))
      title = "ROC Curve"
    else
      title=title
    p = ggplot(df, aes(x=Specificity,
                       y=Sensitivity,
                       group = Outcome,
                       colour = Outcome)) +
      xlab("100 - Specificity (%)") +
      ylab("Sensitivity (%)") +
      geom_line(size = 0.75) +
      scale_x_continuous(breaks=seq(0, 100, by = 10)) +
      scale_y_continuous(breaks=seq(0, 100, by = 10)) +
      theme_bw() 
    
    p = p + geom_abline(intercept = 1) +
      theme(legend.key.size = unit(1, "cm"),
            plot.title = element_text(lineheight=.8, face="bold"),
            legend.title = element_text(size=14, face="bold")) +
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5),
                             panel.grid.minor = element_blank(),
                             panel.grid.major = element_blank()) + coord_fixed()
    
    if (is.null(vec_col)) {
      p = p + scale_color_brewer(palette="Spectral")  
    } else {
      p = p + scale_color_manual(values=vec_col) 
    }
    
    plot(p)
  } else {
    p=NULL
  }
  return(list(ann_text,graph=p))
  
}




my_statauc_multiple <- function(data = NULL, plot = FALSE, title = NULL, vec_col = NULL){
  
  list_res.predict = data$data; outcome = data$outcome
  print(names(list_res.predict))
  ann_text = matrix(ncol=2,nrow=nlevels(outcome))
  colnames(ann_text) = c("AUC", "p-value")
  
  list_df <- NULL
  
  for (dn in seq_along(list_res.predict)) {
    
    name_predict <- names(list_res.predict)[dn]
    res.predict <- list_res.predict[[dn]]
    
    df = NULL; seauc = NULL; zauc = NULL
    for (i in 1 : nlevels(outcome)){
      tempout = outcome
      levels(tempout)[-i] = "Other(s)"
      tempout = factor(tempout, c(levels(outcome)[i], "Other(s)"))
      temp = my_roc.default(response = tempout,
                            predictor = as.matrix(res.predict[, i],ncol=1))
      
      temp_plot <- pROC::roc(
        response = tempout, #outcome, # two_class_example$truth,
        predictor = res.predict[, i],#res.predict[,1], # two_class_example$Class1,
        #levels =rev(levels(outcome)),# rev(levels(two_class_example$truth)),
        #direction = "<",
        smooth = TRUE
      )
      
      # print(coords(temp,x="local maximas" ,ret="threshold"))
      # print(coords(temp,x="best" ,ret="threshold",best.method = "youden"))
      # print(coords(temp,x="best" ,ret="threshold",best.method = "closest.topleft"))
      
      seauc = sqrt((0.25 + (sum(table(tempout)) -2) * 0.083333)/(prod(table(tempout))))
      zauc = (temp$auc-0.5)/seauc
      zauc[zauc < 0] = - zauc[zauc < 0]
      for (k in unique(temp_plot$specificities)){
        temp_plot$sensitivities[which(temp_plot$specificities == k)] = temp_plot$sensitivities[rev(which(temp_plot$specificities == k))]
      }
      if (nlevels(outcome) == 2){
        ann_text = matrix(ncol=2,nrow=1)
        colnames(ann_text) = c("AUC", "p-value")
        df = rbind(df, cbind(temp_plot$specificities, temp_plot$sensitivities, 
                             paste(gsub("_", " ", name_predict), "\n", #paste(levels(outcome)[1], levels(outcome)[2], sep = " vs "),"\n",
                                   "AUC: ", signif(temp$auc*100, 3), "%")))
        #ann_text[i , 1] =
        ann_text[i , 1] = signif(temp$auc, 3)
        ann_text[i , 2] = signif((1 - pnorm(zauc, 0, 1))*2, 3)
        break
      } else {
        df = rbind(df, cbind(temp_plot$specificities, temp_plot$sensitivities, 
                             paste(gsub("_", " ", name_predict), "\n", # levels(outcome)[i], "vs Other(s)\n",
                                   "AUC: ", signif(temp$auc*100, 3), "%")))
        #ann_text[i , 1] =
        ann_text[i , 1] = as.numeric(signif(temp$auc, 3))
        ann_text[i , 2] = as.numeric(signif((1 - pnorm(zauc, 0, 1))*2, 3))
      }
    }
    #define rownames for ann_text
    if (nlevels(outcome) == 2){
      rownames(ann_text) = paste(levels(outcome)[1], levels(outcome)[2], sep = " vs ")
    }else{
      rownames(ann_text) = paste(levels(outcome), "vs Other(s)")
    }
    
    df = data.frame(df, stringsAsFactors = FALSE)
    names(df) = c("Specificity", "Sensitivity", "Outcome")
    df$Specificity = 100 - 100*as.numeric(df$Specificity)
    df$Sensitivity = 100*as.numeric(df$Sensitivity)
    
    list_df <- append(list_df, list(df))
    
    print(name_predict)
    print(unique(df$Outcome))
  }
  names(list_df) <- names(list_res.predict)
  
  if(plot)
  {
    # print(list_df[[1]]$Outcome) 
    # p=NULL
    # Sensitivity = Specificity = Outcome = NULL #R check
    if(is.null(title))
      title = "ROC Curve"
    else
      title=title
    p = ggplot(list_df[[1]], aes(x=Specificity,
                                 y=Sensitivity,
                                 group = Outcome,
                                 colour = "white")) +
      xlab("100 - Specificity (%)") +
      ylab("Sensitivity (%)") +
      geom_line(size = 0.75) +
      scale_x_continuous(breaks=seq(0, 100, by = 10)) +
      scale_y_continuous(breaks=seq(0, 100, by = 10)) +
      theme_bw()  #+ scale_color_manual(values=vec_col[1])
    
    for (ii in 1:length(list_df)) {
      # print(names(list_df)[ii])
      if (is.null(vec_col)) {
        p = p + geom_line(list_df[[ii]], mapping = aes(x=Specificity,
                                                       y=Sensitivity,
                                                       group = Outcome,
                                                       colour = Outcome), size = 0.75) 
        
      } else {
        p = p + geom_line(list_df[[ii]], mapping = aes(x=Specificity,
                                                       y=Sensitivity,
                                                       group = Outcome,
                                                       colour = Outcome), size = 0.75, col = vec_col[ii])
      }
    }
    if (is.null(vec_col)) {
      p = p + scale_color_brewer(palette = "Spectral") #scale_color_brewer(palette="Spectral")
    } else {
      # print(vec_col)
      p = p + scale_colour_manual(values=vec_col, labels = unique(as.vector(sapply(list_df, "[[", "Outcome"))))
    }
    
    p = p + geom_abline(intercept = 1) +
      theme(legend.key.size = unit(1, "cm"),
            plot.title = element_text(lineheight=.8, face="bold"),
            legend.title = element_text(size=14, face="bold")) +
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5),
                             panel.grid.minor = element_blank(),
                             panel.grid.major = element_blank()) + coord_fixed() 
    # p <- p + scale_fill_discrete(breaks = names(list_df))
    # print(names(list_df))
    
    
    
    plot(p)
  } else {
    p=NULL
  }
  return(list(ann_text,graph=p))
  
}




my_roc.utils.perfs <- function(threshold, controls, cases, direction) {
  if (direction == '>') {
    tp <- sum(cases <= threshold)
    tn <- sum(controls > threshold)
  }
  else if (direction == '<') {
    tp <- sum(cases >= threshold)
    tn <- sum(controls < threshold)
  }
  # return(c(sp, se))
  return(c(sp=tn/length(controls), se=tp/length(cases)))
}

my_roc.default <- function(response, predictor,
                           
                           auc=TRUE,
                           
                           levels=base::levels(as.factor(response))
                           
                           
) {
  
  
  # Response / Predictor
  original.predictor <- predictor # store a copy of the original predictor (before converting ordered to numeric and removing NA)
  original.response <- response # store a copy of the original predictor (before converting ordered to numeric)
  
  # remove NAs if requested
  nas <- is.na(response) | is.na(predictor)
  if (any(nas)) {
    na.action <- grep(TRUE, nas)
    class(na.action) <- "omit"
    response <- response[!nas]
    attr(response, "na.action") <- na.action
    predictor <- predictor[!nas]
    attr(predictor, "na.action") <- na.action
  }
  
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(levels[1])]]
  cases <- splitted[[as.character(levels[2])]]
  
  # Remove patients not in levels
  patients.in.levels <- response %in% levels
  if (!all(patients.in.levels)) {
    response <- response[patients.in.levels]
    predictor <- predictor[patients.in.levels]
  }
  
  # update 13/01/17: first level as control to force directionality:
  # > : if the predictor values for the control group are higher than the values of
  #the case group (controls > t >= cases)
  
  #if (median(controls) <= median(cases))
  #direction <- "<"
  #else if (median(controls) > median(cases))
  direction <- ">"
  
  
  # create the roc object
  roc <- list()
  class(roc) <- "roc"
  
  # compute SE / SP
  thresholds <-((c(-Inf, sort(unique(c(controls, cases)))) + c(sort(unique(c(controls, cases))), +Inf))/2)
  perf.matrix <- sapply(thresholds, my_roc.utils.perfs, controls=controls, cases=cases, direction=direction)
  perfs <- list(se=perf.matrix[2,], sp=perf.matrix[1,])
  
  se <- perfs$se
  sp <- perfs$sp
  
  
  # store the computations in the roc object
  roc$sensitivities <- se
  roc$specificities <- sp
  roc$thresholds <- thresholds
  roc <- my_sort.roc(roc)
  roc$direction <- direction
  roc$cases <- cases
  roc$controls <- controls
  
  # compute AUC
  if (auc)
    roc$auc <- my_auc_roc(roc)
  
  roc$call <- match.call()
  roc$original.predictor <- original.predictor
  roc$original.response <- original.response
  roc$predictor <- predictor
  roc$response <- response
  roc$levels <- levels
  
  
  # return roc
  return(roc)
}

my_sort.roc <- function(roc) {
  if (is.unsorted(roc$specificities)) {
    roc$sensitivities <- rev(roc$sensitivities)
    roc$specificities <- rev(roc$specificities)
    roc$thresholds <- rev(roc$thresholds)
  }
  return(roc)
}

my_auc_roc <- function(roc,
                       # Partial auc definition
                       partial.auc=FALSE, # false (consider total area) or numeric length 2: boundaries of the AUC to consider, between 0 and 1, or 0 and 100 if percent is TRUE
                       partial.auc.focus=c("specificity", "sensitivity"), # if partial.auc is not FALSE: do the boundaries
                       partial.auc.correct=FALSE,
                       allow.invalid.partial.auc.correct = FALSE,
                       ... # unused required to allow roc passing arguments to plot or ci.
) {
  if (!identical(partial.auc, FALSE)) {
    partial.auc.focus <- match.arg(partial.auc.focus)
  }
  
  percent <- FALSE
  
  # Validate partial.auc
  if (! identical(partial.auc, FALSE) & !(is.numeric(partial.auc) && length(partial.auc)==2))
    stop("partial.auc must be either FALSE or a numeric vector of length 2")
  
  # Ensure partial.auc is sorted with partial.auc[1] >= partial.auc[2]
  partial.auc <- sort(partial.auc, decreasing=TRUE)
  # Get and sort the sensitivities and specificities
  roc <- my_sort.roc(roc)
  se <- roc$se
  sp <- roc$sp
  
  # Full area if partial.auc is FALSE
  if (identical(partial.auc, FALSE)) {
    if (is(roc, "smooth.roc") && ! is.null(roc$smoothing.args) && roc$smoothing.args$method == "binormal") {
      coefs <- coefficients(roc$model)
      auc <- unname(pnorm(coefs[1] / sqrt(1+coefs[2]^2)) * ifelse(percent, 100^2, 1))
    }
    else {
      diffs.x <- sp[-1] - sp[-length(sp)]
      means.vert <- (se[-1] + se[-length(se)])/2
      auc <- sum(means.vert * diffs.x)
    }
  }
  # Partial area
  else {
    if (partial.auc.focus == "sensitivity") {
      # if we focus on SE, just swap and invert x and y and the computations for SP will work
      x <- rev(se)
      y <- rev(sp)
    }
    else {
      x <- sp
      y <- se
    }
    
    # find the SEs and SPs in the interval
    x.inc <- x[x <= partial.auc[1] & x >= partial.auc[2]]
    y.inc <- y[x <= partial.auc[1] & x >= partial.auc[2]]
    # compute the AUC strictly in the interval
    diffs.x <- x.inc[-1] - x.inc[-length(x.inc)]
    means.vert <- (y.inc[-1] + y.inc[-length(y.inc)])/2
    auc <- sum(means.vert * diffs.x)
    # add the borders:
    if (length(x.inc) == 0) { # special case: the whole AUC is between 2 se/sp points. Need to interpolate from both
      diff.horiz <- partial.auc[1] - partial.auc[2]
      # determine indices
      idx.hi <- match(FALSE, x < partial.auc[1])
      idx.lo <- idx.hi - 1
      # proportions
      proportion.hi <- (x[idx.hi] - partial.auc[1]) / (x[idx.hi] - x[idx.lo])
      proportion.lo <- (partial.auc[2] - x[idx.lo]) / (x[idx.hi] - x[idx.lo])
      # interpolated y's
      y.hi <- y[idx.hi] + proportion.hi * (y[idx.lo] - y[idx.hi])
      y.lo <- y[idx.lo] - proportion.lo * (y[idx.lo] - y[idx.hi])
      # compute AUC
      mean.vert <- (y.hi + y.lo)/2
      auc <- mean.vert*diff.horiz
    }
    else { # if the upper limit is not exactly present in SPs, interpolate
      if (!(partial.auc[1] %in% x.inc)) {
        # find the limit indices
        idx.out <- match(FALSE, x < partial.auc[1])
        idx.in <- idx.out - 1
        # interpolate y
        proportion <- (partial.auc[1] - x[idx.out]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.out] + proportion * (y[idx.in] - y[idx.out])
        # add to AUC
        auc <- auc + (partial.auc[1] - x[idx.in]) * (y[idx.in] + y.interpolated)/2
      }
      if (!(partial.auc[2] %in% x.inc)) { # if the lower limit is not exactly present in SPs, interpolate
        # find the limit indices in and out
        #idx.out <- length(x) - match(TRUE, rev(x) < partial.auc[2]) + 1
        idx.out <- match(TRUE, x > partial.auc[2]) - 1
        idx.in <- idx.out + 1
        # interpolate y
        proportion <- (x[idx.in] - partial.auc[2]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.in] + proportion * (y[idx.out] - y[idx.in])
        # add to AUC
        auc <- auc + (x[idx.in] - partial.auc[2]) * (y[idx.in] + y.interpolated)/2
      }
    }
  }
  
  # In percent, we have 100*100 = 10,000 as maximum area, so we need to divide by a factor 100
  if (percent)
    auc <- auc/100
  
  # Correction according to McClish DC, 1989
  if (all(!identical(partial.auc, FALSE), partial.auc.correct)) { # only for pAUC
    min <- roc.utils.min.partial.auc(partial.auc, percent)
    max <- roc.utils.max.partial.auc(partial.auc, percent)
    # The correction is defined only when auc >= min
    if (!allow.invalid.partial.auc.correct && auc < min) {
      warning("Partial AUC correction not defined for ROC curves below the diagonal.")
      auc <- NA
    }
    else if (percent) {
      auc <- (100+((auc-min)*100/(max-min)))/2 # McClish formula adapted for %
    }
    else {
      auc <- (1+((auc-min)/(max-min)))/2 # original formula by McClish
    }
  }
  # Prepare the AUC to return with attributes
  auc <- as.vector(auc) # remove potential pre-existing attributes
  attr(auc, "partial.auc") <- partial.auc
  attr(auc, "percent") <- percent
  attr(auc, "roc") <- roc
  if (!identical(partial.auc, FALSE)) {
    attr(auc, "partial.auc.focus") <- partial.auc.focus
    attr(auc, "partial.auc.correct") <- partial.auc.correct
  }
  class(auc) <- c("auc", class(auc))
  return(auc)
}

