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
    
    formula_colnames <- paste0("X", gsub("\\+", ".", gsub(":", ".",  gsub("\\(", ".", gsub("\\)", ".", gsub(" ", ".", gsub("-", ".", var_names))))))) # colnames starting with a digit cause an error
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
    if (isTRUE(all.equal(sort(selected_severity_groups), sort(c("HC", "A", "B", "C", "D", "E"))))) {
      title <- "COVID+ vs HC"
    } else if (isTRUE(all.equal(sort(selected_severity_groups), sort(c("HC", "B", "C", "D", "E"))))) {
      title <- "COVID+ symptomatics vs HC"
    } else if (isTRUE(all.equal(sort(selected_severity_groups), sort(c("HC", "A", "B", "C", "D"))))) {
      title <- "Non severe vs HC"
    } else if (isTRUE(all.equal(sort(selected_severity_groups), sort(c("HC", "A"))))) {
      title <- "Asymptomatics vs HC" 
    } else if (isTRUE(all.equal(sort(selected_severity_groups), sort(c("HC", "E"))))) {
      title <- "Severe vs HC" 
    } else if (isTRUE(all.equal(sort(selected_severity_groups), sort(c("HC", "A", "B"))))) {
      title <- "Mild vs HC"
    } else if (isTRUE(all.equal(sort(selected_severity_groups), sort(c("HC", "C", "D"))))) {
      title <- "Moderate vs HC"
    } else {
      stop("Invalid comparison.")
    }
  } else if (test_covariate == "log_CRP") {
    title <- paste0("CRP (", 
                    paste0(intersect(c("A", "B", "C", "D", "E"), selected_severity_groups), 
                           collapse = ", "), ")")
  } else if (test_covariate %in% c("bmi", "log_bmi", vec_cytokines, vec_cplt, vec_ig, vec_glyc)) {
    title <- paste0(test_covariate, " (", 
                    paste0(intersect(c("A", "B", "C", "D", "E"), selected_severity_groups), 
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
                               title = "Volcano plot", subtitle = "Bioconductor package EnhancedVolcano", 
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


