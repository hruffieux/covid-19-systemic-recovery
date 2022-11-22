rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

main_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/scripts/")
data_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/data/preprocessed_data/")
out_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/output/")

setwd(main_dir)

require(nlme) # directly provides p-values, unlike lme4
require(splines)
require(dplyr)

source("fun_utils.R")

my_seed <- 123
set.seed(my_seed)

load(file.path(data_dir, "metabo_data.RData"))
load(file.path(data_dir, "cell_subset_data.RData"))

selected_severity_groups <- c("B", "C", "D", "E") 

trunc_days <- 50

score_types <- "CRP"
nb_gr <- 3
group_type <- "clusters_mclust"

bool_save <- T

all_data <- c("NMR_LP_abundances",
              "MS",
              "cell_types",
              "cytokines",
              "glyc",
              "log_ratios")

for (data_id in seq_along(all_data)) {
  
  data_name <- all_data[data_id]
  
  df_comb <- list(df_nmr_lp_comb, 
                  df_ms_comb, 
                  df_ct_comb,
                  df_cytokine_comb,
                  df_glyc_comb, 
                  df_all_ratios_comb)[[data_id]]
  
  df <- list(df_nmr_lp, 
             df_ms,
             df_ct,
             df_cytokine,
             df_glyc,
             df_all_ratios)[[data_id]]
  
  resp_var <- names(list(df_nmr_lp, 
                         df_ms, 
                         df_ct, 
                         df_cytokine, 
                         df_glyc,
                         df_all_ratios)[[data_id]])
  
  if (data_name == "log_ratios") {
    for (resp in resp_var) {
      df_comb[, resp] <- log2(df_comb[, resp] + 1)
      df[, resp] <- log2(df[, resp] + 1)
    }
  }  
  
  bool_covariates <- FALSE 
  
  bool_log10_display <- FALSE
  
  bool_disp_pval <- TRUE  # to display the baseline and interaction p-values on the plots
  bool_load_pval <- FALSE # to load the adjusted p-values and display them on 
  # the plots (the code need first to be run with bool_save_pval)
  bool_save_pval <- FALSE # to compute the adjusted p-values and save them in 
  # order to then be able to run with bool_load_pval <- T
  
  
  if (grepl("clust", group_type)) {
    
    days_thres <- 50
    load(file.path(out_dir, "all_scores.RData"))  # new, includes with mclust
    
    df_comb <- left_join(df_comb, subset(df_scores, 
                                         select = unique(c("sample_id_v2", 
                                                           "subject_id", 
                                                           "days_from_sx_or_swab_imputed", 
                                                           "severity", 
                                                           group_type))))
    df_comb <- df_comb[!is.na(df_comb[, group_type]) | df_comb$severity == "HC", ]
    
    
    df_comb[,group_type] <- as.factor(df_comb[,group_type])
    
    
    if (nb_gr == 3) {
      vec_col_gr <- vec_col_gr_3
    } else if (nb_gr == 2) {
      vec_col_gr <- vec_col_gr_2
    }
    
  }
  
} else {
  
  df_comb <- df_comb[!is.na(df_comb[, group_type]) | df_comb$severity == "HC", ]
  vec_col_gr <- vec_col[setdiff(selected_severity_groups, "HC")]
  
}   


# Split controls and cases
# ------------------------
#

# Cases:
#
df_comb_cases <- df_comb[df_comb$severity != "HC" &
                           df_comb$days_from_sx_or_swab_imputed < trunc_days,, 
                         drop = FALSE]

df_comb_cases <- choose_severity_groups(df_comb_cases, selected_severity_groups)

table(df_comb_cases$severity)
df_comb_cases$severity <- droplevels(df_comb_cases$severity) # drop HC level


# Controls:
#
df_comb_controls <- df_comb[df_comb$severity == "HC",, drop = FALSE]
nrow(df_comb_controls)

if(bool_save){
  save_path <- paste0(out_dir, "longitudinal", mess_version, "_clin_v", 
                      version_info, "_", 
                      paste0(selected_severity_groups, collapse = "-"),
                      "_data_name_", data_name, "_group_", group_type,
                      "_", trunc_days, "_days", 
                      ifelse(bool_covariates, "_age_gender",  ""), "/")
  dir.create(save_path, showWarnings = TRUE)
} else {
  save_path <- NULL
}


# For dot labels
#
df_comb_cases$status <- df_comb_cases$hospital_outcome
df_comb_cases$status[df_comb_cases$hospital_outcome %in% c("hospital", 
                                                           "other_facility", 
                                                           "other_hospital")] <- "hospital or facility"
df_comb_cases$status <- as.factor(df_comb_cases$status)
label_var <- "status"

# With confidence bands
bool_cb <- TRUE

if (bool_load_pval) { # pvalues and FDR
  
  ff <- paste0(save_path, "/pvalues.RData")
  if (file.exists(ff)) {
    load(ff)
  } else {
    stop("No existing p-value file.")
  }
  
} else {
  
  vec_pval_base <- vec_pval_int <- NULL
  vec_fdr_base <- vec_fdr_int <- NULL
  
}

df_comb_cases$group_type <- df_comb_cases[,group_type]

excluded_var <- NULL

for (ii in seq_along(resp_var)) {
  
  dep_var <- resp_var[ii]
  
  dep_var_controls <- df_comb_controls[, dep_var]
  
  df_comb_cases_dep <- df_comb_cases[!is.na(df_comb_cases[, dep_var]),, drop = FALSE]
  dep_var_cases <- df_comb_cases_dep[, dep_var]
  
  excluded_var <- c(excluded_var, tryCatch({
    
    if (bool_covariates) {
      out_lmm <- nlme::lme(dep_var_cases ~ bs(days_from_sx_or_swab_imputed, degree = 2) * group_type + age + gender,
                           random = ~ 1 | subject_id,
                           data = df_comb_cases_dep,
                           na.action = na.exclude,
                           method = "ML")
    } else {
      out_lmm <- nlme::lme(dep_var_cases ~ bs(days_from_sx_or_swab_imputed, degree = 2) * group_type,
                           random = ~ 1 | subject_id,
                           data = df_comb_cases_dep,
                           na.action = na.exclude,
                           method = "ML")
      
    }
    
    NULL
    
  }, error = function(e) {
    print(paste0("No fit for ", dep_var, " as overparametrised (singular convergence)."))
    return(dep_var)
  }))
  
}
resp_var <- resp_var[!(resp_var %in% excluded_var)]


for (ii in seq_along(resp_var)) {
  
  dep_var <- resp_var[ii]
  
  print(paste0("Dependent variable: ", dep_var))
  
  dep_var_controls <- df_comb_controls[, dep_var]
  
  df_comb_cases_dep <- df_comb_cases[!is.na(df_comb_cases[, dep_var]),, drop = FALSE]
  dep_var_cases <- df_comb_cases_dep[, dep_var]
  
  
  if (bool_covariates) {
    out_lmm <- nlme::lme(dep_var_cases ~ bs(days_from_sx_or_swab_imputed, degree = 2) * group_type + age + gender,
                         random = ~ 1 | subject_id,
                         data = df_comb_cases_dep,
                         na.action = na.exclude,
                         method = "ML")
  } else {
    out_lmm <- nlme::lme(dep_var_cases ~ bs(days_from_sx_or_swab_imputed, degree = 2) * group_type,
                         random = ~ 1 | subject_id,
                         data = df_comb_cases_dep,
                         na.action = na.exclude,
                         method = "ML")
  }
  
  
  if (bool_load_pval) {
    
    pval_base <- vec_pval_base[dep_var]
    fdr_base <- vec_fdr_base[dep_var]
    
    pval_int <- vec_pval_int[dep_var]
    fdr_int <- vec_fdr_int[dep_var]
    
  } else {
    
    if (bool_covariates) {
      out_lmm_f0_base <- nlme::lme(dep_var_cases ~ bs(days_from_sx_or_swab_imputed, degree = 2) + age + gender +
                                     bs(days_from_sx_or_swab_imputed, degree = 2):group_type,
                                   random = ~ 1 | subject_id,
                                   data = df_comb_cases_dep,
                                   na.action = na.exclude,
                                   method = "ML")
    } else {
      out_lmm_f0_base <- nlme::lme(dep_var_cases ~ bs(days_from_sx_or_swab_imputed, degree = 2) +
                                     bs(days_from_sx_or_swab_imputed, degree = 2):group_type,
                                   random = ~ 1 | subject_id,
                                   data = df_comb_cases_dep,
                                   na.action = na.exclude,
                                   method = "ML")
    }
    
    # likelihood ratio test for baseline term
    #
    pval_base <- anova(out_lmm, out_lmm_f0_base)$`p-value`[2] 
    fdr_base <- NULL
    
    if (bool_covariates) {
      out_lmm_f0_int <- nlme::lme(dep_var_cases ~ bs(days_from_sx_or_swab_imputed, degree = 2) + age + gender + group_type,
                                  random = ~ 1 | subject_id,
                                  data = df_comb_cases_dep,
                                  na.action = na.exclude,
                                  method = "ML")
    } else {
      out_lmm_f0_int <- nlme::lme(dep_var_cases ~ bs(days_from_sx_or_swab_imputed, degree = 2) + group_type,
                                  random = ~ 1 | subject_id,
                                  data = df_comb_cases_dep,
                                  na.action = na.exclude,
                                  method = "ML")
    }
    
    
    # likelihood ratio test for interaction term
    #
    pval_int <- anova(out_lmm, out_lmm_f0_int)$`p-value`[2] 
    fdr_int <- NULL
  }
  
  
  pval_base_label <- add_signif_label(pval_base)
  pval_int_label <- add_signif_label(pval_int)
  
  if (!is.null(fdr_base)) {
    fdr_base_label <- add_signif_label(fdr_base)
  }
  
  if (!is.null(fdr_int)) {
    fdr_int_label <- add_signif_label(fdr_int)
  }
  
  pred <- predict(out_lmm, level = 0) 
  
  if(bool_cb) { # add confidence bands
    des <- model.matrix(formula(out_lmm)[-2], df_comb_cases_dep)
    predvar <- diag(des %*% vcov(out_lmm) %*% t(des))
    df_comb_cases_dep$lower <- with(df_comb_cases_dep, pred - 1.96 * sqrt(predvar))
    df_comb_cases_dep$upper <- with(df_comb_cases_dep, pred + 1.96 * sqrt(predvar))
  }
  
  require(ggplot2)
  
  ymin <- min(dep_var_cases, na.rm = TRUE)
  ymax <- max(dep_var_cases, na.rm = TRUE)
  
  
  p <- ggplot(data = df_comb_cases_dep,
              aes(x = days_from_sx_or_swab_imputed,
                  y = dep_var_cases,
                  colour = group_type))
  
  if (any(!is.na(dep_var_controls))) {
    
    # HC stripe
    #
    lo <- quantile(dep_var_controls, probs = 0.25, na.rm = TRUE)
    up <- quantile(dep_var_controls, probs = 0.75, na.rm = TRUE)
    if (lo == up) { # for the case the upper quantile = lower quantile, so we see the line
      lo <- lo * 0.8
      up <- up * 1.2
    }
    
    p <- p + annotate("rect",
                      ymin = lo,
                      ymax = up,
                      xmin = -Inf, xmax = Inf, fill="black", # "pink",
                      alpha = 0.075)
    
  }
  
  
  # Confidence bands
  #
  if (bool_cb) {
    p <- p + geom_ribbon(data = df_comb_cases_dep,
                         aes(y = NULL,
                             ymin = lower,
                             ymax = upper,
                             color = NULL,
                             fill = group_type),
                         alpha = 0.1)
    p <- p + scale_fill_manual(values=vec_col_gr)
  }
  
  # Sample points
  #
  p <- p + geom_jitter(width = 0,
                       height = 0,
                       size = 1.5,
                       alpha=0.5)
  
  # Labels for points
  #
  if (!is.null(label_var)) {
    
    df_comb_cases_dep[, label_var] <- as.factor(df_comb_cases_dep[, label_var])
    
    p <- p + geom_point(aes_string(shape = label_var),
                        color="black", alpha=0.5, size = 1.5)
    p <- p + scale_shape_manual(values = c(13, 3, 1:20)[seq_along(levels(df_comb_cases_dep[, label_var]))],
                                na.translate = FALSE)
  }
  
  # Added fitted value lines
  #
  p <- p + geom_line(data = df_comb_cases_dep, aes(y = pred), size = 0.9)
  
  # Plot parameters
  #
  
  p <- p + theme_light()
  p <- p + ggtitle(gsub("_", " ", dep_var))
  p <- p + ylab(paste0("log2(.+1) ", gsub("_", " ", dep_var)))
  p <- p + theme(plot.title = element_text(size = 15, face = "bold"))
  p <- p + xlab("Days post symptom onset")
  p <- p + scale_color_manual(values = vec_col_gr)
  
  if (bool_log10_display) {
    
    rymin <- round_power_down(ifelse(bool_cb, min(df_comb_cases_dep$lower, ymin), ymin))
    rymax <- round_power_up(ifelse(bool_cb, max(df_comb_cases_dep$upper, ymax), ymax))
    
    p <- p + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                                breaks = c(0, 0.1, 1, 10, 100, 1000, 10000, 100000),
                                labels = c(0, 0.1, 1, 10, 100, 1000, 10000, 100000),
                                limits = c(rymin, ifelse((rymin < 1 & rymax > 1) |
                                                           (rymin < 10 & rymax > 10) |
                                                           (rymin < 100 & rymax > 100)|
                                                           (rymin < 1000 & rymax > 1000), ymax, rymax))
    )
    
  }
  
  if (bool_disp_pval) {
    if (is.null(fdr_base) & is.null(fdr_int)) {
      p <- p + labs(subtitle = paste0("LR test: baseline p = ",
                                      ifelse(pval_base < 0.0001, 
                                             format(pval_base, scientific = TRUE, digits = 2), 
                                             format(pval_base, digits = 2)),
                                      ifelse(pval_base_label != "",  paste0(" (", pval_base_label, ")"), ""),
                                      ", interaction p = ",
                                      ifelse(pval_int < 0.0001, 
                                             format(pval_int, scientific = TRUE, digits = 2), 
                                             format(pval_int, digits = 2)),
                                      ifelse(pval_int_label != "",  paste0(" (", pval_int_label, ")"), "")
      ))
    } else if (!is.null(fdr_base) & !is.null(fdr_int)) {
      p <- p + labs(subtitle = paste0("baseline p = ",
                                      ifelse(pval_base < 0.0001, 
                                             format(pval_base, scientific = TRUE, digits = 2), 
                                             format(pval_base, digits = 2)),
                                      ", adj p = ",
                                      ifelse(fdr_base < 0.0001, 
                                             format(fdr_base, scientific = TRUE, digits = 2), 
                                             format(fdr_base, digits = 2)),
                                      ifelse(fdr_base_label != "",  
                                             paste0(" (", fdr_base_label, ")"), ""),
                                      ", interaction p = ",
                                      ifelse(pval_int < 0.0001, 
                                             format(pval_int, scientific = TRUE, digits = 2), 
                                             format(pval_int, digits = 2)),
                                      ", adj p = ",
                                      ifelse(fdr_int < 0.0001, 
                                             format(fdr_int, scientific = TRUE, digits = 2), 
                                             format(fdr_int, digits = 2)),
                                      ifelse(fdr_int_label != "",  
                                             paste0(" (", fdr_int_label, ")"), "")
      ))
    } else {
      stop("fdr_base and fdr_int must be both NULL or both non-NULL.")
    }
    p <- p + theme(plot.subtitle=element_text(size = 9, color = "black"))
  }
  if (!is.null(save_path)) {
    pdf(paste0(save_path, "/longitudinal_",
               gsub(",", "", gsub(" ", "_", gsub("/", "_", dep_var))),
               ifelse(bool_cb, "_with_cb", ""),
               ifelse(bool_log10_display, "_log10_ylab", ""), ".pdf"),
        width = 5.75, height = 4.5, paper='special')
  }
  print(p)
  if (!is.null(save_path)) {
    dev.off()
  }
  
  if (!bool_load_pval) {
    vec_pval_base <- c(vec_pval_base, pval_base)
    vec_pval_int <- c(vec_pval_int, pval_int)
  }
  
}

if (!bool_load_pval) {
  names(vec_pval_base) <- names(vec_pval_int) <- resp_var
}

if (bool_save_pval) {
  vec_fdr_base <- p.adjust(vec_pval_base, metho = "BH") # FDR
  names(vec_fdr_base) <- names(vec_pval_base)
  
  vec_fdr_int <- p.adjust(vec_pval_int, metho = "BH") # FDR
  names(vec_fdr_int) <- names(vec_pval_int)
  
  save(vec_fdr_base, vec_pval_base, vec_fdr_int, vec_pval_int,
       file = paste0(save_path, "/pvalues.RData"))
  df_pval <- data.frame(
    "p-values baseline" = vec_pval_base, 
    "adjusted p-values baseline" = vec_fdr_base,
    "p-values interaction" = vec_pval_int, 
    "adjusted p-values interaction" = vec_fdr_int, check.names = FALSE)
  df_pval <- df_pval[order(df_pval$`adjusted p-values interaction`),]
  df_pval <- df_pval[order(df_pval$`adjusted p-values baseline`),]
  rownames(df_pval) <- paste0("'", rownames(df_pval), "'")
  write.csv(df_pval,
            file = paste0(save_path, "/table_longitudinal_evidence.csv"),
            row.names = TRUE)
  
  if (!is.null(save_path)) {
    pdf(paste0(save_path, "/plot_minlog10_adjpval_baseline_interaction.pdf"),
        width = 5.5, height = 5.5, paper='special')
  }
  
  fdr_thres <- 0.05
  vec_signif_col <- rep("black", nrow(df_pval))
  vec_signif_col[df_pval$`adjusted p-values baseline` <= fdr_thres & 
                   df_pval$`adjusted p-values interaction` <= fdr_thres] <- "springgreen4"
  vec_signif_col[df_pval$`adjusted p-values baseline` <= fdr_thres & 
                   df_pval$`adjusted p-values interaction` > fdr_thres] <- adjustcolor( "springgreen4", alpha.f = 0.7)
  vec_signif_col[df_pval$`adjusted p-values baseline` > fdr_thres & 
                   df_pval$`adjusted p-values interaction` <= fdr_thres] <- adjustcolor( "springgreen4", alpha.f = 0.7)
  
  names(vec_signif_col) <- rownames(df_pval)
  par(pty="s")
  plot(-log10(df_pval$`adjusted p-values baseline`),
       -log10(df_pval$`adjusted p-values interaction`), pch = 20, col = vec_signif_col,
       xlab = "Baseline effect",
       ylab = "Interaction effect",
       main = "Association with recovery groups\n (-log10 adj. p-value)")
  
  if (any(vec_signif_col != "black")) {
    text(-log10(df_pval$`adjusted p-values baseline`[vec_signif_col != "black"]),
         -log10(df_pval$`adjusted p-values interaction`[vec_signif_col != "black"]),
         labels = names(vec_signif_col)[vec_signif_col != "black"],
         col = vec_signif_col[vec_signif_col != "black"],
         pos = 2, cex = 0.6)
  }
  abline(v = -log10(fdr_thres), col =  adjustcolor( "springgreen4", alpha.f = 0.7), lwd = 1, lty = 2)
  abline(h = -log10(fdr_thres), col = adjustcolor( "springgreen4", alpha.f = 0.7), lwd = 1, lty = 2)
  
  if (!is.null(save_path)) {
    dev.off()
  }
  
}
