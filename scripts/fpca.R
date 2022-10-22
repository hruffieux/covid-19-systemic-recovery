rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

main_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/scripts/")
data_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/data/preprocessed_data/")
out_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/output/")

setwd(main_dir)

require(mfaces)    

source("fun_utils.R")

my_seed <- 123
set.seed(my_seed)

load(file.path(data_dir, "metabo_data.RData"))
load(file.path(data_dir, "cell_subset_data.RData"))

data_id <- 1

data_name <- c("CRP", 
               "cytokines",
               "glyc",
               "KP_metabolites", 
               "lipoproteins",
               "lymphocytes")[data_id]

df_comb <- list(df_info,
                df_cytokine_comb,
                df_glyc_comb, 
                df_ms_comb, 
                df_nmr_lp_comb, 
                df_flow_main_comb)[[data_id]]

all_var <- names(list(df_info,
                      df_cytokine, 
                      df_glyc,
                      df_ms,
                      df_nmr_lp, 
                      df_flow_main)[[data_id]])


bool_save <- T

days_thres <- 50
tnew <- 0:(days_thres-1) # timepoints to evaluate the prediction at

min_timepoints <- 2 # keep individuals with at least two samples 

if (data_name == "CRP") {
  
  vec_var <- "log_CRP"
  
  vec_flip <- c(1, -1) # swap sign for eigenfunctions (for interpretation purposes only)
  var_disp_names <- "CRP"
  
} else if (data_name == "cytokines") {
  
  vec_var <- all_var
  
  vec_flip <- c(-1, -1) 
  var_disp_names <- NULL
  
} else if (data_name == "glyc") { 
  
  vec_var <- c("GlycA", "GlycB")
  
  vec_flip <- c(-1, -1) 
  var_disp_names <- NULL
  
} else if (data_name == "KP_metabolites") {
  
  vec_var <- c("3-hydroxykynurenine", 
               "Kynurenine", 
               "Quinolinic acid",
               "Tryptophan")
  
  vec_flip <- c(-1, 1)
  var_disp_names <- NULL
  
} else if (data_name == "lipoproteins") {
  
  vec_var <- c("HDA1", "HDA2", "VLAB")
  
  vec_flip <- c(1, 1)
  var_disp_names <- NULL
  
} else if (data_name == "lymphocytes") {
  
  df_comb <- inner_join(df_comb, df_flow_comb)
  
  list_vec_var <- list(c("CD4 Naive T", 
                         "CD4 EMRA T", 
                         "CD4 Non-naive HLA-DR+CD38+ T",
                         "CD8 Naive T", 
                         "CD8 EMRA T", 
                         "CD8 Non-naive HLA-DR+CD38+ T",
                         "CD19+", 
                         "NK", 
                         "NKT", 
                         "MAIT", 
                         "Plasmablasts", 
                         "gd T"))
  
  vec_flip <- c(1, -1)
  var_disp_names <- NULL
  
} 

selected_severity_groups <- c("B", "C", "D", "E") # only symptomatic patients

df_comb_no_HC <- choose_severity_groups(df_comb, selected_severity_groups)
df_comb_no_HC <- df_comb_no_HC[df_comb_no_HC$days_from_sx_or_swab_imputed < days_thres, ]

vec_var <- sort(vec_var)


if (bool_save) {
  res_dir <- paste0(out_dir, "/FPCA_", 
                    paste0(selected_severity_groups, collapse = "-"),
                    "_days_thres_", days_thres, "_",
                    ifelse(length(vec_var)>10, 
                           paste0(substring(vec_var, 1, 7), collapse = "-"), 
                           paste0(vec_var, collapse = "-")), "/")
  dir.create(res_dir)
}

ff <- file.path(res_dir, "output.RData")

if (file.exists(ff)) {
  load(ff)
} else {
  data <- get_mface_data(df_comb_no_HC, vec_var, min_nb_timepoints = min_timepoints)
  
  if (length(vec_var) > 1) {
    
    # multivariate
    #
    fit <- mface.sparse(data, argvals.new = tnew, calculate.scores = TRUE)
    
  } else {
    
    # univariate:
    #
    fit <- face::face.sparse(data, argvals.new = tnew, calculate.scores = TRUE) 
   
  }
  
}

scores <- get_mface_scores(fit, data, tnew, df_comb_no_HC, vec_col = vec_col)
eigenfunctions <- fit$eigenfunctions
eigenvalues <- fit$eigenvalues

if ("EF2" %in% names(scores)) {
  
  print("/!\\ if flip not all 1, sign from object fit won't correspond to sign of eigenfunctions and scores.")
  
  eigenfunctions[,1:2] <- sweep(fit$eigenfunctions[,1:2], 2, vec_flip, "*")
  scores[,1:2] <- sweep(scores[,1:2], 2, vec_flip, "*")
  
  if (bool_save) {
    pdf(paste0(res_dir, "/scatterplot_scores_col_by_severity.pdf"),
        width = 5, height = 5, paper='special')
  }
  
  scatterplot_scores(fit, scores, color = scores$color_severity)
  
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/scatterplot_scores_col_by_severity_grey.pdf"),
        width = 5, height = 5, paper='special')
  }
  
  col_sev_grey <- scores$color_severity
  col_sev_grey[scores$severity == "A"] <- "grey90"
  col_sev_grey[scores$severity == "B"] <- "grey80"
  col_sev_grey[scores$severity == "C"] <- "grey65"
  col_sev_grey[scores$severity == "D"] <- "grey45"
  col_sev_grey[scores$severity == "E"] <- "grey5"
  scatterplot_scores(fit, scores, color = col_sev_grey)
  
  if (bool_save) {
    dev.off()
  } 
  
  nb_to_show <- 6
  for (show_extreme_subjects in c(TRUE, FALSE)) {
    if (bool_save) {
      pdf(paste0(res_dir, "/predicted_trajectories", 
                 ifelse(show_extreme_subjects, "_extreme_subj", ""), 
                 "_trunc.pdf"), 
          width = ifelse(data_name == "KP_metabolites", 3.25, 3)*length(vec_var), 
          height = ifelse(data_name == "KP_metabolites", 3, 2.75)*nb_to_show, paper='special')
    }
    
    subjects_to_show <- get_subjects_to_show(data, scores, 
                                             nb_to_show = nb_to_show, 
                                             show_extreme_subjects = show_extreme_subjects,
                                             min_timepoints = min_timepoints) 
    
    
    ylim_offset <- 3.2
    plot_subject_trajectories(fit, scores, tnew, data, 
                              subjects_to_show, 
                              vec_var, df_comb,
                              ylim_offset = ylim_offset,
                              var_disp_names = var_disp_names,
                              bool_trunc = TRUE)
    
    if (bool_save) {
      dev.off()
    } 
  }
  
  
  if (length(vec_var) > 1) {
    all_subjects <- sort(unique(unlist(lapply(data, function(ll) ll$subj))))
  } else {
    all_subjects <- sort(unique(data$subj))
  }
  
  for (subjects_to_show in split(all_subjects, ceiling(seq_along(all_subjects)/nb_to_show))) {
    
    if (bool_save) {
      pdf(paste0(res_dir, "/predicted_trajectories_subj_", 
                 paste0(subjects_to_show, collapse = "-"),
                 "_trunc.pdf"), 
          width = ifelse(data_name == "KP_metabolites", 3.25, 3)*length(vec_var), 
          height = ifelse(data_name == "KP_metabolites", 3, 2.75)*nb_to_show, paper='special')
    }
    
    plot_subject_trajectories(fit, scores, tnew, data, subjects_to_show, 
                              vec_var, df_comb, ylim_offset = ylim_offset,
                              var_disp_names = var_disp_names, bool_trunc = TRUE)
    
    if (bool_save) {
      dev.off()
    } 
  }
  
  
  if (bool_save) {
    pdf(paste0(res_dir, "/eigenfunctions_bw.pdf"), 
        width = ifelse(length(vec_var) == 1, 4.65, 2.75*length(vec_var)), 
        height = ifelse(length(vec_var) == 1, 5.1, 3), paper='special')
  }
  
  plot_eigenfunctions(eigenfunctions, eigenvalues, data, vec_var, days_thres)
  
  if (bool_save) {
    dev.off()
  } 
  
}
if (bool_save) {
  pdf(paste0(res_dir, "/variances.pdf"), 
      width = ifelse(length(vec_var) == 1, 5.5, 2.75*length(vec_var)), 
      height = ifelse(length(vec_var) == 1, 5, 2.75), paper='special')
} 

plot_variance(fit, data, tnew, vec_var)

if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/correlation.pdf"), 
      width = ifelse(length(vec_var) == 1, 6, 9), 
      height = ifelse(length(vec_var) == 1, 5.5, 8.5), paper='special')
} 
nb_spaces <- 2*length(vec_var) 
nb_ticks <- 4
plot_correlation(fit, data, tnew, vec_var, nb_spaces, nb_ticks)
if (bool_save) {
  dev.off()
}

if (bool_save) {
  save(data, fit, eigenfunctions, scores, tnew, vec_var, 
       file = file.path(res_dir, "output.RData"))
}


