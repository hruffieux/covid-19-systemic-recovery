rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

main_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/scripts/")
data_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/data/preprocessed_data/")
out_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/output/")

setwd(main_dir)

require(dplyr)
require(metafor)

source("fun_utils.R")

my_seed <- 123
set.seed(my_seed)

load(file.path(data_dir, "metabo_data.RData"))
load(file.path(data_dir, "cell_subset_data.RData"))

bool_save <- T

sample_type <- c("all", "early", "late")[2]

if (sample_type == "early") {
  days_thres <- 49 # from 0 to days_thres days (incl.) after symptom onset of first swab
  mess_samples <- paste0("_samples_before_", days_thres, "_days")
  bool_early <- TRUE
} else if (sample_type == "late") {
  days_thres <- 28 # from day days_thres to end of measurements.
  bool_early <- FALSE
  mess_samples <- paste0("_samples_after_", days_thres, "_days")
} else {
  mess_samples <- "" # all samples
}

nuis_covariates <- c("gender", "age", "ethnicity", "bmi", "log_CRP")[c(1,2)] 

fdr_thres <- 0.05

vec_test_covariates <- c("covid_status", "log_CRP")

selected_severity_groups <- c("HC", "A", "B", "C", "D", "E") 

vec_test_covariates <- filter_invalid_test_covariates(vec_test_covariates, 
                                                      selected_severity_groups)

all_data_types <- c("MS", "cell_types", "log_ratios", "glprot")

for(data_id in seq_along(all_data_types)) {
  
  data_name <- all_data_types[data_id]
  
  if (data_name == "glprot") {
    var_names <- c(names(df_nmr_lp), c("GlycA", "GlycB", "SPC"))
  } else {
    var_names <- names(list(df_ms, df_ct, df_all_ratios)[[data_id]]) 
  }
  
  df_comb <- list(df_ms_comb, df_ct_comb, df_all_ratios_comb, df_nmr_lp_comb)[[data_id]]
  
  if (bool_save) {
    res_dir <- paste0(out_dir, "DE_adj_", 
                      paste0(nuis_covariates, collapse = "-"),
                      "_severity_groups_", 
                      paste0(selected_severity_groups, collapse = "-"), 
                      mess_samples, "/")
    dir.create(res_dir)
  } else {
    res_dir <- NULL
  }
  
  df_comb <- choose_severity_groups(df_comb, selected_severity_groups)
  
  if (sample_type!="all") {
    list_fs <- get_early_or_late_samples(df_comb, vec_var = var_names, 
                                         bool_early = bool_early,
                                         trunc_days = days_thres, 
                                         single_sample_per_subject = FALSE,
                                         enforce_subj_w_samples_avail_after_trunc_days = FALSE)
    df_comb <- list_fs$df_comb
  }
  
  for(test_covariate in vec_test_covariates) {
    
    print(test_covariate)
    title <- set_analysis_title(test_covariate, selected_severity_groups)
    
    list_out <- run_mixed_models(df_comb, 
                                 test_covariate = test_covariate, 
                                 nuis_covariates = nuis_covariates, 
                                 var_names = var_names, 
                                 fdr_thres = fdr_thres)
    
    df_res <- list_out$df_res
    df_comb <- list_out$df_comb
    var_names <- list_out$var_names
    
    if (!is.null(df_res)) {
      
      if (test_covariate %in% c("log_CRP", "bmi")){
        
        lab_x <- "estimate"
        
      } else {
        
        lab_x <- "log2 fold change"
        
      }
      
      if (!is.null(res_dir)) {
        
        if (grepl("glprot", data_name)) {
          height <- 8 
          width <- 7.8
        } else {
          height <- 7
          width <- 8
        }
        pdf(paste0(res_dir, "/volcano_plot_", test_covariate, "_", data_name, ".pdf"),
            width = width, height = height, paper='special')
      }
      
     if (data_name == "glprot") {
       
        classes_glprot <- rbind(classes_all_nmr_lp, 
                                c("GlycA", "", "Glycoproteins-SPC", ""), 
                                c("GlycB", "", "Glycoproteins-SPC", ""), 
                                c("SPC", "", "Glycoproteins-SPC", ""))
        
        vec_classes <- classes_glprot$Compound[match(df_res$var_name, classes_glprot$Key)]
        
        palette <- c("#6388B4FF", "#FFAE34FF", "#EF6F6AFF", "#8CC2CAFF", 
                     "#55AD89FF", "#C3BC3FFF", "#BB7693FF", "#BAA094FF", 
                     "#A9B5AEFF", "#767676FF")
        
        col_classes <- palette[as.numeric(as.factor(vec_classes))] 
        names(col_classes) <- classes_glprot$Key[match(df_res$var_name, classes_glprot$Key)]
        col_by_classes <- col_classes
        names(col_by_classes) <- vec_classes
        stopifnot(all.equal(names(col_classes), df_res$var_name))
        legend <- col_by_classes
        
      } else {
        
        col_classes <- col_by_classes <- NULL
        legend <- NULL
        
      }
      
      stopifnot(all.equal(rownames(df_res), var_names))
      print(myEnhancedVolcano(df_res,
                              lab = gsub("_", " ", var_names),
                              subtitle = gsub("_", " ", data_name),
                              title = title,
                              x = "estimate",
                              y = "adj_pval",
                              legend = legend,
                              xlim = c(-1*max(abs(df_res$estimate), na.rm = TRUE), 
                                       1*max(abs(df_res$estimate), na.rm = TRUE)),
                              ylim = c(0, max(-log10(df_res$adj_pval))),
                              transcriptPointSize = 2.5,
                              transcriptLabSize = 4.75,
                              drawConnectors = TRUE,
                              widthConnectors = 0.2,
                              typeConnectors = "open",
                              endsConnectors = "last",
                              transcriptLabhjust = 0, 
                              transcriptLabvjust = 0.5, 
                              colCustom = col_by_classes,
                              gridlines.major = TRUE, 
                              gridlines.minor = FALSE,
                              xlab = lab_x,
                              colAlpha = 1,
                              ylab = bquote(~-Log[10] ~ "adj. p"),
                              FCcutoff = 0, 
                              pCutoff = fdr_thres,
                              pLabellingCutoff = fdr_thres,
                              FCLabellingCutoff = 0))
      
      if (!is.null(res_dir)) {
        dev.off()
      }
      
      if (test_covariate %in% c("log_CRP", "bmi")) {
        
        if (data_name == "glprot") {
          
          vec_cat <- c("Apolipoprotein", "Cholesterol", "Phospholipids", 
                       "Triglycerides", "Glycoproteins")
          
          for (cat in vec_cat) {
            
            df_res_sub <- df_res[df_res$var_name %in% classes_glprot$Key[grepl(cat, classes_glprot$Compound)],, drop = FALSE]
            
            if (!is.null(res_dir)) {
              pdf(paste0(res_dir, "/forest_plot_", test_covariate, "_", cat, "_", data_name, ".pdf"),
                  width = 6.8, height = 1.5 + 0.2*nrow(df_res_sub), paper='special')
            }
            
            vec_col <- rep("black", nrow(df_res_sub))
            vec_col[df_res_sub$sign == "Upregulated"] <- "red"
            vec_col[df_res_sub$sign == "Downregulated"] <- "blue"
            
            forest(df_res_sub$estimate, 
                            sei = df_res_sub$se,
                            slab = df_res_sub$var_name, 
                            main = paste0("Associations with ", 
                                          test_covariate, #ifelse(test_covariate == "log_CRP", "CRP levels", "severity"), 
                                          "\n ", gsub(" abundances", "", gsub("_", " ", data_name)), ", ",
                                          cat),
                            cex.main = 2,
                            col = vec_col)
            
            if (!is.null(res_dir)) {
              dev.off()
            }
            
          }
          
      
          sub_var <- c("GlycA", "GlycB", "SPC", 
                       sort(setdiff(trimws(gsub("\\w*[0-9]+\\w*\\s*", "", colnames(df_nmr_lp))), ""))) 
          
          df_res_sub <- df_res[match(sub_var, df_res$var_name),, drop = FALSE]
          
          if (!is.null(res_dir)) {
            pdf(paste0(res_dir, "/forest_plot_", test_covariate, "_subset_glyco_lipo_", data_name, ".pdf"),
                width = 6.8, height = 1.5 + 0.2*length(sub_var), paper='special')
          }
          
          vec_col <- rep("black", nrow(df_res_sub))
          vec_col[df_res_sub$sign == "Upregulated"] <- "red"
          vec_col[df_res_sub$sign == "Downregulated"] <- "blue"
          
          forest(df_res_sub$estimate, 
                          sei = df_res_sub$se,
                          slab = df_res_sub$var_name, 
                          main = paste0("Associations with ", 
                                        test_covariate, 
                                        "\n ", gsub(" abundances", "", gsub("_", " ", data_name)), ", ",
                                        "subset"),
                          cex.main = 2,
                          col = vec_col)
          
          if (!is.null(res_dir)) {
            dev.off()
          }
          
          
        } else {
          
          vec_col <- rep("black", nrow(df_res))
          vec_col[df_res$sign == "Upregulated"] <- "red"
          vec_col[df_res$sign == "Downregulated"] <- "blue"
          
          
          if (grepl("ratios", data_name)) {
            
            if (!is.null(res_dir)) {
              pdf(paste0(res_dir, "/forest_plot_", test_covariate, "_", data_name, ".pdf"),
                  width = 8, height = 1.5 + 0.2*nrow(df_res), paper='special')
            }
            
            forest(df_res$estimate, 
                            sei = df_res$se,
                            slab = gsub("_", " ", df_res$var_name), 
                            xlim = c(-max(abs(df_res$estimate)+df_res$se)*2.3,
                                     max(abs(df_res$estimate)+df_res$se)*2.3),
                            main = paste0("Associations with ",
                                          test_covariate,
                                          "\n ", gsub(" abundances", "", gsub("_", " ", data_name))),
                            cex.main = 2,
                            col = vec_col)
          } else {
            
            if (!is.null(res_dir)) {
              pdf(paste0(res_dir, "/forest_plot_", test_covariate, "_", data_name, ".pdf"),
                  width = 9.25, height = 8.25, paper='special')
            }
            
            
            forest(df_res$estimate, 
                            sei = df_res$se,
                            slab = df_res$var_name, 
                            main = paste0("Associations with ",
                                          test_covariate,
                                          "\n ", gsub(" abundances", "", gsub("_", " ", data_name))),
                            cex.main = 2,
                            col = vec_col)
          }
          
          if (!is.null(res_dir)) {
            dev.off()
          }
        }
      }
      
      
      list_col_annot <- list(severity = list_col$severity, gender = list_col$gender)
      
      split_rows <- df_res$sign
      names(split_rows) <- rownames(df_res)
      
      if (!is.null(res_dir)) {
        save(df_res, file = paste0(res_dir, test_covariate, "_", data_name, ".RData"))
      }
      
    }
    
  }
}

