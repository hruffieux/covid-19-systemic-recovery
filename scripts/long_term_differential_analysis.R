rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")
CORE_DIR_ICLOUD <- Sys.getenv("CORE_DIR_ICLOUD")

main_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/scripts/")
data_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/data/preprocessed_data/")
out_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/output/")

setwd(main_dir)

require(dplyr)

source("fun_utils.R")

my_seed <- 123
set.seed(my_seed)

bool_save <- FALSE

load(file.path(data_dir, "metabo_data.RData"))
load(file.path(data_dir, "cell_subset_data.RData"))

list_windows <- list(c(0, 21), # 3 weeks
                     c(21, 49), # 3-7
                     c(49, 84), # 7-12
                     c(84, 189), # 12-27
                     c(189, 365)) # 27 - 1 year 

n_windows <- length(list_windows)

data_id <- 1

data_name <- c("NMR_LP_abundances",
               "MS",
               "cell_types",
               "cytokines",
               "glyc",
               "log_ratios")[data_id]

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

all_var <- names(df)

group_names <- "clusters_mclust"

selected_severity_groups <- c("HC", "B", "C", "D", "E") 

trunc_days <- 50

load(file.path(out_dir, "all_scores.RData")) 

df_comb <- left_join(df_comb, subset(df_scores, 
                                     select = c("sample_id_v2", 
                                                "subject_id", 
                                                "days_from_sx_or_swab_imputed", 
                                                "severity",
                                                "clusters_mclust")))

df_comb <- choose_severity_groups(df_comb, selected_severity_groups)

classes <- c("i", "ii", "iii")


if (bool_save) {
  res_dir <- paste0(out_dir, "recovery_vs_HC_days_", 
                    paste0(unique(unlist(list_windows)), collapse="-"),
                    group_names, "/")
  dir.create(res_dir)
} else {
  res_dir <- NULL
}

list_res_window <- NULL

for (w_id in 1:n_windows) {
  
  print(paste0("==== ", w_id, " ===="))
  
  single_sample_per_subject <- F
  
  if (w_id == 1) {
    list_samples <- get_early_or_late_samples(df_comb, 
                                              all_var, 
                                              bool_early = TRUE, 
                                              drop_HC = FALSE,
                                              trunc_days = list_windows[[w_id]][2], 
                                              single_sample_per_subject = single_sample_per_subject) 
    
  } else {
    
    list_samples <- get_early_or_late_samples(df_comb, all_var, bool_early = TRUE, 
                                              drop_HC = FALSE,
                                              trunc_days = list_windows[[w_id]][2], 
                                              single_sample_per_subject = single_sample_per_subject) 
    df_comb_window <- list_samples$df_comb
    list_samples <- get_early_or_late_samples(df_comb_window, 
                                              all_var, 
                                              bool_early = FALSE, 
                                              drop_HC = FALSE,
                                              trunc_days = list_windows[[w_id]][1], 
                                              single_sample_per_subject = single_sample_per_subject) 
    
  } 
  
  df_comb_window <- list_samples$df_comb
  table(df_comb_window$severity)
  
  list_res <- lapply(group_names, function(group_type) {
    
    df_comb_window[df_comb_window$severity == "HC", group_type] <- "HC"
    df_comb_window[, group_type] <- factor(df_comb_window[, group_type], 
                                           levels = c("HC", classes))
    print(table(df_comb_window[, group_type]))
    print(table(df_comb_window[!duplicated(df_comb_window$subject_id), group_type]))
    
    mat_log_fc <- t(sapply(all_var, function(vv) {
      med_HC <- median(df_comb_window[df_comb_window[, group_type] == "HC", vv], na.rm = T)
      vec_log_fc <- sapply(classes, function(cc) {
        log2(median(df_comb_window[df_comb_window[, group_type] == cc, vv], na.rm = T) / med_HC)
      })
      vec_log_fc
    }))
    
    if (single_sample_per_subject) {
      list_lm <- lapply(all_var, function(vv) {
        out_lm <- lm(df_comb_window[, vv] ~ df_comb_window[, group_type])
        ss <- summary(out_lm)
        
        tval <- sapply(classes, function(cc) { 
          ss$coefficients[paste0(group_type, cc), 3]})
        pval <- sapply(classes, function(cc) { 
          ss$coefficients[paste0(group_type, cc), 4]})
        
        create_named_list(tval, pval)
      })
    } else {
      require(lmerTest)
      list_lm <- lapply(all_var, function(vv) {
        
        vv_new <- paste0("X", gsub("\\+", ".", 
                                   gsub(":", ".",  
                                        gsub("\\(", ".", 
                                             gsub("\\)", ".", 
                                                  gsub(" ", ".", 
                                                       gsub("-", ".", vv)))))))
        colnames(df_comb_window)[colnames(df_comb_window) == vv] <- vv_new
        formula_string <- paste0(vv_new, " ~ ", group_type, " + (1 | subject_id)")
        
        tryCatch({ 
          out_lm <- lmer(formula_string, # dep_var ~ covid_status + gender + age + (1|subject_id),
                         data=df_comb_window,
                         REML=F)
          ss <- summary(out_lm)
          
          
          tval <- sapply(classes, function(cc) { 
            ss$coefficients[paste0(group_type, cc), 4]})
          pval <- sapply(classes, function(cc) { 
            ss$coefficients[paste0(group_type, cc), 5]})
          
          create_named_list(tval, pval)
          
        }, error=function(e){
          tval <- pval <- rep(NA, nb_gr)
          create_named_list(tval, pval)
        })
        
        
        
        
      })
    }
    
    names(list_lm) <- all_var
    mat_tval <- t(sapply(list_lm, "[[", "tval"))
    mat_pval <- t(sapply(list_lm, "[[", "pval"))
    colnames(mat_log_fc) <- colnames(mat_tval) <- colnames(mat_pval) <- classes
    
    mat_adj_pval <- apply(mat_pval, 2, 
                          function(vec_pval) 
                            p.adjust(vec_pval, method = "BH"))
    mat_signif <- apply(mat_adj_pval, 2, 
                        function(vec_adj_pval) 
                          sapply(vec_adj_pval, 
                                 function(pp) add_signif_label(pp)))
    
    create_named_list(mat_log_fc, mat_tval, mat_pval, mat_adj_pval, 
                      mat_adj_pval, mat_signif)
  })
  names(list_res) <- group_names
  
  list_res_window <- append(list_res_window, list(list_res))
}
names(list_res_window) <- paste0("window_", 1:n_windows)



list_group <- lapply(group_names, function(group_name) 
  lapply(list_res_window, "[[", group_name)) 
names(list_group) <- group_names

if (n_windows == 4) {
  
  list_res <- lapply(list_group, function(ll) mapply(cbind, 
                                                     ll$window_1, 
                                                     ll$window_2, 
                                                     ll$window_3, 
                                                     ll$window_4, 
                                                     SIMPLIFY = FALSE))
  list_res <- lapply(list_res, 
                     function(ll) 
                       lapply(ll, function(mat) { 
                         colnames(mat) <- paste0(c(rep("window_1_", nb_gr), 
                                                   rep("window_2_", nb_gr), 
                                                   rep("window_3_", nb_gr), 
                                                   rep("window_4_", nb_gr)), 
                                                 classes)
                         mat}))
  
} else if (n_windows == 5) {
  
  list_res <- lapply(list_group, 
                     function(ll) mapply(cbind, 
                                         ll$window_1, 
                                         ll$window_2, 
                                         ll$window_3, 
                                         ll$window_4, 
                                         ll$window_5, 
                                         SIMPLIFY = FALSE))
  
  list_res <- lapply(list_res, 
                     function(ll) 
                       lapply(ll, function(mat) { 
                         colnames(mat) <- paste0(c(rep("window_1_", nb_gr), 
                                                   rep("window_2_", nb_gr), 
                                                   rep("window_3_", nb_gr), 
                                                   rep("window_4_", nb_gr), 
                                                   rep("window_5_", nb_gr)), 
                                                 classes)
                         mat})) 
  
} else if (n_windows == 6) {
  
  list_res <- lapply(list_group, 
                     function(ll) mapply(cbind, 
                                         ll$window_1, 
                                         ll$window_2, 
                                         ll$window_3, 
                                         ll$window_4, 
                                         ll$window_5, 
                                         ll$window_6, 
                                         SIMPLIFY=FALSE))
  
  list_res <- lapply(list_res, 
                     function(ll) lapply(ll, function(mat) { 
                       colnames(mat) <- paste0(c(rep("window_1_", nb_gr), 
                                                 rep("window_2_", nb_gr), 
                                                 rep("window_3_", nb_gr), 
                                                 rep("window_4_", nb_gr), 
                                                 rep("window_5_", nb_gr), 
                                                 rep("window_6_", nb_gr)), 
                                               classes)
                       mat}))
  
} else {
  
  stop("Number of windows not implemented.")
  
}

stopifnot(length(unique(lapply(list_res, rownames))) == 1) # check that same rownames

mat_log_fc <- do.call(cbind, lapply(list_res, "[[", "mat_log_fc"))
mat_tval <- do.call(cbind, lapply(list_res, "[[", "mat_tval"))
mat_pval <- do.call(cbind, lapply(list_res, "[[", "mat_pval"))
mat_adj_pval <- do.call(cbind, lapply(list_res, "[[", "mat_adj_pval"))
mat_signif <- do.call(cbind, lapply(list_res, "[[", "mat_signif"))


colnames(mat_log_fc) <- colnames(mat_tval) <- 
  colnames(mat_pval) <- colnames(mat_adj_pval) <- 
  colnames(mat_signif) <- 
  paste0(rep(names(list_res),
             each = n_windows*as.numeric(gsub("b", "", nb_gr))), 
         "_", colnames(mat_pval))


require(gplots)
require(RColorBrewer)

reord <- c(1, 4, 7, 10, 13, 2, 5, 8, 11, 14, 3, 6, 9, 12, 15)
col_sep <- n_windows*c(1, 2)

long_path <- paste0(out_dir, "longitudinal_",  
                    paste0(setdiff(selected_severity_groups, "HC"), collapse = "-"),
                    "_data_name_", data_name, "_group_", group_names,
                    "_", trunc_days, "_days", "/")

load(file.path(long_path, "pvalues.RData"))


sub_vec_fdr_base <- vec_fdr_base[match(rownames(mat_adj_pval), names(vec_fdr_base))]
sub_vec_fdr_int <- vec_fdr_int[match(rownames(mat_adj_pval), names(vec_fdr_int))]

stopifnot(all.equal(rownames(mat_adj_pval), names(sub_vec_fdr_base)))

mat_reord <- as.data.frame(cbind(mat_tval[, reord], 
                                 -log10(sub_vec_fdr_base), 
                                 -log10(sub_vec_fdr_int)))
mat_signif_reord <- cbind(mat_signif[, reord], 
                          sapply(sub_vec_fdr_base, add_signif_label),
                          sapply(sub_vec_fdr_int, add_signif_label))

names(mat_reord) <- c(colnames(mat_tval)[reord],  "Baseline", "Interaction")

row_min <- apply(mat_adj_pval, 1, min)
row_ord <- order(row_min)

mat_reord <- mat_reord[row_ord, ]
mat_reord <- mat_reord[order(mat_reord$Interaction, decreasing = T), ]
mat_reord <- mat_reord[order(mat_reord$Baseline, decreasing = T), ]

mat_signif_reord <- mat_signif_reord[match(rownames(mat_reord), 
                                           rownames(mat_signif_reord)), ]

pal <- colorRampPalette(c("blue", "white", "red"))(1000)
br <- seq(-max(abs(mat_reord), na.rm = T), max(abs(mat_reord), na.rm = T), 
          length.out=1001)

window_labels <- rep(sapply(list_windows, 
                            function(window) 
                              paste0("(", window[1], ", ", window[2], "]")),
                     as.numeric(gsub("b", "", nb_gr)))


if (bool_save) {
  pdf(paste0(res_dir, "heatmap_", data_name, ".pdf"),
      width = 12, height = ifelse(grepl("ratios", data_name), 7, 9.25), paper='special')
}
par(cex.main=1)

main_hp <- heatmap.2(as.matrix(mat_reord),
                     dendrogram='none', 
                     Rowv=TRUE, Colv=FALSE, 
                     trace='none', scale="none",
                     labCol = c(window_labels, "Baseline", "Interaction"),
                     labRow = gsub("_", " ", rownames(mat_reord)),
                     xlab = "",
                     ylab = "",
                     cexCol = 1.2,
                     cexRow = 1,
                     key=FALSE, keysize=1,
                     margins = c(10, 20),
                     key.title="",
                     sepcolor = "white",
                     colsep = c(col_sep, col_sep[length(col_sep)]+n_windows, 
                                col_sep[length(col_sep)]+n_windows),
                     sepwidth = 0.1,
                     cellnote = mat_signif_reord,
                     notecex = 1.25,
                     notecol="black",
                     main = ifelse(as.numeric(gsub("b", "", nb_gr)) == 3,
                                   ifelse(grepl("clust", group_names),
                                          c("\n\n\n \n\n\n \n\n\n \n\n\n\n\n            i                           ii                          iii                               "),
                                          c("\n\n\n \n\n\n \n\n\n \n\n\n\n\n Deteriorates          Stable           Recovers                               ")),
                                   ifelse(grepl("clust", group_names),
                                          c("-               +                "),
                                          c("Deteriorates  Recovers    "))),
                     key.xlab="", key.ylab="", col=pal, density.info="none",
                     breaks=br,na.color="grey90")

if (bool_save) {
  dev.off()
  save(mat_tval, mat_signif, 
       vec_pval_base, vec_fdr_base, vec_pval_int, vec_fdr_int, 
       file = paste0(res_dir, "output_", data_name, ".RData"))
}
