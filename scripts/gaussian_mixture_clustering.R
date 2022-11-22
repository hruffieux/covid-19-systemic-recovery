rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")
CORE_DIR_ICLOUD <- Sys.getenv("CORE_DIR_ICLOUD")

main_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/scripts/")
data_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/data/preprocessed_data/")
out_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/output/")

setwd(main_dir)

require(mclust)

res_dir_fpca_crp <- file.path(out_dir, 
                              "../FPCA_clin_B-C-D-E_days_thres_50_log_CRP/") 
load(file.path(res_dir_fpca_crp, "output.RData"))

bool_save <- FALSE

if (bool_save) {
  res_dir <- paste0(out_dir, "gaussian_mixture_clusters/")
  dir.create(res_dir)
} else {
  res_dir <- NULL
}

X <- scores[, c("EF1", "EF2")]

gmm <- Mclust(X, 
               G = 3, 
               modelNames = "VVV", # ellipsoidal, varying volume, shape, and orientation
               control = emControl(tol = 1e-5))

scores$clusters_mclust <- gmm$classification
scores$clusters_mclust[scores$clusters_mclust == 1] <- "iii"
scores$clusters_mclust[scores$clusters_mclust == 2] <- "ii"
scores$clusters_mclust[scores$clusters_mclust == 3] <- "i"

scores$clusters_mclust_2gr <- scores$clusters_mclust
scores$clusters_mclust_2gr[scores$clusters_mclust %in% c("i", "ii")] <- "i+ii"

scores$mclust_prob_assignment <- apply(gmm$z, 1, function(rr) max(rr))

scores_mclust <- subset(scores, 
                        select = c("subject_id", 
                                   "clusters_mclust", 
                                   "clusters_mclust_2gr", 
                                   "mclust_prob_assignment"))
df_scores <- dplyr::left_join(df_info, scores_mclust)

if (bool_save) {
  save(df_scores, file = file.path(out_dir, "all_scores.RData"))
}
