rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

main_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/scripts/")
data_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/data/preprocessed_data/")
out_dir <- file.path(CORE_DIR, "covid-19-systemic-recovery/output/")

setwd(main_dir)

require(dplyr)
require(survival)
require(survminer)

source("fun_utils.R")

my_seed <- 123
set.seed(my_seed)

trunc_days <- 50
load(file.path(out_dir, "all_scores.RData")) # new, includes with mclust

group_type <- "clusters_mclust"


bool_save <- FALSE
if (bool_save) {
  res_dir <- paste0(out_dir, "survival_", group_type, "/")
  dir.create(res_dir)
} else {
  res_dir <- NULL
}

df_surv <- df_scores[!is.na(df_scores[, group_type]), ]

df_surv <- get_early_or_late_samples(df_surv, 
                                     vec_var = NULL, 
                                     bool_early = FALSE, 
                                     single_sample_per_subject = TRUE)$df_comb
df_surv$time <- df_surv$hospital_outcome_days_from_start_sx_or_first_pos_swab
df_surv$time[!(df_surv$hospital_outcome %in% "dead")] <- df_surv$last_observation_days_from_start_sx_or_first_pos_swab[!(df_surv$hospital_outcome %in% "dead")]

df_surv$hospital_outcome[df_surv$severity %in% c("HC", "A", "B")] <- "alive"

df_surv$status <- ifelse(df_surv$hospital_outcome == "dead", 2, 1) # 1 = cencored, 2 = dead

fit <- survfit(Surv(time, status) ~ df_surv[,group_type], data = df_surv)
print(fit)
summary(fit)
summary(fit)$table

d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower)

if (bool_save) {
  pdf(paste0(res_dir, "survival.pdf"), width = 4, height = 5.5, paper='special')
}
ggsurvplot(fit,
           pval = T, conf.int = TRUE,
           pval.coord = c(5, 0.025),
           pval.size = 4,
           risk.table = T, #TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           color = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           xlab="Days from symptom onset",
           ylab="Survival probability",
           risk.table.y.text = FALSE,
           conf.int.alpha = 0.15,
           fontsize = 3.5,
           palette = as.vector(vec_col_gr_3)) 
if (bool_save) {
  dev.off()
}
