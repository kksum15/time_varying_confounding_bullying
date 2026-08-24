##########################################################################################
# Date script created: 26 August 2025
# Estimating causal effects using IPW: Modified treatment policy 
# Bullying-internalizing symptoms
# Dataset: file.path(data_dir,"processed_data/imp_merged_female_072025.RData") #8808 records; 52 variables
#          file.path(data_dir,"processed_data/imp_merged_male_072025.RData")
#############################################################################################

####################
# Load packages
####################

library(mice)
library(VIM)
library(sandwich)
library(tidyverse)
library(summarytools)
library(janitor)
library(patchwork)
library(parallel)
library(writexl)

####################
# Load datasets
####################

load(file = file.path(data_dir,"processed_data/imp_merged_female_3int_072025.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_male_3int_072025.RData"))

load(file = file.path(data_dir,"processed_data/ipw_female.RData")) #list containing all imputed data with stabilized weights 
load(file = file.path(data_dir,"processed_data/ipw_male.RData"))

########################
# Source functions 
#######################

source("/Volumes/005/working/scripts/bullying_ob/msm_lmtp_functions.R")

##########################
# Define vectors and lists
##########################
fac_vars <- c("e0_victimization2", "e1_victimization2", "e2_victimization2", "depression_icd_10", "anxiety_cisr")
time_fixed_covs <- c("c0_bw_g", "c0_m_parity", "c0_m_age_delivery", "c0_medu", "c0_fedu", "c0_hhincome", "c0_m_marital", "c0_n_siblings3")
time_varying_t0 <- c("c0_sdq_total", "c0_m_mh1", "c0_f_mh1", "c0_child_maltx1", "c0_sib_aggression_bin", "c0_sib_victimization_bin", 
                     "c0_dv1", "c0_sd_bmi")
time_varying_t1 <- c("c1_sdq_total", "c1_m_mh", "c1_f_mh", "c1_child_maltx","c1_sib_aggression", "c1_sib_victimization", "c1_dv", 
                     "c1_sd_bmi")
time_varying_t2 <- c("c2_sdq_total", "c2_m_mh", "c2_f_mh", "c2_child_maltx","c2_sib_aggression", "c2_sib_victimization", "c2_dv2", 
                     "c2_sd_bmi")
pubertal_covs_t0 <- list(
  female = c("c0_pub_f_breast", "c0_pub_f_hair"),
  male = c("c0_pub_m_genitalia", "c0_pub_m_hair")
)
pubertal_covs_t1 <- list(
  female = c("c1_pub_f_breast", "c1_pub_f_hair"),
  male = c("c1_pub_m_genitalia", "c1_pub_m_hair")
)
pubertal_covs_t2 <- list(
  female = c("c2_pub_f_breast", "c2_pub_f_hair"),
  male = c("c2_pub_m_genitalia", "c2_pub_m_hair")
)
outcomes <- c("dep_anx", "depression_icd_10", "anxiety_cisr")
sexes <- c("female", "male")
datasets <- list(female=ipw_female, male=ipw_male)  # 100 weighted imputed datasets

datasets$female <- convert_list(datasets$female, fac_vars)
datasets$male   <- convert_list(datasets$male, fac_vars)

##############################################
# 1. Define the intervention (reduce by 15%, 25%, 40%)
# 2. Construct the dataset with both the original and shifted data (redefine intervention)
# 3. Run regression with lambda as the outcome conditional on At and Ht; will have three density ratios for each time point
# 4. Predict the model to estimate the conditional probability that lambda = 1 
# 5. Obtain the density ratios = prob lambda / 1 - probability lambda 
# 6. Obtain the weights = product of three density ratios 
# 7. Calculate the marginal mean outcome 

# 1. Define stochastic intervention functions
##############################################
set.seed(78263)
policies <- list( 
  d15 = make_d(0.15), 
  d25 = make_d(0.25), 
  d40 = make_d(0.40))

##############################################
# 2. Generate datasets with weights
##############################################
datasets_with_cdr <- list()

for (sex in sexes) {
  datasets_sex <- datasets[[sex]]
  datasets_with_cdr[[sex]] <- list()
  
  for (i in seq_along(datasets_sex)) {
    dat <- datasets_sex[[i]]
    dat <- estimate_policies_cdr(dat, sex)  # compute all policies at once
    dat$imp <- i
    datasets_with_cdr[[sex]][[i]] <- dat
  }
}

save(datasets_with_cdr, file = file.path(data_dir,"processed_data/datasets_with_cdr.RData"))

cdr_female <- datasets_with_cdr$female 
cdr_female_df <- bind_rows(datasets_with_cdr$female, .id = "imp") %>%
  mutate(imp = as.numeric(imp))

cdr_male <- datasets_with_cdr$male 
cdr_male_df <- bind_rows(datasets_with_cdr$male, .id = "imp") %>%
  mutate(imp = as.numeric(imp))

cdr_diag <- function(df, sex_label){
  
  trunc_vars <- grep("^cdr_.*_trunc$", names(df), value = TRUE)
  diag_list <- lapply(trunc_vars, function(v){
    raw <- sub("_trunc$", "", v)
    truncated_flag <- df[[raw]] > df[[v]]
    
    data.frame(
      sex = sex_label,
      policy = raw,
      n = nrow(df),
      
      n_truncated = sum(truncated_flag, na.rm = TRUE),
      pct_truncated = mean(truncated_flag, na.rm = TRUE) * 100,
      
      raw_mean = mean(df[[raw]], na.rm = TRUE),
      trunc_mean = mean(df[[v]], na.rm = TRUE),
      
      raw_min = min(df[[raw]], na.rm = TRUE),
      raw_max = max(df[[raw]], na.rm = TRUE),
      
      trunc_min = min(df[[v]], na.rm = TRUE),
      trunc_max = max(df[[v]], na.rm = TRUE)
    )
  })
  
  bind_rows(diag_list)
}

cdr_diag_by_imp <- function(df, sex_label){
  
  trunc_vars <- grep("^cdr_.*_trunc$", names(df), value = TRUE)
  diag_list <- lapply(trunc_vars, function(v){
    raw <- sub("_trunc$", "", v)
    
    df %>%
      group_by(imp) %>%
      summarise(
        sex = sex_label,
        policy = raw,
        
        n = n(),
        n_truncated = sum(.data[[raw]] > .data[[v]], na.rm = TRUE),
        pct_truncated = mean(.data[[raw]] > .data[[v]], na.rm = TRUE) * 100,
        
        raw_mean = mean(.data[[raw]], na.rm = TRUE),
        trunc_mean = mean(.data[[v]], na.rm = TRUE),
        
        raw_max = max(.data[[raw]], na.rm = TRUE),
        trunc_max = max(.data[[v]], na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  bind_rows(diag_list)
}

diag_female <- cdr_diag(cdr_female_df, "female")
diag_male   <- cdr_diag(cdr_male_df, "male")
diag <- bind_rows(diag_female, diag_male)

diag_imp_female <- cdr_diag_by_imp(cdr_female_df, "female")
diag_imp_male   <- cdr_diag_by_imp(cdr_male_df, "male")
diag_imp <- bind_rows(diag_imp_female, diag_imp_male)

write_xlsx(
  list(
    cdr_diag_pooled = diag,
    cdr_diag_by_imp = diag_imp
  ),
  path = file.path(output_dir, "msm/cdr_diagnostics.xlsx")
)


# check CDR dist 
cdr_female_df1 <- cdr_female_df %>%
  select(imp, cdr_d15_trunc, cdr_d25_trunc, cdr_d40_trunc) %>%
  mutate(sex="Female")

cdr_male_df1 <- cdr_male_df %>%
  select(imp, cdr_d15_trunc, cdr_d25_trunc, cdr_d40_trunc) %>%
  mutate(sex="Male")

cdr_overall_df <- bind_rows(cdr_female_df1, cdr_male_df1 ) %>%
  filter(imp <10) %>%
  pivot_longer(
    cols = ends_with("_trunc"),  
    names_to = "int",           
    values_to = "cdr"            
  ) %>%
  mutate(int = sub("cdr_", "", int))

cdr_overall_df %>%
  ggplot(aes(x = factor(imp), y = cdr, col = sex)) +
  geom_boxplot() +
  facet_grid(int ~ .) +
  labs(title = "CDR distribution", y = "Density Ratio", x = "Imputation") +
  scale_y_continuous(limits = c(0, NA)) + 
  scale_fill_brewer(palette = "Set2") +
  theme (
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    panel.background = element_rect(fill="white")
  )

ggsave(file.path(output_dir, "msm/cdr_dist_overall.png"), width = 8, height = 6, dpi = 300)

cdr_female_df %>%
  select(imp, cdr_d15_trunc, cdr_d25_trunc, cdr_d40_trunc) %>%
  filter(imp < 10) %>%
  pivot_longer(
    cols = ends_with("_trunc"),  
    names_to = "int",           
    values_to = "cdr"            
  ) %>%
  mutate(int = sub("cdr_", "", int)) %>%
  
  ggplot(aes(x = factor(imp), y = cdr)) +
  geom_boxplot() +
  facet_grid(int ~ .) +
  labs(title = "CDR distribution (female)", y = "Density Ratio", x = "Imputation") +
  scale_y_continuous(limits = c(0, NA)) + 
  scale_fill_brewer(palette = "Set2") +
  theme (
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    panel.background = element_rect(fill="white")
  )

ggsave(file.path(output_dir, "msm/cdr_dist_female.png"), width = 8, height = 6, dpi = 300)

cdr_male_df %>%
  select(imp, cdr_d15_trunc, cdr_d25_trunc, cdr_d40_trunc) %>%
  filter(imp < 10) %>%
  pivot_longer(
    cols = ends_with("_trunc"),  
    names_to = "int",           
    values_to = "cdr"            
  ) %>%
  mutate(int = sub("cdr_", "", int)) %>%
  
  ggplot(aes(x = factor(imp), y = cdr)) +
  geom_boxplot() +
  facet_grid(int ~ .) +
  labs(title = "CDR distribution (male)", y = "Density Ratio", x = "Imputation") +
  scale_y_continuous(limits = c(0, NA)) + 
  scale_fill_brewer(palette = "Set2") +
  theme (
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    panel.background = element_rect(fill="white")
  )

ggsave(file.path(output_dir, "msm/cdr_dist_male.png"), width = 8, height = 6, dpi = 300)

##############################################
# 3. Estimate risks from weighted datasets
##############################################

df_risks <- purrr::map_dfr(sexes, function(sex) {
  purrr::map_dfr(seq_along(datasets_with_cdr[[sex]]), function(i) {
    dat <- datasets_with_cdr[[sex]][[i]]
    purrr::map_dfr(outcomes, function(out) {
      risk_none <- mean(dat[[out]])
      risks <- estimate_risks_from_cdr(dat, out)
      tibble::tibble(sex, 
                     outcome = out, 
                     imp = i,
                     policy = c("cdr_d0", names(risks)), 
                     risk = c(risk_none, risks)
      )
    })
  })
})

##############################################
# 4. Compute RRs and RDs per imputation, then pool
##############################################
df_effects <- df_risks %>%
  pivot_wider(names_from = policy, values_from = risk) %>%
  mutate(
    rr_d15_vs_d0 = cdr_d15 / cdr_d0,
    rr_d25_vs_d0 = cdr_d25 / cdr_d0,
    rr_d40_vs_d0 = cdr_d40 / cdr_d0,
    rr_d25_vs_d15 = cdr_d25 / cdr_d15,
    rr_d40_vs_d15 = cdr_d40 / cdr_d15,
    rd_d15_vs_d0 = cdr_d15 - cdr_d0,
    rd_d25_vs_d0 = cdr_d25 - cdr_d0,
    rd_d40_vs_d0 = cdr_d40 - cdr_d0,
    rd_d25_vs_d15  = cdr_d25 - cdr_d15,
    rd_d40_vs_d15  = cdr_d40 - cdr_d15
  )

# Pool across imputations (mean)
df_pooled <- df_effects %>%
  group_by(sex, outcome) %>%
  summarise(across(-imp, list(est = mean, between_var = var), 
                   .names = "{.col}_{.fn}"),
            .groups = "drop")

df_pooled <- df_pooled %>%
  pivot_longer(
    cols = -c(sex, outcome),  # keep sex & outcome fixed
    names_to = c("param", ".value"),  # split into param + value columns
    names_pattern = "(.*)_(est|between_var)" # regex to capture base + suffix
  ) %>%
  mutate(param = sub("^cdr_", "risk_cdr_", param))

saveRDS(df_pooled, file.path(output_dir, "msm/lmtp/df_pooled_rr_rd.RDS"))

##############################################
# Bootstrapping to get SE 

# 1. Bootstrap within each imputation set to get the SE (bootstrap SD)

# 2. Repeat bootstrap for 100 imputations 

# 3. Rubin's rule to get the within imputation variance 
# withinimp_var = mean(SE)^2
# betweenimp_var = var_rr and var_rd from combined_summary_df 
# get the total_var =  withinimp_var + (1 + 1/m) * betweenimp_var 

# 4. Loop for each sex and outcome 
##############################################

m = 100
B = 500

cl <- makeCluster(detectCores() - 1)
clusterSetRNGStream(cl, 198676)
clusterExport(cl, varlist = c("datasets", "estimate_policies_cdr", "bootstrap_lmtp", 
                              "make_d", "policies", "m", "B", "time_fixed_covs", 
                              "time_varying_t0", "pubertal_covs_t0",
                              "time_varying_t1", "pubertal_covs_t1", 
                              "time_varying_t2", "pubertal_covs_t2"))
clusterEvalQ(cl, {
  library(dplyr)
  library(purrr)
})

lmtp_results <- list()

start <- Sys.time()
for (sex in sexes) {
  for (outcome in outcomes) {
    message("Running: ", sex, " × ", outcome)
    lmtp_se_df <- estimate_lmtp_par(sex, outcome, m)
    lmtp_within_var_df <- combine_bootstrap_variance_lmtp(lmtp_se_df, sex, outcome)
    lmtp_results[[paste(sex, outcome, sep = "_")]] <- lmtp_within_var_df
  }
}

# Bind everything into one tibble
lmtp_results <- dplyr::bind_rows(lmtp_results)

stopCluster(cl)
end <- Sys.time()
end - start 

saveRDS(lmtp_results, file.path(output_dir, "msm/lmtp/final_bootstrap_var.RDS"))

overall_results_lmtp <- lmtp_results %>% 
  left_join(df_pooled, by = c("sex", "outcome", "param")) %>% 
  mutate(total_var = within_var + ((1 + 1/m) * between_var),
         se = sqrt(total_var),
         lci = est + qnorm(0.025) * se,
         uci = est + qnorm(0.975) * se
  ) %>%
  separate(param, into = c("measure", "comparison"), sep = "_", extra = "merge")


overall_results_lmtp$formatted <- ifelse(
  overall_results_lmtp$measure == "rd",
  paste0(
    sprintf("%.3f", overall_results_lmtp$est), " [",
    sprintf("%.3f", overall_results_lmtp$lci), ", ",
    sprintf("%.3f", overall_results_lmtp$uci), "]"
  ),
  paste0(
    sprintf("%.2f", overall_results_lmtp$est), " [",
    sprintf("%.2f", overall_results_lmtp$lci), ", ",
    sprintf("%.2f", overall_results_lmtp$uci), "]"
  )
)

overall_results_lmtp <- overall_results_lmtp %>% 
  mutate(outcome = factor(outcome, 
                          levels = c("dep_anx", "anxiety_cisr", "depression_icd_10"),
                          labels = c("Depression or Anxiety", "Anxiety", "Depression")),
         sex = factor(sex,
                      levels = c("female", "male"),
                      labels = c("Female", "Male")), 
         comparison = factor(comparison, 
                             levels = c("d15_vs_d0", "d25_vs_d0", "d40_vs_d0", "d25_vs_d15", "d40_vs_d15")))

saveRDS(overall_results_lmtp, file.path(output_dir, "msm/lmtp/overall_results_lmtp.RDS"))

ggplot(filter(overall_results_lmtp, measure == "rr", str_detect(comparison, "_d15")), aes(x = est, y = comparison, color = sex)) +
  geom_point(position = position_dodge(width = 0.3), size = 1.5) +
  geom_errorbarh(aes(xmin = lci, xmax = uci),
                 position = position_dodge(width = 0.3), height = 0.1) +
  geom_vline(xintercept = 1, linetype = "dashed") +  # reference line at RR=1
  geom_text(
    aes(
      label = formatted,
      x = 1.05, 
      group = sex                      
    ),
    position = position_dodge(width = 0.3),
    hjust = 0,
    size = 2.7,
    color = "black",
    show.legend = FALSE
  ) +
  facet_grid(. ~ outcome, space = "free_x") +  # horizontal facets, shared y
  theme (
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    panel.background = element_rect(fill="white"),
    axis.text.y = element_text(size = 9),
    panel.spacing = unit(0.1, "lines"),
    legend.position = "bottom"
  ) +
  scale_x_continuous(limits = c(0.8, 1.2, na.rm = TRUE),
                     breaks = c(0.9, 1, 1.1)) +
  scale_color_brewer(palette ="Set2") +
  labs(x = "Risk ratio", y = NULL, color = "Group") 

ggsave(file.path(output_dir, "msm/lmtp/lmtp_rr.png"), width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "msm/lmtp/lmtp_rr_d0.png"), width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "msm/lmtp/lmtp_rr_d15.png"), width = 8, height = 6, dpi = 300)

ggplot(filter(overall_results_lmtp, measure == "rd", str_detect(comparison, "d0")), aes(x = est, y = comparison, color = sex)) +
  geom_point(position = position_dodge(width = 0.3), size = 1.5) +
  geom_errorbarh(aes(xmin = lci, xmax = uci),
                 position = position_dodge(width = 0.3), height = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed") +  # reference line at RD=0
  geom_text(
    aes(
      label = formatted,
      x = 0.007, 
      group = sex                      
    ),
    position = position_dodge(width = 0.3),
    hjust = 0,
    size = 2.7,
    color = "black",
    show.legend = FALSE
  ) +
  facet_grid(. ~ outcome, space = "free_x") +  # horizontal facets, shared y
  theme (
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    panel.background = element_rect(fill="white"),
    axis.text.y = element_text(size = 9),
    panel.spacing = unit(0.1, "lines"),
    legend.position = "bottom"
  ) +
  scale_x_continuous(limits = c(-0.03, 0.04, na.rm = TRUE)) +
  scale_color_brewer(palette ="Set2") +
  labs(x = "Risk difference", y = NULL, color = "Group") 

ggsave(file.path(output_dir, "msm/lmtp/lmtp_rd.png"), width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "msm/lmtp/lmtp_rd_d0.png"), width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "msm/lmtp/lmtp_rd_d15.png"), width = 8, height = 6, dpi = 300)




