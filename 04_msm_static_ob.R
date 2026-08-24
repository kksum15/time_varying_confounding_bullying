##########################################################################################
# Date script created: 28 Oct 2025
# MSM: Static interventions 
# Bullying-overweight
# Dataset: file.path(data_dir,"processed_data/imp_merged_female_ob.RData") 
#          file.path(data_dir,"processed_data/imp_merged_male_ob.RData") 
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

####################
# Load datasets
####################

load(file = file.path(data_dir,"processed_data/imp_merged_female_ob.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_male_ob.RData"))

load(file = file.path(data_dir,"processed_data/ipw_female_ob.RData")) #list containing all imputed data with stabilized weights 
load(file = file.path(data_dir,"processed_data/ipw_male_ob.RData"))

exp_output_dir <- file.path(output_dir, "msm", "exposure_regime")

##############################################
# 1. Use SW computed from msm_joint_effects (ipw_female, ipw_male )
# 2. Fit a model for the outcome with interaction term
# 3. Estimate marginal risks
# 3a. Create 4 new datatsets containing these patterns: (YYY, YYN, YNN, NNN); change A1, A2, A3 to the respective patterns
# 3b. Estimate the marginal risks for each dataset 
# 4. Calculate the risk ratios and risk differences with always bullied as the reference (YYY)
##############################################

time_fixed_covs <- c("c0_bw_g", "c0_m_parity", "c0_m_age_delivery", "c0_medu", "c0_fedu", "c0_hhincome", "c0_m_marital", "c0_n_siblings3")

time_varying_t0 <- c("c0_sdq_total", "c0_m_mh1", "c0_f_mh1", "c0_child_maltx1", "c0_sib_aggression_bin", "c0_sib_victimization_bin", 
                     "c0_dv1", "c0_bmi_z", "c0_food_insecurity")

time_varying_t1 <- c("c1_sdq_total", "c1_m_mh", "c1_f_mh", "c1_child_maltx","c1_sib_aggression", "c1_sib_victimization", "c1_dv", 
                     "c1_bmi_z", "c1_food_insecurity")

time_varying_t2 <- c("c2_sdq_total", "c2_m_mh", "c2_f_mh", "c2_child_maltx","c2_sib_aggression", "c2_sib_victimization", "c2_dv2", 
                     "c2_bmi_z", "c2_food_insecurity")

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

outcome <- "overweight"
sexes <- c("female", "male")
datasets <- list(female=ipw_female, male=ipw_male)  # 100 weighted imputed datasets

rr_all <- list()
rd_all <- list()
risk_all <- list()
summary_all <- list()

for (sex in sexes) {
    cat("Processing:", sex, outcome, "\n")
    
    # Initialize for storing across imputations
    rr_list <- list()
    rd_list <- list()
    risk_list <- list()
    
    for (i in 1:100) {
      data <- datasets[[sex]][[i]]
      
      # Build formula
      formula_current <- as.formula(
        paste0("as.factor(", outcome, ") ~ e0_victimization2 * e1_victimization2 * e2_victimization2 + ",
               paste(c(time_fixed_covs, time_varying_t0, pubertal_covs_t0[[sex]]), collapse = " + "))
      )
      
      # Fit the outcome model with interaction term 
      model <- glm(
        formula = formula_current,
        data = data,
        family = binomial(link = "logit"),
        weights = sw_overall_trunc
      )
      
      # Create 4 new datatsets containing these patterns: (YYY, YYN, YNN, NNN)
      pattern_data <- list(
        YYY = data %>% mutate(e0_victimization2 = "Yes", e1_victimization2 = "Yes", e2_victimization2 = "Yes"),
        YYN = data %>% mutate(e0_victimization2 = "Yes", e1_victimization2 = "Yes", e2_victimization2 = "No"),
        YNY = data %>% mutate(e0_victimization2 = "Yes", e1_victimization2 = "No", e2_victimization2 = "Yes"),
        YNN = data %>% mutate(e0_victimization2 = "Yes", e1_victimization2 = "No",  e2_victimization2 = "No"),
        NYY = data %>% mutate(e0_victimization2 = "No",  e1_victimization2 = "Yes",  e2_victimization2 = "Yes"),
        NYN = data %>% mutate(e0_victimization2 = "No",  e1_victimization2 = "Yes",  e2_victimization2 = "No"),
        NNY = data %>% mutate(e0_victimization2 = "No",  e1_victimization2 = "No",  e2_victimization2 = "Yes"),
        NNN = data %>% mutate(e0_victimization2 = "No",  e1_victimization2 = "No",  e2_victimization2 = "No")
      )
      
      # Estimate the marginal risks and average over all individuals for each dataset  
      mean_risks <- sapply(pattern_data, function(df) {
        mean(predict(model, newdata = df, type = "response"))
      })
      
      # Store risks
      risk_list[[i]] <- data.frame(
        sex = sex,
        outcome = outcome,
        pattern = names(mean_risks),
        imputation = i,
        mean_risk = as.numeric(mean_risks)
      )
      
      # Compute RR and RD relative to YYY
      rr_i <- mean_risks[c("YYN", "YNY", "YNN", "NYY", "NYN", "NNY", "NNN")] / mean_risks["YYY"]
      rd_i <- mean_risks[c("YYN", "YNY", "YNN", "NYY", "NYN", "NNY", "NNN")] - mean_risks["YYY"]
      
      # Store RR and RD
      rr_list[[i]] <- data.frame(
        sex = sex,
        outcome = outcome,
        comparison = names(rr_i),
        imputation = i,
        risk_ratio = as.numeric(rr_i)
      )
      
      rd_list[[i]] <- data.frame(
        sex = sex,
        outcome = outcome,
        comparison = names(rd_i),
        imputation = i,
        risk_difference = as.numeric(rd_i)
      )
    }
    
    # Bind across 100 imputations
    risk_df <- do.call(rbind, risk_list)
    rr_df <- do.call(rbind, rr_list)
    rd_df <- do.call(rbind, rd_list)
    
    # Store full data
    risk_all[[paste0(sex, "_", outcome)]] <- risk_df
    rr_all[[paste0(sex, "_", outcome)]] <- rr_df
    rd_all[[paste0(sex, "_", outcome)]] <- rd_df
    
    # Calculate mean RR and RD across imputations
    rr_summary <- rr_df %>%
      group_by(sex, outcome, comparison) %>%
      summarise(rr = mean(risk_ratio), 
                var_rr = var(risk_ratio),
                .groups = "drop")
    
    rd_summary <- rd_df %>%
      group_by(sex, outcome, comparison) %>%
      summarise(rd = mean(risk_difference),
                var_rd = var(risk_difference),
                .groups = "drop")
    
    # Merge and store summary
    summary_df <- left_join(rr_summary, rd_summary, by = c("sex", "outcome", "comparison"))
    summary_all[[paste0(sex, "_", outcome)]] <- summary_df
}

# Final combined outputs
combined_rr_ob_df <- do.call(rbind, rr_all)
combined_rd_ob_df <- do.call(rbind, rd_all)
combined_risk_ob_df <- do.call(rbind, risk_all)
combined_rr_var_ob_df <- do.call(rbind, summary_all)

risk_summary_ob <- combined_risk_ob_df %>%
  group_by(sex, outcome, pattern) %>%
  summarise(risk = mean(mean_risk),
            .groups = "drop")

risk_summary_ob %>% filter(pattern == "YYY")

save(combined_rr_var_ob_df, file = file.path(exp_output_dir, "summary_level/combined_rr_var_ob_df.RDS"))
save(risk_summary_ob, file = file.path(exp_output_dir, "summary_level/risk_summary_ob.RDS"))

load(file.path(exp_output_dir, "summary_level/combined_rr_var_ob_df.RDS"))

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

# Source functions 
source("/Volumes/005/working/scripts/bullying_ob/msm_exposure_regime_functions.R")

# Set up the cluster
cl <- makeCluster(detectCores() - 1)
clusterSetRNGStream(cl, 198798)
clusterExport(cl, varlist = c("bootstrap_margins", "datasets", "m", "B",
                              "time_fixed_covs", "time_varying_t0", "pubertal_covs_t0",
                              "time_varying_t1", "pubertal_covs_t1", 
                              "time_varying_t2", "pubertal_covs_t2", 
                              "exp_output_dir"))
clusterEvalQ(cl, {
  library(dplyr)
})

bootstrap_var_list <- list()

start <- Sys.time()
for (sex in sexes) {
    message("Running: ", sex, " × ", outcome)
    boot_se_df <- get_bootstrap_SEs_par(sex, outcome, m)
    combined_var_df <- combine_bootstrap_variance(boot_se_df, sex, outcome)
    bootstrap_var_list[[paste0(sex, "_", outcome)]] <- combined_var_df
}
final_bootstrap_var_ob <- do.call(rbind, bootstrap_var_list)
stopCluster(cl)
end <- Sys.time()
end - start # 1.286675 hours

saveRDS(final_bootstrap_var_ob, file.path(exp_output_dir, "summary_level/final_bootstrap_var_ob.RDS"))

final_bootstrap_var_ob
rownames(final_bootstrap_var_ob) <- NULL

overall_results_ob <- final_bootstrap_var_ob %>% 
  left_join(combined_rr_var_ob_df, by = c("sex", "outcome", "comparison")) %>%
  rename(between_var_rr = var_rr,
         between_var_rd = var_rd) %>%
  mutate(total_var_rr = within_var_rr + ((1 + 1/m) * between_var_rr),
         total_var_rd = within_var_rd + ((1 + 1/m) * between_var_rd),
         rr_se = sqrt(total_var_rr),
         rd_se = sqrt(total_var_rd),
         lci_rr = rr + qnorm(0.025) * rr_se,
         uci_rr = rr + qnorm(0.975) * rr_se,
         lci_rd = rd + qnorm(0.025) * rd_se,
         uci_rd = rd + qnorm(0.975) * rd_se)

overall_results_ob$formatted_rr <- paste0(
  sprintf("%.2f", overall_results_ob$rr), " [",
  sprintf("%.2f", overall_results_ob$lci_rr), ", ", 
  sprintf("%.2f", overall_results_ob$uci_rr), "]"
)

overall_results_ob$formatted_rd <- paste0(
  sprintf("%.3f", overall_results_ob$rd), " [",
  sprintf("%.3f", overall_results_ob$lci_rd), ", ", 
  sprintf("%.3f", overall_results_ob$uci_rd), "]"
)

overall_results_ob <- overall_results_ob %>% 
  mutate(sex = factor(sex,
                      levels = c("female", "male"),
                      labels = c("Female", "Male")))

saveRDS(overall_results_ob, file.path(exp_output_dir, "summary_level/final_summary_df_ob.RDS"))

overall_results_is <- readRDS(file.path(exp_output_dir, "summary_level/final_summary_df.RDS"))

overall_results_comb <- rbind(overall_results_ob, overall_results_is) %>%
  mutate(outcome = factor(outcome, 
                          levels = c("Depression or Anxiety", "Anxiety", "Depression", "overweight"),
                          labels = c("Depression or Anxiety", "Anxiety", "Depression", "Overweight")))

saveRDS(overall_results_comb, file.path(exp_output_dir, "summary_level/final_summary_df_comb.RDS"))

ggplot(overall_results_comb, aes(x = rr, y = comparison, color = sex)) +
  geom_point(position = position_dodge(width = 0.3), size = 1.5) +
  geom_errorbarh(aes(xmin = lci_rr, xmax = uci_rr),
                 position = position_dodge(width = 0.3), height = 0.1) +
  geom_vline(xintercept = 1, linetype = "dashed") +  # reference line at RR=1
  geom_text(
    aes(
      label = formatted_rr,
      x = 3, 
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
  scale_x_continuous(limits = c(-0.5, 5, na.rm = TRUE)) +
  scale_color_brewer(palette ="Set2") +
  labs(x = "Risk ratio", y = NULL, color = "Group") 

ggsave(file.path(exp_output_dir, "exp_regime_rr_comb.png"), width = 9, height = 6, dpi = 300)

ggplot(overall_results_comb, aes(x = rd, y = comparison, color = sex)) +
  geom_point(position = position_dodge(width = 0.3), size = 1.5) +
  geom_errorbarh(aes(xmin = lci_rd, xmax = uci_rd),
                 position = position_dodge(width = 0.3), height = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed") +  # reference line at RD=0
  geom_text(
    aes(
      label = formatted_rd,
      x = 0.19, 
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
  scale_x_continuous(limits = c(-0.25, 0.75, na.rm = TRUE)) +
  scale_color_brewer(palette ="Set2") +
  labs(x = "Risk difference", y = NULL, color = "Group") 

ggsave(file.path(exp_output_dir, "exp_regime_rd_comb.png"), width = 9, height = 6, dpi = 300)


