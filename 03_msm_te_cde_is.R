##########################################################################################
# Date script created: 14 July 2025
# MSM: Total effects and controlled direct effects 
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
library(cobalt)
#library(devtools)
#devtools::install_github("jwjackson/confoundr")
#library(confoundr)

####################
# Load dataset
####################

load(file = file.path(data_dir,"processed_data/imp_merged_female_3int_072025.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_male_3int_072025.RData"))

########################
# Source functions 
#######################

source("/Volumes/005/working/scripts/bullying_ob/msm_joint_effects_functions.R")

####################################################################################################################
# Compute joint effects of bullying (total = 50 variables)
# outcome (3): dep_anx, depression_icd_10, anxiety_cisr

# exposure (3): e0_victimization2, e1_victimization2, e2_victimization2

# time-fixed variables (8): c0_bw_g, c0_m_parity, c0_m_age_delivery, c0_medu, c0_fedu, c0_hhincome, 
#                       c0_m_marital, c0_n_siblings3

# time-varying variables: 
# T0 (12): c0_sdq_total, c0_m_mh1, c0_f_mh1, c0_child_maltx1, c0_sib_aggression_bin, c0_sib_victimization_bin, 
#        c0_dv1, c0_sd_bmi, c0_pub_f_breast, c0_pub_f_hair, c0_pub_m_genitalia, c0_pub_m_hair

# T1 (12): c1_sdq_total, c1_m_mh, c1_f_mh, c1_child_maltx, c1_sib_aggression, c1_sib_victimization, 
#          c1_dv, c1_sd_bmi, c1_pub_f_breast, c1_pub_f_hair, c1_pub_m_genitalia, c1_pub_m_hair

# T2 (12): c2_sdq_total, c2_m_mh, c2_f_mh, c2_child_maltx, c2_sib_aggression, c2_sib_victimization, 
#          c2_dv2, c2_sd_bmi, c2_pub_f_breast, c2_pub_f_hair, c2_pub_m_genitalia, c2_pub_m_hair
###################################################################################################################

# Step 1 - estimate the propensity score in each imputed dataset
########################################################################

m <- 100 #number of imputations 

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

# Model each individual's probability of receiving the treatment from their baseline covariates
# Technically weight for exposure at time 1 = 1 because numerator and denominator predicted probabilities cancel out 

ps_female <- create_ps_models("female", imp_merged_female)
ps_male <- create_ps_models("male", imp_merged_male)
ipw_female <- process_sex_ipw("female", imp_merged_female, ps_female, 100)
ipw_male <- process_sex_ipw("male", imp_merged_male, ps_male, 100)

ipw_female_1 <- ipw_female[[1]]
summary(ipw_female_1$sw_a0)
summary(ipw_female_1$sw_a1)
summary(ipw_female_1$sw_a2)
summary(ipw_female_1$sw_overall)
summary(ipw_female_1$sw_overall_trunc)
ggplot(ipw_female_1, aes(x=sw_overall, y=sw_overall_trunc)) +
  geom_point()
sum(ipw_female_1$sw_overall_trunc == quantile(ipw_female_1$sw_overall, 0.99, na.rm = TRUE), na.rm = TRUE) #45 

ipw_male_1 <- ipw_male[[1]]
summary(ipw_male_1$sw_a0)
summary(ipw_male_1$sw_a1)
summary(ipw_male_1$sw_a2)
summary(ipw_male_1$sw_overall)
ggplot(ipw_male_1, aes(x=sw_overall, y=sw_overall_trunc)) +
  geom_point()
sum(ipw_male_1$sw_overall_trunc == quantile(ipw_male_1$sw_overall, 0.99, na.rm = TRUE), na.rm = TRUE) #44

#save(ipw_female, file = file.path(data_dir,"processed_data/ipw_female.RData"))
#save(ipw_male, file = file.path(data_dir,"processed_data/ipw_male.RData"))

load(file = file.path(data_dir,"processed_data/ipw_female.RData")) #list containing all imputed data with stabilized weights 
load(file = file.path(data_dir,"processed_data/ipw_male.RData"))

# Step 2: estimate the treatment effect via weighted analysis in each imputed dataset
####################################################################################
outcomes <- c("dep_anx", "depression_icd_10", "anxiety_cisr")
datasets <- list(ipw_female, ipw_male)

results_female_depanx <- outcome_model("female", ipw_female, "dep_anx")
results_female_dep <- outcome_model("female", ipw_female, "depression_icd_10")
results_female_anx <- outcome_model("female", ipw_female, "anxiety_cisr")
results_male_depanx <- outcome_model("male", ipw_male, "dep_anx")
results_male_dep <- outcome_model("male", ipw_male, "depression_icd_10")
results_male_anx <- outcome_model("male", ipw_male, "anxiety_cisr")

results <- rbind(results_female_depanx, results_female_dep, results_female_anx, 
                 results_male_depanx, results_male_dep, results_male_anx) %>%
  mutate(effects = case_when(model == "e0_only" ~ "total",
                             model == "e0_e1" ~ "total",
                             model == "e0_e1_e2" & variable == "e0_victimization2" ~ "direct",
                             model == "e0_e1_e2" & variable == "e1_victimization2" ~ "direct", 
                             model == "e0_e1_e2" & variable == "e2_victimization2" ~ "total"))


# Manually combine the results using Rubin's rules 
pooled_results_is <- results %>%
  group_by(sex, outcome, variable, effects) %>%
  summarize(
    pooled_estimate = mean(logOR),
    withinimp_var = mean(stderr^2),
    betweenimp_var = var(logOR),       # var() already divides by (m-1)
    total_var = withinimp_var + (1 + 1/m) * betweenimp_var,
    lci = pooled_estimate + qnorm(0.025) * sqrt(total_var),
    uci = pooled_estimate + qnorm(0.975) * sqrt(total_var),
    .groups = "drop"
  )

rownames(pooled_results_is) <- NULL
pooled_results_is$formatted <- paste0(
  sprintf("%.2f", exp(pooled_results_is$pooled_estimate)), " [",
  sprintf("%.2f", exp(pooled_results_is$lci)), ", ", 
  sprintf("%.2f", exp(pooled_results_is$uci)), "]"
)

variable_labels <- c(
  "e0_victimization2" = "Victimization at 8 years",
  "e1_victimization2" = "Victimization at 10 years",
  "e2_victimization2" = "Victimization at 12 years"
)
pooled_results_is$variable_label <- variable_labels[pooled_results_is$variable]
pooled_results_is$variable_label <- factor(pooled_results_is$variable_label, levels = variable_labels)

pooled_results_is <- pooled_results_is %>% 
  mutate(outcome = factor(outcome, 
                          levels = c("dep_anx", "anxiety_cisr", "depression_icd_10"),
                          labels = c("Depression or Anxiety", "Anxiety", "Depression")),
         sex = factor(sex,
                      levels = c("female", "male"),
                      labels = c("Female", "Male")))

save(pooled_results_is, file = file.path(output_dir,"msm/pooled_results_joint_is.RData"))

# Plot results 
ggplot(filter(pooled_results, effects == "total"), aes(x = exp(pooled_estimate), y = variable_label, color = sex)) +
  geom_point(position = position_dodge(width = 0.3), size = 1.5) +
  geom_errorbarh(aes(xmin = exp(lci), xmax = exp(uci)),
                 position = position_dodge(width = 0.3), height = 0.1) +
  geom_vline(xintercept = 1, linetype = "dashed") +  # reference line at OR=1
  geom_text(
    aes(
      label = formatted,
      x = max(exp(uci)), 
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
  scale_color_brewer(palette ="Set2") + 
  scale_x_continuous(trans = "log10", 
                     limits = c(0.5, max(exp(pooled_results$uci), na.rm = TRUE) * 4),
                     breaks = c(0.5, 1, 2, 4, 8),
                     labels = scales::comma) +
  labs(x = "Odds ratio", y = NULL, color = "Group", title = "Total Effects") 

ggsave(file.path(output_dir, "msm/total_effects.png"), width = 8, height = 6, dpi = 300)

ggplot(filter(pooled_results, effects == "direct"), aes(x = exp(pooled_estimate), y = variable_label, color = sex)) +
  geom_point(position = position_dodge(width = 0.3), size = 1.5) +
  geom_errorbarh(aes(xmin = exp(lci), xmax = exp(uci)),
                 position = position_dodge(width = 0.3), height = 0.1) +
  geom_vline(xintercept = 1, linetype = "dashed") +  # reference line at OR=1
  geom_text(
    aes(
      label = formatted,
      x = max(exp(uci)), 
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
  scale_color_brewer(palette ="Set2") + 
  scale_x_continuous(trans = "log10", 
                     limits = c(0.5, max(exp(pooled_results$uci), na.rm = TRUE) * 4),
                     breaks = c(0.5, 1, 2, 4, 8),
                     labels = scales::comma) +
  labs(x = "Odds ratio", y = NULL, color = "Group", title = "Direct Effects") 

ggsave(file.path(output_dir, "msm/direct_effects.png"), width = 8, height = 6, dpi = 300)


####################################################################################################
# Diagnostics
###################

ipw_female_df <- do.call(
  rbind,
  lapply(seq_along(ipw_female), function(i) {
    # Add the 'imp' column to each data frame
    df <- ipw_female[[i]]
    df$imp <- i
    df
  })
)

ipw_male_df <- do.call(
  rbind,
  lapply(seq_along(ipw_male), function(i) {
    # Add the 'imp' column to each data frame
    df <- ipw_male[[i]]
    df$imp <- i
    df
  })
)

saveRDS(ipw_female_df, file = file.path(data_dir,"processed_data/ipw_female_df.RDS"))
saveRDS(ipw_male_df, file = file.path(data_dir,"processed_data/ipw_male_df.RDS"))

# IPW diagnostics 
####################

ipw_diag<- function(df, weight_raw, weight_trunc, sex_label){
  
  truncated_flag <- df[[weight_raw]] > df[[weight_trunc]]
  
  data.frame(
    sex = sex_label,
    n = nrow(df),
    
    n_truncated = sum(truncated_flag, na.rm = TRUE),
    pct_truncated = mean(truncated_flag, na.rm = TRUE) * 100,
    
    raw_mean = mean(df[[weight_raw]], na.rm = TRUE),
    trunc_mean = mean(df[[weight_trunc]], na.rm = TRUE),
    
    raw_min = min(df[[weight_raw]], na.rm = TRUE),
    raw_max = max(df[[weight_raw]], na.rm = TRUE),
    
    trunc_min = min(df[[weight_trunc]], na.rm = TRUE),
    trunc_max = max(df[[weight_trunc]], na.rm = TRUE)
  )
}


ipw_diag_by_imp <- function(df, weight_raw, weight_trunc, sex_label){
  
  df %>%
    group_by(imp) %>%
    summarise(
      sex = sex_label,
      n = n(),
      
      n_truncated = sum(.data[[weight_raw]] > .data[[weight_trunc]], na.rm = TRUE),
      pct_truncated = mean(.data[[weight_raw]] > .data[[weight_trunc]], na.rm = TRUE) * 100,
      
      raw_mean = mean(.data[[weight_raw]], na.rm = TRUE),
      trunc_mean = mean(.data[[weight_trunc]], na.rm = TRUE),
      
      raw_max = max(.data[[weight_raw]], na.rm = TRUE),
      trunc_max = max(.data[[weight_trunc]], na.rm = TRUE),
      .groups = "drop"
    )
}

ipw_diag_female <- ipw_diag(ipw_female_df, "sw_overall", "sw_overall_trunc", "female")
ipw_diag_male <- ipw_diag(ipw_male_df, "sw_overall", "sw_overall_trunc", "male")
ipw_diag_pooled <- bind_rows(ipw_diag_female, ipw_diag_male)

ipw_diag_imp_female <- ipw_diag_by_imp(ipw_female_df, "sw_overall", "sw_overall_trunc", "female")
ipw_diag_imp_male <- ipw_diag_by_imp(ipw_male_df, "sw_overall", "sw_overall_trunc", "male")
ipw_diag_imp <- bind_rows(ipw_diag_imp_female, ipw_diag_imp_male)

write_xlsx(
  list(
    ipw_diag_pooled = ipw_diag_pooled,
    ipw_diag_imp = ipw_diag_imp
  ),
  path = file.path(output_dir, "msm/ipw_diagnostics.xlsx")
)

# PS overlap plot
####################

ps_female_avg <-  ipw_female_df %>%
  select(imp, ps_a0_denom, ps_a1_denom, ps_a2_denom, e0_victimization2, e1_victimization2, e2_victimization2) %>%
  mutate(rowid = rep(1:nrow(ipw_female[[1]]), length(ipw_female))) %>%
  group_by(rowid) %>%
  summarise(across(starts_with("ps_"), mean, na.rm = TRUE),
            across(starts_with("e"), first),
            .groups = "drop") %>%
  rename_with(~ str_replace(., "ps_a([0-9])_denom", "ps_\\1"), starts_with("ps")) %>%
  rename_with(~ str_replace(., "e([0-9])_victimization2", "victimization_\\1"), starts_with("e")) %>%
  pivot_longer(
    cols = -rowid,
    names_to = c(".value", "timepoint"),
    names_sep = "_"
  ) %>%
  mutate(sex = "Female")

ps_male_avg <-  ipw_male_df %>%
  select(imp, ps_a0_denom, ps_a1_denom, ps_a2_denom, e0_victimization2, e1_victimization2, e2_victimization2) %>%
  mutate(rowid = rep(1:nrow(ipw_male[[1]]), length(ipw_male))) %>%
  group_by(rowid) %>%
  summarise(across(starts_with("ps_"), mean, na.rm = TRUE),
            across(starts_with("e"), first),
            .groups = "drop") %>%
  rename_with(~ str_replace(., "ps_a([0-9])_denom", "ps_\\1"), starts_with("ps")) %>%
  rename_with(~ str_replace(., "e([0-9])_victimization2", "victimization_\\1"), starts_with("e")) %>%
  pivot_longer(
    cols = -rowid,
    names_to = c(".value", "timepoint"),
    names_sep = "_"
  ) %>%
  mutate(sex = "Male")

ps_overall_avg <- bind_rows(ps_female_avg, ps_male_avg)

ps_overall_avg %>%
  ggplot(., aes(x = ps, fill = victimization)) +
  geom_density(alpha = 0.4) +
  facet_grid(sex~timepoint) +
  labs(title = "PS Overlap",
       x = "Propensity Score", y = "Density") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") 

ggsave(file.path(output_dir, "msm/ps_overlap_avg.png"), width = 8, height = 6, dpi = 300)

ipw_male[[1]] %>%
  #ipw_female[[1]] %>%
  select(ps_a0_denom, ps_a1_denom, ps_a2_denom, e0_victimization2, e1_victimization2, e2_victimization2) %>%
  rename_with(~ str_replace(., "ps_a([0-9])_denom", "ps_\\1"), starts_with("ps")) %>%
  rename_with(~ str_replace(., "e([0-9])_victimization2", "victimization_\\1"), starts_with("e")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "timepoint"),
    names_sep = "_" 
  ) %>%
  
  ggplot(., aes(x = ps, fill = victimization)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~timepoint) +
  labs(title = "PS Overlap (Male)",
       x = "Propensity Score", y = "Density") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

ggsave(file.path(output_dir, "msm/ps_overlap_female.png"), width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "msm/ps_overlap_male.png"), width = 8, height = 6, dpi = 300)

ipw_female_df %>%
  filter(imp < 10) %>%
  ggplot(aes(x = sw_overall_trunc)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  facet_wrap(.~imp) +
  labs(title = "IPW distribution (female)", y = "Weight", x = "Count") +
  scale_x_continuous(limits = c(0, NA)) + 
  theme (
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    panel.background = element_rect(fill="white")
  )

ipw_female_sw_df <- ipw_female_df %>%
  select(imp, sw_overall_trunc) %>%
  mutate(sex="Female")

ipw_male_sw_df <- ipw_male_df %>%
  select(imp, sw_overall_trunc) %>%
  mutate(sex="Male")

ipw_overall_sw_df <- bind_rows(ipw_female_sw_df, ipw_male_sw_df) %>%
  filter(imp <10) 

ipw_overall_sw_df %>%
  ggplot(aes(x = factor(imp), y = sw_overall_trunc, col = sex)) +
  geom_boxplot() +
  labs(title = "IPW distribution", y = "Weight", x = "Imp") +
  scale_y_continuous(limits = c(0, NA)) + 
  scale_fill_brewer(palette = "Set1") +
  theme (
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    panel.background = element_rect(fill="white")
  ) 

ggsave(file.path(output_dir, "msm/weight_dist_overall.png"), width = 8, height = 6, dpi = 300)

ipw_female_df %>%
  filter(imp < 10) %>%
  ggplot(aes(x = factor(imp), y = sw_overall_trunc)) +
  geom_boxplot() +
  labs(title = "IPW distribution (female)", y = "Weight", x = "Imp") +
  scale_y_continuous(limits = c(0, NA)) + 
  scale_fill_brewer(palette = "Set2") +
  theme (
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    panel.background = element_rect(fill="white")
  )

ggsave(file.path(output_dir, "msm/weight_dist_female.png"), width = 8, height = 6, dpi = 300)

ipw_male_df %>%
  filter(imp < 10) %>%
  ggplot(aes(x = sw_overall_trunc)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  facet_wrap(.~imp) +
  labs(title = "IPW distribution (male)", y = "Weight", x = "Count") +
  scale_x_continuous(limits = c(0, NA)) + 
  theme_minimal()

ipw_male_df %>%
  filter(imp < 10) %>%
  ggplot(aes(x = factor(imp), y = sw_overall_trunc)) +
  geom_boxplot() +
  labs(title = "IPW distribution (male)", y = "Weight", x = "Imp") +
  scale_y_continuous(limits = c(0, NA)) + 
  scale_fill_brewer(palette = "Set2") +
  theme (
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    panel.background = element_rect(fill="white")
  )

ggsave(file.path(output_dir, "msm/weight_dist_male.png"), width = 8, height = 6, dpi = 300)

# Examine balance across all timepoints: using cobalt package 
####################

for (sex in c("male", "female")) {
  assign(paste0("formula_a0_denom_", sex),
         paste("e0_victimization2 ~",
               paste(c(time_fixed_covs, time_varying_t0, pubertal_covs_t0[[sex]]), collapse = " + ")))
  
  assign(paste0("formula_a1_denom_", sex),
         paste("e1_victimization2 ~ e0_victimization2 +",
               paste(c(time_fixed_covs, time_varying_t0, pubertal_covs_t0[[sex]], time_varying_t1, 
                       pubertal_covs_t1[[sex]]), collapse = " + ")))
  
  assign(paste0("formula_a2_denom_", sex),
         paste("e2_victimization2 ~ e0_victimization2 + e1_victimization2 +",
               paste(c(time_fixed_covs, time_varying_t0, pubertal_covs_t0[[sex]], time_varying_t1, 
                       pubertal_covs_t1[[sex]], time_varying_t2, pubertal_covs_t2[[sex]]), collapse = " + ")))
}


bal.tab(list(as.formula(formula_a0_denom_female),
             as.formula(formula_a1_denom_female),
             as.formula(formula_a2_denom_female)),
        data=ipw_female_df,
        imp = "imp",
        weights = "sw_overall_trunc")


love_plot_female <- love.plot(list(as.formula(formula_a0_denom_female),
                                   as.formula(formula_a1_denom_female),
                                   as.formula(formula_a2_denom_female)),
                              imp = "imp",
                              data = ipw_female_df, 
                              weights = "sw_overall_trunc", 
                              s.d.denom ="pooled",
                              standardize = T,
                              thresholds = c(m = .1),
                              #stars="std", 
                              binary = "std") +
  theme (
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white"),
    axis.text.y = element_text(size = 9),
    panel.spacing = unit(0.1, "lines")
  )

love_plot_female
ggsave(file.path(output_dir, "msm/love_plot_female.png"), width = 8.5, height = 7, dpi = 300)
ggsave(file.path(output_dir, "msm/love_plot_female_smd.png"), width = 8.5, height = 7, dpi = 300)

love_plot_male <- love.plot(list(as.formula(formula_a0_denom_male),
                                 as.formula(formula_a1_denom_male),
                                 as.formula(formula_a2_denom_male)),
                            imp = "imp",
                            data = ipw_male_df, 
                            weights = "sw_overall_trunc", 
                            thresholds = c(m = .1),
                            #stars="std",
                            binary ="std") +
  theme (
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color="black"),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white"),
    axis.text.y = element_text(size = 9),
    panel.spacing = unit(0.1, "lines"),
  )

love_plot_male
ggsave(file.path(output_dir, "msm/love_plot_male.png"), width = 8.5, height = 7, dpi = 300)
ggsave(file.path(output_dir, "msm/love_plot_male_smd.png"), width = 8.5, height = 7, dpi = 300)

#ipw_male[[1]] %>%
ipw_female[[1]] %>%
  select(c0_sdq_total, c1_sdq_total, c2_sdq_total, e0_victimization2, e1_victimization2, e2_victimization2) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("timepoint", ".value"),
    names_pattern = "^.([0-2])_(.*)"  
  ) %>%
  
  ggplot(., aes(x = sdq_total, fill = victimization2)) +
  geom_density(alpha = 0.3) +   # alpha for transparency
  facet_wrap(~timepoint, scales = "free") + 
  labs(
    title = "Density of SDQ Total Scores by Pattern (Female)",
    x = "SDQ Total Score",
    y = "Density",
    fill = "Victimization"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "msm/sdq_male.png"), width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "msm/sdq_female.png"), width = 8, height = 6, dpi = 300)
