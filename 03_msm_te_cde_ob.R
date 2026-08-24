##########################################################################################
# Date script created: 27 Oct 2025
# MSM: Total effects and controlled direct effects  
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
library(cobalt)

####################
# Load dataset
####################

load(file = file.path(data_dir,"processed_data/imp_merged_female_ob.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_male_ob.RData"))

########################
# Source functions 
#######################

source("/Volumes/005/working/scripts/bullying_ob/msm_joint_effects_functions.R")

####################################################################################################################
# Compute joint effects of bullying (total = 51 variables)
# outcome (1): overweight

# exposure (3): e0_victimization2, e1_victimization2, e2_victimization2

# time-fixed variables (8): c0_bw_g, c0_m_parity, c0_m_age_delivery, c0_medu, c0_fedu, c0_hhincome, 
#                       c0_m_marital, c0_n_siblings3

# time-varying variables: 
# T0 (13): c0_sdq_total, c0_m_mh1, c0_f_mh1, c0_child_maltx1, c0_sib_aggression_bin, c0_sib_victimization_bin, 
#        c0_dv1, c0_bmi_z, c0_pub_f_breast, c0_pub_f_hair, c0_pub_m_genitalia, c0_pub_m_hair, c0_food_insecurity

# T1 (13): c1_sdq_total, c1_m_mh, c1_f_mh, c1_child_maltx, c1_sib_aggression, c1_sib_victimization, 
#          c1_dv, c1_bmi_z, c1_pub_f_breast, c1_pub_f_hair, c1_pub_m_genitalia, c1_pub_m_hair, c1_food_insecurity 

# T2 (13): c2_sdq_total, c2_m_mh, c2_f_mh, c2_child_maltx, c2_sib_aggression, c2_sib_victimization, 
#          c2_dv2, c2_bmi_z, c2_pub_f_breast, c2_pub_f_hair, c2_pub_m_genitalia, c2_pub_m_hair, c2_food_insecurity 
###################################################################################################################

# Step 1 - estimate the propensity score in each imputed dataset
########################################################################

m <- 100 #number of imputations 

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

# Model each individual's probability of receiving the treatment from their baseline covariates
# Technically weight for exposure at time 1 = 1 because numerator and denominator predicted probabilities cancel out 
ps_female <- create_ps_models("female", imp_merged_female_ob)
ps_male <- create_ps_models("male", imp_merged_male_ob)
ipw_female <- process_sex_ipw("female", imp_merged_female_ob, ps_female, 100)
ipw_male <- process_sex_ipw("male", imp_merged_male_ob, ps_male, 100)

ipw_female_1 <- ipw_female[[1]]
summary(ipw_female_1$sw_a0) # should be all 1 
summary(ipw_female_1$sw_a1) # mean should center around 1
summary(ipw_female_1$sw_a2) 
summary(ipw_female_1$sw_overall) 
summary(ipw_female_1$sw_overall_trunc)
ggplot(ipw_female_1, aes(x=sw_overall, y=sw_overall_trunc)) +
  geom_point()
sum(ipw_female_1$sw_a2 >= quantile(ipw_female_1$sw_overall, 0.99, na.rm = TRUE), na.rm = TRUE) #23 

ipw_male_1 <- ipw_male[[1]]
summary(ipw_male_1$sw_a0)
summary(ipw_male_1$sw_a1)
summary(ipw_male_1$sw_a2)
summary(ipw_male_1$sw_overall)
ggplot(ipw_male_1, aes(x=sw_overall, y=sw_overall_trunc)) +
  geom_point()
sum(ipw_male_1$sw_a2 >= quantile(ipw_male_1$sw_overall, 0.99, na.rm = TRUE), na.rm = TRUE) #22

save(ipw_female, file = file.path(data_dir,"processed_data/ipw_female_ob.RData"))
save(ipw_male, file = file.path(data_dir,"processed_data/ipw_male_ob.RData"))

# Step 2: estimate the treatment effect via weighted analysis in each imputed dataset
####################################################################################
datasets <- list(ipw_female, ipw_male)

results_female_ow <- outcome_model("female", ipw_female, "overweight")
results_male_ow <- outcome_model("male", ipw_male, "overweight")

results <- rbind(results_female_ow, results_male_ow) %>%
  mutate(effects = case_when(model == "e0_only" ~ "total",
                             model == "e0_e1" ~ "total",
                             model == "e0_e1_e2" & variable == "e0_victimization2" ~ "direct",
                             model == "e0_e1_e2" & variable == "e1_victimization2" ~ "direct", 
                             model == "e0_e1_e2" & variable == "e2_victimization2" ~ "total"))


# Manually combine the results using Rubin's rules 
pooled_results_ob <- results %>%
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

rownames(pooled_results_ob) <- NULL
pooled_results_ob$formatted <- paste0(
  sprintf("%.2f", exp(pooled_results_ob$pooled_estimate)), " [",
  sprintf("%.2f", exp(pooled_results_ob$lci)), ", ", 
  sprintf("%.2f", exp(pooled_results_ob$uci)), "]"
)

variable_labels <- c(
  "e0_victimization2" = "Victimization at 8 years",
  "e1_victimization2" = "Victimization at 10 years",
  "e2_victimization2" = "Victimization at 12 years"
)
pooled_results_ob$variable_label <- variable_labels[pooled_results_ob$variable]
pooled_results_ob$variable_label <- factor(pooled_results_ob$variable_label, levels = variable_labels)

pooled_results_ob <- pooled_results_ob %>% 
  mutate(sex = factor(sex,
                      levels = c("female", "male"),
                      labels = c("Female", "Male")))

#save(pooled_results_ob, file = file.path(output_dir,"/msm/pooled_results_joint_ob.RData"))

load(file = file.path(output_dir,"msm/pooled_results_joint_is.RData"))

pooled_results_overall <- rbind(pooled_results_ob, pooled_results_is) %>%
  mutate(outcome = factor(outcome, 
                          levels = c("Depression or Anxiety", "Anxiety", "Depression", "overweight"),
                          labels = c("Depression or Anxiety", "Anxiety", "Depression", "Overweight")))

save(pooled_results_overall, file = file.path(output_dir,"/msm/pooled_results_joint_overall.RData"))

# Plot results 
ggplot(filter(pooled_results_overall, effects == "total"), aes(x = exp(pooled_estimate), y = variable_label, color = sex)) +
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
                     limits = c(0.5, max(exp(pooled_results_overall$uci), na.rm = TRUE) * 4),
                     breaks = c(0.5, 1, 2, 4, 8),
                     labels = scales::comma) +
  labs(x = "Odds ratio", y = NULL, color = "Group", title = "Total Effects") 

ggsave(file.path(output_dir, "msm/total_effects_overall.png"), width = 9, height = 6, dpi = 300)

ggplot(filter(pooled_results_overall, effects == "direct"), aes(x = exp(pooled_estimate), y = variable_label, color = sex)) +
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
                     limits = c(0.5, max(exp(pooled_results_overall$uci), na.rm = TRUE) * 4),
                     breaks = c(0.5, 1, 2, 4, 8),
                     labels = scales::comma) +
  labs(x = "Odds ratio", y = NULL, color = "Group", title = "Direct Effects") 

ggsave(file.path(output_dir, "msm/direct_effects_overall.png"), width = 9, height = 6, dpi = 300)

####################################################################################################
# Diagnostics
####################

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

# Plot weight distribution 
###########################
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

ggsave(file.path(output_dir, "msm/weight_dist_female_ob.png"), width = 8, height = 6, dpi = 300)

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

ggsave(file.path(output_dir, "msm/weight_dist_male_ob.png"), width = 8, height = 6, dpi = 300)

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


love_plot_female <- love.plot(list(as.formula(formula_a0_denom_female),
                                   as.formula(formula_a1_denom_female),
                                   as.formula(formula_a2_denom_female)),
                              imp = "imp",
                              data = ipw_female_df, 
                              weights = "sw_overall_trunc", 
                              s.d.denom ="pooled",
                              standardize = T,
                              thresholds = c(m = .1),
                              stars="std") +
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
ggsave(file.path(output_dir, "msm/love_plot_female_ob.png"), width = 8.5, height = 7, dpi = 300)

love_plot_male <- love.plot(list(as.formula(formula_a0_denom_male),
                                 as.formula(formula_a1_denom_male),
                                 as.formula(formula_a2_denom_male)),
                            imp = "imp",
                            data = ipw_male_df, 
                            weights = "sw_overall_trunc", 
                            thresholds = c(m = .1),
                            stars="std") +
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
ggsave(file.path(output_dir, "msm/love_plot_male_ob.png"), width = 8.5, height = 7, dpi = 300)

