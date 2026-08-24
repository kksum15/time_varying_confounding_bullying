##########################################################################################
# Date script created: Fri Apr 18th 2025
# Descriptives
# Dataset: file.path(data_dir, "raw_data/alspac_derived_bullying_ob_haven.RDS") #15236 aln; 22 variables
#############################################################################################

library(summarytools)
library(gtsummary)
library(readxl)
library(writexl)
library(tidyverse)
library(tableone)

comb_data_final <- readRDS(file = file.path(data_dir, "processed_data/comb_data_final.RDS"))
names(comb_data_final_sub)

final_vars_is <- read_xlsx(file.path(data_dir, "raw_data/data_dic_bullying_is_short.xlsx"), sheet="data_final1") %>%
  pull(variables)
comb_vars_bmi <- c("c0_food_insecurity", "c2_food_insecurity", "c0_bmi_z", "c1_bmi_z", "c2_bmi_z", "y17_bmi_z", "overweight")
final_vars_comb <- setdiff(c(final_vars_is, comb_vars_bmi), c("c0_sd_bmi", "c1_sd_bmi", "c2_sd_bmi"))
var_bin <- c("e0_victimization2", "e1_victimization2", "e2_victimization2", "depression_icd_10", "anxiety_cisr", "dep_anx", "overweight", 
             "c0_m_marital", "c0_m_mh1", "c1_m_mh", "c2_m_mh", "c0_f_mh1", "c1_f_mh", "c2_f_mh", "c0_child_maltx1", "c1_child_maltx", 
             "c2_child_maltx", "c0_sib_aggression_bin", "c0_sib_victimization_bin", "c2_sib_aggression", "c2_sib_victimization", 
             "c0_dv1", "c1_dv", "c2_dv2", "c0_food_insecurity", "c2_food_insecurity")
var_cont <- c("c0_bw_g", "c0_m_age_delivery", "c0_sdq_total", "c1_sdq_total", "c2_sdq_total", "c0_n_siblings3", "c0_pub_f_breast", 
              "c0_pub_f_hair", "c1_pub_f_breast", "c1_pub_f_hair", "c2_pub_f_breast", "c2_pub_f_hair", "c0_pub_m_genitalia", "c0_pub_m_hair",
              "c1_pub_m_genitalia", "c1_pub_m_hair", "c2_pub_m_genitalia", "c2_pub_m_hair", "c0_bmi_z", "c1_bmi_z", "c2_bmi_z")
var_cat <- c("c0_m_parity", "c0_medu", "c0_fedu", "c0_hhincome")
var_nnorm <- c("c0_sdq_total", "c1_sdq_total", "c2_sdq_total", "c0_n_siblings3", "c0_pub_f_breast", 
               "c0_pub_f_hair", "c1_pub_f_breast", "c1_pub_f_hair", "c2_pub_f_breast", "c2_pub_f_hair", "c0_pub_m_genitalia", "c0_pub_m_hair",
               "c1_pub_m_genitalia", "c1_pub_m_hair", "c2_pub_m_genitalia", "c2_pub_m_hair")


############################################
# OBSERVED SAMPLE (8808)
#################################

comb_data_final_sub <- comb_data_final %>% 
  filter(!is.na(f8003b) | !is.na(fd003b) | !is.na(ff0011b)) %>%
  select(aln, qlet, alnqlet, final_vars_comb, smfq_13_sum_prorated_scores, smfq_16_sum_prorated_scores) 

tab_obs_comb <- CreateTableOne(vars = final_vars_comb, # set descriptive variables
                               data = comb_data_final_sub, # baseline
                               strata = "c0_female",
                               addOverall = T,
                               #includeNA = T,
                               factorVars = c(var_bin, var_cat)) # define categorical variables

tab_obs_comb <- print(tab_obs_comb,
                      nonnormal = var_nnorm,
                      formatOptions = list(big.mark = ","),
                      showAllLevels = TRUE,
                      test = FALSE,
                      quote = FALSE, 
                      noSpaces = TRUE, 
                      printToggle = FALSE, 
                      missing = T)

write.csv(tab_obs_comb, file = file.path(output_dir, "desc/tab_obs_comb.csv"))

tab_obs_comb_female <- CreateTableOne(vars = final_vars_comb, # set descriptive variables
                               data = comb_data_final_sub[comb_data_final_sub$c0_female == "Female",],
                               strata = "e0_victimization2",
                               #includeNA = T,
                               factorVars = c(var_bin, var_cat)) # define categorical variables

tab_obs_comb_female <- print(tab_obs_comb_female,
                      nonnormal = var_nnorm,
                      formatOptions = list(big.mark = ","),
                      showAllLevels = TRUE,
                      test = FALSE,
                      quote = FALSE, 
                      noSpaces = TRUE, 
                      printToggle = FALSE, 
                      missing = T)

write.csv(tab_obs_comb_female, file = file.path(output_dir, "desc/tab_obs_comb_female_nomiss.csv"))

tab_obs_comb_male <- CreateTableOne(vars = final_vars_comb, # set descriptive variables
                                      data = comb_data_final_sub[comb_data_final_sub$c0_female == "Male",],
                                      strata = "e0_victimization2",
                                      #includeNA = T,
                                      factorVars = c(var_bin, var_cat)) 

tab_obs_comb_male <- print(tab_obs_comb_male,
                             nonnormal = var_nnorm,
                             formatOptions = list(big.mark = ","),
                             showAllLevels = TRUE,
                             test = FALSE,
                             quote = FALSE, 
                             noSpaces = TRUE, 
                             printToggle = FALSE, 
                             missing = T)

write.csv(tab_obs_comb_male, file = file.path(output_dir, "desc/tab_obs_comb_male_nomiss.csv"))

vars_to_summarise <- setdiff(final_vars_comb, "c0_female")

# Overall missing
obs_miss_overall <- comb_data_final_sub %>%
  select(all_of(vars_to_summarise)) %>%
  summarise(across(everything(), 
                   list(n = ~sum(is.na(.)),
                        pct = ~round(mean(is.na(.)) * 100, 1)))) %>%
  pivot_longer(everything(),
               names_to = c("var", ".value"),
               names_sep = "_(?=n$|pct$)") %>%
  mutate(overall_missing = paste0(n, " (", pct, ")")) %>%
  select(var, overall_missing)

# Stratified missing by c0_female
obs_miss_strat <- comb_data_final_sub %>%
  group_by(c0_female) %>%
  summarise(across(all_of(vars_to_summarise), 
                   list(n = ~sum(is.na(.)),
                        pct = ~round(mean(is.na(.)) * 100, 1)))) %>%
  pivot_longer(-c0_female,
               names_to = c("var", ".value"),
               names_sep = "_(?=n$|pct$)") %>%
  mutate(missing_str = paste0(n, " (", pct, ")"),
         group = paste0("missing_", c0_female)) %>%
  select(var, group, missing_str) %>%
  pivot_wider(names_from = group, values_from = missing_str)

# Combine
obs_miss_all <- obs_miss_overall %>%
  left_join(obs_miss_strat, by = "var")

write_xlsx(obs_miss_all, file.path(output_dir, "desc/obs_miss_all.xlsx"))

############################################
# ALIVE BY ONE YEAR SAMPLE (14859)
#################################

comb_data_alive <- comb_data_final %>% 
  filter(kz011b == 1)

tab_alive_comb <- CreateTableOne(vars = final_vars_comb, 
                                 data = comb_data_alive, 
                                 strata = "c0_female",
                                 addOverall = T,
                                 #includeNA = T,
                                 factorVars = c(var_bin, var_cat)) 

tab_alive_comb <- print(tab_alive_comb,
                        nonnormal = var_nnorm,
                        formatOptions = list(big.mark = ","),
                        showAllLevels = TRUE,
                        test = FALSE,
                        quote = FALSE, 
                        noSpaces = TRUE, 
                        printToggle = FALSE)

write.csv(tab_alive_comb, file = file.path(output_dir, "desc/tab_alive_comb.csv"))

# Overall missing
alive_miss_overall <- comb_data_alive %>%
  select(all_of(vars_to_summarise )) %>%
  summarise(across(everything(), 
                   list(n = ~sum(is.na(.)),
                        pct = ~round(mean(is.na(.)) * 100, 1)))) %>%
  pivot_longer(everything(),
               names_to = c("var", ".value"),
               names_sep = "_(?=n$|pct$)") %>%
  mutate(overall_missing = paste0(n, " (", pct, ")")) %>%
  select(var, overall_missing)

# Stratified missing by c0_female
alive_miss_strat <- comb_data_alive %>%
  group_by(c0_female) %>%
  summarise(across(all_of(vars_to_summarise), 
                   list(n = ~sum(is.na(.)),
                        pct = ~round(mean(is.na(.)) * 100, 1)))) %>%
  pivot_longer(-c0_female,
               names_to = c("var", ".value"),
               names_sep = "_(?=n$|pct$)") %>%
  mutate(missing_str = paste0(n, " (", pct, ")"),
         group = paste0("missing_", c0_female)) %>%
  select(var, group, missing_str) %>%
  pivot_wider(names_from = group, values_from = missing_str)

# Combine
alive_miss_all <- alive_miss_overall %>%
  left_join(alive_miss_strat, by = "var")

write_xlsx(alive_miss_all, file.path(output_dir, "desc/alive_miss_all.xlsx"))
