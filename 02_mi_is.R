##########################################################################################
# Date script created: 23 Apr 2025
# Account for missing data: Multiple imputation 
# Bullying-internalizing symptoms
# Dataset: file.path(data_dir, "processed_data/my_data_final1.RDS") #8808 records; 53 variables
#############################################################################################

####################
# Load packages
####################

library(haven)
library(tidyverse)
library(summarytools)
library(mice)
library(VIM)
library(ggmice)
library(readxl)
library(naniar)
library(parallel)
library(gtsummary)
library(writexl)
library(janitor)
library(lsr)
library(patchwork)
library(tableone)

####################
# File directory
####################

# data_dir <- "/Volumes/005/working/data"
# output_dir <- "/Volumes/005/working/results

####################
# Import dataset 
####################

my_data_final1 <- readRDS(file.path(data_dir, "processed_data/my_data_final1.RDS")) %>%
  mutate(c1_sib_aggression = c0_sib_aggression_bin,
         c1_sib_victimization = c0_sib_victimization_bin)

#########################################
# Explore the patterns of missing data 
########################################

final_vars <- read_xlsx(file.path(data_dir, "raw_data/data_dic_bullying_is_short.xlsx"), sheet="data_final1") %>%
  pull(variables)
final_vars <- c(final_vars, "c1_sib_aggression", "c1_sib_victimization")

final_vars_female <- final_vars[!grepl("pub_m", final_vars) & final_vars != "c0_female"]
final_vars_male <- final_vars[!grepl("pub_f", final_vars) & final_vars != "c0_female"]
final_vars_female_aux <- c(final_vars_female, "smfq_13_sum_prorated_scores", "smfq_16_sum_prorated_scores")
final_vars_male_aux <- c(final_vars_male, "smfq_13_sum_prorated_scores", "smfq_16_sum_prorated_scores")                                             

female <- my_data_final1 %>% 
  filter(c0_female == "Female") %>%
  select("aln", "qlet", final_vars_female_aux) %>%
  mutate(qlet = as.factor(qlet),
         e0xe1 = NA,
         e1xe2 = NA, 
         e0xe2 = NA,
         e0xe1xe2 = NA) #4408

male <- my_data_final1 %>% 
  filter(c0_female == "Male") %>%
  select("aln", "qlet",final_vars_male_aux) %>%
  mutate(qlet = as.factor(qlet),
         e0xe1 = NA,
         e1xe2 = NA, 
         e0xe2 = NA,
         e0xe1xe2 = NA) #4400

#########################################
# Explore the mice function
#########################################
# continuous values: predictive mean matching (pmm)
# binary outcomes (factors with two levels): logistic regression imputation (logreg)
# unordered categorical data (factors with more than two levels): multinomial logistic regression imputation (polyreg) 
# ordered categorical data (ordered factors with more than two levels): ordered logistic regression (polr)
names(female)
testimp_female <- mice(data=female, maxit=0, printFlag = F)
testimp_female
testimp_female$formulas
testimp_female$predictorMatrix #rows = variables to be imputed; columns = variables used in the imputation; 0 = variables excluded from imputation model

testimp_male <- mice(data=male, maxit=0, printFlag = F)
testimp_male

# Change the methods for 
# (1) normally distributed variables from pmm to norm (linear regression)
# (2) ordered variables from polyreg to polr 
# (3) dep_anx, depression_icd_10, anxiety_cisr to logreg
# (4) passive imputation for dep_anx and c1_sib_aggression, c1_sib_victimization, and all the intreaction terms 
method_female <- testimp_female$method
method_female[c("c0_bw_g", "c0_m_age_delivery", "c0_sd_bmi", "c1_sd_bmi", "c2_sd_bmi")] <- "norm"
method_female[c("c0_m_parity", "c0_medu", "c0_fedu", "c0_hhincome")] <- "polr"
method_female["dep_anx"] <- "~ I(depression_icd_10 == 'Yes' | anxiety_cisr == 'Yes')" #passive imputation
method_female["c1_sib_aggression"] <- "~ c0_sib_aggression_bin"
method_female["c1_sib_victimization"] <- "~ c0_sib_victimization_bin"
method_female["e0xe1"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e1_victimization2))"
method_female["e1xe2"] <- "~ I(as.numeric(e1_victimization2)*as.numeric(e2_victimization2))"
method_female["e0xe2"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e2_victimization2))"
method_female["e0xe1xe2"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e1_victimization2)*as.numeric(e2_victimization2))"
print(method_female)

method_male <- testimp_male$method
method_male[c("c0_bw_g", "c0_m_age_delivery", "c0_sd_bmi", "c1_sd_bmi", "c2_sd_bmi")] <- "norm"
method_male[c("c0_m_parity", "c0_medu", "c0_fedu", "c0_hhincome")] <- "polr"
method_male["dep_anx"] <- "~ I(depression_icd_10 == 'Yes' | anxiety_cisr == 'Yes')" #passive imputation
method_male["c1_sib_aggression"] <- "~ c0_sib_aggression_bin"
method_male["c1_sib_victimization"] <- "~ c0_sib_victimization_bin"
method_male["e0xe1"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e1_victimization2))"
method_male["e1xe2"] <- "~ I(as.numeric(e1_victimization2)*as.numeric(e2_victimization2))"
method_male["e0xe2"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e2_victimization2))"
method_male["e0xe1xe2"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e1_victimization2)*as.numeric(e2_victimization2))"
print(method_male)

# nothing predicts aln, qlet, dep_anx, c1_sib_aggression, and c1_sib_victimization (remove all)
predictor_female <- testimp_female$predictorMatrix
predictor_male <- testimp_male$predictorMatrix

predictor_female["aln",] <- 0  
predictor_female["qlet",] <- 0
predictor_male["aln",] <- 0  
predictor_male["qlet",] <- 0
predictor_female["c1_sib_aggression",] <- 0  # done automatically 
predictor_female["c1_sib_victimization",] <- 0
predictor_male["c1_sib_aggression",] <- 0  
predictor_male["c1_sib_victimization",] <- 0
predictor_female["dep_anx",] <- 0
predictor_male["dep_anx",] <- 0

# omit c1_sib_aggression, and c1_sib_victimization from predicting other variables to avoid collinearity;
predictor_female[, "c1_sib_aggression" ] <- 0  
predictor_female[, "c1_sib_victimization"] <- 0
predictor_male[, "c1_sib_aggression"] <- 0  
predictor_male[, "c1_sib_victimization"] <- 0
predictor_female[, "aln" ] <- 0  
predictor_female[, "qlet"] <- 0
predictor_male[, "aln"] <- 0  
predictor_male[, "qlet"] <- 0

# omit dep_anx from predicting depression_icd_10, anxiety_cisr
predictor_female[c("depression_icd_10", "anxiety_cisr"), "dep_anx"] <- 0
predictor_male[c("depression_icd_10", "anxiety_cisr"), "dep_anx"] <- 0

# allow to predict other variables but omit e1_victimization2 and e2_victimization2 from predicting e1_e2
predictor_female[, "e0xe1"] <- 1
predictor_female[, "e1xe2"] <- 1
predictor_female[, "e0xe2"] <- 1
predictor_female[, "e0xe1xe2"] <- 1

predictor_male[, "e0xe1"] <- 1
predictor_male[, "e1xe2"] <- 1
predictor_male[, "e0xe2"] <- 1
predictor_male[, "e0xe1xe2"] <- 1

predictor_female[c("e0_victimization2", "e1_victimization2"), "e0xe1"] <- 0
predictor_female[c("e1_victimization2", "e2_victimization2"), "e1xe2"] <- 0
predictor_female[c("e0_victimization2", "e2_victimization2"), "e0xe2"] <- 0
predictor_female[c("e0_victimization2", "e1_victimization2", "e2_victimization2"), "e0xe1xe2"] <- 0

predictor_male[c("e0_victimization2", "e1_victimization2"), "e0xe1"] <- 0
predictor_male[c("e1_victimization2", "e2_victimization2"), "e1xe2"] <- 0
predictor_male[c("e0_victimization2", "e2_victimization2"), "e0xe2"] <- 0
predictor_male[c("e0_victimization2", "e1_victimization2", "e2_victimization2"), "e0xe1xe2"] <- 0

predictor_female["e0xe1", ] <- 0
predictor_female["e1xe2",] <- 0
predictor_female["e0xe2",] <- 0
predictor_female["e0xe1xe2",] <- 0
predictor_male["e0xe1", ] <- 0
predictor_male["e1xe2",] <- 0
predictor_male["e0xe2",] <- 0
predictor_male["e0xe1xe2",] <- 0

# omit from imputation model for c0_sd_bmi (only allows c2_sd_bmi to impute for c0_sd_bmi)
predictor_female["c0_sd_bmi","c1_sd_bmi"] <- 0
predictor_male["c0_sd_bmi","c1_sd_bmi"] <- 0

# omit from imputation model for c1_sd_bmi (only c0_sd_bmi to impute c1_sd_bmi)
predictor_female["c1_sd_bmi","c2_sd_bmi"] <- 0
predictor_male["c1_sd_bmi","c2_sd_bmi"] <- 0

# omit from imputation model for c2_sd_bmi (only c1_sd_bmi to impute c2_sd_bmi)
predictor_female["c2_sd_bmi","c0_sd_bmi"] <- 0
predictor_male["c2_sd_bmi","c0_sd_bmi"] <- 0

options(max.print=10000)
print(predictor_female)
print(predictor_male)

####################
# Run imputation 
###################

# 1. Create one imputation set first with 100 iterations - assess the number of iterations required for convergence (mean and SD becomes stable) using
# trace plots for each partially observed variable

# final_model 
############################
start <- Sys.time()
imp_m1_female_single_m100_3int_072025 <-  mice(data=female,
                                               method = method_female, 
                                               predictorMatrix = predictor_female,
                                               m = 1, maxit = 100, printFlag = FALSE, seed = 1045) 
end <- Sys.time()
end - start #7.671074 mins
save(imp_m1_female_single_m100_3int_072025 , file = file.path(data_dir,"processed_data/imp_m1_female_single_m100_3int_072025.RData"))
load(file = file.path(data_dir,"processed_data/imp_m1_female_single_m100_3int_072025.RData"))

start <- Sys.time()
imp_m1_male_single_m100_3int_072025 <-  mice(data=male,
                                             method = method_male, 
                                             predictorMatrix = predictor_male,
                                             m = 1, maxit = 100, printFlag = FALSE, seed = 1045)
end <- Sys.time()
end - start #7.604345 mins
save(imp_m1_male_single_m100_3int_072025, file = file.path(data_dir,"processed_data/imp_m1_male_single_m100_3int_072025.RData"))
load(file = file.path(data_dir,"processed_data/multiple_imputation/test2/imp_m1_male_single_m100_3int_072025.RData")) 

##################
# 200 iterations
###################

start <- Sys.time()
imp_m1_female_single_m200_3int_072025 <-  mice(data=female,
                                               method = method_female, 
                                               predictorMatrix = predictor_female,
                                               m = 1, maxit = 200, printFlag = FALSE, seed = 1045) 
end <- Sys.time()
end - start #15.34936 mins
save(imp_m1_female_single_m200_3int_072025 , file = file.path(data_dir,"processed_data/imp_m1_female_single_m200_3int_072025.RData"))
load(file = file.path(data_dir,"processed_data/imp_m1_female_single_m200_3int_072025.RData"))

start <- Sys.time()
imp_m1_male_single_m200_3int_072025 <-  mice(data=male,
                                             method = method_male, 
                                             predictorMatrix = predictor_male,
                                             m = 1, maxit = 200, printFlag = FALSE, seed = 1045)
end <- Sys.time()
end - start #15.16666 mins
save(imp_m1_male_single_m200_3int_072025, file = file.path(data_dir,"processed_data/imp_m1_male_single_m200_3int_072025.RData"))
load(file = file.path(data_dir,"processed_data/imp_m1_male_single_m200_3int_072025.RData")) 

############################
# Produce a trace plot for each partially observed variable 
############################
#vars_imputed_female <- final_vars_female[!final_vars_female %in% c("depression_icd_10", "anxiety_cisr", "dep_anx", "c1_n_siblings3",
#"c2_n_siblings3", "c0_adiposity", "c1_m_marital", "c2_m_marital")]
#vars_imputed_female <- c(vars_imputed_female, "f7_adiposity")

vars_imputed_female <- c(final_vars_female_aux, "e0xe1", "e1xe2", "e0xe2", "e0xe1xe2")

for (var in vars_imputed_female) {
  formula <- as.formula(paste(var, "~ .it | .ms"))
  #file_name <- file.path(output_dir, "multiple_imputation/test/imputation_m1_100_3int_072025", paste0("female_plot_", var, ".png"))
  file_name <- file.path(output_dir, "multiple_imputation/test/imputation_m1_200_3int_072025", paste0("female_plot_", var, ".png"))
  png(filename = file_name, width = 800, height = 600)
  #plot <- plot(imp_m1_female_single_m100_3int_072025, formula, layout = c(2, 1), main = var)
  plot <- plot(imp_m1_female_single_m200_3int_072025, formula, layout = c(2, 1), main = var)
  print(plot)
  dev.off()
}

vars_imputed_male <- c(final_vars_male_aux, "e0xe1", "e1xe2", "e0xe2", "e0xe1xe2")

for (var in vars_imputed_male) {
  formula <- as.formula(paste(var, "~ .it | .ms"))
  #file_name <- file.path(output_dir, "multiple_imputation/test/imputation_m1_100_3int_072025", paste0("male_plot_", var, ".png"))
  file_name <- file.path(output_dir, "multiple_imputation/test/imputation_m1_200_3int_072025", paste0("male_plot_", var, ".png"))
  png(filename = file_name, width = 800, height = 600)
  #plot <- plot(imp_m1_male_single_m100_3int_072025, formula, layout = c(2, 1), main = var)
  plot <- plot(imp_m1_male_single_m200_3int_072025, formula, layout = c(2, 1), main = var)
  print(plot)
  dev.off()
}

# 2. Run final imputations with the number of iterations determined with 100 imputation sets
# Only 6.96% of females and 4.66% of males are complete cases, so technically minimum imputation sets should be around 95; 
# rule of thumb: number of imputation sets should be at least as many as % of missing cases but Jon usually does 100 minimum 
# try to run imputation in parallel using multiple cores - speed-up mainly affects the number of imputations not iterations 
##################################################################################################################

detectCores() #10

# Female 
###########

# Set up parameters
total_imputations <- 100 # Total number of imputations you want
cores_2_use <- detectCores() - 2 # Leave 2 cores free
imps_per_core <- ceiling(total_imputations / cores_2_use) # Compute how many imputations each core should do

# Create cluster
cl <- makeCluster(cores_2_use)
clusterSetRNGStream(cl, 9700) # Ensure reproducibility by setting seed
clusterExport(cl, c("female", "method_female", "predictor_female", "imps_per_core")) # Export required objects to each worker
clusterEvalQ(cl, library(mice)) # Load required package on each worker

start <- Sys.time()
# Run imputations in parallel
imp_pars_female <- 
  parLapply(cl = cl, X = 1:cores_2_use, fun = function(no){
    mice(female, 
         method = method_female, 
         predictorMatrix = predictor_female,
         m = imps_per_core, maxit = 50,
         printFlag = FALSE)
  })

end <- Sys.time()
end - start #55.70214 mins

# Stop the cluster
stopCluster(cl)

# Combine mids objects into a single object
imp_merged_female <- imp_pars_female[[1]]
for (n in 2:length(imp_pars_female)){
  imp_merged_female <- 
    ibind(imp_merged_female,
          imp_pars_female[[n]])
}

# Trim to exactly the number of imputations desired
if (imp_merged_female$m > total_imputations) {
  imp_merged_female$imp <- lapply(imp_merged_female$imp, function(x) x[, 1:total_imputations])
  imp_merged_female$m <- total_imputations
}

# Confirm
print(paste("Final number of imputations:", imp_merged_female$m))

summary(imp_merged_female)

save(imp_merged_female, file = file.path(data_dir,"processed_data/imp_merged_female_3int_072025.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_female_3int_072025.RData"))

# Male 
###########
# Create cluster

cl <- makeCluster(cores_2_use)
clusterSetRNGStream(cl, 9896) # Ensure reproducibility by setting seed
clusterExport(cl, c("male", "method_male", "predictor_male", "imps_per_core")) # Export required objects to each worker
clusterEvalQ(cl, library(mice)) # Load required package on each worker

start <- Sys.time()
# Run imputations in parallel
imp_pars_male <- 
  parLapply(cl = cl, X = 1:cores_2_use, fun = function(no){
    mice(male, 
         method = method_male, 
         predictorMatrix = predictor_male,
         m = imps_per_core, maxit = 50,
         printFlag = FALSE)
  })
end <- Sys.time()
end - start # 1.533272 hours

# Stop the cluster
stopCluster(cl)

# Combine mids objects into a single object
imp_merged_male <- imp_pars_male[[1]]
for (n in 2:length(imp_pars_male)){
  imp_merged_male <- 
    ibind(imp_merged_male,
          imp_pars_male[[n]])
}

# Trim to exactly the number of imputations desired
if (imp_merged_male$m > total_imputations) {
  imp_merged_male$imp <- lapply(imp_merged_male$imp, function(x) x[, 1:total_imputations])
  imp_merged_male$m <- total_imputations
}

# Confirm
print(paste("Final number of imputations:", imp_merged_male$m))

summary(imp_merged_male)

save(imp_merged_male, file = file.path(data_dir,"processed_data/imp_merged_male_3int_072025.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_male_3int_072025.RData"))

# Combine both imputed dataset and add sibling victimization and aggression at t1 to the dataset 
#############################################################################################

# Set number of imputations (assuming both have same m)
m <- imp_merged_female$m

# Get all variable names across both datasets
vars_female <- names(complete(imp_merged_female, 1))
vars_male <- names(complete(imp_merged_male, 1))
all_vars <- union(vars_female, vars_male)

# Add the new indicator column name
all_vars <- union(all_vars, "female")

# Function to align columns and add binary female indicator
align_vars <- function(data, all_vars, female_value) {
  data$female <- female_value  # Add binary indicator (1 = female, 0 = male)
  for (var in setdiff(all_vars, names(data))) {
    data[[var]] <- NA  # Add missing columns as NA
  }
  data <- data[ , all_vars]  # Reorder columns
  return(data)
}

# ---------------------
# 1. Create the `.imp == 0` data (original, un-imputed)
# ---------------------
original_female <- align_vars(imp_merged_female$data, all_vars, 1)
original_male   <- align_vars(imp_merged_male$data, all_vars, 0)

original_combined <- rbind(original_female, original_male)
original_combined$.imp <- 0
original_combined$.id <- 1:nrow(original_combined)

# ---------------------
# 2. Create completed datasets for each imputation
# ---------------------

long_list <- list()
long_list[[1]] <- original_combined  # first entry is always .imp == 0

for (i in 1:m) {
  complete_female <- align_vars(complete(imp_merged_female, i), all_vars, 1)
  complete_male   <- align_vars(complete(imp_merged_male, i), all_vars, 0)
  
  complete_combined <- rbind(complete_female, complete_male)
  complete_combined$.imp <- i
  complete_combined$.id <- 1:nrow(complete_combined)
  
  long_list[[i + 1]] <- complete_combined
}

# ---------------------
# 3. Combine to long-format data frame and convert
# ---------------------
long_data <- do.call(rbind, long_list)

# Convert to mids
imp_merged <- as.mids(long_data)

save(imp_merged, file = file.path(data_dir,"processed_data/imp_merged_3int_072025.RData"))

################################################################
# Diagnostics 

load(file = file.path(data_dir,"processed_data/imp_merged.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_female.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_male.RData"))

# To access the first imputed dataset
imp_merged_female_m1 <- complete(imp_merged_female,1)
imp_merged_male_m1 <- complete(imp_merged_male,1)
imp_merged_m1 <- complete(imp_merged,1)

# Listing the first few rows
head(imp_merged_female_m1)
head(imp_merged_male_m1)
head(imp_merged_m1)

# summary
summary(imp_merged_female_m1)
summary(imp_merged_male_m1)
summary(imp_merged_m1)

imp_merged_m1 %>%
  dplyr::count(depression_icd_10, anxiety_cisr, dep_anx)

imp_merged_m1 %>%
  dplyr::count(c0_sib_aggression_bin, c1_sib_aggression)

imp_merged_m1 %>%
  dplyr::count(c0_sib_victimization_bin, c1_sib_victimization)

imp_merged_m1 %>%
  dplyr::count(e0_victimization2, e1_victimization2, e0xe1)

imp_merged_m1 %>%
  dplyr::count(e1_victimization2, e2_victimization2, e1xe2)

imp_merged_m1 %>%
  dplyr::count(e0_victimization2, e2_victimization2, e0xe2)

imp_merged_m1 %>%
  dplyr::count(e0_victimization2, e1_victimization2, e2_victimization2, e0xe1xe2)

var_bin <- c("e0_victimization2", "e1_victimization2", "e2_victimization2", "depression_icd_10", "anxiety_cisr", "dep_anx", 
             "c0_m_marital", "c0_m_mh1", "c1_m_mh", "c2_m_mh", "c0_f_mh1", "c1_f_mh", "c2_f_mh", "c0_child_maltx1", "c1_child_maltx", 
             "c2_child_maltx", "c0_sib_aggression_bin", "c0_sib_victimization_bin", "c2_sib_aggression", "c2_sib_victimization", 
             "c0_dv1", "c1_dv", "c2_dv2")

var_cont <- c("c0_bw_g", "c0_m_age_delivery", "c0_sdq_total", "c1_sdq_total", "c2_sdq_total", "c0_n_siblings3", "c0_pub_f_breast", 
              "c0_pub_f_hair", "c1_pub_f_breast", "c1_pub_f_hair", "c2_pub_f_breast", "c2_pub_f_hair", "c0_pub_m_genitalia", "c0_pub_m_hair",
              "c1_pub_m_genitalia", "c1_pub_m_hair", "c2_pub_m_genitalia", "c2_pub_m_hair", "c0_sd_bmi", "c1_sd_bmi", "c2_sd_bmi")

var_cat <- c("c0_m_parity", "c0_medu", "c0_fedu", "c0_hhincome")

var_nnorm <- c("c0_sdq_total", "c1_sdq_total", "c2_sdq_total", "c0_n_siblings3", "c0_pub_f_breast", 
               "c0_pub_f_hair", "c1_pub_f_breast", "c1_pub_f_hair", "c2_pub_f_breast", "c2_pub_f_hair", "c0_pub_m_genitalia", "c0_pub_m_hair",
               "c1_pub_m_genitalia", "c1_pub_m_hair", "c2_pub_m_genitalia", "c2_pub_m_hair")

# boxplot for continuous 
########################

boxplot_dir <- file.path(output_dir, "multiple_imputation/mi_diagnostics/boxplot_3int_072025")

for (var in var_cont) {
  p <- ggmice(imp_merged, aes(x = .imp, y = .data[[var]])) +
    geom_boxplot() +
    scale_x_discrete(drop = FALSE) +
    labs(
      x = "Imputation number",
      y = var,
      title = paste("Boxplot of", var, "by Imputation")
    )
  
  ggsave(
    filename = file.path(boxplot_dir, paste0("boxplot_", var, ".png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# stripp plot continuous 
########################

stripp_dir <- file.path(output_dir, "multiple_imputation/mi_diagnostics/stripp_plot_3int_072025")

for (var in var_cont) {
  p <- ggmice(imp_merged, aes(x = .imp, y = .data[[var]])) +
    geom_jitter() +
    labs(
      x = "Imputation number",
      y = var,
      title = paste("Stripp plot of", var, "by Imputation")
    )
  
  ggsave(
    filename = file.path(stripp_dir, paste0("stripp_plot_", var, ".png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# density plot for continuous 
########################

density_dir <- file.path(output_dir, "multiple_imputation/mi_diagnostics/density_plot_3int_072025")

for (var in var_cont) {
  p <- ggmice(imp_merged, aes(x = .data[[var]], group = .imp)) +
    geom_density() +
    labs(x = var, 
         title = paste("Density plot of", var, "by Imputation")
    )
  
  ggsave(
    filename = file.path(density_dir, paste0("density_plot_", var, ".png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# bar plots for categorical / binary 
source("/Volumes/005/working/scripts/bullying_is/function_for_prop_plots.R")

save_prop_plots(imp_merged, output_dir = output_dir)
save_prop_plots(imp_merged_female, output_dir = output_dir, gender_label ="female")
save_prop_plots(imp_merged_male, output_dir = output_dir, gender_label ="male")

# Create summary statistics for observed and imputed values 
# Observed 
tab_obs <- CreateTableOne(vars = final_vars, # set descriptive variables
                          data = my_data_final1, # baseline
                          factorVars = c(var_bin, var_cat)) # define categorical variables

tab_obs <- print(tab_obs,
                 nonnormal = var_nnorm,
                 formatOptions = list(big.mark = ","),
                 showAllLevels = TRUE,
                 test = FALSE,
                 quote = FALSE, 
                 noSpaces = TRUE, 
                 printToggle = FALSE)

write.csv(tab_obs, file = file.path(output_dir, "desc/tab_obs.csv"))

tab_obs_sex <- CreateTableOne(vars = final_vars, # set descriptive variables
                          data = my_data_final1, # baseline
                          strata = "c0_female",
                          includeNA = T,
                          factorVars = c(var_bin, var_cat)) # define categorical variables

tab_obs_sex <- print(tab_obs_sex,
                 nonnormal = var_nnorm,
                 formatOptions = list(big.mark = ","),
                 showAllLevels = TRUE,
                 test = FALSE,
                 quote = FALSE, 
                 noSpaces = TRUE, 
                 printToggle = FALSE)

write.csv(tab_obs_sex, file = file.path(output_dir, "desc/tab_obs_sex.csv"))

# Sex specific stratified by victimization status 
tab_obs_female <- CreateTableOne(vars = final_vars_female, # set descriptive variables
                                 strata = "e0_victimization2",
                                 data = female, 
                                 #includeNA = T,
                                 factorVars = c(var_bin, var_cat)) # define categorical variables

tab_obs_female <- print(tab_obs_female,
                 nonnormal = var_nnorm,
                 formatOptions = list(big.mark = ","),
                 showAllLevels = TRUE,
                 test = FALSE,
                 quote = FALSE, 
                 noSpaces = TRUE, 
                 printToggle = FALSE)

summary(tab_obs_female)

write.csv(tab_obs_female, file = file.path(output_dir, "desc/tab_obs_female.csv"))
write.csv(tab_obs_female, file = file.path(output_dir, "desc/tab_obs_female_nomiss.csv"))

tab_obs_male <- CreateTableOne(vars = final_vars_male, # set descriptive variables
                               strata = "e0_victimization2",
                               data = male, 
                               #includeNA = T,
                               factorVars = c(var_bin, var_cat)) # define categorical variables

tab_obs_male <- print(tab_obs_male,
                      nonnormal = var_nnorm,
                      formatOptions = list(big.mark = ","),
                      showAllLevels = TRUE,
                      test = FALSE,
                      quote = FALSE, 
                      noSpaces = TRUE, 
                      printToggle = FALSE)

summary(tab_obs_male)

write.csv(tab_obs_male, file = file.path(output_dir, "desc/tab_obs_male.csv"))
write.csv(tab_obs_male, file = file.path(output_dir, "desc/tab_obs_male_nomiss.csv"))

# count missing for continuous variables 
female_cont_miss <- female %>%
  group_by(e0_victimization2) %>%
  summarise(across(all_of(var_cont),
                   list(
                     missing_n = ~ sum(is.na(.)),
                     missing_pct = ~ round(100 * mean(is.na(.)), 1)
                   ),
                   .names = "{.col}_{.fn}"
  ))

male_cont_miss <- male %>%
  group_by(e0_victimization2) %>%
  summarise(across(all_of(var_cont),
                   list(
                     missing_n = ~ sum(is.na(.)),
                     missing_pct = ~ round(100 * mean(is.na(.)), 1)
                   ),
                   .names = "{.col}_{.fn}"
  ))

# Imputed 
imp_merged_long <- complete(imp_merged,"long",include = F)
tab_imputed <- CreateTableOne(vars = final_vars, 
                              data = imp_merged_long , 
                              factorVars = c(var_bin, var_cat))

tab_imputed <- print(tab_imputed,
                     nonnormal = var_nnorm,
                     formatOptions = list(big.mark = ","),
                     showAllLevels = TRUE,
                     test = FALSE,
                     quote = FALSE, 
                     noSpaces = TRUE, 
                     printToggle = FALSE)

write.csv(tab_imputed, file = file.path(output_dir, "desc/tab_imputed_3int_072025.csv"))

#imp_merged_long <- complete(imp_merged,"long",include = F)
tab_imputed_sex <- CreateTableOne(vars = final_vars, 
                              data = imp_merged_m1,
                              strata = "female",
                              addOverall = T,
                              factorVars = c(var_bin, var_cat))

tab_imputed_sex <- print(tab_imputed_sex,
                     nonnormal = var_nnorm,
                     formatOptions = list(big.mark = ","),
                     showAllLevels = TRUE,
                     test = FALSE,
                     quote = FALSE, 
                     noSpaces = TRUE, 
                     printToggle = FALSE)

write.csv(tab_imputed_sex, file = file.path(output_dir, "desc/tab_imputed_sex_3int_072025.csv"))

tab_imputed_female <- CreateTableOne(vars = final_vars, 
                                  data = imp_merged_female_m1,
                                  strata = "e0_victimization2",
                                  factorVars = c(var_bin, var_cat))

tab_imputed_female <- print(tab_imputed_female,
                         nonnormal = var_nnorm,
                         formatOptions = list(big.mark = ","),
                         showAllLevels = TRUE,
                         test = FALSE,
                         quote = FALSE, 
                         noSpaces = TRUE, 
                         printToggle = FALSE)

write.csv(tab_imputed_female, file = file.path(output_dir, "desc/tab_imputed_female.csv"))

tab_imputed_male <- CreateTableOne(vars = final_vars, 
                                     data = imp_merged_male_m1,
                                     strata = "e0_victimization2",
                                     factorVars = c(var_bin, var_cat))

tab_imputed_male <- print(tab_imputed_male,
                            nonnormal = var_nnorm,
                            formatOptions = list(big.mark = ","),
                            showAllLevels = TRUE,
                            test = FALSE,
                            quote = FALSE, 
                            noSpaces = TRUE, 
                            printToggle = FALSE)

write.csv(tab_imputed_male, file = file.path(output_dir, "desc/tab_imputed_male.csv"))
