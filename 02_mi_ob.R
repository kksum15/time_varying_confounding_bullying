##########################################################################################
# Date script created Wed Oct 22 15:39:18 2025 ------------------------------
# Account for missing data: Multiple imputation 
# Bullying-overweight 
# Dataset: file.path(data_dir, "processed_data/comb_data_final_sub.RDS") #8808 records; 60 variables
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
# Import dataset 
####################

mi_data <- readRDS(file.path(data_dir, "processed_data/comb_data_final_sub.RDS")) %>%
  mutate(c1_sib_aggression = c0_sib_aggression_bin,
         c1_sib_victimization = c0_sib_victimization_bin, 
         c1_food_insecurity = c0_food_insecurity)

female <- mi_data %>%
  filter(c0_female == "Female") %>%
  select(-matches("pub_m"), -c0_female) %>%
  mutate(qlet = as.factor(qlet),
         e0xe1 = NA,
         e1xe2 = NA, 
         e0xe2 = NA,
         e0xe1xe2 = NA)

male <- mi_data %>%
  filter(c0_female == "Male") %>%
  select(-matches("pub_f"), -c0_female) %>%
  mutate(qlet = as.factor(qlet),
         e0xe1 = NA,
         e1xe2 = NA, 
         e0xe2 = NA,
         e0xe1xe2 = NA)

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
# (4) passive imputation for dep_anx and c1_sib_aggression and c1_sib_victimization
method_female <- testimp_female$method
method_female[c("c0_bw_g", "c0_m_age_delivery", "c0_bmi_z", "c1_bmi_z", "c2_bmi_z", "y17_bmi_z")] <- "norm"
method_female[c("c0_m_parity", "c0_medu", "c0_fedu", "c0_hhincome")] <- "polr"
method_female["overweight"] <- "~ as.numeric(y17_bmi_z >= 1.04)" #passive imputation
method_female["c1_sib_aggression"] <- "~ c0_sib_aggression_bin"
method_female["c1_sib_victimization"] <- "~ c0_sib_victimization_bin"
method_female["c1_food_insecurity"] <- "~ c0_food_insecurity"
method_female["e0xe1"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e1_victimization2))"
method_female["e1xe2"] <- "~ I(as.numeric(e1_victimization2)*as.numeric(e2_victimization2))"
method_female["e0xe2"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e2_victimization2))"
method_female["e0xe1xe2"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e1_victimization2)*as.numeric(e2_victimization2))"
print(method_female)

method_male <- testimp_male$method
method_male[c("c0_bw_g", "c0_m_age_delivery", "c0_bmi_z", "c1_bmi_z", "c2_bmi_z", "y17_bmi_z")] <- "norm"
method_male[c("c0_m_parity", "c0_medu", "c0_fedu", "c0_hhincome")] <- "polr"
method_male["overweight"] <- "~ as.numeric(y17_bmi_z >= 1.04)" #passive imputation
method_male["c1_sib_aggression"] <- "~ c0_sib_aggression_bin"
method_male["c1_sib_victimization"] <- "~ c0_sib_victimization_bin"
method_male["c1_food_insecurity"] <- "~ c0_food_insecurity"
method_male["e0xe1"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e1_victimization2))"
method_male["e1xe2"] <- "~ I(as.numeric(e1_victimization2)*as.numeric(e2_victimization2))"
method_male["e0xe2"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e2_victimization2))"
method_male["e0xe1xe2"] <- "~ I(as.numeric(e0_victimization2)*as.numeric(e1_victimization2)*as.numeric(e2_victimization2))"
print(method_male)

# nothing predicts aln, qlet, overweight, c1_sib_aggression, c1_sib_victimization, c1_food_insecurity (remove all)
predictor_female <- testimp_female$predictorMatrix
predictor_male <- testimp_male$predictorMatrix

predictor_female["aln",] <- 0  
predictor_female["qlet",] <- 0
predictor_female["alnqlet",] <- 0
predictor_male["aln",] <- 0  
predictor_male["qlet",] <- 0
predictor_male["alnqlet",] <- 0
predictor_female["c1_sib_aggression",] <- 0  # done automatically 
predictor_female["c1_sib_victimization",] <- 0
predictor_female["c1_food_insecurity",] <- 0
predictor_male["c1_sib_aggression",] <- 0  
predictor_male["c1_sib_victimization",] <- 0
predictor_male["c1_food_insecurity",] <- 0
predictor_female["overweight",] <- 0
predictor_male["overweight",] <- 0

# omit c1_sib_aggression, c1_sib_victimization, c1_food_insecurity from predicting other variables to avoid collinearity;
predictor_female[, "c1_sib_aggression" ] <- 0  
predictor_female[, "c1_sib_victimization"] <- 0
predictor_female[, "c1_food_insecurity"] <- 0
predictor_male[, "c1_sib_aggression"] <- 0  
predictor_male[, "c1_sib_victimization"] <- 0
predictor_male[, "c1_food_insecurity"] <- 0
predictor_female[, "aln" ] <- 0  
predictor_female[, "qlet"] <- 0
predictor_female[, "alnqlet"] <- 0
predictor_male[, "aln"] <- 0  
predictor_male[, "qlet"] <- 0
predictor_male[, "alnqlet"] <- 0

# omit overweight from y17_bmi_z
predictor_female[c("y17_bmi_z"), "overweight"] <- 0
predictor_male[c("y17_bmi_z"), "overweight"] <- 0

# allow to predict other variables but omit e0_victimization2, e1_victimization2, e2_victimization2 from predicting all interaction terms
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

# omit from imputation model for c0_bmi_z (only allows c2_bmi_z to impute for c0_bmi_z)
predictor_female["c0_bmi_z","c1_bmi_z"] <- 0
predictor_male["c0_bmi_z","c1_bmi_z"] <- 0

# omit from imputation model for c1_bmi_z (only c0_bmi_z to impute c1_bmi_z)
predictor_female["c1_bmi_z","c2_bmi_z"] <- 0
predictor_male["c1_bmi_z","c2_bmi_z"] <- 0

# omit from imputation model for c2_bmi_z (only c1_bmi_z to impute c2_bmi_z)
predictor_female["c2_bmi_z","c0_bmi_z"] <- 0
predictor_male["c2_bmi_z","c0_bmi_z"] <- 0

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
imp_m1_female_single_m100_ob <-  mice(data=female,
                                      method = method_female, 
                                      predictorMatrix = predictor_female,
                                      m = 1, maxit = 100, printFlag = FALSE, seed = 1066) 
end <- Sys.time()
end - start #7.306182 mins
save(imp_m1_female_single_m100_ob, file = file.path(data_dir,"processed_data/imp_m1_female_single_m100_ob.RData"))
load(file = file.path(data_dir,"processed_data/imp_m1_female_single_m100_ob.RData"))

start <- Sys.time()
imp_m1_male_single_m100_ob <-  mice(data=male,
                                    method = method_male, 
                                    predictorMatrix = predictor_male,
                                    m = 1, maxit = 100, printFlag = FALSE, seed = 1066)
end <- Sys.time() # 7.350649 mins
end - start 
save(imp_m1_male_single_m100_ob, file = file.path(data_dir,"processed_data/imp_m1_male_single_m100_ob.RData"))
load(file = file.path(data_dir,"processed_data/multiple_imputation/test2/imp_m1_male_single_m100_ob.RData")) 

##################
# 200 iterations
###################

start <- Sys.time()
imp_m1_female_single_m200_ob<-  mice(data=female,
                                     method = method_female, 
                                     predictorMatrix = predictor_female,
                                     m = 1, maxit = 200, printFlag = FALSE, seed = 1066) 
end <- Sys.time()
end - start #15.91298 mins
save(imp_m1_female_single_m200_ob, file = file.path(data_dir,"processed_data/imp_m1_female_single_m200_ob.RData"))
load(file = file.path(data_dir,"processed_data/imp_m1_female_single_m200_ob.RData"))

start <- Sys.time()
imp_m1_male_single_m200_ob <-  mice(data=male,
                                    method = method_male, 
                                    predictorMatrix = predictor_male,
                                    m = 1, maxit = 200, printFlag = FALSE, seed = 1045)
end <- Sys.time()
end - start #14.53456 mins
save(imp_m1_male_single_m200_ob, file = file.path(data_dir,"processed_data/imp_m1_male_single_m200_ob.RData"))
load(file = file.path(data_dir,"processed_data/imp_m1_male_single_m200_ob.RData")) 

############################
# Produce a trace plot for each partially observed variable 
############################

final_vars_female <- names(mi_data)[
  !grepl("pub_m", names(mi_data)) &
    !(names(mi_data) %in% c("c0_female", "aln", "qlet", "alnqlet"))
]

final_vars_male <- names(mi_data)[
  !grepl("pub_f", names(mi_data)) &
    !(names(mi_data) %in% c("c0_female", "aln", "qlet", "alnqlet"))
]

vars_imputed_female <- c(final_vars_female, "e0xe1", "e1xe2", "e0xe2", "e0xe1xe2")
vars_imputed_male <- c(final_vars_male, "e0xe1", "e1xe2", "e0xe2", "e0xe1xe2")

for (var in vars_imputed_female) {
  formula <- as.formula(paste(var, "~ .it | .ms"))
  #file_name <- file.path(output_dir, "multiple_imputation/test/imputation_m1_100_ob", paste0("female_plot_", var, ".png"))
  file_name <- file.path(output_dir, "multiple_imputation/test/imputation_m1_200_ob", paste0("female_plot_", var, ".png"))
  png(filename = file_name, width = 800, height = 600)
  #plot <- plot(imp_m1_female_single_m100_ob, formula, layout = c(2, 1), main = var)
  plot <- plot(imp_m1_female_single_m200_ob, formula, layout = c(2, 1), main = var)
  print(plot)
  dev.off()
}

for (var in vars_imputed_male) {
  formula <- as.formula(paste(var, "~ .it | .ms"))
  #file_name <- file.path(output_dir, "multiple_imputation/test/imputation_m1_100_ob", paste0("male_plot_", var, ".png"))
  file_name <- file.path(output_dir, "multiple_imputation/test/imputation_m1_200_ob", paste0("male_plot_", var, ".png"))
  png(filename = file_name, width = 800, height = 600)
  #plot <- plot(imp_m1_male_single_m100_ob, formula, layout = c(2, 1), main = var)
  plot <- plot(imp_m1_male_single_m200_ob, formula, layout = c(2, 1), main = var)
  print(plot)
  dev.off()
}

# check all variables that are passively imputed 
imp_female_test1 <- complete(imp_m1_female_single_m100_ob ,1) %>% 
  select(c0_sib_victimization_bin, c1_sib_victimization, c0_sib_aggression_bin, c1_sib_aggression, c0_food_insecurity, c1_food_insecurity,
         y17_bmi_z, overweight, e0_victimization2, e1_victimization2, e2_victimization2, e0xe1, e1xe2, e0xe1xe2)

imp_male_test1 <- complete(imp_m1_male_single_m100_ob ,1) %>%
  select(c0_sib_victimization_bin, c1_sib_victimization, c0_sib_aggression_bin, c1_sib_aggression, c0_food_insecurity, c1_food_insecurity,
         y17_bmi_z, overweight, e0_victimization2, e1_victimization2, e2_victimization2, e0xe1, e1xe2, e0xe1xe2)


# 2. Run final imputations with the number of iterations determined with 100 imputation sets
# Only 6.96% of females and 4.66% of males are complete cases, so technically minimum imputation sets should be around 95; 
# rule of thumb: number of imputation sets should be at least as many as % of missing cases but run 100 minimum 
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
end - start #1.088272 hours

# Stop the cluster
stopCluster(cl)

# Combine mids objects into a single object
imp_merged_female_ob <- imp_pars_female[[1]]
for (n in 2:length(imp_pars_female)){
  imp_merged_female_ob <- 
    ibind(imp_merged_female_ob,
          imp_pars_female[[n]])
}

# Trim to exactly the number of imputations desired
if (imp_merged_female_ob$m > total_imputations) {
  imp_merged_female_ob$imp <- lapply(imp_merged_female_ob$imp, function(x) x[, 1:total_imputations])
  imp_merged_female_ob$m <- total_imputations
}

# Confirm
print(paste("Final number of imputations:", imp_merged_female_ob$m))

save(imp_merged_female_ob, file = file.path(data_dir,"processed_data/imp_merged_female_ob.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_female_ob.RData"))

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
end - start # 

# Stop the cluster
stopCluster(cl)

# Combine mids objects into a single object
imp_merged_male_ob <- imp_pars_male[[1]]
for (n in 2:length(imp_pars_male)){
  imp_merged_male_ob <- 
    ibind(imp_merged_male_ob,
          imp_pars_male[[n]])
}

# Trim to exactly the number of imputations desired
if (imp_merged_male_ob$m > total_imputations) {
  imp_merged_male_ob$imp <- lapply(imp_merged_male_ob$imp, function(x) x[, 1:total_imputations])
  imp_merged_male_ob$m <- total_imputations
}

# Confirm
print(paste("Final number of imputations:", imp_merged_male_ob$m))

save(imp_merged_male_ob, file = file.path(data_dir,"processed_data/imp_merged_male_ob.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_male_ob.RData"))

# Combine both imputed dataset 
#############################################################################################

# Set number of imputations (assuming both have same m)
m <- imp_merged_female_ob$m

# Get all variable names across both datasets
vars_female <- names(complete(imp_merged_female_ob, 1))
vars_male <- names(complete(imp_merged_male_ob, 1))
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
original_female <- align_vars(imp_merged_female_ob$data, all_vars, 1)
original_male   <- align_vars(imp_merged_male_ob$data, all_vars, 0)

original_combined <- rbind(original_female, original_male)
original_combined$.imp <- 0
original_combined$.id <- 1:nrow(original_combined)

# ---------------------
# 2. Create completed datasets for each imputation
# ---------------------

long_list <- list()
long_list[[1]] <- original_combined  # first entry is always .imp == 0

for (i in 1:m) {
  complete_female <- align_vars(complete(imp_merged_female_ob, i), all_vars, 1)
  complete_male   <- align_vars(complete(imp_merged_male_ob, i), all_vars, 0)
  
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
imp_merged_ob <- as.mids(long_data)
save(imp_merged_ob, file = file.path(data_dir,"processed_data/imp_merged_ob.RData"))

################################################################
# Diagnostics 

load(file = file.path(data_dir,"processed_data/imp_merged_ob.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_female_ob.RData"))
load(file = file.path(data_dir,"processed_data/imp_merged_male_ob.RData"))

# To access the first imputed dataset
imp_merged_female_ob_m1 <- complete(imp_merged_female_ob,1)
imp_merged_male_ob_m1 <- complete(imp_merged_male_ob,1)
imp_merged_ob_m1 <-  complete(imp_merged_ob,1)

# Listing the first few rows
head(imp_merged_female_ob_m1)
head(imp_merged_male_ob_m1)

# summary
summary(imp_merged_female_ob_m1)
summary(imp_merged_male_ob_m1)

x <- imp_merged_female_ob_m1 %>%
  filter(y17_bmi_z >=1.04) %>%
  select(y17_bmi_z, overweight)

imp_merged_female_ob_m1 %>%
  dplyr::count(c0_sib_aggression_bin, c1_sib_aggression)

imp_merged_male_ob_m1 %>%
  dplyr::count(c0_sib_victimization_bin, c1_sib_victimization)

imp_merged_female_ob_m1 %>%
  dplyr::count(e0_victimization2, e1_victimization2, e0xe1)

imp_merged_male_ob_m1 %>%
  dplyr::count(e1_victimization2, e2_victimization2, e1xe2)

imp_merged_female_ob_m1 %>%
  dplyr::count(e0_victimization2, e2_victimization2, e0xe2)

imp_merged_male_ob_m1 %>%
  dplyr::count(e0_victimization2, e1_victimization2, e2_victimization2, e0xe1xe2)

var_bin <- c("e0_victimization2", "e1_victimization2", "e2_victimization2", "overweight", "c0_food_insecurity", "c2_food_insecurity",
             "c0_m_marital", "c0_m_mh1", "c1_m_mh", "c2_m_mh", "c0_f_mh1", "c1_f_mh", "c2_f_mh", "c0_child_maltx1", "c1_child_maltx", 
             "c2_child_maltx", "c0_sib_aggression_bin", "c0_sib_victimization_bin", "c2_sib_aggression", "c2_sib_victimization", 
             "c0_dv1", "c1_dv", "c2_dv2")

var_cont_female <- c("c0_bw_g", "c0_m_age_delivery", "c0_sdq_total", "c1_sdq_total", "c2_sdq_total", "c0_n_siblings3", "c0_pub_f_breast", 
              "c0_pub_f_hair", "c1_pub_f_breast", "c1_pub_f_hair", "c2_pub_f_breast", "c2_pub_f_hair", "c0_bmi_z", "c1_bmi_z", "c2_bmi_z")

var_cont_male <- c("c0_bw_g", "c0_m_age_delivery", "c0_sdq_total", "c1_sdq_total", "c2_sdq_total", "c0_n_siblings3", "c0_pub_m_genitalia", 
                   "c0_pub_m_hair", "c1_pub_m_genitalia", "c1_pub_m_hair", "c2_pub_m_genitalia", "c2_pub_m_hair", "c0_bmi_z", "c1_bmi_z", 
                   "c2_bmi_z")

var_cat <- c("c0_m_parity", "c0_medu", "c0_fedu", "c0_hhincome")

var_nnorm <- c("c0_sdq_total", "c1_sdq_total", "c2_sdq_total", "c0_n_siblings3", "c0_pub_f_breast", 
               "c0_pub_f_hair", "c1_pub_f_breast", "c1_pub_f_hair", "c2_pub_f_breast", "c2_pub_f_hair", "c0_pub_m_genitalia", "c0_pub_m_hair",
               "c1_pub_m_genitalia", "c1_pub_m_hair", "c2_pub_m_genitalia", "c2_pub_m_hair")

# boxplot for continuous 
########################

boxplot_dir <- file.path(output_dir, "multiple_imputation/mi_diagnostics/boxplot_ob")

for (var in var_cont_female) {
  p <- ggmice(imp_merged_female_ob, aes(x = .imp, y = .data[[var]])) +
    geom_boxplot() +
    scale_x_discrete(drop = FALSE) +
    labs(
      x = "Imputation number",
      y = var,
      title = paste("Boxplot of", var, "by Imputation")
    )
  
  ggsave(
    filename = file.path(boxplot_dir, paste0("boxplot_", var, "_female.png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}


for (var in var_cont_male) {
  p <- ggmice(imp_merged_male_ob, aes(x = .imp, y = .data[[var]])) +
    geom_boxplot() +
    scale_x_discrete(drop = FALSE) +
    labs(
      x = "Imputation number",
      y = var,
      title = paste("Boxplot of", var, "by Imputation")
    )
  
  ggsave(
    filename = file.path(boxplot_dir, paste0("boxplot_", var, "_male.png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# stripp plot continuous 
########################

stripp_dir <- file.path(output_dir, "multiple_imputation/mi_diagnostics/stripp_plot_ob")

for (var in var_cont_female) {
  p <- ggmice(imp_merged_female_ob, aes(x = .imp, y = .data[[var]])) +
    geom_jitter() +
    labs(
      x = "Imputation number",
      y = var,
      title = paste("Stripp plot of", var, "by Imputation")
    )
  
  ggsave(
    filename = file.path(stripp_dir, paste0("stripp_plot_", var, "_female.png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

for (var in var_cont_male) {
  p <- ggmice(imp_merged_male_ob, aes(x = .imp, y = .data[[var]])) +
    geom_jitter() +
    labs(
      x = "Imputation number",
      y = var,
      title = paste("Stripp plot of", var, "by Imputation")
    )
  
  ggsave(
    filename = file.path(stripp_dir, paste0("stripp_plot_", var, "_male.png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# density plot for continuous 
########################

density_dir <- file.path(output_dir, "multiple_imputation/mi_diagnostics/density_plot_ob")

for (var in var_cont_female) {
  p <- ggmice(imp_merged_female_ob, aes(x = .data[[var]], group = .imp)) +
    geom_density() +
    labs(x = var, 
         title = paste("Density plot of", var, "by Imputation")
    )
  
  ggsave(
    filename = file.path(density_dir, paste0("density_plot_", var, "_female.png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

for (var in var_cont_male) {
  p <- ggmice(imp_merged_male_ob, aes(x = .data[[var]], group = .imp)) +
    geom_density() +
    labs(x = var, 
         title = paste("Density plot of", var, "by Imputation")
    )
  
  ggsave(
    filename = file.path(density_dir, paste0("density_plot_", var, "_male.png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# bar plots for categorical / binary 
source("/Volumes/005/working/scripts/bullying_is/function_for_prop_plots.R")

save_prop_plots(imp_merged_female_ob, output_dir = output_dir, gender_label ="female")
save_prop_plots(imp_merged_male_ob, output_dir = output_dir, gender_label ="male")

final_vars_is <- read_xlsx(file.path(data_dir, "raw_data/data_dic_bullying_is_short.xlsx"), sheet="data_final1") %>%
  pull(variables)
comb_vars_bmi <- c("c0_food_insecurity", "c2_food_insecurity", "c0_bmi_z", "c1_bmi_z", "c2_bmi_z", "y17_bmi_z", "overweight")
final_vars <- setdiff(c(final_vars_is, comb_vars_bmi), c("c0_sd_bmi", "c1_sd_bmi", "c2_sd_bmi", "depression_icd_10", "anxiety_cisr", "dep_anx"))

# Create summary statistics for observed and imputed values 
# Observed 
tab_obs_ob <- CreateTableOne(vars = final_vars, # set descriptive variables
                          data = comb_data_final_sub, # baseline
                          factorVars = c(var_bin, var_cat)) # define categorical variables

tab_obs_ob <- print(tab_obs_ob,
                 nonnormal = var_nnorm,
                 formatOptions = list(big.mark = ","),
                 showAllLevels = TRUE,
                 test = FALSE,
                 quote = FALSE, 
                 noSpaces = TRUE, 
                 printToggle = FALSE)

summary(tab_obs_ob)
write.csv(tab_obs_ob, file = file.path(output_dir, "desc/tab_obs_ob.csv"))

tab_obs_ob_sex <- CreateTableOne(vars = final_vars, # set descriptive variables
                                    strata = "c0_female",
                                    data = comb_data_final_sub, 
                                    includeNA = T,
                                    factorVars = c(var_bin, var_cat)) # define categorical variables

tab_obs_ob_sex <- print(tab_obs_ob_sex,
                           nonnormal = var_nnorm,
                           formatOptions = list(big.mark = ","),
                           showAllLevels = TRUE,
                           test = FALSE,
                           quote = FALSE, 
                           noSpaces = TRUE, 
                           printToggle = FALSE)

write.csv(tab_obs_ob_sex, file = file.path(output_dir, "desc/tab_obs_ob_sex.csv"))

# Sex specific stratified by victimization status 
tab_obs_ob_female <- CreateTableOne(vars = final_vars_female, # set descriptive variables
                                 strata = "e0_victimization2",
                                 data = female, 
                                 includeNA = T,
                                 factorVars = c(var_bin, var_cat)) # define categorical variables

tab_obs_ob_female <- print(tab_obs_ob_female,
                        nonnormal = var_nnorm,
                        formatOptions = list(big.mark = ","),
                        showAllLevels = TRUE,
                        test = FALSE,
                        quote = FALSE, 
                        noSpaces = TRUE, 
                        printToggle = FALSE)

write.csv(tab_obs_ob_female, file = file.path(output_dir, "desc/tab_obs_ob_female.csv"))

tab_obs_ob_male <- CreateTableOne(vars = final_vars_male, # set descriptive variables
                               strata = "e0_victimization2",
                               data = male, 
                               includeNA = T,
                               factorVars = c(var_bin, var_cat)) # define categorical variables

tab_obs_ob_male <- print(tab_obs_ob_male,
                      nonnormal = var_nnorm,
                      formatOptions = list(big.mark = ","),
                      showAllLevels = TRUE,
                      test = FALSE,
                      quote = FALSE, 
                      noSpaces = TRUE, 
                      printToggle = FALSE)
write.csv(tab_obs_ob_male, file = file.path(output_dir, "desc/tab_obs_ob_male.csv"))

# count missing for continuous variables 
female_cont_miss <- female %>%
  #group_by(e0_victimization2) %>%
  summarise(across(all_of(var_cont_female),
                   list(
                     missing_n = ~ sum(is.na(.)),
                     missing_pct = ~ round(100 * mean(is.na(.)), 1)
                   ),
                   .names = "{.col}_{.fn}"
  ))

male_cont_miss <- male %>%
  #group_by(e0_victimization2) %>%
  summarise(across(all_of(var_cont_male),
                   list(
                     missing_n = ~ sum(is.na(.)),
                     missing_pct = ~ round(100 * mean(is.na(.)), 1)
                   ),
                   .names = "{.col}_{.fn}"
  ))

# Imputed 
imp_merged_ob_long <- complete(imp_merged_ob,"long",include = F)
tab_imputed_ob <- CreateTableOne(vars = final_vars, 
                              data = imp_merged_ob_long, 
                              factorVars = c(var_bin, var_cat))

tab_imputed_ob <- print(tab_imputed_ob,
                     nonnormal = var_nnorm,
                     formatOptions = list(big.mark = ","),
                     showAllLevels = TRUE,
                     test = FALSE,
                     quote = FALSE, 
                     noSpaces = TRUE, 
                     printToggle = FALSE)

write.csv(tab_imputed_ob, file = file.path(output_dir, "desc/tab_imputed_ob.csv"))

tab_imputed_ob_sex <- CreateTableOne(vars = final_vars, 
                                 data = imp_merged_ob_m1, 
                                 strata = "female",
                                 addOverall = T,
                                 factorVars = c(var_bin, var_cat))

tab_imputed_ob_sex <- print(tab_imputed_ob_sex,
                        nonnormal = var_nnorm,
                        formatOptions = list(big.mark = ","),
                        showAllLevels = TRUE,
                        test = FALSE,
                        quote = FALSE, 
                        noSpaces = TRUE, 
                        printToggle = FALSE)

write.csv(tab_imputed_ob_sex, file = file.path(output_dir, "desc/tab_imputed_ob_sex.csv"))

tab_imputed_ob_female <- CreateTableOne(vars = final_vars_female, 
                                     data = imp_merged_female_ob_m1,
                                     strata = "e0_victimization2",
                                     factorVars = c(var_bin, var_cat))

tab_imputed_ob_female <- print(tab_imputed_ob_female,
                            nonnormal = var_nnorm,
                            formatOptions = list(big.mark = ","),
                            showAllLevels = TRUE,
                            test = FALSE,
                            quote = FALSE, 
                            noSpaces = TRUE, 
                            printToggle = FALSE)

write.csv(tab_imputed_ob_female, file = file.path(output_dir, "desc/tab_imputed_ob_female.csv"))

tab_imputed_ob_male <- CreateTableOne(vars = final_vars_male, 
                                   data = imp_merged_male_ob_m1,
                                   strata = "e0_victimization2",
                                   factorVars = c(var_bin, var_cat))

tab_imputed_ob_male <- print(tab_imputed_ob_male,
                          nonnormal = var_nnorm,
                          formatOptions = list(big.mark = ","),
                          showAllLevels = TRUE,
                          test = FALSE,
                          quote = FALSE, 
                          noSpaces = TRUE, 
                          printToggle = FALSE)

write.csv(tab_imputed_ob_male, file = file.path(output_dir, "desc/tab_imputed_ob_male.csv"))

