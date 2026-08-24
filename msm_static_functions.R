##############################################
# MSM exposure regime functions 
##############################################

# 1. Define a bootstrap function for a single dataset:

bootstrap_margins <- function(data, outcome, sex, B) {
  
  pattern_names <- c("YYN", "YNY", "YNN", "NYY", "NYN", "NNY", "NNN")  # comparisons to YYY
  all_patterns <- c("YYY", pattern_names)
  
  boot_rr <- matrix(NA, nrow = B, ncol = length(pattern_names))
  boot_rd <- matrix(NA, nrow = B, ncol = length(pattern_names))
  boot_risks <- matrix(NA, nrow = B, ncol = length(all_patterns))
  
  for (b in 1:B) {
    # 1: Bootstrap sample 
    index <- sample(1:nrow(data), replace = TRUE)  # resample indices
    boot_data <- data[index, ]   
    
    # 2: Propensity score models 
    pubertal_t0 <- pubertal_covs_t0[[sex]]
    pubertal_t1 <- pubertal_covs_t1[[sex]]
    pubertal_t2 <- pubertal_covs_t2[[sex]]
    
    formula_a0_denom <- paste("e0_victimization2 ~",
                              paste(c(time_fixed_covs, time_varying_t0, pubertal_t0), collapse = " + "))
    formula_a1_denom <- paste("e1_victimization2 ~ e0_victimization2 +",
                              paste(c(time_fixed_covs, time_varying_t0, pubertal_t0, time_varying_t1, pubertal_t1), collapse = " + "))
    formula_a2_denom <- paste("e2_victimization2 ~ e0_victimization2 + e1_victimization2 +",
                              paste(c(time_fixed_covs, time_varying_t0, pubertal_t0, time_varying_t1, pubertal_t1, time_varying_t2, pubertal_t2), collapse = " + "))
    formula_a0_num <- formula_a0_denom
    formula_a1_num <- paste("e1_victimization2 ~ e0_victimization2 +",
                            paste(c(time_fixed_covs, time_varying_t0, pubertal_t0), collapse = " + "))
    formula_a2_num <- paste("e2_victimization2 ~ e0_victimization2  + e1_victimization2 +",
                            paste(c(time_fixed_covs, time_varying_t0, pubertal_t0), collapse = " + "))
    
    # 3: Predict probabilities
    boot_data$ps_a0_denom <- predict(glm(formula_a0_denom, data = boot_data, family = binomial), type = "response")
    boot_data$ps_a1_denom <- predict(glm(formula_a1_denom, data = boot_data, family = binomial), type = "response")
    boot_data$ps_a2_denom <- predict(glm(formula_a2_denom, data = boot_data, family = binomial), type = "response")
    
    boot_data$ps_a0_num <- predict(glm(formula_a0_num, data = boot_data, family = binomial), type = "response")
    boot_data$ps_a1_num <- predict(glm(formula_a1_num, data = boot_data, family = binomial), type = "response")
    boot_data$ps_a2_num <- predict(glm(formula_a2_num, data = boot_data, family = binomial), type = "response")
    
    # 4: Compute stabilized weights
    boot_data$sw_a0 <- ifelse(boot_data$e0_victimization2 == "Yes",
                              boot_data$ps_a0_num / boot_data$ps_a0_denom,
                              (1 - boot_data$ps_a0_num) / (1 - boot_data$ps_a0_denom))
    boot_data$sw_a1 <- ifelse(boot_data$e1_victimization2 == "Yes",
                              boot_data$ps_a1_num / boot_data$ps_a1_denom,
                              (1 - boot_data$ps_a1_num) / (1 - boot_data$ps_a1_denom))
    boot_data$sw_a2 <- ifelse(boot_data$e2_victimization2 == "Yes",
                              boot_data$ps_a2_num / boot_data$ps_a2_denom,
                              (1 - boot_data$ps_a2_num) / (1 - boot_data$ps_a2_denom))
    boot_data$sw_overall <- boot_data$sw_a0 * boot_data$sw_a1 * boot_data$sw_a2
    
    # Truncate weights
    upper_bound <- quantile(boot_data$sw_overall, 0.99, na.rm = TRUE)
    boot_data$sw_overall_trunc <- boot_data$sw_overall
    boot_data$sw_overall_trunc[boot_data$sw_overall_trunc >= upper_bound] <- upper_bound
    
    # 5: Fit outcome model with interaction
    formula_outcome <- as.formula(paste0("as.factor(", outcome, ") ~ e0_victimization2 * e1_victimization2 * e2_victimization2 + ",
                                         paste(c(time_fixed_covs, time_varying_t0, pubertal_t0), collapse = " + ")))
    
    model <- glm(
      formula = formula_outcome,
      data = boot_data,
      family = binomial(link = "logit"),
      weights = sw_overall_trunc
    )
    
    # 6: Estimate marginal risks 
    pattern_data <- list(
      YYY = boot_data %>% mutate(e0_victimization2 = "Yes", e1_victimization2 = "Yes", e2_victimization2 = "Yes"),
      YYN = boot_data %>% mutate(e0_victimization2 = "Yes", e1_victimization2 = "Yes", e2_victimization2 = "No"),
      YNY = boot_data %>% mutate(e0_victimization2 = "Yes", e1_victimization2 = "No", e2_victimization2 = "Yes"),
      YNN = boot_data %>% mutate(e0_victimization2 = "Yes", e1_victimization2 = "No",  e2_victimization2 = "No"),
      NYY = boot_data %>% mutate(e0_victimization2 = "No",  e1_victimization2 = "Yes",  e2_victimization2 = "Yes"),
      NYN = boot_data %>% mutate(e0_victimization2 = "No",  e1_victimization2 = "Yes",  e2_victimization2 = "No"),
      NNY = boot_data %>% mutate(e0_victimization2 = "No",  e1_victimization2 = "No",  e2_victimization2 = "Yes"),
      NNN = boot_data %>% mutate(e0_victimization2 = "No",  e1_victimization2 = "No",  e2_victimization2 = "No")
    )
    
    mean_risks <- sapply(pattern_data, function(df) {
      mean(predict(model, newdata = df, type = "response"))
    })
    
    boot_risks[b, ] <- mean_risks  
    rr_b <- mean_risks[pattern_names] / mean_risks["YYY"]
    rd_b <- mean_risks[pattern_names] - mean_risks["YYY"]
    
    boot_rr[b, ] <- rr_b
    boot_rd[b, ] <- rd_b
  }
  
  rr_se <- apply(boot_rr, 2, sd, na.rm = TRUE)
  rd_se <- apply(boot_rd, 2, sd, na.rm = TRUE)
  
  boot_rr_df <- as.data.frame(boot_rr)
  colnames(boot_rr_df) <- pattern_names
  boot_rr_df$bootstrap <- 1:B
  
  boot_rd_df <- as.data.frame(boot_rd)
  colnames(boot_rd_df) <- pattern_names
  boot_rd_df$bootstrap <- 1:B
  
  boot_risks_df <- as.data.frame(boot_risks)
  colnames(boot_risks_df) <- all_patterns
  boot_risks_df$bootstrap <- 1:B
  
  se_df <- data.frame(
    comparison = pattern_names,
    rr_se = rr_se,
    rd_se = rd_se
  )
  
  list(
    boot_rr = boot_rr_df,
    boot_rd = boot_rd_df,
    boot_risks = boot_risks_df,
    se_df = se_df
  )
}

# 2. Loop through imputations and apply the bootstrap 

# Parallelize over imputations 
get_bootstrap_SEs_par <- function(sex, outcome, m) {
  dataset_list <- datasets[[sex]]
  
  all_boot_SEs <- parLapply(cl, 1:m, function(i, outcome, sex) {
    df <- dataset_list[[i]]
    boot_se_i <- bootstrap_margins(df, outcome = outcome, sex = sex, B = B)
    
    # Save full bootstrap result
    #saveRDS(boot_se_i, file = file.path(exp_output_dir, "bootstrap_level", paste0(sex, "_", outcome, "_imp", i, "_boot.RDS")))
    
    se_df <- boot_se_i[["se_df"]]
    se_df$imputation <- i
    se_df$sex <- sex
    se_df$outcome <- outcome
    
    return(se_df)
  }, outcome = outcome, sex = sex)
  
  all_se_df <- do.call(rbind, all_boot_SEs)
  #saveRDS(all_se_df, file = file.path(exp_output_dir, "imputation_level", paste0(sex, "_", outcome, "_all_imp_SEs.RDS")))
  saveRDS(all_se_df, file = file.path("/Users/zz23061/Library/CloudStorage/OneDrive-UniversityofBristol/PhD/paper1/results", paste0(sex, "_", outcome, "_all_imp_SEs.RDS")))
  
  all_se_df
}


# 3. Combine SEs using Rubin's rule 

combine_bootstrap_variance <- function(boot_se_df, sex, outcome) {
  pattern_names <- c("YYN", "YNY", "YNN", "NYY", "NYN", "NNY", "NNN")
  final_rows <- list()
  
  for (comp in pattern_names) {
    subset_boot <- boot_se_df %>% filter(sex == sex, outcome == outcome, comparison == comp)
    var_rr <-(subset_boot$rr_se)^2
    var_rd <- (subset_boot$rd_se)^2
    
    within_var_rr <- mean(var_rr)
    within_var_rd <- mean(var_rd)
    
    final_rows[[comp]] <- data.frame(
      sex = sex,
      outcome = outcome,
      comparison = comp,
      within_var_rr = within_var_rr,
      within_var_rd = within_var_rd
    )
  }
  
  do.call(rbind, final_rows)
  
}
