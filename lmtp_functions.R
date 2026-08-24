#######################
# MSM LMTP functions 
#######################

convert_list <- function(datalist, vars) {
  lapply(datalist, function(df) {
    for (v in vars) {
      if (is.factor(df[[v]])) {
        df[[v]] <- as.numeric(df[[v]] == "Yes")
      }
    }
    df
  })
}

make_d <- function(p) {
  function(a) {
    epsilon <- runif(length(a))
    ifelse(epsilon < p & a == 1, 0, a)
  }
}

estimate_policies_cdr <- function(ori, sex) {
  pubertal_t0 <- pubertal_covs_t0[[sex]]
  pubertal_t1 <- pubertal_covs_t1[[sex]]
  pubertal_t2 <- pubertal_covs_t2[[sex]]
  
  for (pol in names(policies)) {
    d_fun <- policies[[pol]]
    n <- nrow(ori)
    density_ratios <- matrix(NA, nrow = n, ncol = 3)
    
    for (t in 1:3) {
      a_name <- c("e0_victimization2", "e1_victimization2", "e2_victimization2")[t]
      shifted <- ori
      shifted[[a_name]] <- d_fun(ori[[a_name]])
      
      ori$lambda <- 0
      shifted$lambda <- 1
      density_ratio_data <- rbind(ori, shifted)
      
      if (t == 1) parents <- c("e0_victimization2", time_fixed_covs, time_varying_t0, pubertal_t0)
      if (t == 2) parents <- c("e0_victimization2", "e1_victimization2", time_fixed_covs, time_varying_t0, pubertal_t0, 
                               time_varying_t1, pubertal_t1)
      if (t == 3) parents <- c("e0_victimization2", "e1_victimization2", "e2_victimization2", time_fixed_covs, time_varying_t0, 
                               pubertal_t0, time_varying_t1, pubertal_t1, time_varying_t2, pubertal_t2)
      
      formula <- as.formula(paste("lambda ~", paste(parents, collapse = "+")))
      model <- glm(formula, data = density_ratio_data, family = binomial)
      p1 <- predict(model, ori, type = "response")
      odds <- p1 / (1 - p1)
      density_ratios[, t] <- odds
    }
    
    cdr <- apply(density_ratios, 1, prod)
    cdr <- pmin(cdr, quantile(cdr, 0.99)) # truncate at 99th percentile 
    ori[[paste0("cdr_", pol)]] <- cdr
  }
  
  return(ori)
}

estimate_risks_from_cdr <- function(dat, outcome) {
  risks <- c()
  for (cdr in c("cdr_d15", "cdr_d25", "cdr_d40")) {
    risks[cdr] <- mean(dat[[cdr]] * dat[[outcome]], na.rm = TRUE)
  }
  return(risks)
}

#1 Bootstrap function for one dataset (one imputation)

bootstrap_lmtp <- function(data, sex, outcome, B) {
  boot_results <- vector("list", B)
  
  for (b in 1:B) {
    #1: Bootstrap sample 
    index <- sample(1:nrow(data), replace = TRUE)  # resample indices
    boot_data <- data[index, ] 
    
    #2: Estimate density ratios for this bootstrap sample
    boot_data <- estimate_policies_cdr(boot_data, sex)
    
    #3: Compute weighted risk 
    risks <- c(
      cdr_d0  = mean(boot_data[[outcome]], na.rm = TRUE), # natural risk
      cdr_d15 = mean(boot_data$cdr_d15 * boot_data[[outcome]]),
      cdr_d25 = mean(boot_data$cdr_d25 * boot_data[[outcome]]),
      cdr_d40 = mean(boot_data$cdr_d40 * boot_data[[outcome]])
    )
    
    #4: Store risk and effect measures 
    boot_results[[b]] <- tibble::tibble(
      risk_cdr_d0  = risks["cdr_d0"],
      risk_cdr_d15 = risks["cdr_d15"],
      risk_cdr_d25 = risks["cdr_d25"],
      risk_cdr_d40 = risks["cdr_d40"],
      rr_d15_vs_d0 = risks["cdr_d15"] / risks["cdr_d0"],
      rr_d25_vs_d0 = risks["cdr_d25"] / risks["cdr_d0"],
      rr_d40_vs_d0 = risks["cdr_d40"] / risks["cdr_d0"],
      rr_d25_vs_d15 = risks["cdr_d25"] / risks["cdr_d15"],
      rr_d40_vs_d15 = risks["cdr_d40"] / risks["cdr_d15"],
      rd_d15_vs_d0 = risks["cdr_d15"] - risks["cdr_d0"],
      rd_d25_vs_d0 = risks["cdr_d25"] - risks["cdr_d0"],
      rd_d40_vs_d0 = risks["cdr_d40"] - risks["cdr_d0"],
      rd_d25_vs_d15 = risks["cdr_d25"] - risks["cdr_d15"],
      rd_d40_vs_d15 = risks["cdr_d40"] - risks["cdr_d15"]
    )
  }
  
  #5: Combine all bootstrap estimates 
  boot_df <- bind_rows(boot_results)
  
  results <- boot_df %>%
    summarise(across(
      everything(), 
      list(est = mean, se = sd),
      na.rm = TRUE, 
      .names = "{.col}_{.fn}"
    ))
  
  return(results)
}


#2 Loop through imputations; get output per imputation

estimate_lmtp_par <- function(sex, outcome, m) {
  imputation_results <- parLapply(cl, 1:m, function(i, sex, outcome) {
    datasets_list <- datasets[[sex]]
    dat <- datasets_list[[i]]
    res <- bootstrap_lmtp(dat, sex, outcome, B = B)
    res$imp <- i
    res
  }, sex = sex, outcome = outcome)
  
  dplyr::bind_rows(imputation_results)
}


#3 Pool across imputations with Rubin’s Rules

combine_bootstrap_variance_lmtp <- function(df_imp, sex, outcome) {
  se_cols <- grep("_se$", names(df_imp), value = TRUE)
  
  final_rows <- purrr::map_dfr(se_cols, function(se_col) {
    var <- (df_imp[[se_col]])^2
    within_var <- mean(var)
    
    tibble(
      sex = sex,
      outcome = outcome,
      param = gsub("_se$", "", se_col),
      within_var = within_var
    )
  })
  
  final_rows
}
