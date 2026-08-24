##############################################
# MSM joint effects functions 
##############################################

robres <- function(mod) {
  cov <- vcovHC(mod, type = "HC0")
  stderr <- sqrt(diag(cov))
  z <- coef(mod) / stderr
  pval <- 2 * pnorm(-1 * abs(z))
  qval <- qnorm(.975)
  results <- round(cbind(
    logOR = coef(mod),
    stderr = stderr,
    LL = coef(mod) - qval * stderr,
    UL = coef(mod) + qval * stderr
  ), 6)
  cbind(results, pval)
}

create_ps_models <- function(sex, imp_data) {
  pubertal_t0 <- pubertal_covs_t0[[sex]]
  pubertal_t1 <- pubertal_covs_t1[[sex]]
  pubertal_t2 <- pubertal_covs_t2[[sex]]
  
  # Denominator A0
  formula_a0_denom <- paste("e0_victimization2 ~",
                            paste(c(time_fixed_covs, time_varying_t0, pubertal_t0), collapse = " + "))
  
  # Denominator A1
  formula_a1_denom <- paste("e1_victimization2 ~ e0_victimization2 +",
                            paste(c(time_fixed_covs, time_varying_t0, pubertal_t0, time_varying_t1, 
                                    pubertal_t1), collapse = " + "))
  
  # Denominator A2
  formula_a2_denom <- paste("e2_victimization2 ~ e0_victimization2 + e1_victimization2 +",
                            paste(c(time_fixed_covs, time_varying_t0, pubertal_t0, time_varying_t1, 
                                    pubertal_t1, time_varying_t2, pubertal_t2), collapse = " + "))
  
  # Numerator A0
  formula_a0_num <- formula_a0_denom
  
  # Numerator A1
  formula_a1_num <- paste("e1_victimization2 ~ e0_victimization2 +",
                          paste(c(time_fixed_covs, time_varying_t0, pubertal_t0), collapse = " + "))
  
  # Numerator A2
  formula_a2_num <- paste("e2_victimization2 ~ e0_victimization2  + e1_victimization2 +",
                          paste(c(time_fixed_covs, time_varying_t0, pubertal_t0), collapse = " + "))
  
  list(
    a0_denom = with(imp_data, glm(family = binomial(link = "logit"), formula = as.formula(formula_a0_denom))),
    a1_denom = with(imp_data, glm(family = binomial(link = "logit"), formula = as.formula(formula_a1_denom))),
    a2_denom = with(imp_data, glm(family = binomial(link = "logit"), formula = as.formula(formula_a2_denom))),
    a0_num = with(imp_data, glm(family = binomial(link = "logit"), formula = as.formula(formula_a0_num))),
    a1_num = with(imp_data, glm(family = binomial(link = "logit"), formula = as.formula(formula_a1_num))),
    a2_num = with(imp_data, glm(family = binomial(link = "logit"), formula = as.formula(formula_a2_num)))
  )
}

compute_sw <- function(df) {
  df$sw_a0 <- ifelse(df$e0_victimization2 == "Yes",
                     df$ps_a0_num / df$ps_a0_denom,
                     (1 - df$ps_a0_num) / (1 - df$ps_a0_denom))
  
  df$sw_a1 <- ifelse(df$e1_victimization2 == "Yes",
                     df$ps_a1_num / df$ps_a1_denom,
                     (1 - df$ps_a1_num) / (1 - df$ps_a1_denom))
  
  df$sw_a2 <- ifelse(df$e2_victimization2 == "Yes",
                     df$ps_a2_num / df$ps_a2_denom,
                     (1 - df$ps_a2_num) / (1 - df$ps_a2_denom))
  
  df$sw_a0_a1 <- df$sw_a0 * df$sw_a1
  
  df$sw_overall <- df$sw_a0 * df$sw_a1 * df$sw_a2
  
  upper_bound1 <- quantile(df$sw_a0_a1, 0.99, na.rm = TRUE)
  upper_bound2 <- quantile(df$sw_overall, 0.99, na.rm = TRUE)
  
  df$sw_a0_a1_trunc <- df$sw_a0_a1
  df$sw_overall_trunc <- df$sw_overall
  
  df$sw_a0_a1_trunc[df$sw_a0_a1_trunc >= upper_bound1] <- upper_bound1
  df$sw_overall_trunc[df$sw_overall_trunc >= upper_bound2] <- upper_bound2
  
  df
}

process_sex_ipw <- function(sex, imp_data, ps_models, m) {
  ipw_list <- vector("list", m)
  
  for (i in 1:m) {
    df <- complete(imp_data, i)
    df$ps_a0_denom <- predict(ps_models$a0_denom$analyses[[i]], newdata = df, type = "response")
    df$ps_a1_denom <- predict(ps_models$a1_denom$analyses[[i]], newdata = df, type = "response")
    df$ps_a2_denom <- predict(ps_models$a2_denom$analyses[[i]], newdata = df, type = "response")
    df$ps_a0_num <- predict(ps_models$a0_num$analyses[[i]], newdata = df, type = "response")
    df$ps_a1_num <- predict(ps_models$a1_num$analyses[[i]], newdata = df, type = "response")
    df$ps_a2_num <- predict(ps_models$a2_num$analyses[[i]], newdata = df, type = "response")
    
    ipw_list[[i]] <- compute_sw(df)
  }
  
  ipw_list
}

outcome_model <- function(sex, ipw_list, outcome) {
  
  results <- list()  # store results here
  weight_lookup <- c(
    e0_only = "sw_a0",
    e0_e1 = "sw_a0_a1_trunc",
    e0_e1_e2 = "sw_overall_trunc"
  )
  
    for (i in seq_along(ipw_list)) {
      
      # define formulas
      formula_list <- list(
        e0_only = as.formula(
          paste0(outcome, " ~ e0_victimization2 + ",
                 paste(c(time_fixed_covs, time_varying_t0, pubertal_covs_t0[[sex]]), collapse = " + "))
        ),
        e0_e1 = as.formula(
          paste0(outcome, " ~ e0_victimization2 + e1_victimization2 + ",
                 paste(c(time_fixed_covs, time_varying_t0, pubertal_covs_t0[[sex]]), collapse = " + "))
        ),
        e0_e1_e2 = as.formula(
          paste0(outcome, " ~ e0_victimization2 + e1_victimization2 + e2_victimization2 + ",
                 paste(c(time_fixed_covs, time_varying_t0, pubertal_covs_t0[[sex]]), collapse = " + "))
        )
      )
      
      for (model_type in names(formula_list)) {
        
        weight_var <- weight_lookup[model_type]
        weights_local <- ipw_list[[i]][[weight_var]]
        
        # fit model
        model <- glm(
          formula = formula_list[[model_type]],
          data = ipw_list[[i]],
          family = binomial(link = "logit"),
          weights = weights_local
        )
        
        res <- robres(model)
        coef_names <- rownames(res)
        
        # select relevant coefficients
        if (model_type == "e0_only") {
          sel <- coef_names[grepl("^e0_victimization2", coef_names)]
        } else if (model_type == "e0_e1") {
          sel <- coef_names[grepl("^e1_victimization2", coef_names)]
        } else {
          sel <- coef_names[grepl("^e[0-2]_victimization2", coef_names)]
        }
        
        results[[length(results) + 1]] <- data.frame(
          imp = i,
          model = model_type,
          variable = sub("Yes$", "", sel),
          logOR = res[sel, "logOR"],
          stderr = res[sel, "stderr"],
          sex = sex,
          outcome = outcome,
          stringsAsFactors = FALSE
        )
      }
    }
  
  do.call(rbind, results)
}

