# Effects of peer victimization on internalizing disorders and overweight: causal analysis accounting for time-varying confounding 

## Overview 
This repository contains R code for a causal analysis of the effects of peer victimization from mid-childhood to early adolescence (ages 8, 10, and 12 years) on internalizing disorders and overweight status in late adolescence (age 17) using data from Avon Longitudinal Study of Parents and Children (ALSPAC). To address time-varying confounding, we implemented **marginal structural models (MSMs)** estimated via stabilised inverse probability weighting to estimate:

- **Total effects (TEs)** of peer victimization at each time point
- **Controlled direct effects (CDEs)** isolating timing-specific effects
- **Joint effects** of fixed peer victimization patterns across all three ages 

We additionally used **longitudinal modified treatment policies (LMTPs)** to estimate the effects of hypothetical population-level reductions in peer victimization of 15%, 25%, and 40% on internalizing disorder risk, providing causal estimates under policy-relevant intervention scenarios. 

## File structure 
### Data preparation
- 01_desc_alive.R - Descriptive statistics and sample characteristics
- 02_mi_is.R - Multiple imputation for internalizing disorders outcome
- 02_mi_ob.R - Multiple imputation for overweight outcome 

### Main analyses 
- 03_msm_te_cde_is.R - Estimation of TEs and CDEs for internalizing disorders
- 03_msm_te_cde_ob.R - Estimation of TEs and CDEs for overweight
- 04_msm_static_is.R - Estimation of joint effects of victimization patterns for internalizing disorders
- 04_msm_static_ob.R - Estimation of joint effects of victimization patterns for overweight
- 05_lmtp_is.R - Stochastic intervention analysis for internalizing disorders using LMTP

### Functions 
- function for prop_plots.R - Functions for diagnostic plots of multiply imputed datasets
- msm_te_cde_functions.R - Functions for the estimation of TEs and CDEs
- msm_static_functions.R - Functions for the estimation of joint effects of victimization patterns
- lmtp_functions.R - Functions for stochastic intervention analysis using LTMP 



