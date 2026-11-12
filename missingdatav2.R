# ==============================================================================
# DrPH CAUSAL EPIDEMIOLOGY: MISSING DATA WORKFLOW
# Part A: Identify, Classify (m-DAG), and Check (MCAR vs. MAR)
# Part B: Handle Missing Data via 4 Approaches:
#         1. Listwise Deletion (Complete Case Analysis)
#         2. Simple Imputation (Mean, Median, Mode)
#         3. Single Imputation (Deterministic Regression, Stochastic Regression, Decision Tree)
#         4. Multiple Imputation (MICE + Causal G-Computation & IPTW + Rubin's Rules)
# Part C: Benchmark All Methods Against Ground Truth
# ==============================================================================

# 0. INSTALL & LOAD PACKAGES ---------------------------------------------------
required_pkgs <- c(
  "tidyverse", "naniar", "dagitty", "rpart", "mice", 
  "WeightIt", "cobalt", "survey", "marginaleffects"
)
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(tidyverse)
library(naniar)           # Missingness identification & Little's MCAR test
library(dagitty)          # Causal & Missingness DAGs (m-DAGs)
library(rpart)            # Decision tree (CART) single imputation
library(mice)             # Multiple Imputation by Chained Equations
library(WeightIt)         # Inverse Probability of Treatment Weighting (IPTW)
library(cobalt)           # Covariate balance diagnostics
library(survey)           # Robust variance for IPTW
library(marginaleffects)  # Marginal causal Risk Difference (ATE) via G-computation

set.seed(2026)

# SIMULATE EPIDEMIOLOGIC COHORT (N = 2,000) ------------------------------------
# Exposure (A): Occupational/Environmental exposure (0/1)
# Outcome  (Y): 1-year Incident Cardiovascular Disease (0/1)
# Confounders:
#   L1 = Age in years (continuous, symmetric, fully observed)
#   L2 = Baseline Hypertension (binary 0/1, ~15% MAR)
#   L3 = Socioeconomic Index (continuous symmetric, ~20% MAR)
#   L4 = Inflammatory Biomarker CRP (continuous right-skewed, ~20% MAR)
# Auxiliary (Z): Clinic Distance / Access Score (predicts missingness & L3/L4)

n <- 2000
cohort_full <- tibble(
  id = 1:n,
  L1 = rnorm(n, mean = 55, sd = 10),
  L2 = rbinom(n, size = 1, prob = 0.40),
  L3 = rnorm(n, mean = 0, sd = 1),
  L4 = rlnorm(n, meanlog = 0.5, sdlog = 0.65), # Right-skewed biomarker
  Z  = 0.60 * L3 - 0.30 * log(L4) + rnorm(n, mean = 0, sd = 0.70)
) %>%
  mutate(
    # Confounded treatment/exposure mechanism
    p_A = plogis(-2.1 + 0.03 * L1 + 0.55 * L2 - 0.60 * L3 + 0.15 * L4),
    A   = rbinom(n, size = 1, prob = p_A),
    
    # True outcome mechanism
    p_Y = plogis(-2.8 + 0.80 * A + 0.025 * L1 + 0.65 * L2 - 0.45 * L3 + 0.18 * L4),
    Y   = rbinom(n, size = 1, prob = p_Y)
  )

# Inject MAR missingness driven by observed A, L1, L2, and Auxiliary Z
cohort <- cohort_full %>%
  mutate(
    R_L2 = rbinom(n, 1, plogis(-2.9 + 0.02 * L1 + 0.60 * Z)),
    R_L3 = rbinom(n, 1, plogis(-2.7 + 0.02 * L1 + 0.90 * Z)),
    R_L4 = rbinom(n, 1, plogis(-2.6 + 0.40 * A  + 0.85 * Z)),
    R_Y  = rbinom(n, 1, plogis(-2.6 + 0.75 * A  + 0.50 * L1/50 + 0.85 * Z)),
    
    L2 = ifelse(R_L2 == 1, NA, L2),
    L3 = ifelse(R_L3 == 1, NA, L3),
    L4 = ifelse(R_L4 == 1, NA, L4),
    Y  = ifelse(R_Y  == 1, NA, Y)
  ) %>%
  select(id, A, Y, L1, L2, L3, L4, Z)


# ==============================================================================
# PART A: IDENTIFY, CLASSIFY (m-DAG), AND CHECK MISSINGNESS
# ==============================================================================

# 1. IDENTIFY: Quantify missingness by variable and case ----------------------
miss_var_summary(cohort %>% select(-id))
pct_miss_case(cohort %>% select(-id))     # % of cohort dropped under listwise deletion
md.pattern(cohort %>% select(-id), rotate.names = TRUE)
vis_miss(cohort %>% select(-id), cluster = TRUE)

# 2. CLASSIFY: Formalize Missingness DAG (m-DAG) -------------------------------
# Arrows into R_Y from A and Z show loss to follow-up depends on exposure and
# auxiliary factors correlated with confounders -> MAR (biases Complete Case Analysis)
m_dag <- dagitty('dag {
  L1 -> A; L1 -> Y; L1 -> R_L2; L1 -> R_L3
  L2 -> A; L2 -> Y
  L3 -> A; L3 -> Y; L3 <-> Z
  L4 -> A; L4 -> Y; L4 <-> Z
  Z  -> R_L2; Z -> R_L3; Z -> R_L4; Z -> R_Y
  A  -> Y; A -> R_L4; A -> R_Y
}')
adjustmentSets(m_dag, exposure = "A", outcome = "Y")
plot(graphLayout(m_dag))

# 3. CHECK: Empirical Diagnostics (MCAR vs. MAR) -------------------------------
# 3A. Little's MCAR Test (p < 0.05 rejects MCAR)
mcar_test(cohort %>% select(-id))

# 3B. Regress missingness indicators (R) on fully observed covariates (A, L1, Z)
cohort_diag <- cohort %>%
  mutate(
    miss_L2 = as.integer(is.na(L2)),
    miss_L3 = as.integer(is.na(L3)),
    miss_L4 = as.integer(is.na(L4)),
    miss_Y  = as.integer(is.na(Y))
  )

summary(glm(miss_L3 ~ A + L1 + Z, data = cohort_diag, family = binomial)) # Z & L1 significant -> MAR
summary(glm(miss_Y  ~ A + L1 + Z, data = cohort_diag, family = binomial)) # A & Z significant -> MAR

# 3C. Covariate imbalance between observed vs. missing outcome (|SMD| > 0.10)
bal.tab(miss_Y ~ A + L1 + Z, data = cohort_diag, stats = "m", s.d.denom = "pooled")


# ==============================================================================
# HELPER FUNCTION: CAUSAL G-COMPUTATION ESTIMATOR (ATE: RISK DIFFERENCE)
# Fits standardized causal outcome model and computes marginal Risk Difference
# ==============================================================================
estimate_ate <- function(df, method_name) {
  mod <- glm(Y ~ A + L1 + L2 + L3 + L4, data = df, family = binomial())
  ate <- avg_comparisons(mod, variables = "A", comparison = "difference")
  
  tibble(
    Method     = method_name,
    N_Analyzed = nrow(df),
    ATE_RD     = ate$estimate,
    SE         = ate$std.error,
    CI_Low     = ate$conf.low,
    CI_High    = ate$conf.high
  )
}

# Benchmark 0: Ground Truth (Full data before missingness was injected)
res_truth <- estimate_ate(cohort_full, "0. Ground Truth (Full Cohort)")


# ==============================================================================
# PART B: HANDLE MISSING DATA — THE FOUR APPROACHES
# ==============================================================================

# ------------------------------------------------------------------------------
# APPROACH 1: LISTWISE DELETION (Complete Case Analysis)
# Drops any participant with >= 1 missing value across analysis variables.
# ------------------------------------------------------------------------------
df_listwise <- cohort %>%
  select(A, Y, L1, L2, L3, L4) %>%
  drop_na()

res_1_listwise <- estimate_ate(df_listwise, "1. Listwise Deletion (CCA)")
res_1_listwise

# ------------------------------------------------------------------------------
# APPROACH 2: SIMPLE IMPUTATION (Unconditional Substitution)
# Replaces missing values with marginal summary statistics (ignores relationships).
# ------------------------------------------------------------------------------

# 2A. Mean Substitution (Arithmetic mean; binary rounded to nearest integer)
df_mean <- cohort %>%
  mutate(
    L2 = ifelse(is.na(L2), round(mean(L2, na.rm = TRUE)), L2),
    L3 = ifelse(is.na(L3), mean(L3, na.rm = TRUE), L3),
    L4 = ifelse(is.na(L4), mean(L4, na.rm = TRUE), L4),
    Y  = ifelse(is.na(Y),  round(mean(Y,  na.rm = TRUE)), Y)
  )
res_2a_mean <- estimate_ate(df_mean, "2a. Simple: Mean Substitution")

# 2B. Median Substitution (Preferred over mean for skewed variables like L4)
df_median <- cohort %>%
  mutate(
    L2 = ifelse(is.na(L2), median(L2, na.rm = TRUE), L2),
    L3 = ifelse(is.na(L3), median(L3, na.rm = TRUE), L3),
    L4 = ifelse(is.na(L4), median(L4, na.rm = TRUE), L4),
    Y  = ifelse(is.na(Y),  median(Y,  na.rm = TRUE), Y)
  )
res_2b_median <- estimate_ate(df_median, "2b. Simple: Median Substitution")

# 2C. Mode Substitution (Most frequent value; strictly suited for categorical data)
get_mode <- function(x) {
  x_clean <- na.omit(x)
  # Bin continuous variables to 1 decimal place to identify modal bin
  if (length(unique(x_clean)) > 10) x_clean <- round(x_clean, 1)
  ux <- unique(x_clean)
  ux[which.max(tabulate(match(x_clean, ux)))]
}

df_mode <- cohort %>%
  mutate(
    L2 = ifelse(is.na(L2), get_mode(L2), L2),
    L3 = ifelse(is.na(L3), get_mode(L3), L3),
    L4 = ifelse(is.na(L4), get_mode(L4), L4),
    Y  = ifelse(is.na(Y),  get_mode(Y),  Y)
  )
res_2c_mode <- estimate_ate(df_mode, "2c. Simple: Mode Substitution")


# ------------------------------------------------------------------------------
# APPROACH 3: SINGLE IMPUTATION (Model-Based, m = 1)
# Uses conditional relationships (including auxiliary Z) to impute once (m = 1).
# Sequential chain ensures predictors are complete at each step.
# ------------------------------------------------------------------------------

# 3A. Deterministic Regression Imputation (Predicted conditional mean, NO noise)
df_det_reg <- cohort

# Step 1: Impute L2 via logistic regression
m_L2 <- glm(L2 ~ A + L1 + Z, data = cohort, family = binomial())
idx_L2 <- is.na(df_det_reg$L2)
df_det_reg$L2[idx_L2] <- as.integer(predict(m_L2, newdata = df_det_reg[idx_L2, ], type = "response") >= 0.5)

# Step 2: Impute L3 via linear regression
m_L3 <- lm(L3 ~ A + L1 + L2 + Z, data = cohort)
idx_L3 <- is.na(df_det_reg$L3)
df_det_reg$L3[idx_L3] <- predict(m_L3, newdata = df_det_reg[idx_L3, ])

# Step 3: Impute L4 via linear regression
m_L4 <- lm(L4 ~ A + L1 + L2 + L3 + Z, data = cohort)
idx_L4 <- is.na(df_det_reg$L4)
df_det_reg$L4[idx_L4] <- predict(m_L4, newdata = df_det_reg[idx_L4, ])

# Step 4: Impute Y via logistic regression
m_Y <- glm(Y ~ A + L1 + L2 + L3 + L4 + Z, data = cohort, family = binomial())
idx_Y <- is.na(df_det_reg$Y)
df_det_reg$Y[idx_Y] <- as.integer(predict(m_Y, newdata = df_det_reg[idx_Y, ], type = "response") >= 0.5)

res_3a_det <- estimate_ate(df_det_reg, "3a. Single: Deterministic Regression")


# 3B. Stochastic Regression Imputation (Predicted mean + random residual draw)
df_stoch_reg <- cohort

# Step 1: Stochastic draw for binary L2
p_L2 <- predict(m_L2, newdata = df_stoch_reg[idx_L2, ], type = "response")
df_stoch_reg$L2[idx_L2] <- rbinom(sum(idx_L2), size = 1, prob = p_L2)

# Step 2: Predicted L3 + random N(0, sigma) residual
pred_L3 <- predict(m_L3, newdata = df_stoch_reg[idx_L3, ])
df_stoch_reg$L3[idx_L3] <- pred_L3 + rnorm(sum(idx_L3), mean = 0, sd = sigma(m_L3))

# Step 3: Predicted L4 + random N(0, sigma) residual
pred_L4 <- predict(m_L4, newdata = df_stoch_reg[idx_L4, ])
df_stoch_reg$L4[idx_L4] <- pred_L4 + rnorm(sum(idx_L4), mean = 0, sd = sigma(m_L4))

# Step 4: Stochastic draw for binary Y
p_Y <- predict(m_Y, newdata = df_stoch_reg[idx_Y, ], type = "response")
df_stoch_reg$Y[idx_Y] <- rbinom(sum(idx_Y), size = 1, prob = p_Y)

res_3b_stoch <- estimate_ate(df_stoch_reg, "3b. Single: Stochastic Regression")


# 3C. Decision Tree Imputation (CART via rpart; captures non-linearities)
df_tree <- cohort

# Step 1: Classification tree for binary L2
t_L2 <- rpart(factor(L2) ~ A + L1 + Z, data = cohort, method = "class")
df_tree$L2[idx_L2] <- as.integer(as.character(predict(t_L2, newdata = df_tree[idx_L2, ], type = "class")))

# Step 2: Regression tree for continuous L3
t_L3 <- rpart(L3 ~ A + L1 + L2 + Z, data = cohort, method = "anova")
df_tree$L3[idx_L3] <- predict(t_L3, newdata = df_tree[idx_L3, ])

# Step 3: Regression tree for skewed continuous L4
t_L4 <- rpart(L4 ~ A + L1 + L2 + L3 + Z, data = cohort, method = "anova")
df_tree$L4[idx_L4] <- predict(t_L4, newdata = df_tree[idx_L4, ])

# Step 4: Classification tree for binary Y
t_Y <- rpart(factor(Y) ~ A + L1 + L2 + L3 + L4 + Z, data = cohort, method = "class")
df_tree$Y[idx_Y] <- as.integer(as.character(predict(t_Y, newdata = df_tree[idx_Y, ], type = "class")))

res_3c_tree <- estimate_ate(df_tree, "3c. Single: Decision Tree (CART)")


# ------------------------------------------------------------------------------
# APPROACH 4: MULTIPLE IMPUTATION (MICE + Causal Estimation + Rubin's Rules)
# Generates m = 30 completed datasets preserving uncertainty, skewness (PMM),
# and auxiliary predictors (Z), then pools causal contrasts via Rubin's Rules.
# ------------------------------------------------------------------------------
imp_setup_data <- cohort %>%
  mutate(
    L2_fac = factor(L2),
    Y_fac  = factor(Y)
  ) %>%
  select(A, Y_fac, L1, L2_fac, L3, L4, Z)

# Configure MICE methods:
# - L2_fac, Y_fac -> "logreg" (Bayesian logistic regression)
# - L3, L4        -> "pmm" (Predictive Mean Matching; preserves skewness of L4)
ini  <- mice(imp_setup_data, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

imp_mice <- mice(
  imp_setup_data,
  m = 30,
  maxit = 20,
  method = meth,
  predictorMatrix = pred,
  seed = 2026,
  printFlag = FALSE
)

# Check MICE convergence (trace lines should mix with no trend)
plot(imp_mice, layout = c(2, 2))

# 4A. Multiple Imputation + Causal G-Computation (Directly comparable to 1-3)
mi_gcomp_list <- lapply(1:imp_mice$m, function(i) {
  dat_i <- complete(imp_mice, action = i) %>%
    mutate(
      L2 = as.numeric(as.character(L2_fac)),
      Y  = as.numeric(as.character(Y_fac))
    )
  mod_i <- glm(Y ~ A + L1 + L2 + L3 + L4, data = dat_i, family = binomial())
  avg_comparisons(mod_i, variables = "A", comparison = "difference")
})

pooled_mi_gcomp <- summary(pool(mi_gcomp_list), conf.int = TRUE)

res_4a_mi_gcomp <- tibble(
  Method     = "4a. Multiple Imputation (MICE + G-Comp)",
  N_Analyzed = nrow(cohort),
  ATE_RD     = pooled_mi_gcomp$estimate,
  SE         = pooled_mi_gcomp$std.error,
  CI_Low     = pooled_mi_gcomp$conf.low,
  CI_High    = pooled_mi_gcomp$conf.high
)

# 4B. Multiple Imputation + Doubly Robust IPTW ("Within" Approach)
mi_iptw_list <- lapply(1:imp_mice$m, function(i) {
  dat_i <- complete(imp_mice, action = i) %>%
    mutate(
      L2 = as.numeric(as.character(L2_fac)),
      Y  = as.numeric(as.character(Y_fac))
    )
  # Estimate propensity score weights inside each imputed dataset
  w_i <- weightit(A ~ L1 + L2 + L3 + L4, data = dat_i, method = "glm", estimand = "ATE")
  dat_i$w <- w_i$weights
  
  # Weighted outcome regression + G-computation
  mod_i <- glm(Y ~ A + L1 + L2 + L3 + L4, data = dat_i, weights = w, family = quasibinomial())
  avg_comparisons(mod_i, variables = "A", comparison = "difference", vcov = "HC3")
})

pooled_mi_iptw <- summary(pool(mi_iptw_list), conf.int = TRUE)

res_4b_mi_iptw <- tibble(
  Method     = "4b. Multiple Imputation (MICE + DR-IPTW)",
  N_Analyzed = nrow(cohort),
  ATE_RD     = pooled_mi_iptw$estimate,
  SE         = pooled_mi_iptw$std.error,
  CI_Low     = pooled_mi_iptw$conf.low,
  CI_High    = pooled_mi_iptw$conf.high
)


# ==============================================================================
# PART C: BENCHMARK COMPARISON AGAINST GROUND TRUTH
# ==============================================================================
true_ate <- res_truth$ATE_RD

comparison_table <- bind_rows(
  res_truth,
  res_1_listwise,
  res_2a_mean,
  res_2b_median,
  res_2c_mode,
  res_3a_det,
  res_3b_stoch,
  res_3c_tree,
  res_4a_mi_gcomp,
  res_4b_mi_iptw
) %>%
  mutate(
    Bias         = ATE_RD - true_ate,
    Rel_Bias_Pct = (Bias / true_ate) * 100,
    Covers_Truth = (CI_Low <= true_ate) & (CI_High >= true_ate)
  ) %>%
  mutate(across(c(ATE_RD, SE, CI_Low, CI_High, Bias), ~ round(.x, 4)),
         Rel_Bias_Pct = round(Rel_Bias_Pct, 1))

print(comparison_table, n = Inf)

