# ==============================================================================
# DrPH Causal Inference: Comprehensive Reliability, Validity & Information Bias
# ==============================================================================

library(dplyr)
library(ggplot2)
install.packages("irr")
library(irr)     # Cohen's kappa, Fleiss' kappa, ICC
library(psych)   # Cronbach's alpha, Exploratory Factor Analysis (EFA)
library(lavaan)  # Confirmatory Factor Analysis (CFA)
library(pROC)    # ROC curves and AUC

set.seed(2026)

# ------------------------------------------------------------------------------
# 1. GENERATE EXAMPLE EPIDEMIOLOGICAL STUDY DATA (N = 400)
# ------------------------------------------------------------------------------
N <- 400

# True Latent Constructs:
# F1: Occupational Stress (primary exposure of interest)
# F2: Workplace Ergonomic Support (distinct construct, weakly correlated r = 0.15)
F1_true <- rnorm(N, mean = 0, sd = 1)
F2_true <- 0.15 * F1_true + sqrt(1 - 0.15^2) * rnorm(N, mean = 0, sd = 1)

# A. Multi-Item Survey Instrument (8 items, 1-7 Likert scale)
# Items S1-S4 measure Occupational Stress (F1); Items E1-E4 measure Support (F2)
survey_df <- data.frame(
  S1 = 4.0 + 1.15 * F1_true + rnorm(N, 0, 0.75),
  S2 = 4.1 + 1.10 * F1_true + rnorm(N, 0, 0.80),
  S3 = 3.9 + 1.20 * F1_true + rnorm(N, 0, 0.70),
  S4 = 4.0 + 1.05 * F1_true + rnorm(N, 0, 0.85),
  E1 = 4.2 + 1.10 * F2_true + rnorm(N, 0, 0.78),
  E2 = 4.0 + 1.15 * F2_true + rnorm(N, 0, 0.75),
  E3 = 4.1 + 1.08 * F2_true + rnorm(N, 0, 0.82),
  E4 = 3.9 + 1.12 * F2_true + rnorm(N, 0, 0.76)
)

# Composite Scale Scores
stress_scale_score  <- rowMeans(survey_df[, c("S1", "S2", "S3", "S4")])
support_scale_score <- rowMeans(survey_df[, c("E1", "E2", "E3", "E4")])

# B. Test-Retest Measurements (Administered at Baseline T1 and 2 Weeks Later T2)
stress_T1 <- 25 + 5.0 * F1_true + rnorm(N, mean = 0, sd = 1.8)
stress_T2 <- 25 + 5.0 * F1_true + rnorm(N, mean = 0, sd = 1.8) # No systematic drift

# Binary high-stress classification at T1 and T2 (Cutoff >= 26.5)
high_stress_T1 <- as.integer(stress_T1 >= 26.5)
high_stress_T2 <- as.integer(stress_T2 >= 26.5)

# C. Criterion Gold Standard & Construct Validity Variables
# True Clinical Burnout Diagnosis (Structured Psychiatric Interview = Gold Standard)
Gold_Standard <- rbinom(N, size = 1, prob = plogis(-1.2 + 2.2 * F1_true))

# Convergent Validity Biomarker: Hair Cortisol Concentration (pg/mg, strongly tied to F1)
hair_cortisol <- 15 + 4.5 * F1_true + rnorm(N, mean = 0, sd = 2.8)

# Discriminant (Divergent) Validity Variable: Dietary Fiber Score (theoretically unrelated)
diet_fiber_score <- 50 + 0.5 * F1_true + rnorm(N, mean = 0, sd = 10)

# D. Causal Outcome: Incident Hypertension (Y) caused by True Occupational Stress
# True causal OR per 1-SD increase in F1_true = exp(0.75) ≈ 2.12
Hypertension_Y <- rbinom(N, size = 1, prob = plogis(-1.5 + 0.75 * F1_true))

cohort <- cbind(
  data.frame(
    id = 1:N, F1_true, stress_scale_score, support_scale_score,
    stress_T1, stress_T2, high_stress_T1, high_stress_T2,
    Gold_Standard, hair_cortisol, diet_fiber_score, Hypertension_Y
  ),
  survey_df
)

# E. Clinical Substudy for Inter-Rater & Intra-Rater Reliability (n = 120)
# 3 Clinicians independently score continuous severity (0-100) and binary burnout (0/1)
n_sub <- 120
F1_sub <- F1_true[1:n_sub]

rater1_occ1_cont <- 50 + 12 * F1_sub + rnorm(n_sub, 0, 3.2) # Rater 1, Occasion 1
rater1_occ2_cont <- 50 + 12 * F1_sub + rnorm(n_sub, 0, 3.2) # Rater 1, Occasion 2 (Intra-rater)
rater2_cont      <- 51 + 12 * F1_sub + rnorm(n_sub, 0, 3.6) # Rater 2 (Inter-rater)
rater3_cont      <- 49 + 12 * F1_sub + rnorm(n_sub, 0, 3.8) # Rater 3 (Inter-rater)

rater_sub <- data.frame(
  rater1_occ1_cont, rater1_occ2_cont, rater2_cont, rater3_cont,
  rater1_occ1_bin = as.integer(rater1_occ1_cont >= 54),
  rater1_occ2_bin = as.integer(rater1_occ2_cont >= 54),
  rater2_bin      = as.integer(rater2_cont >= 54),
  rater3_bin      = as.integer(rater3_cont >= 54)
)


# ==============================================================================
# PART A: RELIABILITY METRICS & BLAND-ALTMAN PLOT
# ==============================================================================

# ------------------------------------------------------------------------------
# 2.1 TEST-RETEST RELIABILITY (Same Subjects at Time 1 vs. Time 2)
# ------------------------------------------------------------------------------
# A. Intraclass Correlation Coefficient (ICC - Two-way mixed, absolute agreement)
icc_test_retest <- icc(cohort[, c("stress_T1", "stress_T2")],
                       model = "twoway", type = "agreement", unit = "single")
icc_test_retest # icc of 0.75-0.9 is good and >0.9 is excellent  


# B. Cohen's Kappa (Binary test-retest classification)
kappa_test_retest <- kappa2(cohort[, c("high_stress_T1", "high_stress_T2")])
kappa_test_retest #Kappa cut off value 0.6-0.8 is substantial, >0.8 is near perfect 

# C. Within-Subject Coefficient of Variation (CV%)
pair_mean <- (cohort$stress_T1 + cohort$stress_T2) / 2
pair_sd   <- abs(cohort$stress_T1 - cohort$stress_T2) / sqrt(2)
within_cv <- (pair_sd / pair_mean) * 100
mean_cv   <- mean(within_cv)
rms_cv    <- sqrt(mean(within_cv^2)) # Root-Mean-Square CV%
rms_cv # less than 10% is acceptable and <5% is excellence 

# D. Bland-Altman Agreement Statistics & Plot
ba_df <- data.frame(
  mean_score = pair_mean,
  diff_score = cohort$stress_T1 - cohort$stress_T2
)
ba_bias  <- mean(ba_df$diff_score)
ba_sd    <- sd(ba_df$diff_score)
loa_low  <- ba_bias - 1.96 * ba_sd
loa_high <- ba_bias + 1.96 * ba_sd

p_bland_altman <- ggplot(ba_df, aes(x = mean_score, y = diff_score)) +
  geom_point(alpha = 0.55, color = "#2b5c8f") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
  geom_hline(yintercept = ba_bias, color = "#1b7837", linewidth = 0.9) +
  geom_hline(yintercept = c(loa_low, loa_high), linetype = "dashed", color = "#d95f02", linewidth = 0.85) +
  annotate("text", x = min(ba_df$mean_score), y = ba_bias + 0.35,
           label = sprintf("Mean Bias = %.2f", ba_bias), hjust = 0, color = "#1b7837", fontface = "bold") +
  annotate("text", x = min(ba_df$mean_score), y = loa_high + 0.35,
           label = sprintf("+1.96 SD = %.2f", loa_high), hjust = 0, color = "#d95f02") +
  annotate("text", x = min(ba_df$mean_score), y = loa_low - 0.40,
           label = sprintf("-1.96 SD = %.2f", loa_low), hjust = 0, color = "#d95f02") +
  labs(
    title = "Bland-Altman Plot: Test-Retest Agreement (T1 vs. T2)",
    subtitle = sprintf("ICC = %.3f | Mean Within-Subject CV = %.2f%% (RMS-CV = %.2f%%)",
                       icc_test_retest$value, mean_cv, rms_cv),
    x = "Mean of T1 and T2 Stress Scores",
    y = "Difference (T1 - T2)"
  ) +
  theme_minimal(base_size = 12)

print(p_bland_altman)


#Four Visual Checks to Interpret Your Plot

##Check for Systematic Bias (Middle Line vs. Zero)
##What to look for: Compare the middle line (d) to 0.
## How to read it: In a pure test retest reliability study, d should sit very close to 0. 
## If the line is shifted notably above or below 0 (especially if its 95% confidence interval does not include 0),
# there is systematic drift for example, subjects consistently reporting lower stress at T2 because they became familiar with the testing environment.

##Judge the Width of the Limits of Agreement (Random Error)
## What to look for: How far apart are the top and bottom dashed lines?
##How to read it: The LoA tells you the worst-case measurement swing you can expect for 95% of individuals. 
## Unlike CV%, there is no universal statistical cutoff for LoA width; you must decide based on your domain whether a random fluctuation of, say, +/- 4 units is small enough to be practically acceptable. Narrower bands mean higher precision.

## Look for a Slope or "Funnel" Shape (Proportional Bias & Heteroscedasticity)
## Sloped cloud of points: If the dots trend upward or downward from left to right,
## the difference between T1 and T2 depends on the subject's stress level (a common sign of regression to the mean among people with extreme T1 scores).
## Funnel / cone shape (wider on the right): If the vertical spread gets wider as pair mean increases,random error grows proportionally with the magnitude of stress. 
## Seeing a funnel shape directly validates why you calculated RMS-CV earlier because percentage error stays constant even when raw unit error expands.

## Check for Outliers
## What to look for: Points falling outside the dashed 95% LoA lines.
## How to read it: Having roughly 5% of your sample (e.g., 1 out of 20 subjects) sit just outside the lines is mathematically normal. 
## However, points sitting far outside the limits flag potential sensor glitches, data entry errors, or individuals who experienced an unusual acute stressor between T1 and T2.

# ------------------------------------------------------------------------------
# 2.2 INTER-RATER RELIABILITY (2 or More Different Clinicians)
# ------------------------------------------------------------------------------
# A. Cohen's Kappa (Pairwise: Rater 1 vs. Rater 2)
kappa_inter_2 <- kappa2(rater_sub[, c("rater1_occ1_bin", "rater2_bin")])

kappa_inter_2

# B. Fleiss' Kappa (Multi-rater >= 3: Rater 1, Rater 2, and Rater 3)
fleiss_inter_3 <- kappam.fleiss(rater_sub[, c("rater1_occ1_bin", "rater2_bin", "rater3_bin")])
fleiss_inter_3

# C. Continuous Inter-Rater ICC (Two-way random effects, absolute agreement)
icc_inter_3 <- icc(rater_sub[, c("rater1_occ1_cont", "rater2_cont", "rater3_cont")],
                   model = "twoway", type = "agreement", unit = "single")
icc_inter_3

# ------------------------------------------------------------------------------
# 2.3 INTRA-RATER RELIABILITY (Same Clinician Across 2 Occasions)
# ------------------------------------------------------------------------------
kappa_intra <- kappa2(rater_sub[, c("rater1_occ1_bin", "rater1_occ2_bin")])
kappa_intra

icc_intra   <- icc(rater_sub[, c("rater1_occ1_cont", "rater1_occ2_cont")],
                   model = "twoway", type = "agreement", unit = "single")
icc_intra

intra_cv    <- mean((abs(rater_sub$rater1_occ1_cont - rater_sub$rater1_occ2_cont) / sqrt(2)) /
                      ((rater_sub$rater1_occ1_cont + rater_sub$rater1_occ2_cont) / 2)) * 100
intra_cv

# ------------------------------------------------------------------------------
# 2.4 INTERNAL CONSISTENCY (Multi-Item Survey Subscales)
# ------------------------------------------------------------------------------
alpha_stress  <- psych::alpha(cohort[, c("S1", "S2", "S3", "S4")])
alpha_support <- psych::alpha(cohort[, c("E1", "E2", "E3", "E4")])

# Print Consolidated Reliability Summary Table
reliability_summary <- data.frame(
  Domain = c("Test-Retest (Continuous)", "Test-Retest (Binary)",
             "Inter-Rater (2 Raters, Binary)", "Inter-Rater (3 Raters, Binary)",
             "Inter-Rater (3 Raters, Cont.)", "Intra-Rater (Binary)",
             "Intra-Rater (Continuous)", "Internal Consistency (Stress)",
             "Internal Consistency (Support)"),
  Statistic = c("ICC(2,1) / Mean CV%", "Cohen's Kappa",
                "Cohen's Kappa", "Fleiss' Kappa",
                "ICC(2,1)", "Cohen's Kappa",
                "ICC(2,1) / Mean CV%", "Cronbach's Alpha",
                "Cronbach's Alpha"),
  Estimate = c(
    sprintf("%.3f (CV: %.2f%%)", icc_test_retest$value, mean_cv),
    sprintf("%.3f", kappa_test_retest$value),
    sprintf("%.3f", kappa_inter_2$value),
    sprintf("%.3f", fleiss_inter_3$value),
    sprintf("%.3f", icc_inter_3$value),
    sprintf("%.3f", kappa_intra$value),
    sprintf("%.3f (CV: %.2f%%)", icc_intra$value, intra_cv),
    sprintf("%.3f", alpha_stress$total$raw_alpha),
    sprintf("%.3f", alpha_support$total$raw_alpha)
  )
)
cat("\n=== SUMMARY OF RELIABILITY ESTIMATES ===\n")
print(reliability_summary, row.names = FALSE) # 0.8 good, 0.9 excellent

# ==============================================================================
# PART B: VALIDITY METRICS, ROC CURVE & FACTOR ANALYSIS (EFA / CFA)
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1 CRITERION VALIDITY: DIAGNOSTIC ACCURACY (2x2 Table vs. Gold Standard)
# ------------------------------------------------------------------------------
tab_diag <- table(Screening = cohort$high_stress_T1, GoldStandard = cohort$Gold_Standard)
TP <- tab_diag["1", "1"]; FP <- tab_diag["1", "0"]
FN <- tab_diag["0", "1"]; TN <- tab_diag["0", "0"]

sens <- TP / (TP + FN)
spec <- TN / (TN + FP)
ppv  <- TP / (TP + FP)
npv  <- TN / (TN + FN)

cat("\n=== CRITERION DIAGNOSTIC ACCURACY (Cutoff >= 26.5) ===\n")
print(tab_diag)
cat(sprintf("Sensitivity: %.1f%% | Specificity: %.1f%% | PPV: %.1f%% | NPV: %.1f%%\n",
            sens * 100, spec * 100, ppv * 100, npv * 100))

# ------------------------------------------------------------------------------
# 3.2 CRITERION VALIDITY: ROC CURVE & AREA UNDER THE CURVE (AUC)
# ------------------------------------------------------------------------------
roc_obj <- roc(response = cohort$Gold_Standard, predictor = cohort$stress_T1, quiet = TRUE)
auc_val <- auc(roc_obj)
auc_ci  <- ci.auc(roc_obj)
opt_cut <- coords(roc_obj, "best", best.method = "youden",
                  ret = c("threshold", "sensitivity", "specificity"))

cat(sprintf("\nROC AUC: %.3f (95%% DeLong CI: %.3f - %.3f)\n", auc_val, auc_ci[1], auc_ci[3]))
cat(sprintf("Optimal Youden Cutoff: %.2f (Sens = %.1f%%, Spec = %.1f%%)\n",
            opt_cut$threshold, opt_cut$sensitivity * 100, opt_cut$specificity * 100))

p_roc <- ggroc(roc_obj, color = "#0570b0", linewidth = 1.1) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.35, y = 0.25,
           label = sprintf("AUC = %.3f\n(95%% CI: %.3f - %.3f)", auc_val, auc_ci[1], auc_ci[3]),
           size = 4.2, fontface = "bold", color = "#0570b0") +
  labs(
    title = "ROC Curve: Continuous Stress Score vs. Gold-Standard Burnout",
    x = "Specificity (True Negative Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_minimal(base_size = 12)

print(p_roc)

# ------------------------------------------------------------------------------
# 3.3 CONSTRUCT VALIDITY: EXPLORATORY (EFA) & CONFIRMATORY FACTOR ANALYSIS (CFA)
# ------------------------------------------------------------------------------
survey_items <- cohort[, c("S1", "S2", "S3", "S4", "E1", "E2", "E3", "E4")]

# A. Exploratory Factor Analysis (2-factor Maximum Likelihood with Promax rotation)
efa_fit <- fa(r = survey_items, nfactors = 2, rotate = "promax", fm = "ml")
cat("\n=== EXPLORATORY FACTOR ANALYSIS (EFA Standardized Loadings) ===\n")
print(efa_fit$loadings, cutoff = 0.30)

# B. Confirmatory Factor Analysis (2-factor structural measurement model)
cfa_model <- '
  Occupational_Stress =~ S1 + S2 + S3 + S4
  Workplace_Support   =~ E1 + E2 + E3 + E4
'
cfa_fit <- cfa(cfa_model, data = cohort, std.lv = TRUE)
fit_idx <- fitMeasures(cfa_fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

cat("\n=== CONFIRMATORY FACTOR ANALYSIS (CFA Global Fit Indices) ===\n")
print(round(fit_idx, 3)) # Target: CFI/TLI > 0.95, RMSEA < 0.06, SRMR < 0.08

# ------------------------------------------------------------------------------
# 3.4 CONVERGENT & DISCRIMINANT (DIVERGENT) VALIDITY
# ------------------------------------------------------------------------------
# A. Convergent Validity: Stress Scale vs. Hair Cortisol (Target: r > 0.70)
conv_pearson  <- cor.test(cohort$stress_scale_score, cohort$hair_cortisol, method = "pearson")
conv_spearman <- cor.test(cohort$stress_scale_score, cohort$hair_cortisol, method = "spearman")

# B. Discriminant Validity: Stress Scale vs. Unrelated Diet Score & Support Scale (Target: |r| < 0.30)
disc_diet    <- cor.test(cohort$stress_scale_score, cohort$diet_fiber_score, method = "pearson")
disc_support <- cor.test(cohort$stress_scale_score, cohort$support_scale_score, method = "pearson")

cat("\n=== CONVERGENT & DISCRIMINANT VALIDITY CORRELATIONS ===\n")
cat(sprintf("Convergent (Stress vs. Hair Cortisol) Pearson r:   %.3f (p = %.2e) [Target > 0.70]\n",
            conv_pearson$estimate, conv_pearson$p.value))
cat(sprintf("Convergent (Stress vs. Hair Cortisol) Spearman rho: %.3f (p = %.2e)\n",
            conv_spearman$estimate, conv_spearman$p.value))
cat(sprintf("Discriminant (Stress vs. Dietary Fiber) Pearson r:  %.3f (p = %.3f)    [Target < 0.30]\n",
            disc_diet$estimate, disc_diet$p.value))
cat(sprintf("Discriminant (Stress vs. Support Scale) Pearson r:  %.3f (p = %.3f)    [Target < 0.30]\n",
            disc_support$estimate, disc_support$p.value))
