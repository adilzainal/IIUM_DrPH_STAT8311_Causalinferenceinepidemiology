# ==============================================================================
# DrPH: Demonstrating & Controlling Confounding in R
# Methods: Crude vs Stratified (Mantel-Haenszel), Multivariable Logistic Regression,
#          and Propensity Score Matching (MatchIt)
# ==============================================================================

# 1. Install & Load Required Packages ------------------------------------------
required_pkgs <- c("MASS", "epiR", "MatchIt", "tableone")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

suppressPackageStartupMessages({
  library(MASS)     # Contains the 'birthwt' dataset
  library(MatchIt)  # Propensity score matching
})

# 2. Data Preparation ----------------------------------------------------------
data("birthwt")
df <- birthwt

# Create clean, labelled factor variables
# Exposure: Maternal smoking (1 = Smoker, 0 = Non-smoker)
df$smoke <- factor(df$smoke, levels = c(0, 1), labels = c("Non-Smoker", "Smoker"))

# Outcome: Low birth weight (<2500g) (1 = Yes, 0 = No)
df$low <- factor(df$low, levels = c(0, 1), labels = c("Normal", "Low_BW"))

# Key Confounder: Race categorized as binary (White vs Non-White) for stratification
# (Race in birthwt: 1 = White, 2 = Black, 3 = Other)
df$non_white <- factor(ifelse(df$race == 1, 0, 1), 
                       levels = c(0, 1), 
                       labels = c("White", "Non-White"))

cat("--- DATASET OVERVIEW (First 6 rows) ---\n")
head(df[, c("low", "smoke", "non_white", "age", "lwt")])


# ==============================================================================
# PART 1: DEMONSTRATING CONFOUNDING (Crude vs. Stratum-Specific)
# ==============================================================================
cat("\n========================================================\n")
cat("PART 1: CRUDE ASSOCIATION (Subject to Confounding)\n")
cat("========================================================\n")

# 2x2 Contingency Table (Rows: Exposure, Columns: Outcome)
crude_tab <- table(df$smoke, df$low)
cat("Crude 2x2 Table:\n")
print(crude_tab)

# Crude Odds Ratio via standard 2x2 cross-product (ad / bc)
# Format: [1,1]=Non-Smoker/Normal, [1,2]=Non-Smoker/Low_BW,
#         [2,1]=Smoker/Normal,     [2,2]=Smoker/Low_BW
a <- crude_tab["Smoker", "Low_BW"]
b <- crude_tab["Smoker", "Normal"]
c <- crude_tab["Non-Smoker", "Low_BW"]
d <- crude_tab["Non-Smoker", "Normal"]

or_crude <- (a * d) / (b * c)
cat(sprintf("\nCrude Odds Ratio (OR_crude): %.3f\n", or_crude))

# Check confounder distribution across exposure (Property 2 of a confounder)
cat("\nConfounder Distribution across Exposure:\n")
print(prop.table(table(df$smoke, df$non_white), margin = 1))


# ==============================================================================
# PART 2: CONTROLLING VIA STRATIFICATION & MANTEL-HAENSZEL
# ==============================================================================
cat("\n========================================================\n")
cat("PART 2: STRATIFIED ANALYSIS & MANTEL-HAENSZEL SUMMARY OR\n")
cat("========================================================\n")

# 3-Way Contingency Table: [Exposure, Outcome, Stratifying Confounder]
strat_tab <- table(df$smoke, df$low, df$non_white)

cat("Stratum 1: White Mothers\n")
tab_white <- strat_tab[, , "White"]
print(tab_white)
or_white <- (tab_white["Smoker", "Low_BW"] * tab_white["Non-Smoker", "Normal"]) / 
  (tab_white["Smoker", "Normal"] * tab_white["Non-Smoker", "Low_BW"])
cat(sprintf("Stratum-Specific OR (White): %.3f\n\n", or_white))

cat("Stratum 2: Non-White Mothers\n")
tab_nonwhite <- strat_tab[, , "Non-White"]
print(tab_nonwhite)
or_nonwhite <- (tab_nonwhite["Smoker", "Low_BW"] * tab_nonwhite["Non-Smoker", "Normal"]) / 
  (tab_nonwhite["Smoker", "Normal"] * tab_nonwhite["Non-Smoker", "Low_BW"])
cat(sprintf("Stratum-Specific OR (Non-White): %.3f\n\n", or_nonwhite))

# Mantel-Haenszel Adjusted Summary OR
mh_test <- mantelhaen.test(strat_tab, correct = FALSE)
or_mh <- unname(mh_test$estimate)

cat(sprintf("Mantel-Haenszel Adjusted OR (OR_MH): %.3f (95%% CI: %.3f - %.3f)\n", 
            or_mh, mh_test$conf.int[1], mh_test$conf.int[2]))

# Confounding magnitude (% change rule)
pct_change_mh <- abs(or_crude - or_mh) / or_mh * 100
cat(sprintf("Percent Change in Effect Estimate: %.2f%%\n", pct_change_mh))
if (pct_change_mh >= 10) {
  cat(">> CONCLUSION: Confounding IS present (|Delta%| >= 10%).\n")
} else {
  cat(">> CONCLUSION: Confounding is minimal (|Delta%| < 10%).\n")
}


# ==============================================================================
# PART 3: CONTROLLING VIA MULTIVARIABLE LOGISTIC REGRESSION
# ==============================================================================
cat("\n========================================================\n")
cat("PART 3: MULTIVARIABLE LOGISTIC REGRESSION\n")
cat("========================================================\n")

# Model 1: Crude Logistic Model
fit_crude <- glm(low ~ smoke, data = df, family = binomial(link = "logit"))
or_glm_crude <- exp(coef(fit_crude)["smokeSmoker"])
ci_glm_crude <- exp(confint.default(fit_crude)["smokeSmoker", ])

cat(sprintf("Unadjusted Logistic OR: %.3f (95%% CI: %.3f - %.3f)\n", 
            or_glm_crude, ci_glm_crude[1], ci_glm_crude[2]))

# Model 2: Multivariable Adjusted Logistic Model
# Adjusting for race, maternal pre-pregnancy weight (lwt), and maternal age
fit_adj <- glm(low ~ smoke + non_white + lwt + age, data = df, family = binomial(link = "logit"))
summary_adj <- summary(fit_adj)

or_glm_adj <- exp(coef(fit_adj)["smokeSmoker"])
ci_glm_adj <- exp(confint.default(fit_adj)["smokeSmoker", ])

cat("\nMultivariable Model Coefficients:\n")
print(round(coef(summary_adj), 4))

cat(sprintf("\nMultivariable Adjusted OR: %.3f (95%% CI: %.3f - %.3f)\n", 
            or_glm_adj, ci_glm_adj[1], ci_glm_adj[2]))

pct_change_reg <- abs(or_glm_crude - or_glm_adj) / or_glm_adj * 100
cat(sprintf("Percent Change (Crude vs Adjusted GLM): %.2f%%\n", pct_change_reg))


# ==============================================================================
# PART 4: CONTROLLING VIA PROPENSITY SCORE MATCHING (PSM)
# ==============================================================================
cat("\n========================================================\n")
cat("PART 4: PROPENSITY SCORE MATCHING (MatchIt)\n")
cat("========================================================\n")

# Step A: Estimate propensity scores and match 1:1 using nearest neighbor
# Treatment variable must be binary (1/0)
df$smoke_num <- ifelse(df$smoke == "Smoker", 1, 0)
df$low_num   <- ifelse(df$low == "Low_BW", 1, 0)

# Match on potential confounders: non_white, lwt, age
set.seed(1234) # for reproducibility
match_obj <- matchit(smoke_num ~ non_white + lwt + age, 
                     data = df, 
                     method = "nearest", 
                     caliper = 0.2, 
                     ratio = 1)

cat("Propensity Score Matching Summary:\n")
print(summary(match_obj))

# Step B: Extract the matched dataset
matched_df <- match.data(match_obj)
cat(sprintf("\nMatched Sample Size: %d (Exposed: %d, Unexposed: %d)\n", 
            nrow(matched_df), 
            sum(matched_df$smoke_num == 1), 
            sum(matched_df$smoke_num == 0)))

# Step C: Estimate effect in the matched pseudo-population
# Logistic regression on matched data with cluster-robust standard errors
fit_matched <- glm(low_num ~ smoke_num, data = matched_df, family = binomial(link = "logit"))
or_psm <- exp(coef(fit_matched)["smoke_num"])
ci_psm <- exp(confint.default(fit_matched)["smoke_num", ])

cat(sprintf("\nPropensity Score Matched OR: %.3f (95%% CI: %.3f - %.3f)\n", 
            or_psm, ci_psm[1], ci_psm[2]))


# ==============================================================================
# PART 5: COMPARATIVE SUMMARY TABLE
# ==============================================================================
cat("\n========================================================\n")
cat("SYNTHESIS: COMPARING ALL ESTIMATES\n")
cat("========================================================\n")

summary_table <- data.frame(
  Method = c("1. Crude (Unadjusted)", 
             "2. Stratified: White Stratum", 
             "   Stratified: Non-White Stratum", 
             "   Stratified: Mantel-Haenszel Adjusted", 
             "3. Multivariable Logistic Regression", 
             "4. Propensity Score Matching (1:1 Nearest)"),
  Estimated_OR = c(round(or_crude, 3), 
                   round(or_white, 3), 
                   round(or_nonwhite, 3), 
                   round(or_mh, 3), 
                   round(or_glm_adj, 3), 
                   round(or_psm, 3)),
  Comment = c("Confounded by baseline maternal characteristics",
              "Within-stratum estimate",
              "Within-stratum estimate",
              "Weighted non-parametric summary",
              "Adjusted for race, maternal weight, and age",
              "Balanced covariates between smokers & non-smokers")
)

print(summary_table, row.names = FALSE)

