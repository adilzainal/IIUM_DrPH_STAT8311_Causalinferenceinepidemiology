#--------------------------------------------------------
# DrPH Causal inference in epidemiology - Missing data 
#--------------------------------------------------------
# Muhammad Adil ZA 
#--------------------------------------------------------
# Missing data
#--------------------------------------------------------
# In causal inference, missing data is a selection bias problem.

# objective
# 1. explore missing data
# 2. classify missing data
# 3. handle missing data

setwd("~/Desktop/Medical Statistic/R DrPH IIUM/STAT8311")
#--------------------------------------------------------

#missing data
library(Amelia)
library(Rcpp)
library(foreign)
packages_to_install <- c("mice", "VIM", "naniar", "tidyverse", "gtsummary")
install.packages(packages_to_install)
library(mice) # mice: for the single and multiple imputation
library(VIM) # VIM: for missing data exploration
library(naniar) # naniar: for missing data exploration
library(tidyverse) # tidyverse: for data wrangling and manipulation
library(gtsummary) # gtsummary: to provide a nice result in a table

data.sav = read.spss("dataedited.sav", to.data.frame = TRUE)
missmap(data.sav, main = "Missing values vs observed")
#You can add the argument na.rm=TRUE to calculate the result while ignoring the missing values.

#--------------------------------------------------------

#Explore dataset
data.sav = read.spss("healthstatus.sav", to.data.frame = TRUE)
data.sav = read_excel("drph_data copy.xlsx", sheet = 1)
missmap(data.sav, main = "Missing values vs observed")
summary(data.sav)
summary(data.sav$age)
attach(data.sav)
sapply(data.sav, mean, na.rm=TRUE) #Exclude missing value
summary(data.sav$age) # non-normal

#--------------------------------------------------------

#missing data classification

# 1.Missing completely at random (MCAR), Missing data is classified as MCAR if the missingness is unrelated to the data. 
# For example, medical records that is lost due to flood and laboratory equipment malfunction,

# 2.Missing at random (MAR), Missing data is said to be MAR if the missingness is related to other values that we completely observed, but not to the variable with missing values itself. 
# For example, an older person is more likely to complete a survey compared to the younger person

# 3.Missing not at random (MNAR).
#Lastly, the missing data is considered MNAR if the missingness is related to the variable with missing values itself 
# and other variables with completely observed values as well. 
#Also, the missingness is considered MNAR if the causes completely unknown to us. 
# In other words, we can not logically deduce that the missing data fit MCAR or MAR types. 
# For example, missing a weight information for an obese individuals as the normal weighing scale may not able to weigh the individuals. 

#--------------------------------------------------------

# explore missing data 
data(sleep, package = "VIM")
sleep <- 
  sleep %>% 
  mutate_if(is.character, as.factor)
summary(sleep) # change all character to factor using mutate 

data.sav <- 
  data.sav %>% 
  mutate_if(is.character, as.factor)
summary(data.sav) # change all character to factor using mutate 

anyNA(sleep) # to see if there is missing data
prop_miss(sleep) # percentage , if less than 10% desirable
miss_var_summary(sleep) # by variable

missmap(sleep, main = "Missing values vs observed") # visualize
aggr(sleep, numbers = TRUE, prop = FALSE) #vim package
md.pattern(sleep, plot = FALSE) # mice package 

#see correlation Additionally, we can assess the correlation:
# Between variables with missing values
# Between variables with missing values and variable with non-missing values

dummyNA <- 
  as.data.frame(abs(is.na(sleep))) %>% 
  select(NonD, Dream, Sleep, Span, Gest) # pick variable with missing values only
head(dummyNA)

#correlation btw var with miss values only
cor(dummyNA) %>% round(digits = 2)
#no strong correlations
#There is no strong correlation between variables with missing values.
#We can conclude that the missing values in one variable is not related to the missing values in another variable.
#btw all variables
cor(sleep %>% mutate_if(is.factor, as.numeric), 
    dummyNA, use = "pairwise.complete.obs") %>% 
  round(digits = 2)

#no strong correlation
#Lastly, we can do a Little’s test to determine if the missing data is MCAR or other types. 
#The null hypothesis is the missing data is MCAR. 
#Thus, a higher p value indicates probability of missing data is MCAR. 
#The test shows that the missingness in our data is either most likely MCAR 
mcar_test(sleep)

#----------------------------------------------------------------------------------------------------------------

#Handle missing data 
#We going to cover four approaches to handling missing data:

#1.Listwise deletion
#2.Simple imputation:
#    Mean substitution
#    Median substitution
#    Mode substitution
#3.Single imputation:
#    Regression imputation
#    Stochastic regression imputation
#    Decision tree imputation
#4.Multiple imputation

?sleep
---------------------------------------------------------------------------------------------------------------
# 1. listwise deletion
lw <- lm(Sleep ~ ., data = sleep)
summary(lw)

#As shown in the information above (at the bottom), 20 rows or observations were excluded due to missingness. 
#Listwise deletion approach able to produce an unbiased result only when the missing data is MCAR 
# The amount of missing data is relatively small.
# Complete Case Analysis (CCA)

---------------------------------------------------------------------------------------------------------------
#2.  A simple imputation includes mean, median and mode substitution is a relatively easy approach.
# In this approach, the missing values is replaced by a mean or median for numerical variable,
# and mode for categorical variable. 
# Mean, Median, Mode substitution
# This simple imputation approach only appropriate if missing data is MCAR 
# and the amount of missing data is relatively small.
# Mean substitution
# replace_na will replace the missing values in NonD with its mean.

mean_sub <- 
  sleep %>% 
  mutate(NonD = replace_na(NonD, mean(NonD, na.rm = T)))

summary(mean_sub)


#Median substitution
#replace_na will replace the missing values in NonD with its median.

med_sub <- 
  sleep %>% 
  mutate(NonD = replace_na(NonD, median(NonD, na.rm = T)))

summary(med_sub)

#Mode substitution
#We need to find the mode, the most frequent level or group in the variable. 
#In ethnicity variable for example

sleep$newfactor <- as.factor(sleep$Pred)
table(sleep$newfactor)

#replace_na will replace the missing values in race with its mode.

mode_sub <- 
  sleep %>% 
  mutate(newfactor = replace_na(newfactor, "2"))

summary(mode_sub)

#---------------------------------------------------------------------------------------------------------------
#3. single imputation 

#In single imputation, missing data is imputed by any methods producing a single set of a complete dataset.
#In facts, the simple imputation approaches (mean, median and mode) are part of single imputation techniques. 
#Single imputation is better compared to the previous techniques that we have covered so far.
#This approach incorporate more information from other variables to impute the missing values. 
#However, this approach produce a result with a small standard error, which reflect a false precision in the result.
#Additionally, a single imputation approach do not take into account uncertainty about the missing data (except for stochastic regression imputation). 
#Additionally, this approach is applicable if the missing data is at least MAR

#Regression imputation
#We  do a linear regression imputation for numerical variables 
# and multinomial or polytomous logistic regression for categorical variable 
# here all linear regression

#We run mice() with maxit = 0 (zero iteration) to get a model specification 
#and  removed id variable if u have as it is not useful for the analysis. 
#Here, we do not actually run the imputation yet as the iteration is zero (maxit = 0).

id <- seq(from = 1, to = 62, by = 1)
sleep$id <- id
ini <- mice(sleep %>% select(-id), m = 1, maxit = 0) 
ini 

#By default, mice() use predictive mean matching (pmm) for a numerical variable, binary logistic regression (logreg) for a categorical variable with two levels and multinomial or polytomous logistic regression (polyreg) for a categorical variable with more than two level. 
#mice() will not assign any method for variable with non-missing values by default (denote by “” in the section of imputation method above).

#We going to change the method for numerical variable (NonD etc) to a linear regression.


meth <- ini$method
meth[c(3,4,5,6,7)] <- "norm.predict"
meth


# mice() function contains a few arguments:
# m: number of imputed sets (will be used in the multiple imputation later)
# method: to specify a method of imputation
# predictorMatrix: to specify predictors for imputation
# printFlag: print history on the R console
# seed: random number for reproducibility

regImp <- mice(sleep %>% select(-id), 
               m = 1, method = meth, printFlag = F, seed = 123)
regImp # multiple imputation 

# Here, we can see the summary of our imputation model:

# Number of multiple imputation
# Imputation methods
# PredictorMatrix
# We can further assess the full predictor matrix from the model.

regImp$predictorMatrix

#The predictor matrix denotes the variable to be imputed at the left side 
#and the predictors at the top of the column. 
#1 indicates a predictor while 0 zero indicates a non-predictor. 
#This predictor matrix can be changed accordingly if needed by changing 1 to 0 (a predictor to a non-predictor) or vice versa. 
#Noted that the diagonal is 0 as the variable is not allowed to impute itself.

#As an example, to impute missing values in NonD, 
#all of the variables are used as predictors except for NonD itself as shown in the predictor matrix. 
#So, the outcome variable in linear regression equation will be the imputed variable.

#Although it seems that we impute all the missing values in each variable simultaneously, 
#in actuality each imputation model will run separately. 
#So, in our data, on top of the imputation model for NonD, 
#we have another imputation models for Dream and sleep span and Gest variables.
#Fortunately, we do not have to concern much about which variable to be used as a predictor
# as mice() will automatically select the useful predictors for each variables with missing values. 
#We can asses the imputed dataset as follows:

sleep_regImp <- complete(regImp, 1)
summary(sleep_regImp)

missmap(sleep_regImp, main = "Missing values vs observed")

#Decision tree imputation
#Decision tree or classification and regression tree (CART) is a popular method in machine learning area. 
#his method can be applied to impute a missing values for both numerical and categorical variables.

#First, we get a model specification.

ini <- mice(sleep %>% select(-id), m = 1, maxit = 0) 
ini

#Then, we specify the imputation methods for all variable with missing values
# (NonD, Dream, Sleep, Span, Gest) to decision tree.

meth <- ini$method
meth[3:7] <- "cart"
meth

#Next, we run the imputation model.

cartImp <- mice(sleep %>% select(-id), 
                m = 1, method = meth, printFlag = F, seed = 123)
cartImp

#Here is the imputed dataset.

sleep_cartImp <- complete(cartImp, 1)
summary(sleep_cartImp)

---------------------------------------------------------------------------------------------------------------
# 4. multiple imputation

#Multiple imputation is an advanced approach to missing data.
#In this approach, several imputed datasets will be generated. 
#Analyses will be run on each imputed datasets.
#Then, all the results will be combined into a pooled result.
#The advantage of multiple imputation is this approach takes into account uncertainty regarding the missing data, 
#by which the single imputation approach may fail to do. 
#Generally, this approach is applicable if the missing data is at least MAR.

miImp <- mice(sleep %>% select(-id), m = 5, printFlag = F, seed = 123)
miImp

#We can see in the result, the number of imputations is 5.
#We can extract the third imputation set as follows:

complete(miImp, 3) %>% 
  summary()

#Next, we need to check for convergence of the algorithm.

plot(miImp)

#The line in the plot should be intermingled and free of any trend. 
#The number of iteration can be further increased to make sure of this.

miImp2 <- mice.mids(miImp, maxit = 35, printFlag = F)
plot(miImp2)

#Once the imputed datasets are obtained, an analysis (for example, a linear regression) can be run as follows:

lr_mi <- with(miImp, glm(Sleep ~ BodyWgt + BrainWgt + NonD + Dream + Span))
pool(lr_mi) %>% 
  summary(conf.int = T)

#use mutate_if() from dplyr package to round up the numbers to two decimal points.

pool(lr_mi) %>% 
  summary(conf.int = T) %>% 
  as.data.frame() %>% 
  mutate_if(is.numeric, round, 2)

#mice package has provided an easy flow to run the analysis for multiple imputation:

# mice(): impute the missing data
# with(): run a statistical analysis
# pool(): pool the results
# Additionally, a model comparison can be done as well.
# There are three methods available for the model comparison:

# D1(): multivariate Wald test
# D2(): pools test statistics from each analysis of the imputed datasets
# D3(): likelihood-ratio test statistics
# D2() is less powerful compared to the other two methods.

lr_mi2 <- with(miImp, glm(Sleep ~ BrainWgt + NonD + Dream + Span))
summary(D1(lr_mi, lr_mi2))

# Multivariate Wald test is not significant. 
# Hence, removing bodyweight from the model does not reduce its predictive power.
# We can safely exclude bodyweight variable to aim for a parsimonious model.


summary(D2(lr_mi, lr_mi2))
summary(D3(lr_mi, lr_mi2))

#Also, we get a similar result from the remaining two methods for the model comparison. 
#However, this may not always be the case.
#Generally, D1() and D3() are preferred and equally good for a sample size more than 200. 
#However, for a small sample size (n < 200), D1() is better. 
#Besides, D2() should be used with cautious especially in large datasets with many missing values as it may produce a false positive estimate.


#--------------------------------------------------------

# present findings 

#gtsummary package can be used for all the analyses run on the approaches of handling missing data that we have covered
#Here is an example to get a nice table using tbl_regression() from a linear regression model run on the multiple imputation approach.

tbl_regression(lr_mi2)

# references : https://bookdown.org/drki_musa/dataanalysis/






















