#################################################################################################################################
###### 
###### Program:     STAT3306/7306 2024 MR practical
###### R packages loaded:	AER
###### Datafiles :	data.txt	      
###### Objective :	STAT3306/7306 2024 MR practical    
######
#####################################################################################################################

# Make sure to start from a clean environment
rm(list=ls())

# Set your working directory in R to this folder and read in the dataset
setwd("/**YOUR DIRECTORY**/")


example <- read.table("data.txt", header=T) 
attach(example)

# Q1. As you’re running the commands below, fill in the graphical representation 
# of the analyses in the word document figure with the appropriate variables and beta-coefficients

# Look at the data 
# Units: SNP (0,1,2), CRP mmol/L, SBP mmHg, Income per $10,000, HDL mmol/L
head(example) 
summary(example)

# Observational Analyses
# Q2. What does the observational linear regression of SBP on CRP show?

# Run observational OLS regression (ordinary least squares) for SBP & CRP
summary(lm(SBP~CRP)) 

# Plot the observational association between SBP and CRP
plot(CRP,SBP)
abline(lm(SBP~CRP),col="red")

# Q3. What does the OLS regression of the CRP SNP rs3091244 on CRP show?

# Observational OLS regression of CRP on CRP SNP
summary(lm(CRP~rs3091244))

# Plot the relationship between CRP and rs3091244
plot(rs3091244, CRP)
abline(lm(CRP~rs3091244),col="red")

# Q4. What do the OLS regressions of potential confounders (income, HDL) show?
# Confounders
summary(lm(SBP~INCOME))
summary(lm(CRP~INCOME))
summary(lm(INCOME~rs3091244))

summary(lm(SBP~HDL))
summary(lm(CRP~HDL))
summary(lm(HDL~rs3091244))

# Q5. What are the implications for these income and HDL associations for the observational CRP-SBP association?

# Q6. Compare the unadjusted and covariate-adjusted OLS observational regressions. What do they show?

# Run a covariate-adjusted model for the association between CRP & BP
summary(lm(SBP~CRP))
summary(lm(SBP~CRP+INCOME+HDL))

# Q7. What could explain this?

# MR/IV Analyses: Wald Estimator

# Q8. Run the necessary OLS regressions to compute a Wald estimator

# OLS regression of CRP on CRP SNP
summary(lm(CRP~rs3091244))
# OLS regression of SBP on CRP SNP
summary(lm(SBP~rs3091244))

# Q9. From the above output, compute the causal effect using the Wald estimator, 
# as well as it’s SE and 95% CI. What do the results show and what do they mean?

# Q10. Rerun the observational OLS of CRP and SBP and compare with the results 
# from the Wald estimator. What do you notice about the Beta and SEs?

# Observational OLS regression
summary(lm(SBP~CRP))

# MR/IV Analyses: TSLS
# Call the AER library to run TSLS (if the AER package has been installed)
library(AER)

# If AER has not been installed, run the command below first:
#install.packages("AER") 

# General format for TSLS command: 
# summary(ivreg(Outcome~Exposure | Instrument))

# Q11. What do the TSLS results show and did it differ to the Wald estimator?

# TSLS regression
summary(ivreg(SBP~CRP | rs3091244))

# Manual TSLS

# Regress the exposure (CRP) on the instrument (rs3091244)
First_Stage <- lm(CRP~rs3091244)

# Create predicted CRP values, from the first-stage regression
Pred_CRP <- predict(First_Stage)

# Have a quick look at these values
table(Pred_CRP)
plot(rs3091244, CRP)
abline(lm(CRP~ rs3091244), col="red")

# Second stage regression
Second_Stage <- lm(SBP~Pred_CRP)

# Look at the results:
summary(Second_Stage)

# Q12. Are they the same as ‘ivreg’ TSLS function?

# Weak instruments bias 

# Q13. Looking at the F-statistic, determine if weak instruments may be an issue

#Look at F-stat from the first-stage linear regression
summary(lm(CRP~rs3091244))

#Look at F-stat from ‘diagnostics’ by AER package
summary(ivreg(SBP~CRP | rs3091244), diagnostics=T)

# Q14. How would having weak instruments change the causal estimate of CRP on SBP,
# in this study (single sample)?
