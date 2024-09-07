# Week 7 Practical QC (Quality Control Data)
# Hugo Burton
# 04/09/2024

# Packages

library(qqman)

# Globals

plink_results_dir <- "./results/"


# Code

d <- read.table(paste0(plink_results_dir, "gwas_pheno_1_pc.assoc.linear"), header = TRUE)
head(d)

# Create a subset of d

d_subset <- d[which(d$TEST == "ADD"),]
head(d_subset)

manhattan(d_subset)

qq(d_subset$P)
lambda_gc <- qchisq(1 - median(d_subset$P), 1) / qchisq(0.5, 1)
lambda_gc

