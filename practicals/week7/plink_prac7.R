# Week 7 Practical
# Hugo Burton
# 04-09-2024

# === Imports ===

library(qqman)

# === Globals ===
data_folder <- "/data/STAT3306/practical/PartD/data"
plink_results_dir <- "./results/"

# === Functions ===
run_plink <- function(plink_args, out_name) {
	plink_cmd <- paste0("plink --bfile ", data_folder, " ", plink_args, " --out ", plink_results_dir, out_name)
	print(paste("Running:", plink_cmd))
	system(plink_cmd, ignore.stdout = T)
}

# === Data Cleaning ===


maf <- 0.05
geno <- 0.1
mind <- 0.01
hwe <- 0.001
mpheno <- 1
# run_plink(paste("--maf", maf, "--geno", geno, "--mind", mind, "--hwe", hwe, "--assoc", "--pheno", paste0(data_folder, ".pheno"), "--mpheno", mpheno), "gwas_pheno_1")

# === QQ Plot ===

d <- read.table(paste0(plink_results_dir, "gwas_pheno_1.qassoc"), header=TRUE)
manhattan(d)
qq(d$P)
lambda_gc <- qchisq(1 - median(d$P), 1) / qchisq(0.5, 1)
lambda_gc

# === Principal Components ===

pcas <- 10
run_plink(paste("--pca", pcas), "pca")

# Add the PCs as covariate with ‘--linear--covar Results/pca.eigenvec’

run_plink(paste("--linear", "--covar", paste0(plink_results_dir, "pca.eigenvec"), "--pheno", paste0(data_folder, ".pheno"), "--mpheno", mpheno), "gwas_pheno_1_pc")

###

