# Week 7 Practical - GWAS Clumping
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


# Clump

clump_p1 <- 0.5
clump_p2 <- 0.5
clump_r2 <- 0.20
clump_kb <- 500

run_plink(paste("--clump", paste0(plink_results_dir, "gwas_pheno_1_pc.assoc.linear"), "--clump-p1", clump_p1, "--clump-p2", clump_p2, "--clump-r2", clump_r2, "--clump-kb", clump_kb), "gwas_pheno_1_clump")

