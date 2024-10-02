# Week 10 Practical
# 02/10/2024
# Hugo Burton
# STAT3306

source("Practical_partF.R")

data_folder <- "/data/STAT3306/practical/PartF"

grmBinPath <- file.path(data_folder, "QIMRX_no_twin.grm")
cat("Data path ", grmBinPath)

grm <- read_GRMBin(grmBinPath)

preview <- 5

grm[1:preview, 1:preview]
