# Course: STAT3306: Statistical Analysis of Genetic Data
# Taken: Semeste 2, 2024
# Project Part 1
# Author: Hugo Burton
# Date: 7th September 2024
# Due Date: 2pm 04 October 2024


phenotype <- "Fasting Glucose"
pheno_ext <- ".phen"

log <- function(...) {
    msg <- paste0(...)
    cat(msg, "\n")
}

space_to_underscore <- function(str) {
    gsub(" ", "_", str)
}

# Analysis
# Analyse your data set and write a detailed report about all
# your analyses. This should include:
# • SNP QC


# • Sample QC (excluding checking ancestry – you can assume this
#   has been done already)


# • Genome-wide association analysis of the three traits


# • Describe the most associated region of the quantitative trait


# • A comparison of the results for these different traits sets
#   demonstrating an understanding of how the results from the
#   three traits relate to each other.

data_folder <- file.path("/data/STAT3306")
project_data <- file.path(data_folder, "Project")
data_path <- file.path(project_data, "Data")
phenotypes <- file.path(project_data, "Phenotypes")


# fam <- read.table(path(prac_folder, "data.fam"))
# bim <- read.table(path(prac_folder, "data.bim"))

pheno_path <- file.path(phenotypes, paste0(space_to_underscore(phenotype), pheno_ext))
log("Reading phenotype path:", pheno_path)
pheno <- read.table(pheno_path)
