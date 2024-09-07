# Course: STAT3306: Statistical Analysis of Genetic Data
# Taken: Semeste 2, 2024
# Project Part 1
# Author: Hugo Burton
# Date: 7th September 2024
# Due Date: 2pm 04 October 2024

install_if_missing <- function(packages) {
    missing_pkgs <- packages[!(packages %in% installed.packages()[, "Package"])]

    if (length(missing_pkgs) > 0) {
        install.packages(missing_pkgs)
    }
    invisible(NULL)
}

# === Packages ===

required_packages <- c("crayon")
install_if_missing(required_packages)
invisible(lapply(required_packages, require, character.only = TRUE))


# === Logging Config ===

default_log_level <- "INFO"
allowed_log_levels <- c("DEBUG", "INFO", "WARN", "ERROR")
app_log_level <- "DEBUG"
level_colours <- list(
    DEBUG = crayon::cyan,
    INFO = crayon::green,
    WARN = crayon::yellow,
    ERROR = crayon::red
)

# === Globals ===

phenotype <- "Fasting Glucose"
pheno_ext <- ".phen"

# Data Paths
data_folder <- file.path("/data/STAT3306")
project_data <- file.path(data_folder, "Project")
data_path <- file.path(project_data, "Data")
phenotypes <- file.path(project_data, "Phenotypes")

# Plink
plink_out_dir <- file.path("./plink_out")
plink_datafile_basename <- "test"

# === Functions ===

cat0 <- function(...) {
    msg <- paste0(...)
    cat(msg, "\n")
}

pad <- function(width, ...) {
    padded <- sprintf(paste0("%-", width, "s"), paste0(...))
    return(padded)
}

logger <- function(log_level = "INFO", ...) {
    if (!(log_level %in% allowed_log_levels)) {
        logger(default_log_level, log_level, ...)
        return(invisible(NULL))
    }

    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    msg <- paste0("[", timestamp, "]", "  ", pad(7, log_level), ...)
    
    colour_func <- level_colours[[log_level]]
    if (is.null(colour_func)) {
        logger("ERROR", "Could not identify log level")
        colour_func <- identity
    }
    
    coloured_msg <- colour_func(msg)
    cat(coloured_msg, "\n")
    invisible(NULL)
}

space_to_underscore <- function(str) {
    gsub(" ", "_", str)
}

shell_call <- function(...) {
    system(..., ignore.stdout = TRUE)
    invisible(NULL)
}

mkdir_if_not_exist <- function(path) {
    if (!dir.exists(path)) {
        logger("INFO", "Path '", path, "' does not exist. Creating...")
        dir.create(path, recursive = TRUE)
        logger("DEBUG", "Path '", path, "' created.")
    } else {
        logger("DEBUG", "Path '", path, "' already exists.")
    }
}


run_plink <- function(plink_args, out_name) {
    mkdir_if_not_exist(plink_out_dir)
    data_files_pattern <- file.path(data_path, plink_datafile_basename)
    logger("INFO", "Plink Data Files Pattern '", data_files_pattern, "',")
    out_path <- file.path(plink_out_dir, out_name)
    plink_cmd <- paste0("plink --bfile ", data_files_pattern, " ", plink_args, " --out ", out_path)
    logger("Running: ", plink_cmd)
    shell_call(plink_cmd)
    logger("Plink results directed to '", out_path, "'.")
}


# === Main ===

pheno_path <- file.path(phenotypes, paste0(space_to_underscore(phenotype), pheno_ext))
logger("Reading phenotype path:", pheno_path)
pheno <- read.table(pheno_path)
logger("table read")

# Count the number of individuals in phenotype data
dim(pheno)

head(pheno)

length(which(!is.na(pheno[,3])))


run_plink("--missing", "missing")


# === From Task Sheet ===
# Analysis
# Analyse your data set and wridffte a detailed report about all
# your analyses. This should include:
# • SNP QC


# • Sample QC (excluding checking ancestry – you can assume this
#   has been done already)


# • Genome-wide association analysis of the three traits

# • Describe the most associated region of the quantitative trait


# • A comparison of the results for these different traits sets
#   demonstrating an understanding of how the results from the
#   three traits relate to each other.

# fam <- read.table(path(prac_folder, "data.fam"))
# bim <- read.table(path(prac_folder, "data.bim"))


