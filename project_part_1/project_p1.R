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

# Plots
plots_out <- file.path("./plots")

# === Functions ===

cat0 <- function(...) {
    msg <- paste0(...)
    cat(msg, "\n")
}

pad <- function(width, ...) {
    padded <- sprintf(paste0("%-", width, "s"), paste0(...))
    return(padded)
}

get_calling_function <- function(ignore_names) {
    call_stack <- sys.calls()

    # Skip this function
    for (i in (length(call_stack)-1):1) {
        current_call <- call_stack[[i]]
        func_name <- as.character(current_call[[1]])

        if (!(func_name %in% ignore_names)) {
            return(func_name)
        }
    }

    return("Global Environment")
}

logger <- function(log_level = "INFO", ...) {
    if (!(log_level %in% allowed_log_levels)) {
        logger(default_log_level, log_level, ...)
        return(invisible(NULL))
    }
    
    logger_func_name <- deparse(sys.call()[[1]])
    parent_call <- get_calling_function(c(logger_func_name))
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    msg <- paste0("[", timestamp, "]", "  ", pad(7, log_level), pad(25, parent_call), ...)
    
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

file_exists <- function(path) {
    exists <- file.exists(path)
    if (exists) {
        logger("DEBUG", "File exists at '", path, "'.")
    } else {
        logger("DEBUG", "File does not exist at '", path, "'.")
    }
    return(exists)
}

delete_file <- function(path) {
    if (file_exists(path)) {
        logger("WARN", "Deleting file at path '", path, "'.")
        file.remove(path)
    }
}

wrap_read_table <- function(path, ...) {
    if (!file_exists(path)) {
        logger("ERROR", "Could not find file to read at '", path, "'.")
        return(NULL)
    }

    read.table(path, header=TRUE, ...)
}

run_plink_orig_data <- function(plink_args, out_name = NULL) {
    data_files_pattern <- file.path(data_path, plink_datafile_basename)
    logger("DEBUG", "Plink Data Files Pattern '", data_files_pattern, "',")

    run_plink(data_files_pattern, plink_args, out_name)
}

run_plink <- function(bfile, plink_args, out_name = NULL) { 
    plink_base_cmd <- paste0("plink --bfile ", bfile, " ", plink_args)

    if (is.null(out_name)) {
        # Non-outputting plink command. Output to console
        plink_cmd <- plink_base_cmd
        std_out <- TRUE
    } else {
        # Outputs to file
        out_path <- file.path(plink_out_dir, out_name)
        logger("DEBUG", "Plink out path, '", out_path, "'.")
        plink_cmd <- paste0(plink_base_cmd, " --out ", out_path)
        std_out <- FALSE
    }

    logger("Running: ", plink_cmd)

    if (std_out) {
        system(plink_cmd)
    } else {
        shell_call(plink_cmd)
        logger("Plink results directed to '", out_path, "'.")
    }
}

wrap_histogram <- function(df, col_name, out_path, width = 600, height = 350) {
    delete_file(out_path)
    png(out_path, width, height)
    hist(df[[col_name]], main = paste("Histogram of", col_name), xlab = col_name)
    dev.off()
}

init <- function() {
    logger("INFO", "Initialising directories...")
    mkdir_if_not_exist(plots_out)
    mkdir_if_not_exist(plink_out_dir)
    logger("INFO", "Initialisation complete!")
}

init()

# === Main ===

pheno_path <- file.path(phenotypes, paste0(space_to_underscore(phenotype), pheno_ext))
logger("Reading phenotype path:", pheno_path)
pheno <- read.table(pheno_path)
logger("table read")

# Count the number of individuals in phenotype data
dim(pheno)

head(pheno)

length(which(!is.na(pheno[,3])))

# ====== Analysis ======

# === SNP QC (Missing) ===

missing_name <- "missing"
run_plink("--missing", missing_name)

missing_ind <- wrap_read_table(file.path(plink_out_dir, paste0(missing_name, ".imiss")))
dim(missing_ind)
head(missing_ind)

genotype_threshold <- 0.05
hist_out_path <- file.path(plots_out, "fmiss.png")
wrap_histogram(missing_ind, "F_MISS", hist_out_path)
sum(missing_ind$"F_MISS" > genotype_threshold)

# Filter our SNPs with missingness about threshold
filtered_snps_name <- "filtered"
run_plink_orig_data(paste("--geno", genotype_threshold, "--make-bed"), filtered_snps_name)

# === Sample QC (excluding checking ancestry) ===

# Filter individuals with high missingness
filtered_path <- file.path(plink_out_dir, filtered_snps_name)
filtered_indvs_name <- "filtered_individuals"
run_plink(filtered_path, paste("--mind", genotype_threshold, "--make-bed"), filtered_indvs_name)

# Check for Duplicate SNPs
filtered_indvs_path <- file.path(plink_out_dir, filtered_indvs_name)
dup_vars_name <- "duplicate_vars"
run_plink(filtered_indvs_path, "--chr X --list-duplicate-vars", dup_vars_name)

# Check sex
check_sex_name <- "check_sex"
run_plink(filtered_indvs_path, "--check-sex", check_sex_name)

# === Genome-wide association analysis of the three traits ===

# Quantitative Trait
quant_trait_res_name <- "quantitative_trait_results"
run_plink_orig_data("--assoc", quant_trait_res_name)

# Binary Trait (top 20%)
binary_20_name <- "binary_trait_20_results"
run_plink_orig_data("--pheno recorded_pheno.txt --assoc", binary_20_name)

# Binary Trait (top and bottom 30%)
binary_30_name <- "binary_trait_30_results"
run_plink_orig_data("--pheno recorded_pheno_30.txt --assoc", binary_30_name)

# === Describe the most associated region of the quantitative trait ===

# Manhattan Plot
quant_traits_res_path <- file.path(plink_out_dir, paste0(quant_trait_res_name, ".assoc"))
gwas_results <- wrap_read_table(quant_traits_res_path)
if (!is.null(gwas_results)) {
    manhattan_plot_path <- file.path(plots_out, "manhattan.png")
    wrap_histogram(gwas_results, "P", manhattan_plot_path)
}

# A comparison of the results for these different traits sets
#   demonstrating an understanding of how the results from the
#   three traits relate to each other.

binary_20_path <- file.path(plink_out_dir, paste0(binary_20_name, ".assoc"))
binary_20_results <- wrap_read_table(binary_20_path)

binary_30_path <- file.path(plink_out_dir, paste0(binary_30_name, ".assoc"))
binary_30_results <- wrap_read_table(binary_30_path)

logger("DONE!")

