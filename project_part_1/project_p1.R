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

required_packages <- c("crayon", "qqman")
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
fam_ind_cols <- c("FID", "IID")

# Thresholds
genotype_threshold <- 0.05
het_threshold <- 0.2

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
out_dir <- file.path("./out")

# === Functions ===

cat0 <- function(...) {
    msg <- paste0(...)
    cat(msg, "\n")
}

pad <- function(width, ...) {
    padded <- sprintf(paste0("%-", width, "s"), paste0(...))
    return(padded)
}

quotes <- function(path) {
    paste0("'", path, "'")
}

named_flag <- function(fg_name) {
    paste0("--", fg_name)
}

ext <- function(ext_name) {
    paste0(".", ext_name)
}

create_object <- function(items, transform_fn) {
    result <- list()

    for (item in items) {
        if (is.list(item)) {
            name <- names(item)[1]
            value <- item[[1]]
            result[[name]] <- transform_fn(value)
        } else {
            result[[item]] <- transform_fn(item)
        }
    }

    result
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
        logger("INFO", "Path ", quotes(path), " does not exist. Creating...")
        dir.create(path, recursive = TRUE)
        logger("DEBUG", "Path ", quotes(path), " created.")
    } else {
        logger("DEBUG", "Path ", quotes(path), " already exists.")
    }
}

file_exists <- function(path) {
    logger("DEBUG", "Checking if file exists at path ", quotes(path), ".")
    if (length(path) == 0) {
        logger("WARN", "Path is empty")
        return(FALSE)
    }
    exists <- file.exists(path)
    if (exists) {
        logger("DEBUG", "File exists at ", quotes(path), ".")
    } else {
        logger("DEBUG", "File does not exist at ", quotes(path), ".")
    }
    return(exists)
}

delete_file <- function(path) {
    if (file_exists(path)) {
        logger("WARN", "Deleting file at path ", quotes(path), ".")
        file.remove(path)
    }
}

wrap_read_table <- function(path, header = TRUE, ...) {
    if (!file_exists(path)) {
        logger("ERROR", "Could not find file to read at ", quotes(path), ".")
        return(NULL)
    }

    logger("DEBUG", "Reading table at ", quotes(path), ".")
    read.table(path, header = header, ...)
}

wrap_write_table <- function(data, path, row.names = FALSE, ...) {
    if (file_exists(path)) {
        logger("WARN", "Overwriting file at path ", quotes(path), ".")
        delete_file(path)
    }

    logger("DEBUG", "Writing table at ", quotes(path), ".")
    write.table(data, path, row.names = row.names, sep = "\t", ...)
}

plink_orig_data <- function(plink_args, out_name = NULL) {
    data_files_pattern <- file.path(data_path, plink_datafile_basename)
    logger("DEBUG", "Plink Data Files Pattern ", quotes(data_files_pattern), ".")

    plink(data_files_pattern, plink_args, out_name)
}

plink <- function(bfile, plink_args, out_name = NULL) { 
    plink_base_cmd <- paste("plink", pl_fgs$bfile, bfile, plink_args)

    if (is.null(out_name)) {
        # Non-outputting plink command. Output to console
        plink_cmd <- plink_base_cmd
        std_out <- TRUE
    } else {
        # Outputs to file
        out_path <- file.path(plink_out_dir, out_name)
        logger("DEBUG", "Plink out path, ", quotes(out_path), ".")
        plink_cmd <- paste(plink_base_cmd, pl_fgs$out, out_path)
        std_out <- FALSE
    }

    logger("Running: ", quotes(plink_cmd))

    if (std_out) {
        system(plink_cmd)
    } else {
        shell_call(plink_cmd)
        logger("Plink results directed to ", quotes(out_path), ".")
    }
}

construct_plink_table_path <- function(name, ext) {
    path <- file.path(plink_out_dir, paste0(name, ext))
    logger("DEBUG", "Constructed plink out path: ", quotes(path), ".")
    return(path)
}

remove_indices <- function(missing_indices, thresh_col_name, threshold, out_col_names, out_name) {
    ind_to_remove <- which(missing_indices[[thresh_col_name]] > threshold)
    file <- missing_indices[ind_to_remove, out_col_names]
    missing_ind_file_path <- file.path(plink_out_dir, out_name)
    wrap_write_table(file, missing_ind_file_path, col.names = FALSE, quote = FALSE)
    return(missing_ind_file_path)
}

check_ext <- function(out_name, expected_ext, add_if_missing = TRUE) {
    if (length(expected_ext) == 0) {
        logger("ERROR", "Expected extension is blank")
        return(NULL)
    }
    logger("DEBUG", "Checking ", quotes(out_name), " for extension ", quotes(expected_ext), ".")
    if (!endsWith(out_name, expected_ext)) {
        logger("WARN", "Out name ", quotes(out_name), " does not end with ", quotes(expected_ext), ".")
        if (add_if_missing) {
            logger("INFO", "Adding ", quotes(expected_ext), " to file name")   
            out_name <- paste0(out_name, expected_ext)
        }
    } else {
        logger("DEBUG", "Out name ", quotes(out_name), " contains extension ", quotes(expected_ext), ".")
    }

    return(out_name)
}

check_png_ext <- function(out_name, add_if_missing = TRUE) {
    check_ext(out_name, exts$png, add_if_missing)
}

check_txt_ext <- function(out_name, add_if_missing = TRUE) {
    check_ext(out_name, exts$txt, add_if_missing)
}

log_indvs_to_remove <- function(num_indvs) {
    if (num_indvs == 0) {
        logger("INFO", "No individuals to remove")
    } else {
        logger("WARN", "There are ", num_indvs, " individuals to remove")
    }
}

wrap_plot <- function(plot_callback, data, out_name, width = 600, height = 350, ...) {
    out_name <- check_png_ext(out_name, TRUE)
    out_path <- file.path(plots_out, out_name)
    delete_file(out_path)
    png(out_path, width, height)
    plot_callback(data, ...)
    dev.off()
}

wrap_histogram <- function(...) {
    wrap_plot(hist, ...)
}

wrap_scatter <- function(abline_h, abline_col, ...) {
    plot_with_abline <- function(...) {
        plot(...)
        abline(abline_h, abline_col)
    }

    wrap_plot(plot_with_abline, ...)
}

init <- function() {
    logger("INFO", "Initialising directories...")
    mkdir_if_not_exist(plots_out)
    mkdir_if_not_exist(plink_out_dir)
    mkdir_if_not_exist(out_dir)

    logger("INFO", "Initialisation complete!")
}

init()

# Plink Flags
pl_fgs <- create_object(list("remove", "missing", list("mb" = "make-bed"), 
                             "hardy", "het", "mind", "pheno", "covar", 
                             list("dup_vars" = "list-duplicate-vars"), "out", 
                             "bfile", "chr"), 
                        named_flag)

exts <- create_object(list("phen", "imiss", "lmiss", "het", "assoc", "hwe", "txt", "png"), ext)

# === Main ===

# ====== Analysis ======

quality_control <- function() {
    logger("Performing Quality Control")

    excess_missing_genotypes <- function(histogram) {
        logger("Checking for missing genotypes")

        missing_name <- "missing"
        plink_orig_data(pl_fgs$missing, missing_name)
        
        logger("DEBUG", "Reading imiss table...")
        missing_out_path <- construct_plink_table_path(missing_name, exts$imiss)
        missing_ind <- wrap_read_table(missing_out_path)
        dim(missing_ind)
        head(missing_ind)
        
        if (histogram) {
            logger("INFO", "Plotting F_MISS histogram")
            wrap_histogram(missing_ind$"F_MISS", "fmiss.png")
            num_indvs_to_remove <- sum(missing_ind$"F_MISS" > genotype_threshold)
            log_indvs_to_remove(num_indvs_to_remove)
       }

        missing_ind_file_path <- remove_indices(missing_ind, "F_MISS", genotype_threshold, fam_ind_cols, "remove.missing.samples.txt")
        return(missing_ind_file_path)
    }

    outlying_homozygosity <- function(plot) {
        logger("Checking for outlying homozygosity values")

        hz_name <- "hz"
        plink_orig_data(pl_fgs$het, hz_name)
        
        logger("DEBUG", "Reading het table...")
        het_out_path <- construct_plink_table_path(hz_name, exts$het)
        het <- wrap_read_table(het_out_path)
        dim(het)
        head(het)

        if (plot) {
            logger("INFO", "Plotting het histogram")
            wrap_histogram(het$"F", "fhet_hist.png")

            logger("INFO", "Plotting het scatterplot")
            wrap_scatter(0.05, "red", abs(het$"F"), "fhet_scatter.png")
        }

        het_ind_file_path <- remove_indices(het, "F", het_threshold, fam_ind_cols, "remove.het.samples.txt")
        return(het_ind_file_path)
    }

    related_samples <- function() {

    }

    ancestry_outliers <- function() {

    }

    combine_remove_files <- function(...) {
        logger("Combining Remove Files...")

        file_paths <- list(...)

        combined_ind <- data.frame(FID = numeric(), IID = numeric())

        for (file_path in file_paths) {
            ind <- wrap_read_table(file_path, FALSE)
            colnames(ind) <- fam_ind_cols
            combined_ind <- rbind(combined_ind, ind)
        }
        
        # Remove Duplicates
        combined_ind <- unique(combined_ind)

        combined_file_out_path <- file.path(plink_out_dir, "remove.combined.samples.txt")
        wrap_write_table(combined_ind, combined_file_out_path, col.names = FALSE, quote = FALSE)
        return(combined_file_out_path)
    }


    remove_bad_individuals <- function(remove_path) {
        logger("Removing bad samples...")
        out_name <- "test_qc"

        plink_flags <- paste(pl_fgs$mb, pl_fgs$remove, remove_path)
        plink_orig_data(plink_flags, out_name)
        
        out_path <- file.path(plink_out_dir, out_name)
        return(out_path)
    }

    missing_ind_path <- excess_missing_genotypes(TRUE)
    het_ind_path <- outlying_homozygosity(FALSE)

    combined_ind_path <- combine_remove_files(missing_ind_path, het_ind_path)
    qc_data_path <- remove_bad_individuals(combined_ind_path)
    return(qc_data_path)
}

sample_qc <- function(qc_data_path) {
    hw_eq <- function() {
        hw_eq_name <- "hw_eq"
        plink_flags <- paste(pl_fgs$hardy)
        plink(qc_data_path, plink_flags, hw_eq_name)
        
        hw_eq_path <- construct_plink_table_path(hw_eq_name, exts$hwe)
        hw <- wrap_read_table(hw_eq_path)
        
        hw_eq_ind_file_path <- remove_indices(hw, c("MAF", "P"), threshold, 2, "remove.SNPs.txt")
    }

    # Filter individuals with high missingness
    filtered_path <- file.path(plink_out_dir, "filtered")
    filtered_indvs_name <- "filtered_individuals"
    plink(filtered_path, paste(pl_fgs$mind, genotype_threshold, pl_fgs$mb), filtered_indvs_name)

    logger("Checking Duplicate SNPs")
    filtered_indvs_path <- file.path(plink_out_dir, filtered_indvs_name)
    dup_vars_name <- "duplicate_vars"
    plink(filtered_indvs_path, paste(pl_fgs$chr, "X", pl_fgs$dup_vars), dup_vars_name)
}

gwas <- function() {
    phenotype_file_prefix <- space_to_underscore(phenotype)

    get_pheno_path <- function(pheno_suffix) {
        pheno_file_name <- paste0(phenotype_file_prefix, pheno_suffix, exts$pheno)
        file.path(phenotypes, pheno_file_name)
    }

    trait_analysis <- function(include_covariates) {
        logger("Performing trait analysis with covariates: ", include_covariates)

        get_plink_args <- function(pheno_path) {
            base_args <- paste(pl_fgs$pheno, pheno_path, pl_assoc)
            if (include_covariates) {
                age_path <- file.path(data_path, "age.txt")
                gender_path <- file.path(data_path, "gender.txt")
                
                base_args <- paste(base_args, pl_fgs$covar, age_path, pl_covar, gender_path)
            }
            return(base_args)
        }

        # Quantitative Trait Analysis
        pheno_cont_path <- get_pheno_path("")
        plink_orig_data(get_plink_args(pheno_cont_path), "quant_trait_results")

        # Binary Trait (top 20%) Analysis
        binary1_path <- get_pheno_path("_binary1")
        plink_orig_data(get_plink_args(binary1_path), "binary_trait_20_results")

        # Binary Trait (top and bottom 30%) analysis
        binary2_path <- get_pheno_path("_binary2")
        plink_orig_data(get_plink_args(binary2_path), "binary_trait_30_results")
    }

    trait_analysis(TRUE)
}

run_analysis <- function() {
    # === SNP QC (Missing) ===
    genotype_threshold <- 0.05
    snp_qc(genotype_threshold)

    # === Sample QC (excluding checking ancestry) ===
    sample_qc(genotype_threshold)

    # === Genome-wide association analysis of the three traits ===
    gwas()
   
    # === Describe the most associated region of the quantitative trait ===

    # Manhattan Plot
    quant_traits_res_path <- construct_plink_table_path(quant_trait_res_name, exts$assoc)
    gwas_results <- wrap_read_table(quant_traits_res_path)
    if (!is.null(gwas_results)) {
        wrap_histogram(gwas_results$"P", "gwas_manhattan.png")
    }

    # A comparison of the results for these different traits sets
    #   demonstrating an understanding of how the results from the
    #   three traits relate to each other.

    binary_20_path <- construct_plink_table_path(binary_20_name, exts$assoc)
    binary_20_results <- wrap_read_table(binary_20_path)

    binary_30_path <- construct_plink_table_path(binary_30_name, exts$assoc)
    binary_30_results <- wrap_read_table(binary_30_path)

}

read_covariates <- function() {
    age_path <- file.path(data_path, "age.txt")
    gender_path <- file.path(data_path, "gender.txt")

    age <- wrap_read_table(age_path, header = FALSE, col.names = c(fam_ind_cols, "Age"))
    gender <- wrap_read_table(gender_path, header = FALSE, col.names = c(fam_ind_cols, "Gender"))

    return(list(age, gender))
}

read_phenotypes <- function() {
    phenotype_file_prefix <- space_to_underscore(phenotype)

    get_pheno_path <- function(pheno_suffix) {
        pheno_file_name <- paste0(phenotype_file_prefix, pheno_suffix, exts$pheno)
        file.path(phenotypes, pheno_file_name)
    }

    pheno_cont_path <- get_pheno_path("")
    pheno <- wrap_read_table(pheno_cont_path, col.names = c(fam_ind_cols, "Glucose"))

    binary1_path <- get_pheno_path("_binary1")
    binary1 <- wrap_read_table(binary1_path, col.names = c(fam_ind_cols, "Binary1"))
    binary2_path <- get_pheno_path("_binary2")
    binary2 <- wrap_read_table(binary2_path, col.names = c(fam_ind_cols, "Binary2"))

    return(list(pheno, binary1, binary2))
}

compare_minor_allele_freq <- function() {
    snp_freq_path <- file.path(data_path, "reference_allele_frequencies.txt")
    ref <- wrap_read_table(snp_freq_path)
    dim(ref)
    head(ref)

    ind <- match(freq$SNP, ref$V1)
    out <- cbind(freq, ref[ind,])
    head(out)

    res <- out$MAF - out$V2

    wrap_histogram(res, "minor_allele.png")
}

#phenotypes <- read_phenotypes()
#covariates <- read_covariates()

#run_analysis()

qc_data_path <- quality_control()
sample_qc(qc_data_path)

logger("DONE!")

