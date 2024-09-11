# Course: STAT3306: Statistical Analysis of Genetic Data
# Taken: Semester 2, 2024
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
allowed_log_levels <- c("TRACE", "DEBUG", "INFO", "WARN", "ERROR")
app_log_level <- "TRACE"
level_colours <- list(
    TRACE = crayon::silver,
    DEBUG = crayon::magenta,
    INFO = crayon::cyan,
    WARN = crayon::yellow,
    ERROR = crayon::bgRed
)

# === Globals ===

phenotype <- "Fasting Glucose"
fam_ind_cols <- c("FID", "IID")

# Thresholds
genotype_threshold <- 0.05
het_threshold <- 0.2
related_threshold <- 0.025  # Plink Default
hwe_threshold <- 0.001
freq_threshold <- 0.01

# Data Paths
data_folder <- file.path("/data/STAT3306")
project_data <- file.path(data_folder, "Project")
data_path <- file.path(project_data, "Data")
phenotypes <- file.path(project_data, "Phenotypes")

# Plink
plink_datafile_basename <- "test"

# Out Paths
plots_out_dir <- file.path("./plots")
plink_out_dir <- file.path("./plink_out")
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
    logger("TRACE", "Checking if file exists at path ", quotes(path), ".")
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

check_any_empty <- function(...) {
    args <- list(...)
    logger("DEBUG", "Checking if empty: ", quotes(args))
    any(sapply(args, function(x) x == ""))
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
    path <- check_txt_ext(path, exts$txt)
    if (file_exists(path)) {
        logger("WARN", "Overwriting file at path ", quotes(path), ".")
        delete_file(path)
    }

    logger("DEBUG", "Writing table at ", quotes(path), ".")
    write.table(data, path, row.names = row.names, sep = "\t", ...)
}

plink_orig_data <- function(plink_args, out_name = NULL) {
    data_files_pattern <- file.path(data_path, plink_datafile_basename)
    logger("TRACE", "Plink Data Files Pattern ", quotes(data_files_pattern), ".")

    plink(data_files_pattern, plink_args, out_name)
}

plink <- function(bfile, plink_args, out_name = NULL) { 
    # Check if bfile is null
    if (is.null(bfile)) {
        logger("DEBUG", "Using original data")
        return(plink_orig_data(plink_args, out_name))
    } else {
        logger("DEBUG", "Using data file: ", quotes(bfile), ".")
    }
    
    plink_base_cmd <- paste("plink", pl_fgs$bfile, bfile, plink_args)
    logger("TRACE", "Plink base command ", quotes(plink_base_cmd), ".")

    if (is.null(out_name)) {
        # Non-outputting plink command. Output to console
        plink_cmd <- plink_base_cmd
        std_out <- TRUE
    } else {
        # Outputs to file
        plink_out_path <- construct_plink_out_path(out_name)
        logger("TRACE", "Plink out path: ", quotes(plink_out_path), ".")
        plink_cmd <- paste(plink_base_cmd, pl_fgs$out, plink_out_path)
        std_out <- FALSE
    }
    
    # Check double space - indicates missing param
    if (grepl("  ", plink_cmd)) {
        logger("ERROR", "Possible missing argument in plink commmand")
        logger("ERROR", "You were trying to run ", quotes(plink_cmd))
    }

    logger("Running: ", quotes(plink_cmd))

    if (std_out) {
        system(plink_cmd)
    } else {
        shell_call(plink_cmd)
        logger("Plink results directed to ", quotes(plink_out_path), ".")
    }

    return(plink_out_path)
}

add_extension <- function(basename, ...) {
    args <- list(...)
    logger("TRACE", "Checking extensions ", ..., ".")
    if (check_any_empty(...)) {
        logger("ERROR", "Some extensions are empty in ", ..., ".")
    }

    logger("TRACE", "Adding extensions: ", ..., " to basename ", quotes(basename), ".")
    path <- paste0(basename, ...)
    logger("TRACE", "Constructed plink out path: ", path, ".")
    return(path)
}

construct_plink_out_path <- function(...) {
    #' Constructs a file path in the plink out directory given a basename
    #' @param ... {string}: The name of the file (potentially split) without
    #'                      the plink out directory prefixed.
    #' @return path {string}: The relative path to the file.

    basename <- paste0(...)
    file.path(plink_out_dir, basename)
}


construct_out_path <- function(basename) {
    #' Constructs a file path in the out directory given a basename
    #' @param basename {string}: The name of the file without the out directory prefixed.
    #' @return path {string}: The relative path to the file.

    file.path(out_dir, basename)
}

save_removed_indices <- function(table, ind_to_remove, out_cols, out_name) {
    logger("INFO", "Saving removed indices to ", quotes(out_name), ".")
    logger("DEBUG", "Selected columns: ", quotes(out_cols), ".")

    file <- table[ind_to_remove, out_cols]
    ind_out_path <- construct_out_path(out_name)
    wrap_write_table(file, ind_out_path, col.names = FALSE, quote = FALSE)
    return(ind_out_path)
}

remove_indices_by_threshold <- function(table, thresh_col_name, threshold, out_cols, out_name) {
    logger("DEBUG", "Removing indices by threshold (", threshold, ") in table from col name: ",
           quotes(thresh_col_name), ".")
    ind_to_remove <- which(table[[thresh_col_name]] > threshold)
    save_removed_indices(table, ind_to_remove, out_cols, out_name)
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
            logger("DEBUG", "Adding ", quotes(expected_ext), " to file name")   
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
    out_path <- file.path(plots_out_dir, out_name)
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

genomic_inflation_factor <- function(d, df = 1) {
    logger("DEBUG", "Computing genomic inflation factor with ", df, " degree of freedom...")
    
    median_val <- median(d$P)
    logger("TRACE", "Median Value: ", median_val)
    
    chisq_med <- qchisq(0.5, 1)
    qchisq(1 - median_val, df) / chisq_med
}

excess_missing_genotypes <- function(data_path, extension, out_cols, histogram) {
    logger("Checking for missing genotypes (", extension, ")...")

    missing_name <- "missing"
    missing_basename <- plink(data_path, pl_fgs$missing, missing_name)
    
    logger("DEBUG", "Reading ", extension, " table...")
    missing_out_path <- add_extension(missing_basename, extension)
    missing <- wrap_read_table(missing_out_path)
    dim(missing)
    head(missing)
    
    if (histogram) {
        logger("INFO", "Plotting F_MISS histogram")
        wrap_histogram(missing$"F_MISS", paste0("fmiss", extension, exts$png))
        num_indvs_to_remove <- sum(missing$"F_MISS" > genotype_threshold)
        log_indvs_to_remove(num_indvs_to_remove)
    }
    
    out_name <- paste0("remove.missing", extension, ".samples.txt")
    missing_file_path <- remove_indices_by_threshold(missing, "F_MISS", genotype_threshold, 
                                                     out_cols, out_name)
    return(missing_file_path)
}

####

init <- function() {
    logger("INFO", "Initialising directories...")
    mkdir_if_not_exist(plots_out_dir)
    mkdir_if_not_exist(plink_out_dir)
    mkdir_if_not_exist(out_dir)

    logger("INFO", "Initialisation complete!")
}

init()

# Plink Flags
pl_fgs <- create_object(list("remove", "missing", list("mb" = "make-bed"), 
                             "hardy", "het", "mind", "pheno", "covar", 
                             list("dup_vars" = "list-duplicate-vars"), "out", 
                             "bfile", "chr", "freq", "exclude", "mpheno",
                             "pca", "linear", "assoc", "clump", list("cp1" = "clump-p1"), 
                             list("cp2" = "clump-p2"), list("cr2" = "clump-r2"), 
                             list("ckb" = "clump-kb"), list("rel_cutoff" = "rel-cutoff"),
                             "keep"), 
                        named_flag)

exts <- create_object(list("phen", "imiss", "lmiss", "het", "assoc", "hwe", 
                           "frq", "txt", "png", "eigenvec", "eigenval",
                           "qassoc", "linear", "clumped", "rel", "id"), 
                      ext)

# === Main ===

# ====== Analysis ======

quality_control <- function() {
    logger("Performing Quality Control")

    find_individual_missing_genotypes <- function(histogram) {
        #' Finds individuals with excess missing genotypes to remove from the dataset
        #' @param histogram {boolean}: if true, will plot a histogram of the frequency of missingness
        #' @return missing_file_path {string}: The file path to the file containing the individuals to remove.

        excess_missing_genotypes(NULL, exts$imiss, fam_ind_cols, histogram)
    }

    find_outlying_homozygosity <- function(plot) {
        #' Finds individuals with outlying homozygosity values to remove from the dataset
        #' @param plot {boolean}: if true, will plot a histogram of the frequency of hz freqs and a scatterplot
        #'                        of their distribution across the geneome.
        #' @return het_ind_file_path {string}: The file path to the file containing the individuals to remove.

        logger("Checking for outlying homozygosity values...")

        hz_name <- "hz"
        het_basename <- plink_orig_data(pl_fgs$het, hz_name)
        
        logger("DEBUG", "Reading het table...")
        het_out_path <- add_extension(het_basename, exts$het)
        het <- wrap_read_table(het_out_path)
        dim(het)
        head(het)

        if (plot) {
            logger("INFO", "Plotting het histogram")
            wrap_histogram(het$"F", "fhet_hist.png")

            logger("INFO", "Plotting het scatterplot")
            wrap_scatter(0.05, "red", abs(het$"F"), "fhet_scatter.png")
        }

        het_ind_file_path <- remove_indices_by_threshold(het, "F", het_threshold, 
                                                         fam_ind_cols, "remove.het.samples.txt")
        return(het_ind_file_path)
    }

    find_related_samples <- function(threshold = related_threshold) {
        #' Finds individuals from each pair of samples with observed genomic relatedness above a given threshold.
        #' NP-Hard problem so will used cached result if exists.
        #' @param threshold {float}: The threshold to exclude related samples by.
        #' @return related_file_path {string}: File path to a file containing individuals to remove.

        basename_with_extension <- function(basename) {
            add_extension(basename, exts$rel, exts$id)
        }

        logger("Checking for related samples...")
        
        related_name <- "related_samples"
        expected_related_path <- construct_plink_out_path(basename_with_extension(related_name))
        logger("DEBUG", "Checking for cached result at path ", quotes(expected_related_path), "...")
        if (file_exists(expected_related_path)) {
            logger("INFO", "Found cached related samples result")
            return(expected_related_path)
        }

        logger("INFO", "Cached related samples not found. Computing...")

        # Note this file is individuals to KEEP (not remove)
        related_basename <- plink_orig_data(pl_fgs$rel_cutoff, related_name)
        related_path <- basename_with_extension(related_basename)
        return(related_path)
    }

    find_ancestry_outliers <- function() {
        # We can ignore as per task sheet
    }

    combine_remove_files <- function(...) {
        #' Combines the individuasl to remove into a single txt file
        #' @param ... {character}: File paths to individual files containing IDs of individuals to remove.
        #'                         These files are expected to have at least two columns: FID and IID.
        #' @return {string}: File path to the combined output file with individuals to remove.

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

        combined_file_out_path <- construct_out_path("remove.combined.samples.txt")
        wrap_write_table(combined_ind, combined_file_out_path, col.names = FALSE, quote = FALSE)
        return(combined_file_out_path)
    }


    remove_bad_individuals <- function(remove_path) {
        #' Removes individuals from the dataset who are specified in the provided file.
        #' @param remove_path {string}: File path to file containing individuals to remove from the dataset.
        #' @return out_path {string}: File path to subset of original dataset with specified individuals removed.

        logger("Removing bad samples...")
        out_name <- "test_indv_subset"

        plink_flags <- paste(pl_fgs$mb, pl_fgs$remove, remove_path)
        plink_orig_data(plink_flags, out_name)
        
        out_path <- construct_plink_out_path(out_name)
        return(out_path)
    }

    keep_related_samples <- function(related_path, data_subset_path) {
        #' Keeps only related samples from the data subset who are listed in the related file
        #' @param related_path {string}: Path to file containing list of individuals to keep.
        #' @param data_subset_path {string}: Path to subset of data with some individuals already removed.
        #' @return out_path {string}: Path to a further subset of the data

        logger("Retaining unrelated samples...")

        plink_args <- paste(pl_fgs$keep, related_path, pl_fgs$mb)
        out_name <- "test_indv_related"
        plink(data_subset_path, plink_args, out_name)

        out_path <- construct_plink_out_path(out_name)
        return(out_path)
    }

    missing_ind_path <- find_individual_missing_genotypes(TRUE)
    het_ind_path <- find_outlying_homozygosity(FALSE) # Graph to fix

    combined_ind_path <- combine_remove_files(missing_ind_path, het_ind_path)
    data_subset_path <- remove_bad_individuals(combined_ind_path)

    # Keep only individuals specified here
    related_path <- find_related_samples()
    data_related_subset_path <- keep_related_samples(related_path, data_subset_path)

    return(data_related_subset_path)
}

sample_qc <- function(data_subset_path) {
    missing_snps <- function(histogram) {
        missing_file_path <- excess_missing_genotypes(data_subset_path, exts$lmiss, 2, histogram)
        lmiss <- wrap_read_table(missing_file_path)
        return(lmiss)
    }

    hw_eq <- function(histogram) {
        logger("Computing Hardy-Weinberg Equilibrium")
        hwe_name <- "hwe"
        hwe_basename <- plink(data_subset_path, pl_fgs$hardy, hwe_name)
       
        hwe_path <- add_extension(hwe_basename, exts$hwe)
        hwe <- wrap_read_table(hwe_path)
        
        if (histogram) {
            wrap_histogram(hwe, "hwe_deviations.png")
        }

        return(hwe) 
    }

    min_allele_freq <- function(histogram) {
        #' Indentify SNPs with low minor allele frequency
        #' @param histogram {boolean}: if true, frequency of MAF is plotted visually.
        #' @return dataframe of low minor allele frequencies

        logger("Computing Minor Allele Frequencies")
        min_allele_name <- "minor_allele_freq"
        min_allele_basename <- plink(data_subset_path, pl_fgs$freq, min_allele_name)

        min_allele_path <- add_extension(min_allele_basename, exts$frq)
        freq <- wrap_read_table(min_allele_path)

        if (histogram) {
            wrap_histogram(freq, "maf_distribution.png")
        }

        return(freq)
    }

    remove_snps <- function(lmiss, hwe, freq) {
        logger("INFO", "Computing SNPs to remove...")
        
        ind_to_remove <- unique(c(
                                    which(lmiss$F_MISS > genotype_threshold,
                                    freq$MAF < freq_threshold), 
                                    which(hwe$P < hwe_threshold)
                                 )
                               )
        snp_id_col <- 2
        remove_snps_path <- save_removed_indices(freq, ind_to_remove, snp_id_col, "remove.SNPs.txt")
        return(remove_snps_path)
    }

    exclude_snps <- function(remove_snps_path) {
        logger("INFO", "Excluding SNPs...")
        out_name <- "test_qc"
        plink_args <- paste(pl_fgs$mb, pl_fgs$exclude, remove_snps_path)
        plink(data_subset_path, plink_args, out_name)
    }

    compare_minor_allele_freqs <- function(freq, do_plot) {
        logger("Comparing Minor Allele Frequencies...")
        snp_ref_path <- file.path(data_path, "reference_allele_frequencies.txt")
        ref <- wrap_read_table(snp_ref_path)
        dim(ref)
        head(ref)

        ind <- match(freq$SNP, ref$V1)
        out <- cbind(freq, ref[ind,])
        print(head(out))
        # Valentine says need to do 1 - MAF here otherwise will get X (cross)
        # See lcys as
        if (do_plot) {
            logger("Plotting allele frequency comparison with reference...")
            wrap_plot(plot, out$MAF ~ out$V2, "min_allele_freq_comparison.png")
        }

        res <- out$MAF - out$V2
        cat0(head(res))
        allele_freq_threshold <- 0.1
        keep <- c(which(abs(res) <= allele_freq_threshold))
        wrap_histogram(res, "minor_allele.png")
    }
    
    lmiss <- missing_snps(TRUE)
    hwe <- hw_eq(TRUE)
    freq <- min_allele_freq(TRUE)
    
    remove_snps_path <- remove_snps(lmiss, hwe, freq)
    qc_data_path <- exclude_snps(remove_snps_path)

    compare_minor_allele_freqs(freq, TRUE)
    
    return(qc_data_path)
}

gwas <- function(qc_data_path) {
    logger("INFO", "Performing GWAS...")

    get_pheno_path <- function(pheno_suffix) {
        logger("DEBUG", "Retrieving pheno path for suffix ", quotes(pheno_suffix), ".")
        phenotype_file_prefix <- space_to_underscore(phenotype)
        pheno_file_name <- paste0(phenotype_file_prefix, pheno_suffix, exts$phen)
        file.path(phenotypes, pheno_file_name)
    }

    gwas_pheno <- function() {
        logger("Performing Pheno on GWAS...")
        
        alt_mpheno <- 1
        plink_args <- paste(pl_fgs$assoc, pl_fgs$pheno, pheno_path, mpheno_args)
        out_name <- "gwas_pheno_1"
        plink(qc_data_path, plink_args, out_name)
    }

    gwas_plots <- function(d, plots, plot_suffix = "") {
        name_plot <- function(plot_type) {
            file_name <- paste0("gwas_", plot_type)
            if (length(plot_type) > 0) {
                if (plot_type[[1]] != "_") {
                   file_name <- paste0(file_name, "_")
                }

                file_name <- paste0(file_name, plot_suffix)
            }
            
            file_name <- paste0(file_name, ".png")
            logger("DEBUG", "Plot Name: ", quotes(file_name), ".")
        }

        if (plots) {
            man_plot_name <- name_plot("manhattan")
            logger("INFO", "Generating Manhattan Plot ", quotes(man_plot_name), " ...")
            wrap_plot(manhattan, d, man_plot_name)

            qq_plot_name <- name_plot(qq)
            logger("INFO", "Generating QQ Plot ", quotes(qq_plot_name), " ...")
            wrap_plot(qq, d$P, qq_plot_name)
        }

        logger("INFO", "Computing Genomic Inflation Factor (", plot_suffix, ") ...")
        lambda_gc <- genomic_inflation_factor(d, 1)
        logger("INFO", "Lambda_{GC} (", plot_suffix, ") = ", lambda_gc)
        return(lambda_gc)
    }

    check_inflation_factor <- function(pheno_basename) { 
        pheno_qassoc_path <- paste0(pheno_basename, exts$qassoc)
        d <- wrap_read_table(pheno_path)
        gwas_plots(d, FALSE)
    }

    compute_principal_comps <- function(num_components) {
        # Check for existing PCA
        out_name <- "pca"
        pca_path <- construct_plink_out_path(out_name)   
        pca_eig_val <- add_extension(pca_path, exts$eigenval)
        pca_eig_vec <- add_extension(pca_path, exts$eigenvec)

        eig_files <- list(pca_eig_val, pca_eig_vec)
        if (all(sapply(eig_files, file_exists))) {
            logger("INFO", "PCA already exists. Skipping.")
            return(pca_eig_vec)
        }

        # No PCA. Compute:
        logger("Performing Principal Component Analysis on Data QC with ", num_components, 
               " components...")

        plink_args <- paste(pl_fgs$pca, num_components)
        plink(qc_data_path, plink_args, out_name)
        return(pca_eiv_vec)
    }

    add_pc_covariates <- function(pheno_path, pc_eigvec_file) {
        # Check for existing Covariates
        out_name <- "gwas_pheno_1_pc"
        cov_out_path <- construct_plink_out_path(out_name)
        pheno_pc_path <- add_extension(cov_out_path, exts$assoc, exts$linear)
        if (file_exists(pheno_pc_path)) {
            logger("INFO", "Covariates already exist. Skipping.")
            return (pheno_pc_path)
        }
        
        # No covariates. Compute:
        logger("Adding principal components to covariates...")
        logger("DEBUG", "Loading eigenvec file: ", quotes(pc_eigvec_file), ".")
        check_ext(pc_eigvec_file, exts$eigenvec)
        plink_args <- paste(pl_fgs$linear, pl_fgs$covar, pc_eigvec_file, pl_fgs$pheno, pheno_path, mpheno_args)
        plink(qc_data_path, plink_args, out_name)
        return(pheno_pc_path)
    }

    check_pc_inflation_factor <- function(pheno_pc_path) {
        d <- wrap_read_table(pheno_pc_path)

        d_sub <- d[which(d$TEST == "ADD"),]
        cat0(head(d_sub))

        gwas_plots(d_sub, FALSE, "pc")
    }

    clumping <- function(pheno_pc_path) {
        logger("Clumping GWAS Results...")
        out_name <- "gwas_pheno_1_clump"
        clump_p1_val <- 0.5
        clump_p2_val <- 0.5
        clump_r2_val <- 0.2
        clump_kb_val <- 500
        plink_args <- paste(pl_fgs$clump, pheno_pc_path, pl_fgs$cp1, clump_p1_val, pl_fgs$cp2,
                            clump_p2_val, pl_fgs$cr2, clump_r2_val, pl_fgs$ckb, clump_kb_val)
        plink(qc_data_path, plink_args, out_name)
    }

    read_clumps <- function(clump_basename) {
        clump_path <- add_extension(clump_basename, exts$clumped)
        clump <- wrap_read_table(clump_path)

        # Exclude column PC2
        logger("Printing Clump Output")
        clump_out <- clump[, 1:11]
        print(head(clump_out))
        logger("DEBUG", "End of clump output")
        
        # Write to file
        out_path <- construct_out_path("clumps.txt")
        wrap_write_table(clump_out, out_path)
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
    
    # Some scoped variables
    pheno_path <- get_pheno_path("")
    alt_mpheno <- 1
    mpheno_args <- paste(pl_fgs$mpheno, alt_mpheno)
    
    # Main
    pheno_basename <- gwas_pheno()
    check_inflation_factor(pheno_basename)
    pc_eigvec_file <- compute_principal_comps(10)
    pheno_1_pc_path <- add_pc_covariates(pheno_path, pc_eigvec_file)
    check_pc_inflation_factor(pheno_1_pc_path)
    clump_path <- clumping(pheno_1_pc_path)
    read_clumps(clump_path)

    
    #trait_analysis(TRUE)
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
        pheno_file_name <- paste0(phenotype_file_prefix, pheno_suffix, exts$phen)
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

#phenotypes <- read_phenotypes()
#covariates <- read_covariates()

#run_analysis()

data_subset_path <- quality_control()
qc_data_path <- sample_qc(data_subset_path)
gwas(qc_data_path)

logger("DONE!")

