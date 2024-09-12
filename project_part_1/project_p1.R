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
maf_threshold <- 0.1

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

# Overrides
overwrite_ext_plots <- FALSE
overwrite_plink_out <- FALSE

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

list_to_str <- function(lst) {
    paste(lst, collapse = ", ")
}

named_flag <- function(fg_name) {
    paste0("--", fg_name)
}

ext <- function(ext_name) {
    paste0(".", ext_name)
}

wrap_dim <- function(df) {
    #' Gets the dimension of a data frame in a string format.
    #' @param df {data.frame}: Data frame to get the dimension of.
    #' @return {string}: Dimension of data frame in string representation.

    paste(dim(df), collapse = " x ")
}

`?` <- function(x, y) {
    #' A little naughty. We can make R use ? as a ternary operator
    #' just like C! It even works with assignment <- !
    #' @param x {any}: the condition (before the ?).
    #' @param y {any}: both outcomes (separated by :).
    #' @return the outcome based on the condition prior to the `?`.
                              
    xs <- as.list(substitute(x))
    logger("TRACE", "? Eval xs: ", quotes(x), ".")

    # Get context of caller
    env <- parent.frame()

    if (xs[[1]] == as.name("<-")) {
        logger("TRACE", "Assignment detected")
        x <- eval(xs[[3]], envir = env)
    }

    # Parse and eval ternary outcome in parent env
    outcomes <- sapply(strsplit(deparse(substitute(y)), ":"), function(e) parse(text = e))
    r <- eval(outcomes[[2 - as.logical(x)]], envir = env)

    # Handle assignment
    if (xs[[1]] == as.name("<-")) {
        xs[[3]] <- r
        eval.parent(as.call(xs))
    } else {
        r
    }
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

log_df <- function(df, name, log_level = "INFO") {
    #' Logs the dimension and head of a dataframe at the specified log level,
    #' @param df {data.frame}: The data frame to preview,
    #' @param name {string}: The name of the data frame.
    #' @param log_level {string}: The level to log at. Defaults to INFO.
    #' @return {NULL}

    logger(log_level, "Dim ", name, ": ", wrap_dim(df))
    print(head(df))
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

file_exists <- function(path, partial_match = FALSE) {
    logger("TRACE", "Checking if file exists at path ", quotes(path), ".")
    if (length(path) == 0) {
        logger("WARN", "Path is empty")
        return(FALSE)
    }

    if (partial_match) {
        logger("DEBUG", "Searching partial match")
        dir_name <- dirname(path)
        file_pattern <- basename(path)

        matching_files <- list.files(path = dir_name, pattern = paste0("^", file_pattern), full.names = TRUE)

        exists <- (length(matching_files) > 0)
    } else {
        logger("DEBUG", "Searching exact match")
        exists <- file.exists(path)
    }

    if (exists) {
        logger("DEBUG", "File exists at ", quotes(path), ".")
    } else {
        logger("DEBUG", "File does not exist at ", quotes(path), ".")
    }
    return(exists)
}

check_any_empty <- function(...) {
    args <- list(...)
    logger("DEBUG", "Checking if empty: ", quotes(list_to_str(args)))
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


    if (std_out) {
        logger("Running: ", quotes(plink_cmd))
        system(plink_cmd)
    } else {
        # Search for partial match
        logger("DEBUG", "Checking for existing match: ", quotes(plink_out_path), ".")
        if (file_exists(plink_out_path, TRUE) && !overwrite_plink_out) {
            logger("Matching file(s) already exists at: ", quotes(plink_out_path), ".")
            return(plink_out_path)
        }

        logger("Running: ", quotes(plink_cmd))
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

construct_plot_path <- function(basename) {
    #' Constructs a file path in the plots directory given a basename
    #' @param basename {string}: The name of the file without the plots directory prefixed.
    #' @return path {string}: The relative path to the file.

    file.path(plots_dir, basename)
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
    
    if (file_exists(out_path) && !overwrite_ext_plots) {
        logger("INFO", "Plot already saved at: ", quotes(out_path), ".")
        logger("INFO", "Skipping plot")
    } else {
        logger("DEBUG", "Plotting at: ", quotes(out_path), " ...")
        png(out_path, width, height)
        plot_callback(data, ...)
        dev.off()
    }

    return(out_path)
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

excess_missing_genotypes <- function(data_path, extension, out_cols, histogram, suffix) {
    if (length(suffix) == 0) {
        logger("ERROR", "Suffix for missing plot name cannot be empty.")
    }

    # Check for existing missing file
    missing_name <- paste0("missing_", suffix)
    logger("Checking for missing genotypes (", extension, ")...")
    missing_basename <- plink(data_path, pl_fgs$missing, missing_name)
    missing_out_path <- add_extension(missing_basename, extension)

    logger("DEBUG", "Reading ", extension, " table...")
    missing <- wrap_read_table(missing_out_path)
    log_df(missing, "missing")
    
    if (histogram) {
        hist_basename <- add_extension("fmiss", extension, exts$png)
        wrap_histogram(missing$"F_MISS", hist_basename)
        num_indvs_to_remove <- sum(missing$"F_MISS" > genotype_threshold)
        log_indvs_to_remove(num_indvs_to_remove)
    }
    
    out_name <- add_extension("remove.missing", extension, ".samples", exts$txt)
    missing_file_path <- construct_out_path(out_name)
    if (file_exists(missing_file_path)) {
        logger("Missing genotypes already saved at ", quotes(missing_file_path), ".")
    } else {
        missing_file_path <- remove_indices_by_threshold(missing, "F_MISS", genotype_threshold, 
                                                         out_cols, out_name)
    }

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
                             "keep", list("miss_pheno_neg_9" = "1")), 
                        named_flag)

exts <- create_object(list("phen", "imiss", "lmiss", "het", "assoc", "hwe", 
                           "frq", "txt", "png", "eigenvec", "eigenval",
                           "qassoc", "linear", "clumped", "rel", "id", "logistic"), 
                      ext)

# === Main ===

# ====== Analysis ======

quality_control <- function() {
    logger("Performing Quality Control")

    find_individual_missing_genotypes <- function(histogram) {
        #' Finds individuals with excess missing genotypes to remove from the dataset
        #' @param histogram {boolean}: if true, will plot a histogram of the frequency of missingness
        #' @return missing_file_path {string}: The file path to the file containing the individuals to remove.

        excess_missing_genotypes(NULL, exts$imiss, fam_ind_cols, histogram, "individual")
    }

    find_outlying_homozygosity <- function(plot) {
        #' Finds individuals with outlying homozygosity values to remove from the dataset
        #' @param plot {boolean}: if true, will plot a histogram of the frequency of hz freqs and a scatterplot
        #'                        of their distribution across the geneome.
        #' @return het_ind_file_path {string}: The file path to the file containing the individuals to remove.

        logger("Checking for outlying homozygosity values...")

        hz_name <- "hz"
        het_basename <- plink_orig_data(pl_fgs$het, hz_name)
        het_out_path <- add_extension(het_basename, exts$het)

        logger("DEBUG", "Reading het table...")
        het <- wrap_read_table(het_out_path)
        log_df(het, "het")

        if (plot) {
            logger("INFO", "Plotting het histogram")
            wrap_histogram(het$"F", "fhet_hist.png")

            logger("INFO", "Plotting het scatterplot")
            wrap_scatter(0.05, "red", abs(het$"F"), "fhet_scatter.png")
        }

        het_samples_name <- add_extension("remove.het.samples", exts$txt)
        het_ind_file_path <- construct_out_path(het_samples_name)
        if (file_exists(het_ind_file_path)) {
            logger("Het indices saved at: ", quotes(het_ind_file_path))
        } else {
            het_ind_file_path <- remove_indices_by_threshold(het, "F", het_threshold, 
                                                             fam_ind_cols, het_samples_name)
        }
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
        
        combined_file_out_path <- construct_out_path("remove.combined.samples.txt")
        if (file_exists(combined_file_out_path)) {
            logger("INFO", "Combined samples exists at: ", quotes(combined_file_out_path), ".")
            return(combined_file_out_path)
        }

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
        missing_file_path <- excess_missing_genotypes(data_subset_path, exts$lmiss, 2, histogram, "SNPs")
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
            wrap_histogram(hwe$P, "hwe_deviations.png")
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
            wrap_histogram(freq$MAF, "maf_distribution.png")
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

    compare_minor_allele_freqs <- function(freq, do_hist, do_plot) {
        logger("Comparing Minor Allele Frequencies...")
        snp_ref_path <- file.path(data_path, "reference_allele_frequencies.txt")
        ref <- wrap_read_table(snp_ref_path, header = FALSE)
        log_df(ref, "Reference allele requencies")
        
        ind <- match(freq$SNP, ref$V1)
        out <- cbind(freq, ref[ind,])
        log_df(out, "MAFs with Reference")
        
        # Take inverse of allele frequencies which are given with respect to allele 2.
        out_cpy <- out
        flip_allele <- 2
        logger("Inverting Allele Frequncies with respect to allele ", flip_allele, ".")
        out_cpy$MAF[out$A1 == flip_allele] <- 1 - out$MAF[out$A1 == flip_allele]

        # Must be calculated regardless of plots
        res <- out$MAF - out$V2
        log_df(res, "Residual minor allele frequencies")

        if (do_hist) {
            logger("Plotting histogram of Minor Allele Frequency (MAF)")
            keep <- c(which(abs(res) <= maf_threshold))
            wrap_histogram(res, "minor_allele.png")
        }

        if (do_plot) {
            logger("Plotting allele frequency comparison with reference...")
            wrap_plot(plot, out_cpy$MAF ~ out_cpy$V2, "min_allele_freq_comparison.png")

            logger("Removing alleles which deviate significantly from reference...")
            accept_snps <- which(abs(res) <= maf_threshold)
            out_accept <- out_cpy[accept_snps,]
            wrap_plot(plot, out_accept$MAF ~ out_accept$V2, "corrected_min_allele_freq_comparison.png")
        }

        remove_snps_indx <- which(abs(res) > maf_threshold)
        remove_snps <- freq[remove_snps_indx, "SNP"]
        remove_out_path <- construct_out_path(add_extension("maf.referenced.removed", exts$txt))
        wrap_write_table(remove_snps, remove_out_path)
        return(remove_out_path)
    }

    exclude_insig_maf <- function(qc_data_path, remove_out_path) {
        #' Creates a subset of the data which excludes insignificant alleles by their MAF
        #' @param qc_data_path {string}: Path to dataset to exclude from
        #' @param remove_out_path {string}: Path to file detailing which alleles to remove
        #' @return qc_data_maf_path {string}: Path to new subset of data

        logger("Excluding insignificant alleles by MAF...")

        out_name <-"test_qc_insif_maf"
        plink_args <- paste(pl_fgs$mb, pl_fgs$exclude, remove_out_path)
        plink(qc_data_path, plink_args, out_name)
    }
    
    lmiss <- missing_snps(TRUE)
    hwe <- hw_eq(TRUE)
    freq <- min_allele_freq(TRUE)
    
    remove_snps_path <- remove_snps(lmiss, hwe, freq)
    qc_data_path <- exclude_snps(remove_snps_path)

    remove_out_path <- compare_minor_allele_freqs(freq, TRUE, TRUE)
    qc_data_maf_path <- exclude_insig_maf(qc_data_path, remove_out_path)
    
    return(qc_data_maf_path)
}

gwas <- function(qc_data_path) {
    logger("INFO", "Performing GWAS...")

    # Some scoped variables
    alt_mpheno <- 1
    num_pc <- 10
    
    get_mpheno_args <- function(suffix) {
        args <- paste(pl_fgs$mpheno, alt_mpheno)
        if (grepl("binary", suffix)) {
            logger("DEBUG", quotes("binary"), " found in suffix ", quotes(suffix))
            logger("INFO", "Adding missing-phenotype flag")
            args <- paste(args, pl_fgs$miss_pheno_neg_9)
        }
        
        return(args)
    }

    get_pheno_path <- function(pheno_suffix) {
        logger("DEBUG", "Retrieving pheno path for suffix ", quotes(pheno_suffix), ".")
        phenotype_file_prefix <- space_to_underscore(phenotype)
        pheno_file_name <- paste0(phenotype_file_prefix, pheno_suffix, exts$phen)
        file.path(phenotypes, pheno_file_name)
    }

    gwas_pheno <- function(pheno_path, pheno_suffix, mpheno_args) {
        #' Performs association analysis based on the phenotype
        #' defined in the specified file.
        #' @param pheno_suffix {string}: Suffix of phenotype.
        #' @param pc {boolean}: Flag for if we are performing assocation analysis 
        #'                      with principal components.
        #' @return {string}: Path to phenotype association analysis output.

        logger("Performing Pheno Association Analysis on GWAS...")
        
        plink_args <- paste(pl_fgs$assoc, pl_fgs$pheno, pheno_path, mpheno_args)
        out_name <- paste0("gwas_pheno", pheno_suffix)
        plink(qc_data_path, plink_args, out_name)
    }

    get_pheno_analysis_full_path <- function(pheno_basename, phenotype_suffix, pc) {
        logger("DEBUG", "Calculating path for pheno analysis full path")
        
        if (phenotype_suffix == "") {
            if (pc) {
                logger("TRACE", "Principal components enabled.")
                ext_to_add <- paste0(exts$assoc, exts$linear)
            } else {
                logger("TRACE", "Principal components disabled.")
                ext_to_add <- exts$qassoc
            }
        } else {
            if (pc) {
                ext_to_add <- paste0(exts$assoc, exts$logistic)
            } else {
                ext_to_add <- exts$assoc
            }
        }
    
        pheno_path <- add_extension(pheno_basename, ext_to_add)
        logger("DEBUG", "Full pheno path: ", quotes(pheno_path), ".")
        return(pheno_path)
    }
 
    gwas_plots <- function(pheno_analysis_path, plot_suffix = "", pc) {
        name_plot <- function(plot_type) {
            logger("DEBUG", "Naming plot type: ", quotes(plot_type), ".")
            file_name <- paste0("gwas_", plot_type, plot_suffix)
            file_name <- add_extension(file_name, exts$png)
            logger("DEBUG", "GWAS Plot Name: ", quotes(file_name), ".")
            return(file_name)
        }

        logger("DEBUG", "Reading pheno analysis ", quotes(plot_suffix), " ",
                        "Principal components = ", pc, "...")
        d <- wrap_read_table(pheno_analysis_path)
        log_df(d, "D before na omit")
        
        logger("INFO", "Removing n/a p-vals...")
        d <- d[!is.na(d$P), ]

        log_df(d, "D after na omit")

        logger("DEBUG", "Plotting ...")
        man_plot_name <- name_plot("manhattan")
        logger("INFO", "Generating Manhattan Plot ", quotes(man_plot_name), " ...")
        wrap_plot(manhattan, d, man_plot_name)

        qq_plot_name <- name_plot("qq")
        logger("INFO", "Generating QQ Plot ", quotes(qq_plot_name), " ...")
        wrap_plot(qq, d$P, qq_plot_name)

        return(d)
    }

    compute_lambda <- function(d, suffix, pc) {
        logger("INFO", "Computing Genomic Inflation Factor (", suffix, ") with PC: ", 
                       pc ? "enabled" : "disabled", "...")

        lambda_gc <- genomic_inflation_factor(d, 1)
        logger("INFO", "Lambda_{GC} (", suffix, ") = ", lambda_gc)
        return(lambda_gc)
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

    add_pc_covariates <- function(pheno_pc_path, suffix, pc_eigvec_file, mpheno_args) {
        # Check for existing covariates
        out_name <- paste0("gwas_pheno", suffix, "_pc")
        cov_out_path <- construct_plink_out_path(out_name)
        
        if (file_exists(pheno_pc_path)) {
            logger("INFO", "Covariates already exist. Skipping.")
            return (pheno_pc_path)
        }
        
        # No covariates. Compute:
        logger("Adding principal components to covariates...")
        logger("DEBUG", "Loading eigenvec file: ", quotes(pc_eigvec_file), ".")
        check_ext(pc_eigvec_file, exts$eigenvec)
        plink_args <- paste(suffix == "" ? pl_fgs$linear : pl_fgs$logistic, pl_fgs$covar, 
                            pc_eigvec_file, pl_fgs$pheno, pheno_pc_path, mpheno_args)
        plink(qc_data_path, plink_args, out_name)
        return(pheno_pc_path)
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

    read_clumps <- function(clump_basename, suffix) {
        clump_path <- add_extension(clump_basename, exts$clumped)
        clump <- wrap_read_table(clump_path)

        # Exclude column PC2
        logger("Printing Clump Output")
        clump_out <- clump[, 1:11]
        log_df(clump_out, paste0("Clump output", suffix))
        
        # Write to file
        out_path <- construct_out_path(add_extension(paste0("clumps", suffix), exts$txt))
        wrap_write_table(clump_out, out_path)
    }
 
    # Main
    phenotype_suffixes <- list("", "_binary1", "_binary2")
    
    for (pca in list(FALSE, TRUE)) {
        on <- pca == TRUE ? "enabled" : "disabled" 
        logger("INFO", "Performing GWAS with Principal Components ", on, ".")

        if (pca) {
            # Compute principal components once
            pc_eigvec_file <- compute_principal_comps(num_pc)
        }
        
        # Perform analysis for each of the phenotypes
        for (suffix in phenotype_suffixes) {
            logger("INFO", "Inspecting phenotype ", quotes(suffix), ".")
            pheno_path <- get_pheno_path(suffix)
            mpheno_args <- get_mpheno_args(suffix)
           
            if (pca) {
                pheno_basename <- add_pc_covariates(pheno_path, suffix, pc_eigvec_file, mpheno_args)
            } else {
                pheno_basename <- gwas_pheno(pheno_path, suffix, mpheno_args)
            }

            pheno_full_path <- get_pheno_analysis_full_path(pheno_basename, suffix, pca)

            if (pca) {
                clump_path <- clumping(pheno_pc_path)
                read_clumps(clump_path, suffix)
            }

            d <- gwas_plots(pheno_full_path, suffix, pca)
            compute_lambda(d, suffix, pca)
        }
    }
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

