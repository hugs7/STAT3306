# Course: STAT3306: Statistical Analysis of Genetic Data
# Taken: Semester 2, 2024
# Project Part 2: Genetic architecture analyses
# Author: Hugo Burton
# Student #: s4698512
# Due Date: 2pm 25 October 2024

install_if_missing <- function(packages) {
    missing_pkgs <- packages[!(packages %in% installed.packages()[, "Package"])]

    if (length(missing_pkgs) > 0) {
        install.packages(missing_pkgs)
    }
    invisible(NULL)
}

# === Packages ===

required_packages <- c("crayon", "xtable")
install_if_missing(required_packages)
invisible(lapply(required_packages, require, character.only = TRUE))

# === Modules ===

source("./grm_bin.R")

# === Logging Config ===

default_log_level <- "INFO"
allowed_log_levels <- c("ERROR", "WARN", "INFO", "DEBUG", "TRACE")
app_log_level <- "DEBUG"
app_log_index <- which(allowed_log_levels == app_log_level)
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

# Data Paths
course_shared_data_path <- file.path("/data/STAT3306")
project_data <- file.path(course_shared_data_path, "Project")
data_path <- file.path(project_data, "Data_QC")
phenotypes_path <- file.path(project_data, "Phenotypes_QC")

# GCTA
gcta_datafile_basename <- "testFiltered"
gcta_default_tnum <- 2

# Out Paths
plots_out_dir <- file.path("./plots")
gcta_out_dir <- file.path("./gcta_out")
out_dir <- file.path("./out")

# Plotting
plot_w <- 500
plot_h <- 300

# Thresholds
grm_relatedness_threshold <- 0.05

# Overrides
overwrite_ext_plots <- FALSE
overwrite_gcta_out <- FALSE

# === Functions ===

cat0 <- function(...) {
    #' Concatenates strings together and ouputs them to the console
    #' @param ... {string}: String args to concatenate
    #' @return {NULL}

    msg <- paste0(...)
    cat(msg, "\n")
}

pad <- function(width, ...) {
    #' Pads string args to a given width
    #' @param width {integer}: width in characters to pad the string to
    #' @param ... {string}: String args to pad
    #' @return padded {string}: The padded string.

    # Do not trace log here as logger relies on pad.
    clipped <- substr(paste0(...), 1, width - 1)
    padded <- sprintf(paste0("%-", width, "s"), clipped)
    return(padded)
}

quotes <- function(...) {
    #' Surrounds a string with double quotes (")
    #' @param ... {string}: String args to surround with quotes
    #' @return {string}: String surrounded with quotes
 
    str <- args_to_string(...)
    paste0("'", str, "'")
}

brackets <- function(...) {
    #' Surrounds a string with round brackets.
    #' @param ... {string}: String args to surround with brackets.
    #' @return {string}: String surrounded with round brackets.

    str <- args_to_string(...)
    paste0("(", str, ")")
}

latex_math <- function(...) {
    #' Surrounds a string with LaTeX math mode $...$.
    #' @param ... {string}: String args to surround with LaTeX
    #'                      math mode.
    #' @return {string}: String surrounded by LaTex math mode.
    
    str <- args_to_string(...)
    paste0("$", str, "$")
}

to_str <- function(x, collapse = ", ") {
    #' Converts a list or vector to a string separated by commas.
    #' @param x {list|vector}: List or vector to convert to string.
    #' @param collapse {string}: Separator to split items in list
    #'                           or vector. Defaults to ', '.
    #' @return {string}: String representation of list.
    
    if (is.list(x) || is.vector(x)) {
        x <- unlist(x)
    }

    paste(x, collapse = collapse)
}

args_to_string <- function(...) {
    #' Converts any argument which is a list to a string
    #' @param ... {string}: Arguments to convert.
    #' @return {list{string}}: String arguments.

    args <- list(...)
    string_args <- sapply(args, function(arg) {
        if (is.list(arg) || is.vector(arg)) {
            to_str(arg)
        } else {
            as.character(arg)
        }
    })

    str <- paste0(string_args)
    return(str)
}

named_flag <- function(fg_name) {
    #' Prefixes a flag name with --
    #' @param fg_name {string}: The flag name to prefix
    #' @return {string}: The flag with '--' prefixed

    logger("TRACE", "Naming flag: ", quotes(fg_name), ".")
    paste0("--", fg_name)
}

ext <- function(ext_name) {
    #' Prefixes an extension name with a period (.)
    #' @param ext_name {string}: The extension to prefix
    #' @return {string} The extension with a period prefixed

    logger("TRACE", "Initialising extension: ", quotes(ext_name), ".")
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
        logger("TRACE", "Assignment detected.")
        x <- eval(xs[[3]], envir = env)
    }

    # Parse and eval ternary outcome in parent env
    outcomes <- sapply(strsplit(deparse(substitute(y)), ":"), function(e) parse(text = e))
    logger("TRACE", "Outcomes: ", quotes(to_str(outcomes)), ".")
    r <- eval(outcomes[[2 - as.logical(x)]], envir = env)

    # Handle assignment
    if (xs[[1]] == as.name("<-")) {
        logger("TRACE", "Assignment detected.")
        xs[[3]] <- r
        eval.parent(as.call(xs))
    } else {
        logger("TRACE", "No assigment in ternary operation.")
        r
    }
}

create_object <- function(items, transform_fn) {
    #' Creates an object from a list of items where items have a pre-transform
    #' @param items {list}: List of items to insert into the object
    #' @param transform_fn {function}: Pretransform for each item
    #' @return result {list}: List of items which have a key lookup (like an object)

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

    return(result)
}

get_func_name <- function(call) {
    #' Given a call from the call stack, returns
    #' the name of the function called.
    #' @param call {object}: The call to get the function name from.
    #' @return func_name {string}: The function name
    
    if (is.call(call) && length(call) > 0) {
        func_name <- as.character(call[[1]])
    } else {
        func_name <- "Unknown"
    }

    return(func_name)
}

get_calling_function <- function(ignore_names) {
    #' Gets the name of the function calling a function from the call stack
    #' @param ignore_names {list[string]}: List of functions to skip in the call stack
    #' @return func_name {string}: Name of function calling.

    call_stack <- sys.calls()

    # Skip this function
    for (i in (length(call_stack) - 1):1) {
        current_call <- call_stack[[i]]
        func_name <- get_func_name(current_call)

        if (!(func_name %in% ignore_names)) {
            return(func_name)
        }
    }

    return("Global Environment")
}

log_stack <- function(log_level = default_log_level) {
    #' Logs the call stack to the console.
    #' @return {NULL}
   
    call_stack <- sys.calls()

    logger(log_level, "Stack Trace")
    indent <- strrep(" ", 4)

    # Skip this function
    for (i in (length(call_stack) - 1):1) {
        current_call <- call_stack[[i]]
        func_name <- get_func_name(current_call)

        code_line <- paste0(deparse(current_call), collapse = " ")
        logger(log_level, indent, pad(30, func_name), " | Code: ", code_line)
    }
}

logger <- function(log_level = default_log_level, ...) {
    #' A mini logger function with ansi coloured log levels. Outputs timestamp,
    #' calling function, log level and message to the console. Logs the stack
    #' trace if log level is ERROR. Will not log if log_level is higher than
    #' the application log level.
    #' @param log_level {string}: Optional. The level to log at. If not provided,
    #'                            will default to default_log_level.
    #' @return {NULL}

    if (!is.character(log_level) || !(log_level %in% allowed_log_levels)) {
        logger(default_log_level, log_level, ...)
        return(invisible(NULL))
    }

    log_index <- which(allowed_log_levels == log_level)
    if (log_index > app_log_index) {
        return(NULL)
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

    if (log_level == "ERROR" && parent_call != "log_stack") {
        log_stack("ERROR")
    }

    invisible(NULL)
}

log_df <- function(df, name, log_level = default_log_level) {
    #' Logs the dimension and head of a dataframe at the specified log level,
    #' @param df {data.frame}: The data frame to preview,
    #' @param name {string}: The name of the data frame.
    #' @param log_level {string}: The level to log at. Defaults to INFO.
    #' @return {NULL}

    if (is.null(df)) {
        logger("ERROR", "Unexpected! Dataframe is NULL.")
        return
    }

    logger(log_level, "Dim ", name, ": ", wrap_dim(df))
    print(head(df))
}

space_to_underscore <- function(str) {
    #' Replaces all instances of spaces in a string with underscores.
    #' @param str {string}: String to replace from.
    #' @return {string}: string with ' 's as '_'s.

    gsub(" ", "_", str)
}

title_case <- function(str) {
    #' Converts a string to titlecase.
    #' @param str {string}: The string to convert to title case.
    #' @return result {string}: The string in titlecase form.

    if (length(str) == 0) {
        logger("WARN", "String is empty")
        return(str)
    }

    logger("TRACE", "Converting ", quotes(str), " to title case...")
    words <- strsplit(str, " ")[[1]]
    titlecased_words <- paste0(toupper(substring(words, 1, 1)), tolower(substring(words, 2)))
    result <- to_str(paste(titlecased_words), " ")
    
    logger("TRACE", "Title case version: ", quotes(result), ".")
    return(result)
}

regex_escape <- function(string) {
    #' Escapes characters in a string so they can be
    #' used in a RegEx string.
    #' @param string {string|list{string}}: The string or list of strings to be escaped.
    #' @return {string|list{string}}: The same string but RegEx escaped.
    
    escape_string <- function(s) {
        #' Escapes a single string.
        #' @param s {string}: Single string to escape.
        #' @return {string}: RegEx escaped string.

        gsub("([\\W])", "\\\\\\1", s, perl = TRUE)
    }

    multiple <- is.list(string)
    string_word <- paste0("string", multiple ? "s" : "")
    logger("TRACE", "Escaping ", string_word, ": ", quotes(string), ".")

    if (multiple) {
        string <- lapply(string, escape_string)
    } else {
        string <- escape_string(string)
    }

    logger("TRACE", "Escaped ", string_word, ": ", quotes(string), ".")
    return(string)
}

shell_call <- function(...) {
    #' Makes a system call with standard out disabled.
    #' @param ... {string}: String args containing system command.
    #' @return {NULL}

    system(..., ignore.stdout = TRUE)
    invisible(NULL)
}

mkdir_if_not_exist <- function(path) {
    #' Makes a directory at the specified path if it doesn't already exist.
    #' @param path {string}: The path to create a directory at.
    #' @return {NULL}

    if (!dir.exists(path)) {
        logger("Path ", quotes(path), " does not exist. Creating...")
        dir.create(path, recursive = TRUE)
        logger("DEBUG", "Path ", quotes(path), " created.")
    } else {
        logger("DEBUG", "Path ", quotes(path), " already exists.")
    }
}

list_files <- function(dir_name, pattern = NULL, full.names = TRUE, ...) {
    #' Wrapper for list.files to set full.names default to TRUE.
    #' @param dir_name {string}: The directory to search in.
    #' @param pattern {string}: Optionally match filenames with a RegEx pattern.
    #' @param full.names {boolean}: If true, the directory path is prepended
    #'                              to the file names
    #' @param ... {any}: Any extra arguemnts for list.files()
    #' @return {charvec}: Character vector containing filenames/filepaths
    #'                    matching criteria.
    
    logger("DEBUG", "Listing files in dir ", quotes(dir_name), ".")
    if (!is.null(pattern)) {
        logger("DEBUG", "Matching pattern ", quotes(pattern), ".")
    }

    matching_files <- list.files(path = dir_name, full.names = full.names, ...)

    if (!is.null(pattern)) {
        logger("TRACE", "Filtering files...")
        basenames <- basename(matching_files)
        filtered_files <- matching_files[grepl(pattern, basenames, perl = TRUE)]
        matching_files <- filtered_files
    }

    logger("DEBUG", "There are ", length(matching_files), " files matching")
    for (file in matching_files) {
        logger("TRACE", "   ", basename(file))
    }

    return(matching_files)
}

file_exists <- function(path, match_pattern = FALSE) {
    #' Checks to see if a file exists at the specified path
    #' @param path {string}: The specified path to check. Can be relative.
    #' @param match_pattern {bool}: If TRUE, the basename of the path will be considered as a
    #'                        regular expression. Defaults to FALSE where the path is
    #'                        considered exactly.
    #' @return exists {boolean}: True if file exists, false otherwise.

    logger("DEBUG", "Checking if file exists at path ", quotes(path), ".")
    if (length(path) == 0) {
        logger("WARN", "Path is empty")
        return(FALSE)
    }

    if (match_pattern) {
        logger("DEBUG", "Searching pattern match...")
        dir_name <- dirname(path)
        file_pattern <- basename(path)
        matching_files <- list_files(dir_name, file_pattern)
        exists <- (length(matching_files) > 0)
    } else {
        logger("DEBUG", "Searching exact match")
        exists <- file.exists(path)
    }

    logger("DEBUG", "File ", exists ? "exists" : "does not exist", " at ", quotes(path), ".")
    return(exists)
}

latex_table <- function(data, out_name, table_align, caption = NULL, col.names = NULL,
                        digits = 2, line_spacing_factor = 1, hide_row_names = FALSE,
                        size = "normalsize") {
    #' Generates a LaTeX table given a data.frame and saves to a file.
    #' @param data {data.frame}: The data.frame to output as a LaTeX table.
    #' @param out_name {string}: Filename to save the output as.
    #' @param table_align {string}: Latex coding for aligning columns.
    #' @param caption {string|NULL}: Optional caption for the table.
    #' @param col.names {vec|NULL}: Optional column names to provide to the table.
    #' @param digits {integer|vector(integer)}: Number of decimal places to display
    #'                                          numbers as. If digits < 0, values
    #'                                          will be converted to scientific
    #'                                          notation with the negative of the
    #'                                          number of digits prior to passint
    #'                                          to xtable.
    #' @param line_spacing_factor {integer}: Line spacing factor for LaTeX table.
    #' @param hide_row_names {bool}: Whether to include row names from the data.frame
    #'                               as the first column in the table. Defaults to FALSE.
    #' @param size {string}: LaTeX size encoding without \. Defaults to normalsize.
    #' @return path {string}: Path to saved LaTeX table.
   
    log_df(data, paste("Latex table", out_name))
    
    xtable_table_align <- gsub(":", "|", table_align)
    align_num_cols <- nchar(gsub("\\|", "", xtable_table_align))
    logger("DEBUG", "Table align num cols: ", align_num_cols)
    df_ncols <- ncol(data)
   
    if (align_num_cols != df_ncols) {
        logger("WARN", "Number of columns in table alignment ", brackets(align_num_cols),
               " does not match number of cols in df ", brackets(df_ncols), ".")
    }
   
    # Add leading dummy align character for index column.
    xtable_table_align <- paste0("l", xtable_table_align)
    logger("DEBUG", "DF num cols: ", df_ncols)
    
    if (hide_row_names) {
        table_align <- sub("A-Z|a-z", "", table_align)
        logger("TRACE", "Table align with first col hidden: ", quotes(table_align), ".")
    }

    if (df_ncols == 1) {
        logger("WARN", "Something seems off. Dataframe has only 1 column.")
    }

    if (length(digits) > 1 || df_ncols == 1) {
        # Prefix with dummy digit for index column
        digits <- c(0, digits)
    }
    
    logger("DEBUG", "Digits: ", quotes(digits), ", length: ", length(digits), ", typeof: ",
           brackets(typeof(digits)), ".")
    
    table <- xtable(data, align = xtable_table_align, caption = caption, digits = digits)

    if (!is.null(col.names)) {
        colnames(table) <- col.names
    }

    latex <- print.xtable(table, print.results = FALSE, table.placement = "H",
                          comment = FALSE, include.rownames = !hide_row_names,
                          sanitize.colnames.function = function(x) {paste0("\\textbf{", x, "}")})

    latex <- gsub("\\begin{tabular}",
                  paste0(size == "normalsize" ? "" : paste0("\\", size, "\n"),
                         "\\renewcommand{\\arraystretch}{",
                         line_spacing_factor, "}\n\\begin{tabular}"
                        ),
                  latex,
                  fixed=TRUE)
    
    latex <- sub(regex_escape(xtable_table_align), table_align, latex)
    
    logger("DEBUG", "Writing latex table...")
    latex_path <- wrap_write(latex, out_name)
    return(latex_path)
}

get_ext_pattern <- function(extensions) {
    #' Generates a RegEx pattern to match any extension provided
    #' from a list of extensions
    #' @param extensions {list{str}}: List of extensions to match
    #'                                to.
    #' @return {string}: RegEx pattern matching provided extensions.

    ext_pattern <- paste0("(", to_str(regex_escape(extensions),
                          collapse = "|"), ")$")
    logger("TRACE", "Extension pattern: ", quotes(ext_pattern), ".")
    return(ext_pattern)
}

ends_with_extension <- function(ext_pattern, path) {
    #' Checks if a path ends with an extension given in a pattern.
    #' @param ext_pattern {string}: RegEx pattern to match path to.
    #' @param path {string}: Path to check extension in.
    #' @return {bool}: TRUE if the path contains extension in
    #'                 pattern, FALSE otherwise.

    ends_with_ext <- grepl(ext_pattern, path)
    logger("TRACE", "Ends with extension: ", quotes(ext_pattern), " ",
           ends_with_ext, ".")
    return(ends_with_ext)
}

match_not_log <- function(path) {
    #' Given an extensionless file path, constructs
    #' a Regex pattern to exclude matching .log files.
    #' @param path {string}: The extensionless pattern to match.
    #' @return {string}: The Regex pattern which includes a negative
    #'                   lookahead to log files.
   
    # Check the path does not already contain an extension.
    ext_pattern <- get_ext_pattern(exts)
    if (ends_with_extension(ext_pattern, path)) {
        logger("ERROR", "Unexpected! Path ends with extension ",
               quotes(ext_pattern), ".")
        return(path)
    }

    pattern <- paste0(path, "(?!\\.log$).*$")
    logger("DEBUG", "Match not log pattern ", quotes(pattern), ".")
    return(pattern)
}

check_any_empty <- function(...) {
    #' Checks if any string args are empty ("")
    #' @param ... {string}: String args to check
    #' @return {boolean}: True if ANY of the string args are empty, false otherwise.

    args <- list(...)
    logger("DEBUG", "Checking if empty: ", quotes(to_str(args)))
    any(sapply(args, function(x) x == ""))
}

delete_file <- function(path) {
    #' Deletes a file at the specified path if it exists.
    #' @param path {string}: The path of the file to delete.
    #' @return success {boolean}: True if the file was deleted, false otherwise.

    exists <- file_exists(path)

    if (exists) {
        logger("WARN", "Deleting file at path ", quotes(path), ".")
        file.remove(path)
    }

    # Assume if the file exists and we didn't crash, we were able to delete.
    return(exists)
}

wrap_write <- function(content, basename) {
    #' Writes contents to a file. Output is always placed in out directory.
    #' Allows writing as .txt or .tex.
    #' @param content {string}: The content to write.
    #' @param basename {string}: Name of the file to write to.
    #' @return path {string}: The full save path where the file was saved.

    ext_pattern <- get_ext_pattern(list(exts$txt, exts$tex))
    logger("DEBUG", "Initial basename: ", quotes(basename), ".")
    if (!ends_with_extension(ext_pattern, basename)) {
        basename <- check_txt_ext(basename)
    }
    
    path <- construct_out_path(basename)

    if (file_exists(path)) {
        logger("WARN", "Overwriting file at path ", quotes(path), ".")
        delete_file(path)
    }

    if (!is.character(content)) {
        logger("ERROR", "Cannot write content. Type is not character.")
        return(path)
    }
    
    logger("DEBUG", "Writing data to ", quotes(path), " ...")

    write(content, path)
    logger("DEBUG", "Writing complete. Cleaning up...")
    rm(content)
    gc()
    logger("DEBUG", "Cleanup complete.")
    return(path)
}

wrap_read_table <- function(path, header = TRUE, ...) {
    #' Wrapper for reading a table from a file.
    #' @param path {string}: The path of the file to read the table from.
    #' @param header {boolean}. Defaults to true. If true, expects header to exist in file.
    #' @return {data.frame}: Dataframe containing the table data. NULL if file does not exist.

    if (!file_exists(path)) {
        logger("ERROR", "Could not find file to read at ", quotes(path), ".")
        return(NULL)
    }

    logger("DEBUG", "Reading table at ", quotes(path), ".")
    read.table(path, header = header, ...)
}

wrap_write_table <- function(data, basename, row.names = FALSE, col.names = TRUE,
                             sep = "\t", quote = FALSE, ...) {
    #' Wrapper for writing a table to a file. Will overwrite file if it exists
    #' at the same path.
    #' @param data {data.frame}: The data to write.
    #' @param basename {string}: The basename excluding the out dir to  write
    #'                           the file at. Ideally should include extension but if
    #'                           it doesn't this function will add .txt by default.
    #' @param row.names {boolean}: Whether to include row names. Disabled by default.
    #' @param col.names {boolean}: Whether to include col names. Enabled by default but
    #'                             if the data contains default columns, will not be
    #'                             written.
    #' @param sep {string}: The separator to delimit between columns in the table.
    #' @param quote {boolean}: Whether to include quotes for strings in the table
    #'                         data. Disabled by default.
    #' @return path {string}: The full save path where the table was saved.

    ext_pattern <- get_ext_pattern(list(exts$txt, exts$phen, exts$cov))
    logger("DEBUG", "Initial basename: ", quotes(basename), ".")
    if (!ends_with_extension(ext_pattern, basename)) {
        basename <- check_txt_ext(basename)
    }

    path <- construct_out_path(basename)

    if (file_exists(path)) {
        logger("WARN", "Overwriting file at path ", quotes(path), ".")
        delete_file(path)
    }

    if (!is.data.frame(data)) {
        logger("ERROR", "Cannot write data. Type is not dataframe")
        return(path)
    }

    log_df(data, "Writing data...")
    num_cols <- ncol(data)
    if (num_cols == 0) {
        logger("ERROR", "No columns found in data. Failed to write table.")
        return(path)
    }

    default_colnames <- paste0("V", seq_len(num_cols))
    data_colnames <- colnames(data)
    logger("DEBUG", "Data colnames: ", to_str(data_colnames), ".")
    if (col.names == TRUE && all(data_colnames == default_colnames)) {
        logger("Data contains default column names (V1, V2, ...).",
               "Supressing column names.")
        col.names <- FALSE
    }

    logger("DEBUG", "Writing table at ", quotes(path), "...")
    write.table(data, path, row.names = row.names, col.names = col.names,
                sep = sep, quote = quote, ...)
    logger("DEBUG", "Writing table complete. Cleaning up...")
    rm(data)
    gc()
    logger("DEBUG", "Cleanup complete.")
    return(path)
}

gcta_qc_data <- function(gcta_args, out_name = NULL) {
    #' Makes a GCTA call with the QC'd dataset and specific arguments.
    #' Serves as a wrapper to gcta function.
    #' @param gcta_args {string}: Arguments provided to GCTA
    #' @param out_name {string}: Basename for GCTA to output to. Should exclude
    #'                           gcta out directory. If NULL, output is piped to
    #'                           the console.
    #' @return gcta_out_path {string}: Relative path from script to gcta output.
    #'                                 Notable, this path does not contain
    #'                                 extension added by GCTA as this differs
    #'                                 depending upon arguments provided to GCTA.

    data_files_pattern <- construct_data_path(gcta_datafile_basename)
    logger("DEBUG", "GCTA data files pattern: ", quotes(data_files_pattern), ".")
    
    gcta_args <- paste(gcta_fgs$bfile, data_files_pattern, gcta_args
    
    gcta(gcta_args, out_name)
}

gcta <- function(gcta_args, out_name = NULL) {
    #' Runs a GCTA call with the given dataset and specified arguments.
    #' @param gcta_args {string}: Arguments provided to GCTA
    #' @param out_name {string}: Basename for GCTA to output to. Should exclude
    #'                           gcta out directory. If NULL, output is piped to
    #'                           the console.
    #' @return gcta_out_path {string}: Relative path from script to gcta output.
    #'                                 Notable, this path does not contain
    #'                                 extension added by GCTA as this differs
    #'                                 depending upon arguments provided to GCTA.

    gcta_base_cmd <- paste("gcta", gcta_args)
    logger("DEBUG", "GCTA base command: ", quotes(gcta_base_cmd), ".")

    
    if (is.null(out_name)) {
        # Output to console.
        gcta_cmd <- gcta_base_cmd
        std_out <- TRUE
    } else {
        # Outputs to file.
        gcta_out_path <- construct_gcta_out_path(out_name)
        logger("TRACE", "GCTA out path: ", quotes(gcta_out_path), ".")
        gcta_cmd <- paste(gcta_base_cmd, gcta_fgs$out, gcta_out_path)
        std_out <- FALSE
    }

    # Check double space - indicates missing param
    if (grepl(" ", gcta_cmd)) {
        logger("ERROR", "Possible missing argument in gcta command.")
        logger("ERROR", "You were trying to run ", quotes(gcta_cmd))
    }

    if (stdout) {
        logger("Running: ", quotes(gcta_cmd))
        system(gcta_cmd)
    } else {
        # Search for partial matching file
        logger("DEBUG", "Checking for existing match: ", quotes(gcta_out_path), ".")
        gcta_pattern <- match_not_log(gcta_out_path)
        if (file_exists(gcta_pattern< TRUE) && !overwrite_gcta_out) {
            logger("Matching file(s) already exist at: ", quotes(gcta_out_path), ".")
            return(gcta_out_path)
        }

        logger("Running: ", quotes(gcta_cmd))
        shell_call(gcta_cmd)
        logger("GCTA results directed to ", quotes(gcta_out_path), ".")
    }

    gc()
    return(gcta_out_path)
}

add_extension <- function(basename, ...) {
    #' Adds extension(s) to a given basename. The basename can by a file
    #' path or just a file basename.
    #' @param basename {string}: File path or basename to add extension to.
    #' @param ... {string}: String args containing additional extension(s)
    #'                      to append.
    #' @return path {string}: Basename or path with extension(s) appended.

    args <- list(...)
    logger("TRACE", "Checking extensions ", ..., ".")
    if (check_any_empty(...)) {
        logger("ERROR", "Some extensions are empty in ", ..., ".")
    }

    logger("TRACE", "Adding extensions: ", ..., " to basename ", quotes(basename), ".")
    path <- paste0(basename, ...)
    logger("TRACE", "Constructed out path: ", path, ".")
    return(path)
}

construct_phenotypes_path <- function(...) {
    #' Constructs a file base to the phenotypes directory given a basename.
    #' @param ... {string}: The name of the file (potentially split) without
    #'                      the parent directory prefixed.
    #' @return path {string}: The relative path to the file.

    basename <- paste0(...)
    file.path(phenotypes_path, basename)
}

construct_data_path <- function(...) {
    #' Constructs a file base to the data directory given a basename.
    #' @param ... {string}: The name of the file (potentially split) without
    #'                      the parent directory prefixed.
    #' @return path {string}: The relative path to the file.

    basename <- paste0(...)
    file.path(data_path, basename)
}


construct_gcta_out_path <- function(...) {
    #' Constructs a file path in the gcta out directory given a basename.
    #' @param ... {string}: The name of the file (potentially split) without
    #'                      the gcta out directory prefixed.
    #' @return path {string}: The relative path to the file.

    basename <- paste0(...)
    file.path(gcta_out_dir, basename)
}


construct_out_path <- function(basename) {
    #' Constructs a file path in the out directory given a basename.
    #' @param basename {string}: The name of the file without the out directory prefixed.
    #' @return path {string}: The relative path to the file.

    file.path(out_dir, basename)
}

construct_plot_path <- function(basename) {
    #' Constructs a file path in the plots directory given a basename
    #' @param basename {string}: The name of the file without the plots directory prefixed.
    #' @return path {string}: The relative path to the file.

    file.path(plots_out_dir, basename)
}

check_ext <- function(out_name, expected_ext, add_if_missing = TRUE) {
    #' Checks the extension of a file name is present and optionally
    #' adds it if it's not present.
    #' @param out_name {string}: File path or name to check.
    #' @param expected_ext {string}: The expected extension (with .).
    #' @param add_if_missing {boolean}: Adds the expected extension if
    #'                                  not present. Defaults to true.
    #' @return out_name {string|NULL}: The (possibly revised) file path
    #'                                 or name. NULL if expected extension
    #'                                 is blank.

    if (length(expected_ext) == 0) {
        logger("ERROR", "Expected extension is blank")
        return(NULL)
    }

    logger("DEBUG", "Checking ", quotes(out_name), " for extension ", quotes(expected_ext), ".")
    if (!endsWith(out_name, expected_ext)) {
        logger("WARN", "Out name ", quotes(out_name), " does not end with ",
               quotes(expected_ext), ".")
        if (add_if_missing) {
            logger("DEBUG", "Adding ", quotes(expected_ext), " to file name")   
            out_name <- paste0(out_name, expected_ext)
        }
    } else {
        logger("DEBUG", "Out name ", quotes(out_name), " contains extension ",
               quotes(expected_ext), ".")
    }

    return(out_name)
}

check_png_ext <- function(out_name, add_if_missing = TRUE) {
    #' Wrapper to check_ext() for png extension.
    #' @param out_name {string}: File path or name to check.
    #' @param add_if_missing {boolean}: Adds the expected extension if
    #'                                  not present. Defaults to true.
    #' @return out_name {string}: The (possibly revised) file path
    #'                            or name.

    check_ext(out_name, exts$png, add_if_missing)
}

check_txt_ext <- function(out_name, add_if_missing = TRUE) {
    #' Wrapper to check_ext() for txt extension.
    #' @param out_name {string}: File path or name to check.
    #' @param add_if_missing {boolean}: Adds the expected extension if
    #'                                  not present. Defaults to true.
    #' @return out_name {string}: The (possibly revised) file path
    #'                            or name.
    
    check_ext(out_name, exts$txt, add_if_missing)
}

wrap_plot <- function(plot_callback, data, out_name, ..., width = plot_w, height = plot_h) {
    #' Plots a data.frame using the specified callback and saves to a file
    #' in the plots directory.
    #' @param plot_callback {function}: Function used to generate plot.
    #' @param data {data.frame}: The data to plot
    #' @param out_name {string}: The file name to output. Should be a png and
    #'                           should not contain the plots directory.
    #' @param ... {any}: Extra arguemnts to provide to plot callback.
    #' @param width {integer}: The width of the plot to save in pixels. Has a default value.
    #' @param height {integer}: The height of the plot to save in pixels. Has a default value.
    #' @return out_path {string}: Path to saved plot file.

    out_name <- check_png_ext(out_name, TRUE)
    out_path <- construct_plot_path(out_name)
    
    if (file_exists(out_path) && !overwrite_ext_plots) {
        logger("Plot already saved at: ", quotes(out_path), ".")
        logger("Skipping plot")
    } else {
        logger("DEBUG", "Plotting at: ", quotes(out_path), " ...")
        png(out_path, width, height)
        plot_callback(data, ...)
        dev.off()
    }

    return(out_path)
}

wrap_histogram <- function(...) {
    #' Wrapper to wrap_plot but for histograms.
    #' @param ... {any}: Arguments to wrap_plot excluding the plot_callback.
    #' @return out_path {string}: Path to saved plot file.

    wrap_plot(hist, ...)
}

wrap_scatter <- function(abline_h, abline_col, abline_name, ...) {
    #' Wrapper to wrap_plot but for scatterplot. Addionally plots a line on the scatterplot.
    #' @param abline_y {integer | NULL}: Height of the abline.
    #' @param abline_col {integer | NULL}: Colour of the abline.
    #' @param abline_name {string}: Label for the abline.
    #' @param ... {any}: Arguments to wrap_plot excluding the plot_callback.
    #' @return out_path {string}: Path to saved plot file.

    plot_with_abline <- function(...) {
        #' Custom plot callback to generate plot with an abline
        #' @param ... {any}: Arguments to provide to plot.
        #' @return {NULL}

        plot(...)
        abline(h = abline_h, col = abline_col, main = abline_name)
    }

    wrap_plot(plot_with_abline, ...)
}

off_diag <- function(df) {
    #' Returns the off diagonal elements of a data.frame. Returns the upper triangle
    #' @param df {data.frame}: A data.frame to get the off diagonal elements from.
    #' @return {data.frame}: Data.frame containing only the off diagonal elements.

    upper_tri <- df[upper.tri(df)]
    return(upper_tri)
}

get_mpheno_args <- function(mpheno = 1) {
    #' Generates mpheno arguments for GCTA given a (n+2) mpheno column number.
    #' @param mpheno {number}: The (n+2)th column to look at in the covariate
    #'                         file. Defaults to 1
    #' @return mpheno_args {string}: GCTA mpheno arguments.
    
    logger("DEBUG", "Generating mpheno number args with mpheno = ", mpheno, ".")
    mpheno_args <- paste(gcta_fgs$mpheno, mpheno)
    return(mpheno_args)
}

get_thread_args <- function(num_threads = gcta_default_tnum) {
    #' Generates thread number arguments for GCTA given a number of threads.
    #' @param num_threads {number}: How many threads for GCTA to use. Defaults
    #'                              to 2.
    #' @return thread_args {string}: GCTA thread num arguments.
    
    logger("DEBUG", "Generating thread number args with ", num_threads, " threads.")
    thread_args <- paste(gcta_fgs$tnum, num_threads)
    return(thread_args)
}

get_trait_name <- function(suffix) {
    #' Maps suffix to phenotype / trait name. E.g. "_binary1 -> "Binary 1".
    #' @param suffix {string}: Suffix of the filename corresponding to the phenotype.
    #' @return trait_name {string}: Name of the trait in a friendly format.

    if (suffix == "") {
        logger("DEBUG", "Quantitative phenotype ", phenotype, ".")
        return(phenotype)
    }

    logger("DEBUG", "Mapping suffix: ", quotes(suffix), " to trait name...")
    binary_pattern <- "^_?([A-Za-z]+)(\\d+)$"
    trait_name <- gsub(binary_pattern, "\\1 \\2", suffix)
    logger("INFO", "Trait name: ", quotes(trait_name), ".")
    return(trait_name)
}

get_pheno_path <- function(pheno_suffix) {
    #' Gets the phenotype file path for a given trait. the possible traits are:
    #'   1. "": the original, quantitative trait.
    #'   2. "_binary1": The binary derivative in which those scoring in the top 20%
    #'                  of the phenotype are scored 1 = case and the remainder
    #'                  0 = control.
    #'   2. "_binary2": The binary derivative in which those scoring in the top 20%
    #'                  of the phenotype are scored 1 = case and those scoring in
    #'                  the bottom 30% of the phenotype are scored 0 = control. The
    #'                  remainder are scored N/A.
    #' @param pheno_suffix {string}: The suffix corresponding to the phenotype file
    #'                               to compute a filepath for.
    #' @return pheno_path {string}: The complete file path of the phenotype file.

    logger("DEBUG", "Retrieving pheno path for suffix: ", quotes(pheno_suffix), ".")
    phenotype_file_prefix <- space_to_underscore(phenotype)
    pheno_file_name <- paste0(phenotype_file_prefix, pheno_suffix, exts$phen)
    pheno_path <- construct_phenotypes_path(pheno_file_name)
    logger("DEBUG", "Phenotype path: ", quotes(pheno_path), ".")
    return(pheno_path)
}

####

init <- function() {
    #' Performs initialisation for the script.
    #' @return {NULL}

    logger("Initialising...")
    logger("Log level set at ", app_log_level, " ", brackets(app_log_index), ".")

    logger("Initialising directories...")
    mkdir_if_not_exist(plots_out_dir)
    mkdir_if_not_exist(gcta_out_dir)
    mkdir_if_not_exist(out_dir)

    logger("Initialisation complete!")
}

init()

# GCTA Flags
gcta_fgs <- create_object(list("bfile", "mgrm", list("mkgrm" = "make-grm"), "out",
                               list("tnum" = "thread-num"), "autosome", "pheno",
                               "mpheno", "reml", "covar"),
                          named_flag)


# File extensions
exts <- create_object(list("phen", "txt", "png", "cov", "eigenvec", "eigenval",
                           "grm", "bed", "bim", "bin", "fam", "gz", "id"),
                      ext)

# === Main ===

grm_build <- function(run_grm_build) {
    #' Builds the GRM with GCTA from the QC dataset.
    #' @param run_grm_build {boolean}: Whether to run the grm build
    #'                                 process. If FALSE, will try to
    #'                                 use existing file.
    #' @return grm_basepath {string}: Base path to grm files generated.
    
    tnum_args <- get_thread_args()
    gcta_args <- paste(gcta_fgs$mkgrm, gcta_fgs$autosome, tnum_args)
    out_name <- "grm"
    
    if (!run_grm_build) {
        expected_grm_path <- add_extension(construct_gcta_out_path(out_name),
                                           exts$grm)
        
        if (!file_exists(expected_grm_path)) {
            logger("ERROR", "GRM path not found.")
            stop()
        }

        return(expected_grm_path)
    }

    logger("Building GRM...")
    grm_basepath <- gcta_qc_data(gcta_args, out_name)
    logger("DEBUG", "GRM Basepath: ", quotes(grm_basepath), ".")
    logger("GRM Build Complete!")
    return(grm_basepath)
}

estimate_greml_var <- function(grm_basepath) {
    #' Estimates the proportion of phenotypic variance due to genome-wide
    #' SNPs using GCTA.
    #' @param grm_basepath {string}: Basepath to the grm file to estimate from.
    
    estimate_phen_var_prop <- function(suffix) {
        #' Given a suffix and mpheno value, estimates the phenotypic variance
        #' explained by additive genome-wide SNPs for the given phenotype.
        #' @param suffix {string}: The suffix of the phenotype file name which in
        #'                         turn, encodes the phenotype variant.
        #' @return hsq_basepath {string}: Base path to the hsq output from GCTA.
        
        trait_name <- get_trait_name(suffix)
        logger("Estimating Phenotyping Variance Proportion for phenotype: ", trait_name)
        pheno_path <- construct_phenotypes_path(suffix)
        mpheno <- 1
        mpheno_args <- get_mpheno_args(mpheno)
        gcta_args <- paste(gcta_fgs$bfile, bfile, gcta_fgs$pheno, pheno_path,
                           mpheno_args, gcta_fgs$reml)
        out_name <- paste0("greml_var", suffix, "_", mpheno)
        hsq_basepath <- gcta(gcta_args, out_name)
        logger("DEBUG", "HSQ Basepath: ", quotes(hsq_basepath), ".")
        return(hsq_basepath)
    }
    
    read_hsq_res <- function(hsq_basepath) {
        #' Reads a greml file containing the hsq result from GCTA.
        #' @param hsq_basepath {string}: Base path to the hsq output from GCTA.
        #' @return hsq {data.frame}: Data.frame containing the hsq result.
        
        hsq_path <- add_extension(hsq_basepath, exts$hsq)
        logger("Reading HSQ result from: ", quotes(hsq_path), ".")

        hsq <- wrap_read_table(hsq_path, fill = TRUE)

        log_df(hsq, "HSQ Results")
        return(hsq)
    }

    for (suffix in phenotype_suffixes) {
        logger("Inspecting phenotype: ", quotes(suffix), ".")
        hsq_basepath <- estimate_phen_var_prop(suffix)
        hsq <- read_hsq_res(hsq_basepath)
    }
}

unrelated_indvs <- function(grm_basepath) {
    #' Investigates properties of the GRM.
    #' @param grm_basepath {string}: Basepath (no extenision) to the GRM file.
    
    remove_relatedness <- function(grm_path) {
        #' Removes relatness between individuals who may have an affect
        #' on the estimate of heritability. This reduces the maximum
        #' relatedness coefficient via GCTA with --grm-cutoff.
        #' @param grm_path {string}: Path to grm file.
        #' @return grm_path_rr {string}: Path to new grm file with individuals
        #'                               not meeting the relatedness threshold
        #'                               removed from the dataset.

        gcta_args <- paste(gcta_fgs$grm, grm_path, gcta_fgs$grm_cutoff,
                           grm_relatedness_threshold, gcta_fgs$mkgrm)
        out_name <- "grm_rel_rmvd"

        grm_path_rr <- gcta_qc_data(gcta_args, out_name)
        return(grm_path_rr)
    }


    read_grml_no_rel <- function() {
        grm.nr <- read_GRMBin(qimrx_nr_path)
    }


    names_grm <- c("IND_1", "IND_2", "SNP_NUM", "REL")
    names(grm) <- names_grm
    log_df(grm, "GRM (named)")
    
    grm.diag <- diag(grm)
    logger("GRM diag length: ", length(grm.diag), ".")

    log_df(grm.diag, "GRM (diag)")

    grm.off.diag <- off_diag(grm)

    hist_name <- add_extension("grm.diag", exts$png)
    wrap_hist(grm.diag, hist_name, breaks = 2500, freq = FALSE, xlab = "GRM Diagonals",
              xlim = c(0.95, 1.2), main = "GRM Diag Distribution")

    hist_name <- add_extension("grm.off.diag", exts$png)
    wrap_hist(grm.off.diag, hist_name, breaks = 200, freq = FALSE, xlab = "GRM Off-Diagonals",
              xlim = c(0.1, 1.1), main = "GRM Off-Diag Distribution")

    threshold <- 0.1
    grm.off.diag.clipped <- grm.off.diag[which(grm.off.diag > threshold)]
    wrap_hist(grm.off.diag.clipped, hist_name, breaks = 200, freq = FALSE,
              xlab = "GRM Off-Diagonals", xlim = c(0.1, 1.1),
              main = "GRM Off-Diag Distribution")
}

partition_variance <- function(qimrx_cleaned_path) {
    maf_snps_path <- file.path()

    get_gcta_args <- function(snps_path) {
        gcta_args <- paste(gcta_fgs$bfile, qimrx_cleaned_path, gcta_fgs$extract,
                           snps_path, gcta_fgs$autosome, gcta_fgs$mkgrm,
                           gcta_fgs$keep, qimrx_nr_id_path, gcta_fgs$tn, 2)
        logger("DEBUG", "GCTA Args: ", quotes(gcta_args), ".")
        return(gcta_args)
    }

    out_name <- "qimrx_nr_maf_snps"

    gcta_args <- get_gcta_args(maf_snps_path)
    gcta(gcta_args, out_name)
    
    gcta_args <- get_gcta_args(top_snps_path)
    gcta(gcta_args, out_name)
}

estimate_phen_var <- function(mgrm_path) {
    logger("Estimating Phenotypic Variance Proportion...")
    mpheno <- 1
    mpheno_args <- get_mpheno_args(mpheno)
    gcta_args <- paste(gcta_fgs$mgrm, mgrm_path, gcta_fgs$pheno,
                       ht_t_x_pheno_path, mpheno_args, gcta_fgs$reml)
    out_name <- "qimrx_multi_nr"
    
    gcta(gcta_args, out_name)
}

# === Main ===

args <- commandArgs(trailingOnly = TRUE)

# Run all by default
run_make_grm <- TRUE
estimate_grm_var <- TRUE

if (length(args) > 0) {
    if (!("grm" %in% args)) {
        run_make_grm <- FALSE
    }

    if (!("grm_estimate" %in% args)) {
        estimate_grm_var <- FALSE
    }
}

    if (!!("related", %in% args)) {
        run_relatedness <- FALSE
    }
}

grm_basepath <- grm_build(run_make_grm)

if (estimate_grm_var) {
    estimate_greml_var()
}

logger("DONE!")
