# Course: STAT3306: Statistical Analysis of Genetic Data
# Taken: Semester 2, 2024
# Project Part 1
# Author: Hugo Burton
# Student #: s4698512
# Due Date: 2pm 04 October 2024

install_if_missing <- function(packages) {
    missing_pkgs <- packages[!(packages %in% installed.packages()[, "Package"])]

    if (length(missing_pkgs) > 0) {
        install.packages(missing_pkgs)
    }
    invisible(NULL)
}

# === Packages ===

required_packages <- c("crayon", "qqman", "xtable")
install_if_missing(required_packages)
invisible(lapply(required_packages, require, character.only = TRUE))

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
find_related <- TRUE

# Thresholds
genotype_threshold <- 0.05
het_threshold <- 0.2
related_threshold <- 0.025  # Plink Default
hwe_threshold <- 0.001
freq_threshold <- 0.01
maf_threshold <- 0.1

# Data Paths
course_shared_data_path <- file.path("/data/STAT3306")
project_data <- file.path(course_shared_data_path, "Project")
data_path <- file.path(project_data, "Data")
phenotypes <- file.path(project_data, "Phenotypes")

# Plink
plink_datafile_basename <- "test"

# Out Paths
plots_out_dir <- file.path("./plots")
plink_out_dir <- file.path("./plink_out")
out_dir <- file.path("./out")

plot_w <- 700
plot_h <- 400

# Overrides
overwrite_ext_plots <- FALSE
overwrite_plink_out <- FALSE

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
    padded <- sprintf(paste0("%-", width, "s"), paste0(...))
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

to_str <- function(x, collapse = ", ") {
    #' Converts a list or vector to a string separated by commas.
    #' @param x {list|vector}: List or vector to convert to string.
    #' @param collapse {string}: Separator to split items in list
    #'                           or vector. Defaults to ', '.
    #' @return {string}: String representation of list.
    
    if (is.list(x)) {
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
        if (is.list(arg)) {
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
        logger("DEBUG", "String is empty")
        return(str)
    }

    logger("TRACE", "Converting ", quotes(str), " to title case...")
    words <- strsplit(str, " ")[[1]]
    titlecased_words <- paste0(toupper(substring(words, 1, 1)), tolower(substring(words, 2)))
    result <- paste(titlecased_words)

    logger("DEBUG", "Title case version: ", quotes(result), ".")
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
                        digits = 2, line_spacing_factor = 1, hide_row_names = FALSE) {
    #' Generates a LaTeX table given a data.frame and saves to a file.
    #' @param data {data.frame}: The data.frame to output as a LaTeX table.
    #' @param out_name {string}: Filename to save the output as.
    #' @param table_align {string}: Latex coding for aligning columns.
    #' @param caption {string|NULL}: Optional caption for the table.
    #' @param col.names {vec|NULL}: Optional column names to provide to the table.
    #' @param digits {integer}: Number of decimal places to display numbers as.
    #' @param line_spacing_factor {integer}: Line spacing factor for LaTeX table.
    #' @param hide_row_names {bool}: Whether to include row names from the data.frame
    #'                               as the first column in the table. Defaults to FALSE.
    #' @return path {string}: Path to saved LaTeX table.
    
    table <- xtable(data, align = table_align, caption = caption, digits = digits)

    if (!is.null(col.names)) {
        colnames(table) <- col.names
    }

    latex <- print.xtable(table, print.results = FALSE, table.placement = "htb",
                          comment = FALSE, include.rownames = !hide_row_names)

    latex <- gsub("\\begin{tabular}",
                  paste0("\\renewcommand{\\arraystretch}{",
                         line_spacing_factor, "}\n\\begin{tabular}"
                        ),
                  latex,
                  fixed=TRUE)
    
    logger("DEBUG", "Writing latex table...")
    latex_path <- wrap_write(latex, out_name)
    return(latex_path)
}


match_not_log <- function(path) {
    #' Given an extensionless file path, constructs
    #' a Regex pattern to exclude matching .log files.
    #' @param path {string}: The extensionless pattern to match.
    #' @return {string}: The Regex pattern which includes a negative
    #'                   lookahead to log files.
   
    # Check the path does not already contain an extension.
    ext_pattern <- paste0("(", to_str(regex_escape(exts), collapse = "|"), ")$")
    logger("TRACE", "Extension pattern: ", quotes(ext_pattern), ".")

    ends_with_extension <- grepl(ext_pattern, path)
    if (ends_with_extension) {
        logger("ERROR", "Unexpected! Path ends with extension ", quotes(ext_pattern), ".")
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
    #' @param content {string}: The content to write.
    #' @param basename {string}: Name of the file to write to.
    #' @return path {string}: The full save path where the file was saved.

    basename <- check_txt_ext(basename, exts$txt)
    path <- construct_out_path(basename)

    if (file_exsits(path)) {
        logger("WARN", "Overwriting file at path ", quotes(path), ".")
        delete_file(path)
    }

    if (!is.character(content)) {
        logger("ERROR", "Cannot write content. Type is not character.")
        return(path)
    }
    
    logger("DEBUG", "Writing data to ", quotes(path), " ...")

    write(data, path)
    logger("DEBUG", "Writing complete. Cleaning up...")
    rm(data)
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
    #' @param data {data.frame}: The data to write
    #' @param basename {string}: The basename excluding the out dir to  write
    #'                           the file at. Ideally should include extension
    #'                           but if it doesn't this function will add it.
    #' @param row.names {boolean}: Whether to include row names. Disabled by default.
    #' @param col.names {boolean}: Whether to include col names. Enabled by default but
    #'                             if the data contains default columns, will not be
    #'                             written.
    #' @param sep {string}: The separator to delimit between columns in the table.
    #' @param quote {boolean}: Whether to include quotes for strings in the table
    #'                         data. Disabled by default.
    #' @return path {string}: The full save path where the table was saved.

    basename <- check_txt_ext(basename, exts$txt)
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

plink_orig_data <- function(plink_args, out_name = NULL) {
    #' Makes a plink call with the original dataset and specified arguments.
    #' Serves as wrapper to plink function.
    #' @param plink_args {string}: Arguments provided to plink.
    #' @param out_name {string}: Basename for plink to output to. Should
    #'                           exclude plink directory. If NULL, output
    #'                           is piped to console.
    #' @return plink_out_path {string}: Relative path from script to plink
    #'                                  output. Notably, this path does not
    #'                                  contain extension added by plink as
    #'                                  this differs depending upon arguments
    #'                                  provided to plink.

    data_files_pattern <- file.path(data_path, plink_datafile_basename)
    logger("TRACE", "Plink Data Files Pattern ", quotes(data_files_pattern), ".")

    plink(data_files_pattern, plink_args, out_name)
}

plink <- function(bfile, plink_args, out_name = NULL) {
    #' Makes a plink call with the given dataset and specified arguments.
    #' @param bfile {string}: Path to binary fileset. Should not include
    #'                        extension.
    #' @param plink_args {string}: Arguments provided to plink.
    #' @param out_name {string}: Basename for plink to output to. Should
    #'                           exclude plink directory. If NULL, output
    #'                           is piped to console.
    #' @return plink_out_path {string}: Relative path from script to plink
    #'                                  output. Notably, this path does not
    #'                                  contain extension added by plink as
    #'                                  this differs depending upon arguments
    #'                                  provided to plink.

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
        plink_pattern <- match_not_log(plink_out_path)
        if (file_exists(plink_pattern, TRUE) && !overwrite_plink_out) {
            logger("Matching file(s) already exists at: ", quotes(plink_out_path), ".")
            return(plink_out_path)
        }

        logger("Running: ", quotes(plink_cmd))
        shell_call(plink_cmd)
        logger("Plink results directed to ", quotes(plink_out_path), ".")
    }

    gc()
    return(plink_out_path)
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
    logger("TRACE", "Constructed plink out path: ", path, ".")
    return(path)
}

construct_plink_out_path <- function(...) {
    #' Constructs a file path in the plink out directory given a basename.
    #' @param ... {string}: The name of the file (potentially split) without
    #'                      the plink out directory prefixed.
    #' @return path {string}: The relative path to the file.

    basename <- paste0(...)
    file.path(plink_out_dir, basename)
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

save_removed_indices <- function(table, ind_to_remove, out_cols, out_name) {
    #' Given a data.frame and another data.frame containing indices to remove
    #' removes those indices from the table and outputs the result to a file.
    #' Output will go to out_dir.
    #' @param table {data.frame}: Data.frame to cull from.
    #' @param ind_to_remove {data.frame}: Data.frame containing indices to remove.
    #' @param out_cols {list}: Columns to output to file
    #' @param out_name {string}: File name to save to.
    #' @return ind_out_path {string}: Path to saved file.
    
    logger("Saving removed indices to ", quotes(to_str(out_name)), ".")
    logger("DEBUG", "Selected columns: ", quotes(to_str(out_cols)), ".")

    file <- table[ind_to_remove, out_cols, drop = FALSE]
    ind_out_path <- wrap_write_table(file, out_name, col.names = FALSE)
    return(ind_out_path)
}

remove_indices_by_threshold <- function(table, thresh_col, threshold, out_cols, out_name) {
    #' Given a data.frame and a threshold, removes rows from the data.frame where a value falls
    #' above a threshold and saves to a file.
    #' @param table {data.frame}: Data.frame to filter from
    #' @param thresh_col {string | integer}: The column containing the threshold value to cull from.
    #' @param out_cols {list}: Column(s) to output to file.
    #' @param out_name {string}: File name to save to.
    #' @return ind_out_path {string}: Path to saved file.

    logger("DEBUG", "Removing indices by threshold (", threshold, ") in table from col: ",
           quotes(thresh_col), ".")
    ind_to_remove <- which(table[[thresh_col]] > threshold)
    save_removed_indices(table, ind_to_remove, out_cols, out_name)
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

log_indvs_to_remove <- function(num_indvs) {
    #' Logs if there are individuals to remove.
    #' @param num_indvs {integer}: Number of individuals to remove
    #'                             (from a data.frame).
    #' @return {NULL}

    if (num_indvs == 0) {
        logger("No individuals to remove")
    } else {
        logger("There are ", num_indvs, " individuals to remove")
    }
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
    out_path <- file.path(plots_out_dir, out_name)
    
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

genomic_inflation_factor <- function(d, df = 1) {
    #' Computes the genomic inflation factor for a regression analysis
    #' @param d {data.frame}. Data containing P column for computing GIF from.
    #' @param df {integer}. Degrees of freedom. Defaults to 1.
    #' @return lambda {float}: The genomic inflation factor.

    logger("DEBUG", "Computing genomic inflation factor with ", df, " degree of freedom...")
    
    median_val <- median(d$P)
    logger("TRACE", "Median Value: ", median_val)
    
    chisq_med <- qchisq(0.5, 1)
    qchisq(1 - median_val, df) / chisq_med
}

excess_missing_genotypes <- function(data_path, extension, out_cols, histogram, suffix) {
    #' Finds excess missingness in the genetic data. Generic method. Can be used
    #' to find misingness in individuals or SNPs. Includes plotting functionality.
    #' @param data_path {string | NULL}: Path to dataset. If NULL, original dataset
    #'                                   will be used.
    #' @param extension {string}: File extension specifiying what type of missingness to
    #'                            exclude for.
    #' @param out_cols {vec(string)}: Vector of columns to include in save file.
    #' @param histogram {boolean}: If true, will plot histogram of missingness.
    #' @param suffix {string}: Indicates context and prevents save file name collisions.
    #' @return missing_file_path {string}: File containing missing indvs / SNPs to remove.

    if (length(suffix) == 0) {
        logger("ERROR", "Suffix for missing plot name cannot be empty.")
    }

    logger("Checking for missing genotypes (", extension, ")...")

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
        suffix_bracketed <- brackets(suffix)
        plot_title <- paste0("Distribution of Missingness ", suffix_bracketed)
        wrap_histogram(missing$F_MISS, hist_basename, xlab = paste0("Missingness ",
                       suffix_bracketed), main = plot_title)
        num_indvs_to_remove <- sum(missing$F_MISS > genotype_threshold)
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
    #' Performs initialisation for the script.
    #' @return {NULL}

    logger("Initialising...")
    logger("Log level set at ", app_log_level, " ", brackets(app_log_index), ".")

    logger("Initialising directories...")
    mkdir_if_not_exist(plots_out_dir)
    mkdir_if_not_exist(plink_out_dir)
    mkdir_if_not_exist(out_dir)

    logger("Initialisation complete!")
}

init()

# Plink flags
pl_fgs <- create_object(list("remove", "missing", list("mb" = "make-bed"),
                             "hardy", "het", "mind", "pheno", "covar",
                             list("dup_vars" = "list-duplicate-vars"), "out",
                             "bfile", "chr", "freq", "exclude", "mpheno",
                             "pca", "linear", "logistic", "assoc", "clump",
                             list("cp1" = "clump-p1"), list("cp2" = "clump-p2"),
                             list("cr2" = "clump-r2"), list("ckb" = "clump-kb"),
                             list("rel_cutoff" = "rel-cutoff"),"keep",
                             list("miss_pheno_neg_9" = "1"),
                             list("covar_name" = "covar-name")),
                        named_flag)

# File extensions
exts <- create_object(list("phen", "imiss", "lmiss", "het", "assoc", "hwe",
                           "frq", "txt", "png", "eigenvec", "eigenval",
                           "qassoc", "linear", "clumped", "rel", "id", "logistic"),
                      ext)

# === Main ===

# ====== Analysis ======

quality_control <- function(perform) {
    #' Performs stage 1 quality control on individuals in the genomic dataset.
    #' @param perform {bool}: Whether to perform quality control or not. If disabled, will
    #'                        attempt to return existing QC dataset.
    #' @return data_subset_related_path {string}: Path to quality controlled dataset.
    
    if (perform) {
        logger("Performing Individual Quality Control")
    }

    find_individual_missing_genotypes <- function(histogram) {
        #' Finds individuals with excess missing genotypes to remove from the dataset
        #' @param histogram {boolean}: If true, will plot a histogram of the frequency
        #'                             of missingness.
        #' @return missing_file_path {string}: The file path to the file containing
        #'                                     the individuals to remove.

        excess_missing_genotypes(NULL, exts$imiss, fam_ind_cols, histogram, "individual")
    }

    find_outlying_homozygosity <- function(plot) {
        #' Finds individuals with outlying homozygosity values to remove from the dataset
        #' @param plot {boolean}: If true, will plot a histogram of the frequency of hz
        #'                        freqs and a scatterplot of their distribution across
        #"                        the geneome.
        #' @return het_ind_file_path {string}: The file path to the file containing
        #"                                     the individuals to remove.

        logger("Checking for outlying homozygosity values...")

        hz_name <- "hz"
        het_basename <- plink_orig_data(pl_fgs$het, hz_name)
        het_out_path <- add_extension(het_basename, exts$het)

        logger("DEBUG", "Reading het table...")
        het <- wrap_read_table(het_out_path)
        log_df(het, "het")

        if (plot) {
            logger("INFO", "Plotting het histogram")
            wrap_histogram(het$F, add_extension("fhet_hist", exts$png), xlab = "Homozygosity",
                           main = "Distribution of Homozygosity")

            logger("INFO", "Plotting het scatterplot")
            wrap_scatter(0.05, "red", "Distribution of abs heterozygosity", abs(het$F),
                         add_extension("fhet_scatter", exts$png))
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
        #' Finds individuals from each pair of samples with observed genomic relatedness
        #' above a given threshold. NP-Hard problem so will used cached result if exists.
        #' @param threshold {float}: The threshold to exclude related samples by.
        #' @return related_file_path {string}: File path to a file containing individuals
        #'                                     to remove.

        basename_with_extension <- function(basename) {
            add_extension(basename, exts$rel, exts$id)
        }

        logger("Checking for related samples...")
        
        related_name <- "related_samples"
        expected_related_path <- construct_plink_out_path(basename_with_extension(related_name))
        logger("DEBUG", "Checking for cached result at path ", quotes(expected_related_path), "...")
        if (file_exists(expected_related_path) || !perform) {
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
        #' @param ... {character}: File paths to individual files containing IDs of
        #'                         individuals to remove. These files are expected
        #'                         to have at least two columns: FID and IID.
        #' @return {string}: File path to the combined output file with individuals
        #"                   to remove.
        
        combined_basename <- add_extension("remove.combined.samples", exts$txt)
        combined_file_out_path <- construct_out_path(combined_basename)
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

        # Write will return the saved path
        wrap_write_table(combined_ind, combined_basename, col.names = FALSE)
    }

    remove_bad_individuals <- function(remove_path) {
        #' Removes individuals from the dataset who are specified in the provided file.
        #' @param remove_path {string}: File path to file containing individuals to remove
        #'                              from the dataset.
        #' @return out_path {string}: File path to subset of original dataset with specified
        #'                            individuals removed.

        out_name <- "test_indv_subset"
        if (!perform) {
            plink_path <- construct_plink_out_path(out_name)
            logger("Using qc path: ", quotes(plink_path), ".")
            return(plink_path)
        }

        logger("Removing bad samples...")
        plink_flags <- paste(pl_fgs$mb, pl_fgs$remove, remove_path)
        plink_orig_data(plink_flags, out_name)
    }

    keep_related_samples <- function(related_path, data_subset_path) {
        #' Keeps only related samples from the data subset who are listed in the related file.
        #' @param related_path {string}: Path to file containing list of individuals to keep.
        #' @param data_subset_path {string}: Path to subset of data with some individuals
        #'                                   already removed.
        #' @return out_path {string}: Path to a further subset of the data

        
        out_name <- "test_indv_related"
        if (!perform) {
            plink_path <- construct_plink_out_path(out_name)
            return(plink_path)
        }

        logger("Retaining unrelated samples...")
        plink_args <- paste(pl_fgs$keep, related_path, pl_fgs$mb)
        plink(data_subset_path, plink_args, out_name)
    }

    if (perform) {
        missing_ind_path <- find_individual_missing_genotypes(TRUE)
        het_ind_path <- find_outlying_homozygosity(TRUE)

        combined_ind_path <- combine_remove_files(missing_ind_path, het_ind_path)
    } else {
        # Unused - provide dummy path
        combined_ind_path <- ""
    }

    data_subset_path <- remove_bad_individuals(combined_ind_path)
    gc()

    # Keep only individuals specified here
    if (find_related) {
        related_path <- find_related_samples()
        data_related_subset_path <- keep_related_samples(related_path, data_subset_path)

        return(data_related_subset_path)
    } else {
        return(data_subset_path)
    }
}

sample_qc <- function(data_subset_path, perform) {
    #' Performs quality control with respect to SNPs in the genomic dataset
    #' @param data_subset_path {string}: Path to subset of data from individual
    #'                                   quality control step.
    #' @param perform {bool}: Whether to perform quality control or not. If disabled, will
    #'                        attempt to return existing QC dataset.
    #' @return qc_data_maf_path {string}: Path to subset of data from this step
    #'                                    of the quality control process.
    
    if (perform) {
        logger("Performing SNP Quality Control")
    }

    missing_snps <- function(histogram) {
        #' Computes missingness of the SNPs in the dataset.
        #' @param histogram {boolean}: Will plot histogram of missingness if true.
        #' @return lmsis {data.frame}: Data.frame containing missingness of SNPs.

        missing_file_path <- excess_missing_genotypes(data_subset_path, exts$lmiss, 2,
                                                      histogram, "SNPs")
        lmiss <- wrap_read_table(missing_file_path)
        return(lmiss)
    }

    hw_eq <- function(histogram) {
        #' Computes a list of genotype counts and Hardy-Weinberg equilibrium
        #' exact test statistics. Plots histogram if enabled.
        #' @param histogram {boolean}: Will plot histogram of hwe if true.
        #' @return hwe {data.frame}: Data.frame containig hwe exact tests of SNPs.

        logger("Computing Hardy-Weinberg Equilibrium")
        hwe_name <- "hwe"
        hwe_basename <- plink(data_subset_path, pl_fgs$hardy, hwe_name)
       
        hwe_path <- add_extension(hwe_basename, exts$hwe)
        hwe <- wrap_read_table(hwe_path)
        
        if (histogram) {
            wrap_histogram(hwe$P, "hwe_deviations.png", xlab = "p-values of HWE test on SNPs",
                           main = "Distribution of HWE p-values")
        }

        return(hwe)
    }

    min_allele_freq <- function(histogram) {
        #' Indentify SNPs with low minor allele frequency.
        #' @param histogram {boolean}: if true, frequency of MAF is plotted visually.
        #' @return {data.frame}: of low minor allele frequencies.

        logger("Computing Minor Allele Frequencies")
        min_allele_name <- "minor_allele_freq"
        min_allele_basename <- plink(data_subset_path, pl_fgs$freq, min_allele_name)

        min_allele_path <- add_extension(min_allele_basename, exts$frq)
        freq <- wrap_read_table(min_allele_path)

        if (histogram) {
            wrap_histogram(freq$MAF, add_extension("maf_distribution", exts$png),
                           xlab = "Minor Allele Frequency (ref = 1)",
                           main = "Distribution of Minor Allele Frequencies")
        }

        return(freq)
    }

    remove_snps <- function(lmiss, hwe, freq) {
        #' Removes snps from dataset which do not meet thresholds for missingness,
        #' hwe, or minor allele frequency.
        #' @param lmiss {data.frame}: Data.frame for missingness of SNPs.
        #' @param lmiss {data.frame}: Data.frame for hardy statistics of SNPs.
        #' @param lmiss {data.frame}: Data.frame for minor allele frequencies of SNPs.
        #' @return remove_snps_path {string}: Path to a txt file containing
        #'                                    combined SNPs to remove.

        logger("INFO", "Computing SNPs to remove...")
        
        ind_to_remove <- unique(c(
                                    which(lmiss$F_MISS > genotype_threshold,
                                    freq$MAF < freq_threshold),
                                    which(hwe$P < hwe_threshold)
                                 )
                               )
        snp_id_col <- 2
        remove_snps_path <- save_removed_indices(freq, ind_to_remove, snp_id_col,
                                                 add_extension("remove.SNPs", exts$txt))
        return(remove_snps_path)
    }

    exclude_snps <- function(remove_snps_path) {
        #' Removes specified SNPs from the dataset and saves a new subset.
        #' @param remove_snps_path {string}: Path to file containing SNPs to exclude.
        #' @return {string}: New data subset basename path (excluding extension).

        logger("INFO", "Excluding SNPs...")
        out_name <- "test_qc"
        plink_args <- paste(pl_fgs$mb, pl_fgs$exclude, remove_snps_path)
        plink(data_subset_path, plink_args, out_name)
    }

    compare_minor_allele_freqs <- function(freq, do_hist, do_plot) {
        #' Compares minor allele frequencies of the data subset to known
        #' reference allele frequencies. Determines SNPs which fall outside
        #' of threshold and saves file containing SNPs to remove.
        #' @param freq {data.frame}: Data.frame containing minor allele frequencies
        #'                           of the current data subset.
        #' @param do_hist {boolean}: If true, will plot a histogram of minor allele
        #'                           frequencies for this dataset.
        #' @param do_plot {boolean}: If true, will save scatterplot of corrected
        #'                           MAFs both with and without threshold cutoff.
        #' @return remove_out_path {string}: Path to file containins SNPs which fall
        #'                                   outside threshold and are to be removed.

        logger("Comparing Minor Allele Frequencies...")
        snp_ref_path <- file.path(data_path, add_extension("reference_allele_frequencies", exts$txt))
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
        log_df(out_cpy, "Corrected Minor Allele Freqs")

        # Must be calculated regardless of plots
        res <- out_cpy$MAF - out_cpy$V2
        log_df(res, "Residual minor allele frequencies")

        if (do_hist) {
            logger("Plotting histogram of Residual Minor Allele Frequency (MAF)")
            wrap_histogram(res, add_extension("residual_minor_allele", exts$png),
                           xlab = "Residual Minor Allele Frequency",
                           main = "Residual Minor Allele Frequency (r-MAF) (ref = 1)")
        }

        if (do_plot) {
            logger("Plotting allele frequency comparison with reference...")
            min_allele_freq_basename <- "min_allele_freq"
            # y: Observed, x: reference.
            xlab <- "Reference MAF"
            ylab <- "Observed MAF"
            plot_title <- "Reference vs Observed Minor Allele Frequencies"
            wrap_plot(plot, out_cpy$MAF ~ out_cpy$V2,
                      add_extension(min_allele_freq_basename, exts$png),
                      xlab = xlab, ylab = ylab, main = plot_title)

            logger("Removing alleles which deviate significantly from reference...")
            accept_snps <- which(abs(res) <= maf_threshold)
            logger("There were ", length(accept_snps), " accepted.")
            out_accept <- out_cpy[accept_snps, ]
            # y: Observed, x: reference.
            wrap_plot(plot, out_accept$MAF ~ out_accept$V2, add_extension(min_allele_freq_basename,
                                                                          "_comparison", exts$png),
                      xlab = xlab, ylab = ylab, main = paste(plot_title, brackets("Corrected for Ref = 1")))
        }

        remove_snps_indx <- which(abs(res) > maf_threshold)
        to_remove_snps <- freq[remove_snps_indx, "SNP", drop = FALSE]
        remove_basename <- add_extension("maf.referenced.removed", exts$txt)
        wrap_write_table(to_remove_snps, remove_basename, col.names = FALSE)
    }

    exclude_insig_maf <- function(qc_data_path, remove_out_path) {
        #' Creates a subset of the data which excludes insignificant alleles by their MAF.
        #' @param qc_data_path {string}: Path to dataset to exclude from.
        #' @param remove_out_path {string}: Path to file detailing which alleles to remove.
        #' @return qc_data_maf_path {string}: Path to new subset of data.

        out_name <-"test_qc_insif_maf"
        if (!perform) {
            plink_path <- construct_plink_out_path(out_name)
            logger("Using SNP qc path: ", quotes(plink_path), ".")
            return(plink_path)
        }
        
        logger("Excluding insignificant alleles by MAF...")
        plink_args <- paste(pl_fgs$mb, pl_fgs$exclude, remove_out_path)
        plink(qc_data_path, plink_args, out_name)
    }
   
    if (perform) {
        lmiss <- missing_snps(TRUE)
        hwe <- hw_eq(TRUE)
        freq <- min_allele_freq(TRUE)
        
        remove_snps_path <- remove_snps(lmiss, hwe, freq)
        qc_data_path <- exclude_snps(remove_snps_path)

        remove_out_path <- compare_minor_allele_freqs(freq, TRUE, TRUE)
    } else {
        # Provide dummy paths
        qc_data_path <- ""
        remove_out_path <- ""
    }

    qc_data_maf_path <- exclude_insig_maf(qc_data_path, remove_out_path)
    
    gc()
    return(qc_data_maf_path)
}

gwas <- function(qc_data_path) {
    #' Performs a Genome-Wide Association study on the data,
    #' specified by the path.
    #' @param qc_data_path {string}: Path to dataset from the output of
    #'                               quality control.
    #' @return {NULL}

    logger("INFO", "Performing GWAS...")

    # Some scoped variables
    alt_mpheno <- 1
    num_pc <- 10
    
    get_mpheno_args <- function(suffix) {
        #' Constructs the mpheno section of the plink args based on suffix
        #' (which in turn can map to trait).
        #' @param suffix {string}: The suffix of the phenotype file (mapping to trait).
        #' @return args {string}: Mpheno args used for plink.

        args <- paste(pl_fgs$mpheno, alt_mpheno)
        if (grepl("binary", suffix)) {
            logger("DEBUG", quotes("binary"), " found in suffix ", quotes(suffix))
            logger("INFO", "Adding missing-phenotype flag")
            args <- paste(args, pl_fgs$miss_pheno_neg_9)
        }
        
        return(args)
    }

    get_trait_name <- function(suffix) {
        #' Maps suffix to phenotype / trait name. E.g. "_binary1" -> "Binary 1".
        #' @param suffix {string}: Suffix of the filename corresponding to phenotype.
        #' @return trait_name {string}: Name of the trait.
        
        if (suffix == "") {
            return(phenotype)
        }

        logger("DEBUG", "Mapping suffix: ", quotes(suffix), " to trait name...")

        binary_pattern <- "^_?([A-Za-z]+)(\\d+)$"
        trait_name <- gsub(binary_pattern, "\\1 \\2", suffix)
        logger("INFO", "Trait name: ", quotes(trait_name), ".")
    
        return(trait_name)
    }

    get_pheno_path <- function(pheno_suffix) {
        #' Gets the phenotype file path for a given trait. The possible traits are:
        #'   1. "": a quantitative trait.
        #'   2. "_binary1": A binary trait in which those scoring in the top 20% of
        #'                  the phenotype are scored 1 = case and the remainder
        #'                  0 = control.
        #'   3. "_binary2": A binary trait in which those scoring in the top 20% of
        #'                  the phenotype are scored 1 = case and those scoring in
        #'                  the bottom 30% of the phenotype are scored 0 = control.
        #' @param pheno_suffix {string}: The suffix corresponding to the phenotype
        #'                               file compute a filepath for.
        #' @return {string}: Complete file path.
 
        logger("DEBUG", "Retrieving pheno path for suffix ", quotes(pheno_suffix), ".")
        phenotype_file_prefix <- space_to_underscore(phenotype)
        pheno_file_name <- paste0(phenotype_file_prefix, pheno_suffix, exts$phen)
        file.path(phenotypes, pheno_file_name)
    }

    gwas_pheno <- function(pheno_path, pheno_suffix, mpheno_args, covar_file_path, pc) {
        #' Performs association analysis based on the phenotype
        #' defined in the specified file.
        #' @param pheno_path {string}: Path to pheno file for specified trait.
        #' @param pheno_suffix {string}: Suffix of phenotype.
        #' @param mpheno_args {string}: Flags relating to mpheno in plink.
        #' @param covar_file_path {string}: File path to combined covariates.
        #' @param pc {boolean}: Whether we are using principal components as covariates.
        #' @return cov_out_path {string}: Path to phenotype association analysis output.

        logger("Performing Pheno Association Analysis on GWAS ", pc ? "with pc" : "", "...")
        
        out_name <- paste0("gwas_pheno", suffix == "" ? "_quant" : suffix, pc ? "_pc" : "")
        cov_out_path <- construct_plink_out_path(out_name)
        cov_out_pattern <- match_not_log(cov_out_path)
        if (file_exists(cov_out_pattern, TRUE)) {
            logger("INFO", "Covariates already exist at ", quotes(cov_out_pattern), ".")
            return(cov_out_path)
        }

        # No covariates. Compute:
        logger("INFO", "Using covariate file: ", quotes(covar_file_path), ".")
        covar_args <- paste(pl_fgs$covar, covar_file_path)
        regression_args <- pheno_suffix == "" ? pl_fgs$linear : pl_fgs$logistic
        plink_args <- paste(regression_args, pl_fgs$pheno, pheno_path, mpheno_args, covar_args)
        plink(qc_data_path, plink_args, out_name)
    }

    get_pheno_analysis_full_path <- function(pheno_basename, phenotype_suffix) {
        #' Adds the appropriate extensions to the phenotype analysis file basename
        #' based on the trait and if principal components is being used.
        #' @param pheno_basename {string}: Path to pheno analysis out without extensions.
        #' @param phenotype_suffix {string}: Suffix mapping to trait.
        #' @return pheno_path {string}: Full path to pheno analysis file.

        logger("DEBUG", "Calculating path for pheno analysis full path from basename ",
                        quotes(pheno_basename), ".")
        
        if (phenotype_suffix == "") {
            # Quantitative phenotype
            ext_to_add <- paste0(exts$assoc, exts$linear)
        } else {
            # Binary phenotype
            ext_to_add <- paste0(exts$assoc, exts$logistic)
        }
    
        pheno_path <- add_extension(pheno_basename, ext_to_add)
        logger("DEBUG", "Full pheno path: ", quotes(pheno_path), ".")
        return(pheno_path)
    }
 
    gwas_plots <- function(pheno_analysis_path, pc, plot_suffix = "") {
        #' Saves Manhattan and QQ plots for the current pheno analysis.
        #' @param pheno_analysis_path {string}: The path to the pheno analysis file
        #' @param pc {boolean}: Whether principal components is being used.
        #' @param plot_suffix {string}: Suffix mapping to trait.
        #' @return d {data.frame}: Data.frame from phenotype analysis.

        trait_name <- get_trait_name(plot_suffix)
        
        name_plot <- function(plot_type) {
            #' Computes a plot's file name given its type.
            #' @param plot_type {string}: The type of the plot.
            #' @return file_name {string}: File name of plot.

            logger("DEBUG", "Naming plot type: ", quotes(plot_type), ".")
            file_name <- paste0("gwas_", plot_type, plot_suffix, pc ? "_pc" : "")
            file_name <- add_extension(file_name, exts$png)
            logger("DEBUG", "GWAS Plot Name: ", quotes(file_name), ".")
            return(file_name)
        }

        title_plot <- function(plot_type, title_case = TRUE) {
            #' Titles a plot given it's type.
            #' @param plot_type {string}: The type of the plot.
            #' @param title_case {bool}: Whether to use title case on the plot type.
            #'                           Defaults to TRUE.
            #' @return plot_title {string}: Title of the plot.
 
            logger("DEBUG", "Titling plot for type: ", quotes(plot_type), ".")

            if (title_case) {
                plot_type <- title_case(plot_type)
            }

            plot_title <- paste0(plot_type, " plot for ", trait_name,
                                 pc ? paste0(" ", brackets("PC")) : "")
            logger("DEBUG", "Plot title: ", quotes(plot_title), ".")
            return(plot_title)
        }

        logger("DEBUG", "Reading pheno analysis ", quotes(plot_suffix), " ",
                        "Principal components = ", pc, "...")
        d <- wrap_read_table(pheno_analysis_path)
        log_df(d, "D before na omit")

        logger("Retaining ADD tests only...")
        d <- d[which(d$TEST == "ADD"), ]
        log_df(d, "D after retaining ADD tests.")
        
        logger("Removing n/a p-vals...")
        d <- d[!is.na(d$P), ]
        log_df(d, "D after n/a p-val omit")

        logger("Removing 0 p-vals...")
        d <- d[d$P != 0, ]
        log_df(d, "D after 0 p-val omit")

        logger("DEBUG", "Plotting ...")
        man_plot_name <- name_plot("manhattan")
        logger("INFO", "Generating Manhattan Plot ", quotes(man_plot_name), " for trait ",
               quotes(trait_name), " ...")
        wrap_plot(manhattan, d, man_plot_name, main = title_plot("manhattan"))

        qq_plot_name <- name_plot("qq")
        logger("INFO", "Generating QQ Plot ", quotes(qq_plot_name), " ...")
        wrap_plot(qq, d$P, qq_plot_name, main = title_plot("QQ", FALSE))

        return(d)
    }

    compute_lambda <- function(d, suffix, pc) {
        #' Computes Genomic Inflation Factor(GIC) for phenotype analysis.
        #' @param d {data.frame}: Data.frame from phenotype analysis.
        #' @param suffix {string}: Suffix mapping to trait.
        #' @param pc {boolean}: Whether principal components is being used.
        #' @return lambda_bc {float}: The GIC (\lambda) value.

        logger("INFO", "Computing Genomic Inflation Factor (", suffix, ") with PC: ",
                       pc ? "enabled" : "disabled", "...")

        lambda_gc <- genomic_inflation_factor(d, 1)
        logger("INFO", "Lambda_{GC} (", suffix, ") = ", lambda_gc)
        return(lambda_gc)
    }

    combine_covariates <- function() {
        #' Combines the covariate files into a single one, to then be provided.
        #' (Without eigenvectors as covariates).
        #' to plink. Combined covariates file contains headers.
        #' @return combined_covar_file {str}: Path to combined covariate file.

        covariate_basenames <- c("age", "gender")
        
        combined_basename <- add_extension("combined_covariates", exts$txt)
        combined_covar_path <- construct_out_path(combined_basename)

        if (file_exists(combined_covar_path)) {
            logger("DEBUG", "Combined covariates already exists.")
        } else {
            logger("INFO", "Combining covariates...")
            
            combined_covariates <- NULL

            for (covariate in covariate_basenames) {
                covariate_path <- file.path(data_path, add_extension(covariate, exts$txt))
                logger("DEBUG", "Reading covariate ", quotes(covariate), " from path ",
                       quotes(covariate_path), ".")
                covar <- wrap_read_table(covariate_path, header = FALSE)
                colnames(covar) <- c(fam_ind_cols, covariate)

                if (is.null(combined_covariates)) {
                    combined_covariates <- covar
                } else {
                    combined_covariates <- merge(combined_covariates, covar, by = fam_ind_cols)
                }
            }

            log_df(combined_covariates, "Combined covariates")
           
            combined_covar_path <- wrap_write_table(combined_covariates, combined_basename,
                                                    col.names = TRUE)
        }

        return(combined_covar_path)
    }
    
    compute_principal_comps <- function(num_components) {
        #' Computes the specified number of principal components for
        #' for the quality controlled genomic dataset.
        #' @param num_components {integer}: The number of principal components
        #'                                  to compute.
        #' @return pca_eig_vec {string}: Path to eigenvector file.

        # Check for existing PCA
        out_name <- "pca"
        pca_path <- construct_plink_out_path(out_name)   
        pca_eig_val <- add_extension(pca_path, exts$eigenval)
        pca_eig_vec <- add_extension(pca_path, exts$eigenvec)

        # Special: We check for all files present and hence don't rely on
        # plink wrapper function to perform the check for us.
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
        return(pca_eig_vec)
    }

    add_pc_eigvecs_to_covars <- function(existing_covariates_file, pc_eigvec_file) {
        #' Combines pca eigenvector covariates with existing covariates in
        #' a single file.
        #' @param existing_covariates_file {string}: Path to existing (combined) covars file.
        #' @param pc_eigvec_file {string}: Path to eigenvector file (from plink).
        #' @return pc_combined_covars_path {string}: Path to combined covars file with
        #'                                           existing and eigvec covars.

        pc_combined_basename <- add_extension("combined_covariates_pc", exts$txt)
        pc_combined_covars_path <- construct_out_path(pc_combined_basename)

        if (file_exists(pc_combined_covars_path)) {
            logger("DEBUG", "Combined PC covariates already exists.")
        } else {
            logger("INFO", "Combining eigenvectors to existing covariates...")

            existing_covars <- wrap_read_table(existing_covariates_file)
            eigenvecs <- wrap_read_table(pc_eigvec_file, header = FALSE)

            # Set colnames on eigenvector data.frame
            n_cols <- ncol(eigenvecs)
            colnames(eigenvecs)[1:2] <- fam_ind_cols
            if (n_cols > 2) {
                  colnames(eigenvecs)[3:n_cols] <- paste0("PC", seq_len(n_cols - 2))
            } else {
                logger("ERROR", "Expected more than 2 columns for eigenvector file.")
            }
            log_df(eigenvecs, "Eigenvectors")
            log_df(existing_covars, "Existing covars")
            
            pc_combined_covars <- merge(existing_covars, eigenvecs, by = fam_ind_cols)
            log_df(pc_combined_covars, "Combined covariates (with principal components)")

            pc_combined_covars_path <- wrap_write_table(pc_combined_covars, pc_combined_basename,
                                                        col.names = TRUE)
        }

        return(pc_combined_covars_path)
    }

    init_lambdas_df <- function() {
        #' Initialises a data.frame to store lambda values in.
        #' @return {data.frame}: Data.frame to store lambda values in.
        
        lambdas <- data.frame(
            Covariate_Used = character(),
            Fasting_Glucose = numeric(),
            FG_Delta = numeric(),
            Binary_1 = numeric(),
            Binary_1_Delta = numeric(),
            Binary_2 = numeric(),
            Binary_2_Delta = numeric(),
            stringAsFactors = FALSE
        )
        return(lambdas)
    }

    clumping <- function(pheno_pc_path) {
        #' Performs clumping to identify most significant SNP in each LD block.
        #' @param pheno_pc_path {string}: Path to phenotype file for the specified trait.
        #' @return {string}: Path to clump out file from plink.

        logger("Clumping GWAS Results...")
        out_name <- "gwas_pheno_clump"
        clump_p1_val <- 0.5
        clump_p2_val <- 0.5
        clump_r2_val <- 0.2
        clump_kb_val <- 500
        plink_args <- paste(pl_fgs$clump, pheno_pc_path, pl_fgs$cp1, clump_p1_val, pl_fgs$cp2,
                            clump_p2_val, pl_fgs$cr2, clump_r2_val, pl_fgs$ckb, clump_kb_val)
        plink(qc_data_path, plink_args, out_name)
    }

    read_clumps <- function(clump_basename, suffix, pc) {
        #' Saves a copy of the clump file without the last column listing the SNP names
        #' for each clump.
        #' @param clump_basename {string}: Basename path to clump file.
        #' @param suffix {string}: Suffix mapping to trait,
        #' @param pc {boolean}: Whether principal components is being used.
        #' @return {string}: Path to clump file that was saved.

        clump_path <- add_extension(clump_basename, exts$clumped)
        clump <- wrap_read_table(clump_path)

        # Exclude column PC2
        logger("Printing Clump Output")
        clump_out <- clump[, 1:11]
        log_df(clump_out, paste0("Clump output", suffix))
        
        # Write to file
        out_basename <- add_extension(paste0("clumps", suffix, pc ? "_pc" : ""), exts$txt)
        wrap_write_table(clump_out, out_basename, col.names = TRUE)
    }
 
    # Main
    phenotype_suffixes <- list("", "_binary1", "_binary2")
       
    # Combine age and sex covariates
    covar_file_path <- combine_covariates()

    # Compute principal components once
    pc_eigvec_file <- compute_principal_comps(num_pc)
    covar_pc_file_path <- add_pc_eigvecs_to_covars(covar_file_path, pc_eigvec_file)

    lambdas <- init_lambdas_df()
    lambda_df_col_map <- list(
        "Fasting glucose", "Fasting_Glucose"
    )
    covariate_combs <- list(
        "Age, Gender" = FALSE,
        "Age, Gender + 10 PCs" = TRUE
    )
   
    for (covariates in names(covariate_combs)) {
        pc <- covariate_combs[[covariates]]
        on <- pc == TRUE ? "enabled" : "disabled"
        logger("INFO", "Performing GWAS with Principal Components ", on, ".")

        lambda_row <- data.frame(Covariates_Used = covariates)

        # Perform analysis for each of the phenotypes
        for (suffix in phenotype_suffixes) {
            logger("INFO", "Inspecting phenotype ", quotes(suffix), ".")
            pheno_path <- get_pheno_path(suffix)
            mpheno_args <- get_mpheno_args(suffix)
           
            if (pc) {
                covar_file <- covar_pc_file_path
            } else {
                covar_file <- covar_file_path
            }

            pheno_basename <- gwas_pheno(pheno_path, suffix, mpheno_args, covar_pc_file_path, pc)
            pheno_full_path <- get_pheno_analysis_full_path(pheno_basename, suffix)

            if (pc) {
                clump_path <- clumping(pheno_full_path)
                read_clumps(clump_path, suffix, pc)
            }

            d <- gwas_plots(pheno_full_path, pc, suffix)
            lambda <- compute_lambda(d, suffix, pc)
            lambda_delta <- lambda - 1.0
            
            trait_name <- get_trait_name(suffix)
            col_name <- lambda_df_col_map[[trait_name]]
            if (is.null(col_name)) {
                logger("ERROR", "Col name not found from trait: ", quotes(trait_name), ".")
            }

            row[[col_name]] <- lambda
            row[[col_name]] <- lambda_delta

            gc()
        }

        lambdas <- rbind(lambdas, row)
    }

    latex_col_align <- paste0("|l|", paste(rep("l:l", collapse = "|")), "|")
    caption <- "Genomic Inflation Values ($\\lambda$) obtained with different covariates"
    col_names <- c("Covariates Used", sapply(phenotype_suffixes, function(suffix) {
        trait_name <- get_trait_name(suffix)
        c(trait_name, paste0(trait_name, " \\Delta")))
    }

    out_name <- add_extension("lambdas.tex", exts$txt)
    latex_table(lambdas, out_name, latex_col_align, caption, col_names,
                digits = 3, line_spacing_factor = 1.0, hide_row_names = TRUE)
}

args <- commandArgs(trailingOnly = TRUE)

# Run all by default
perform_qc <- TRUE
perform_gwas <- TRUE

if (length(args) > 0) {
    if (!("qc" %in% args)) {
        perform_qc <- FALSE
    }

    if (!("gwas" %in% args)) {
        perform_gwas <- FALSE
    }
}

data_subset_path <- quality_control(perform_qc)
qc_data_path <- sample_qc(data_subset_path, perform_qc)

if (perform_gwas) {
    gwas(qc_data_path)
}

logger("DONE!")

