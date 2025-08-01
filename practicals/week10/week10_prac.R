# Week 10 Practical
# 02/10/2024
# Hugo Burton
# STAT3306

source("Practical_partF.R")

# ==== Globals ====

data_folder <- "/data/STAT3306/practical/PartF"
out_dir <- "./out"
plots_dir <- "./plots"

plot_w <- 430
plot_h <- 300

# ==== Functions ====

shell_call <- function(...) {
    system(..., ignore.stdout = TRUE)
    invisible(NULL)
}

mkdir_if_not_exist <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
}

construct_data_path <- function(basename) {
    file.path(data_folder, basename)
}

construct_out_path <- function(basename) {
    file.path(out_dir, basename)
}

gcta <- function(args, out_name) {
    full_path <- construct_out_path(out_name)
    full_args <- paste(args, "--out", full_path)
    command <- paste("gcta", full_args)
    cat(paste0("Executing '", command, "'."))
    shell_call(command)

    return(full_path)
}

file_exists <- function(path) {
    if (length(path) == 0) {
        return(FALSE)
    }

    exists <- file.exists(path)
    return(exists)
}

wrap_read_table <- function(path, header = TRUE, ...) {
    if (!file_exists(path)) {
        return(NULL)
    }

    read.table(path, header = header, ...)
}

wrap_plot <- function(plot_callback, data, out_name, ..., width = plot_w, height = plot_h) {
    out_path <- file.path(plots_dir, out_name)

    if (file_exists(out_path)) {
        delete_file(out_path)
    }

    png(out_path, width, height)
    plot_callback(data, ...)
    dev.off()

    return(out_path)
}

wrap_histogram <- function(...) {
    wrap_plot(hist, ...)
}

# === Init ===

mkdir_if_not_exist(out_dir)
mkdir_if_not_exist(plots_dir)

# ==== Main ====

grmBinPath <- construct_data_path("QIMRX_no_twin.grm")
cat("Data path ", grmBinPath)

grm <- read_GRMBin(grmBinPath)

preview <- 5

grm[1:preview, 1:preview]

qimrx_no_twin_path <- construct_data_path("QIMRX_no_twin")
ht_t_x_path <- construct_data_path("HT_T_X.pheno")

qimrx_1 <- gcta(paste("--grm", qimrx_no_twin_path, "--pheno", ht_t_x_path, "--mpheno 1", 
                      "--reml"), "QIMRX_1")

qimrx_2 <- gcta(paste("--grm", qimrx_no_twin_path, "--pheno", ht_t_x_path, "--mpheno 2",
                      "--reml"), "QIMRX_2")

qimrx_1 <- paste0(qimrx_1, ".hsq")
qimrx_2 <- paste0(qimrx_2, ".hsq")

hsq.1 <- wrap_read_table(qimrx_1, fill = TRUE)
hsq.2 <- wrap_read_table(qimrx_2, fill = TRUE)

head(hsq.1)
head(hsq.2)

# Unrelated individuals

names(grm) <- c("IND_1", "IND_2", "SNP_NUM", "REL")
dim(grm)

grm.diag <- diag(grm)
length(grm.diag)

head(grm.diag)

grm.off.diag <- grm[upper.tri(grm)]

wrap_histogram(grm.diag, "grm.diag.png", breaks = 2500, freq = FALSE,
    xlab = "GRM diagonals", xlim = c(0.95, 1.2), main = "GRM Diag")

wrap_histogram(grm.off.diag, "grm.off.diag.png", breaks = 2500,
               xlab = "GRM off-diagonals", main = "")

wrap_histogram(grm.off.diag[which(grm.off.diag > 0.1)], 
               "grm.off.diag.0.01.png", breaks = 200,
               freq  FALSE, xlab = "GRM off-diagonals",
               xlim = c(0.1, 1.1), main = "")
               

