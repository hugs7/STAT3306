# Week 6 Practical
# Hugo Burton
# 28-08-2024

# === Globals ===
data_folder <- "/data/STAT3306/practical/PartD/data"
plink_results_dir <- "./results/"

# === Functions ===
run_plink <- function(plink_args, out_name) {
	system(paste0("plink --bfile ", data_folder, " ", plink_args, " --out ", plink_results_dir, out_name), ignore.stdout = T)
}

run_plink("--missing", "missing")

# === Missing ===

missing_path <- paste0(plink_results_dir, "missing.lmiss")
missing <- read.table(missing_path, header=T)
	
print(dim(missing))
print(head(missing))

png("missing_lmiss.png", width=600, height=350)
hist(missing$"F_MISS")
dev.off()

# === Het ===
run_plink("--het", "het")

het_path <- paste0(plink_results_dir, "het.het")
het <- read.table(het_path, header=T)
	
print(dim(het))
print(head(het))


png("het.png", width=600, height=350)
hist(het$"F")
dev.off() 

plot(abs(het$"F"))
abline(h=0.05,col="red")

# Sex

run_plink("--check-sex", "check_sex")

sex <- read.table(paste0(plink_results_dir, "check_sex.sexcheck"), header=T)
summary(sex)
which(sex[,"STATUS"]=="PROBLEM")
