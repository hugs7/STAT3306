# Week 6 Practical
# Date: 28-08-2024
# Author: Hugo Burton

data_folder <- "/data/STAT3306/"
prac_folder <- paste0(data_folder, "practical/PartD/")

fam <- read.table(paste0(prac_folder, "data.fam"))
bim <- read.table(paste0(prac_folder, "data.bim"))
pheno <- read.table(paste0(prac_folder, "data.pheno"))


#count number of individuals
dim(fam)
head(fam)

#count number of genotyped SNPs
dim(bim)
head(bim)


#Count the number of individuals with height and transferrin measurements
dim(pheno)
head(pheno)
#column 3 = height, column 4 = transferrin
length(which(!is.na(pheno[,3])))

#plot distributions- do they look suitable for a GWAS?


# Load missing data from plink
print("-------------")
plink_results_dir <- "./results/"
plink_missing_lmiss <- paste0(plink_results_dir, "missing.lmiss")
plink_missing_imiss <- paste0(plink_results_dir, "missing.imiss")

lmiss <- read.table(plink_missing_lmiss)
imiss <- read.table(plink_missing_imiss)

dim(lmiss)
dim(imiss)

head(lmiss)

png("lmiss.png", width=600, height=350)
hist(lmiss$"F_MISS")
dev.off()
