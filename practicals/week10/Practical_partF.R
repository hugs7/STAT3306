### R code from vignette source 'Sweave/PartE - REML.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: PartE - REML.Rnw:105-143
###################################################
read_GRMBin <- function(prefix, size = 4){
    sum_i <- function(i){
        return(sum(1:i))
    }
    
    ## open file connections and read in data
    BinFileName <- paste(prefix,".bin",sep="")
    NFileName <- paste(prefix,".N.bin",sep="")
    IDFileName <- paste(prefix,".id",sep="")
    id <- read.table(IDFileName)
    n <- dim(id)[1]
    BinFile <- file(BinFileName, "rb")
    grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
    NFile <- file(NFileName, "rb")
    
    ## read in the number of SNPs used to calculate the GRM (does not appear to work)
    N <- readBin(NFile, n=1, what=numeric(0), size=size)
    i <- sapply(1:n, sum_i)
    
    ## clean up file connections
    close(BinFile)
    close(NFile)
    
    ## pull apart diagonal and lower triagular elements
    diag.elem <- grm[i]
    off.diag.elem <- grm[-i]
    
    ## create the full symmetric correlation matrix
    X <- diag(diag.elem)
    X[ upper.tri(X, diag = FALSE) ] <- off.diag.elem
    X <- X + t(X) - diag(diag(X))
    
    ## add sample IDs to rownames and colnames
    rownames(X) <- id$V2
    colnames(X) <- id$V2
    
    return(X)
}


###################################################
### code chunk number 2: PartE - REML.Rnw:146-150
###################################################
# Read in the gzipped GRM file
# will need to use readGRM functions
grm <- read_GRMBin("Data/PartF/QIMRX_no_twin.grm")
grm[1:5,1:5]


###################################################
### code chunk number 3: PartE - REML.Rnw:180-185
###################################################
# Read in GREML result files

hsq.1 <- read.table("Results/QIMRX_1.hsq", header = TRUE, fill = TRUE)
hsq.2 <- read.table("Results/QIMRX_2.hsq", header = TRUE, fill = TRUE)
head(hsq.1)


###################################################
### code chunk number 4: PartE - REML.Rnw:235-248
###################################################
# Name the columns of the GRM
names(grm) <- c("IND_1", "IND_2", "SNP_NUM", "REL")
dim(grm)
# Take out the diagonal elements
grm.diag <- diag(grm)
length(grm.diag)
head(grm.diag)

# Take out the GRM off-diagonal elements
grm.off.diag <- grm[upper.tri(grm)]
# Make a histogram of the diagonals
hist(grm.diag, breaks = 2500, freq = F, 
     xlab = "GRM diagonals", xlim = c(0.95, 1.2), main = "")


###################################################
### code chunk number 5: PartE - REML.Rnw:251-258
###################################################
# Make a histogram of the GRM off-diagonal relatedness estimates
par(mfrow = c(1, 1))
hist(grm.off.diag, breaks = 2500, freq = F, 
      xlab = "GRM off-diagonals", main = "")
hist(grm.off.diag[which(grm.off.diag > 0.1)], 
     breaks = 200, freq = F, 
     xlab = "GRM off-diagonals", xlim = c(0.1, 1.1), main = "") 


###################################################
### code chunk number 6: PartE - REML.Rnw:325-336
###################################################
# Read in GREML result files without relatedness. 
grm.nr <- read_GRMBin("Results/QIMRX_nr.grm")
# Create the same figures as above
# Note that these are only example file extensions
# and may change depending on what you want to call the files
hsq.1.nr <- read.table("Results/QIMRX_1_nr.hsq", 
                       header = TRUE, fill = TRUE)
hsq.2.nr <- read.table("Results/QIMRX_2_nr.hsq", 
                       header = TRUE, fill = TRUE)
hsq.1.nr
hsq.2.nr


###################################################
### code chunk number 7: PartE - REML.Rnw:410-412
###################################################
# the end of that chapter. Command not to be run by students
Stangle('Sweave/PartE - REML.Rnw', encoding = 'utf8', output = "Additional Files For Students/Rcode/PartE - REML.R")


