read_GRMBin <- function(prefix, size = 4) {
    #' Reads a GRM binary file
    #' @param prefix {string}: File basepath.
    #' @param size {integer}: BLock size in the binary file.
    #'                        Defaults to 4.
    #' @return {data.frame}: Data frame from the binary file.

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

