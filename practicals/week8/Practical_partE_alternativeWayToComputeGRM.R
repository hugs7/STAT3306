### R code from vignette source 'Sweave/PartD - BLUP.Rnw'
### Encoding: UTF-8

setwd("/Users/uqjzeng1/Work/Teach/STAT7306/2024/practicals/week\ 7")

###################################################
### code chunk number 1: code
###################################################
nmarkers <- 10      #number of markers
nrecords <- 325     #number of records
lambda    <- 10      #value for lambda


###################################################
### code chunk number 2: code
###################################################
X <- matrix(scan("PartE/xvec_day4.inp"), ncol = nmarkers, byrow = TRUE)
y <- matrix(scan("PartE/yvec_day4.inp"), byrow = TRUE)


###################################################
### code chunk number 3: code
###################################################
ones <- array(1, c(nrecords))
ident_mat <-diag(nmarkers)


###################################################
### code chunk number 4: code
###################################################
coeff <- array(0, c(nmarkers + 1, nmarkers + 1))
coeff[1:1, 1:1] <- t(ones) %*% ones
coeff[1:1, 2:(nmarkers+1)] <- t(ones) %*% X


###################################################
### code chunk number 5: code
###################################################
coeff[2: 2:(nmarkers+1), 1] <- t(X) %*% ones
coeff[2:(nmarkers+1), 2:(nmarkers+1)] <- t(X) %*% X + lambda * ident_mat

rhs = rbind(t(ones) %*% y, t(X) %*% y)


###################################################
### code chunk number 6: code (eval = FALSE)
###################################################
solution_vec <- solve(coeff, rhs)


###################################################
### code chunk number 7: code (eval = FALSE)
###################################################
Xprog <- matrix(scan("PartE/xvec_prog.inp"), ncol = nmarkers, byrow = TRUE)
GEBV = Xprog %*% solution_vec[-1]
  
yprog <- matrix(scan("PartE/yvec_prog.inp"), ncol = 1, byrow = TRUE)
cor(GEBV, yprog)^2


###################################################
### code chunk number 8: code (eval = FALSE)
###################################################
Xall <- rbind(X, Xprog)
nanims = nrow(Xall)


###################################################
### code chunk number 9: code (eval = FALSE)
###################################################
## 
## # p_j for all SNPS
p = colMeans(Xall)/2
## 


###################################################
### code chunk number 10: code (eval = FALSE)
###################################################
W = matrix(0, nrow = nanims, ncol = nmarkers)
for(j in 1:nmarkers)
    W[,j] <- (Xall[,j] - 2*p[j])/sqrt(2*p[j]*(1-p[j]))
## 


###################################################
### code chunk number 11: code (eval = FALSE)
###################################################
G = W%*%t(W)/nmarkers
# The next line adds a small amount to the diagonal of G, 
# otherwise G is not invertable in this small example!
G <- G + diag(nanims)*0.01    
Ginv <- solve(G)


###################################################
### code chunk number 12: code (eval = FALSE)
###################################################
Z1 <-diag(nrecords)
Z2 <-matrix(0, 325, 31)
Z <- cbind(Z1, Z2)


###################################################
### code chunk number 13: code (eval = FALSE)
###################################################
## # coeff
coeff <- array(0, c(nanims + 1, nanims + 1))
coeff[1:1, 1:1] <- t(ones) %*% ones
coeff[1:1, 2:(nanims+1)] <- t(ones) %*% Z
coeff[2: 2:(nanims+1), 1] <- t(Z) %*% ones
coeff[2:(nanims+1), 2:(nanims+1)] <- t(Z) %*% Z + Ginv
## 
rhs = rbind(t(ones) %*% y, t(Z) %*% y)


###################################################
### code chunk number 14: code (eval = FALSE)
###################################################
gblup <- solve(coeff, rhs)


###################################################
### code chunk number 15: code (eval = FALSE)
###################################################
## # the genomic prediction for the 31 selection candidates is
yprog_pred = gblup[-c(1:326)]


###################################################
### code chunk number 16: code (eval = FALSE)
###################################################
## # the accuracy is
cor(yprog, yprog_pred)^2

###################################################
### code chunk number 17: code (eval = FALSE)
###################################################
plot(GEBV,yprog_pred, xlab = "SNP-BLUP", ylab = "GBLUP")
lm(yprog_pred ~ GEBV)
abline(lm(yprog_pred ~ GEBV))


###################################################
### code chunk number 18: code (eval = FALSE)
###################################################
X = W[1:325,]
coeff <- array(0, c(nmarkers + 1, nmarkers + 1))
coeff[1:1, 1:1] <- t(ones) %*% ones
coeff[1:1, 2:(nmarkers+1)] <- t(ones) %*% X
coeff[2: 2:(nmarkers+1), 1] <- t(X) %*% ones
coeff[2:(nmarkers+1), 2:(nmarkers+1)] <- t(X) %*% X + lambda * ident_mat
rhs = rbind(t(ones) %*% y, t(X) %*% y)
solution_vec <- solve(coeff, rhs)
Wprog = W[-c(1:325),] #only W for the progeny
blup_W = Wprog %*% solution_vec[-1]
cor(yprog, blup_W)^2
plot(blup_W, yprog_pred, xlab = "SNP-BLUP with W", ylab = "GBLUP")
lm(yprog_pred ~ blup_W)
abline(lm(yprog_pred ~ blup_W))

###################################################
### code chunk number 19: PartD - BLUP.Rnw:243-245
###################################################
# the end of that chapter. Command not to be run by students
Stangle('Sweave/PartD - BLUP.Rnw', encoding = 'utf8', output = "Additional Files For Students/Rcode/PartD - BLUP.R")


