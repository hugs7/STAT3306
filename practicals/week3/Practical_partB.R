### R code from vignette source 'Sweave/PartB - Quantitative Genetics.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: PartB - Quantitative Genetics.Rnw:31-46
###################################################
data_path = "practicals\\week3\\data\\"
load(paste(data_path,"1000Genomes-MHC.RData", sep=""))
dim(genotypes)
dim(pops)

#look at the data
#genotypes[1:5,1:10]

# Remove first 5 rows as this tells you information about the chromosones (aleles)
geno = genotypes[,-c(1:5)]
rownames(geno) = genotypes$SNP

geno=t(geno) # to put the individuals in rows
#geno[1:5,1:5]

n = ncol(geno)
p = nrow(geno)


###################################################
### code chunk number 2: PartB - Quantitative Genetics.Rnw:53-71 (eval = FALSE)
###################################################
## # observed genotypes
## 
## # q = maf = mean /2
## # p = 1-q
## 
## # calculation expected values if hardy weinberg equilibrium
## # aa = q*q = maf^2
## # AA = p*p = (1-maf)^2
## # Aa = 2*p*q
## 
## # compared to observed values:
## # aa 
## # AA
## # Aa
## 
## 
## # chi square of the actual numbers (not pvalues, but numbers = pvalue*n)
## # (O-E)^2/E


###################################################
### code chunk number 3: PartB - Quantitative Genetics.Rnw:75-79
###################################################
my_pop = "CLM"
my_SNP = geno[which(pops$Population.code == my_pop),]
my_p = ncol(my_SNP) #number of SNPs
my_n = nrow(my_SNP) #number of samples


###################################################
### code chunk number 4: PartB - Quantitative Genetics.Rnw:83-84
###################################################
# Diploid inviduals, so we expect 2 alleles per individual
af_est = colMeans(my_SNP)/2


###################################################
### code chunk number 5: PartB - Quantitative Genetics.Rnw:91-92
###################################################
hist(af_est)


###################################################
### code chunk number 6: code
###################################################
SNP_j = my_SNP[,162]
table(SNP_j) #observed counts


###################################################
### code chunk number 7: code
###################################################
pa=af_est[162]
pa


###################################################
### code chunk number 8: code
###################################################
library(xtable)
Observed = table(SNP_j)
Expected = c(38.30,43.40,12.30)
T=rbind(Observed,Expected)
print(xtable(T,digits=3,caption="Observed and expected counts", label="contOE"))



###################################################
### code chunk number 9: code
###################################################
chi2 = (43-38.30)^2/38.30+(34-43.40)^2/43.40+(17-12.300)^2/12.300
chi2


###################################################
### code chunk number 10: code
###################################################
1-pchisq(chi2, df =  1)


###################################################
### code chunk number 11: PartB - Quantitative Genetics.Rnw:167-171
###################################################
# expected counts if at the Hardy-Weimber equilibrium
E_n_Aa = 2 * af_est * (1-af_est) *my_n # 2p(1-p) *N= 2pq *N
E_n_AA = (1-af_est)^2 *my_n             # p^2 *N
E_n_aa = af_est^2 *my_n                 # q^2 *N


###################################################
### code chunk number 12: PartB - Quantitative Genetics.Rnw:174-178
###################################################
# observed counts
O_n_Aa = apply(my_SNP,2,function(x){sum(x == 1)}) # 2pq *N
O_n_AA = apply(my_SNP,2,function(x){sum(x == 0)}) # p^2 *N
O_n_aa = apply(my_SNP,2,function(x){sum(x == 2)}) # q^2 *N


###################################################
### code chunk number 13: PartB - Quantitative Genetics.Rnw:184-189
###################################################
chi2 = 
  (O_n_AA - E_n_AA)^2/(E_n_AA) + 
  (O_n_Aa - E_n_Aa)^2/(E_n_Aa) + 
  (O_n_aa- E_n_aa)^2/(E_n_aa)
hist(chi2)


###################################################
### code chunk number 14: PartB - Quantitative Genetics.Rnw:193-194
###################################################
sum(chi2 > qchisq(0.99,1), na.rm = T)


###################################################
### code chunk number 15: PartB - Quantitative Genetics.Rnw:198-201
###################################################
pvalues=1-pchisq(chi2,1)
sum(pvalues<0.01,na.rm=TRUE)
sum(pvalues<1e-6,na.rm=TRUE)


###################################################
### code chunk number 16: PartB - Quantitative Genetics.Rnw:229-232
###################################################
set.seed(6155) # for reproducibility
m = 500                                 # number of SNPs
maf = runif(m, 0, .5)                    # random MAF for each SNP


###################################################
### code chunk number 17: PartB - Quantitative Genetics.Rnw:236-237
###################################################
x012 = rbinom(m, 2, maf)


###################################################
### code chunk number 18: PartB - Quantitative Genetics.Rnw:241-243
###################################################
N = 400                                    # number of individuals
x012 = t(replicate(N, rbinom(m, 2, maf)))  # n x m genotype matrix


###################################################
### code chunk number 19: PartB - Quantitative Genetics.Rnw:248-250
###################################################
polymorphic = apply(x012, 2, var) == 0
sum(polymorphic)


###################################################
### code chunk number 20: PartB - Quantitative Genetics.Rnw:254-255
###################################################
x012[1:5, 1:10]


###################################################
### code chunk number 21: PartB - Quantitative Genetics.Rnw:264-299
###################################################

# if dad is AA and mom is AA => dad gives A, mom gives A =>kid is AA
# if dad is aa and mom is aa => dad gives a, mom gives a =>kid is aa
# if dad (or mom) is Aa, then probability 0.5 of giving a


# Given a parent individual, get one of their Alleles in the gamete
getGamete <- function(indiv) {
  indiv_gives = indiv/2 # creates a vector of 0, or 1 . 0:a  and 1:A
    
  ind_1 = which(indiv_gives == 0.5) #we isolate the 1s (Aa) in the parent
  #if Parent is Aa, the gamete is one binomial trial with probability 0.5
  if (length(ind_1) > 0)
    indiv_gives[ind_1] = rbinom(length(ind_1), size=1, 0.5) 
    # gives 0 for a or 1 for A
      
  return(indiv_gives) # vector of 0s and 1s
}
  
# Two parental Gametes combine to form a zygote
combineGametes <- function(momGam, dadGam) {
  kid = (momGam + dadGam) 
  # 0(a) + 0(a) = 0 (aa)
  # 0(a) + 1(A) = 1 (aA)
  # 1(A) + 1(A) = 2 (AA)

  kid
}
  
# Given two individuals, get an offspring for next generation
getOffspring <- function(mom, dad) {
  momGam <- getGamete(mom)
  dadGam <- getGamete(dad)
  return(combineGametes(momGam, dadGam))  
}


###################################################
### code chunk number 22: PartB - Quantitative Genetics.Rnw:303-304
###################################################
kid012 = getOffspring(x012[1, ], x012[2, ])


###################################################
### code chunk number 23: PartB - Quantitative Genetics.Rnw:308-322
###################################################
kid012 = NULL # where we will store the offsprings
pedigree = matrix(0, nrow = N, ncol = N) # we want to record who are the parents
for(i in 1:N){
  # randomly pick two parents
  parents = sample (1:N,2)
  # generate offspring from these two parents
  kid = getOffspring(x012[parents[1], ], x012[parents[2], ])
  
  # 
  pedigree[i, parents] = 1
  
  # store the new offspring with the other ones
  kid012 = rbind(kid012, kid)  
}


###################################################
### code chunk number 24: PartB - Quantitative Genetics.Rnw:327-334
###################################################
kid012[1:5, 1:10]

# who are the parents of offspring 241?
which(pedigree[241,] == 1)

# check that each offspring has only 2 parents
sum(apply(pedigree,1,sum)==2) # should be N


###################################################
### code chunk number 25: PartB - Quantitative Genetics.Rnw:344-349
###################################################
scale_x012 = scale(x012)
scale_x012[which(is.na(scale_x012))] = 0

scale_kid012 = scale(kid012)
scale_kid012[which(is.na(scale_kid012))] = 0


###################################################
### code chunk number 26: PartB - Quantitative Genetics.Rnw:353-362
###################################################
h2 = 0.7
beta = rnorm(m, 0, sqrt(h2/m))
g_parent = scale_x012 %*% beta
e_parent = rnorm(N, 0, sqrt(1-h2))

# phenotype parent
y_parent = g_parent + e_parent

var(y_parent) # should be around 1


###################################################
### code chunk number 27: PartB - Quantitative Genetics.Rnw:366-373
###################################################
g_kid = scale_kid012 %*% beta
e_kid = rnorm(N, 0, sqrt(1-h2))

# phenotype for kid
y_kid = g_kid + e_kid

var(y_kid) # should be around 1


###################################################
### code chunk number 28: PartB - Quantitative Genetics.Rnw:382-388
###################################################
y_midparent = vector(length=N)
for (i in 1:N){
  parents = which(pedigree[i,] == 1)
  y_midparent[i] = mean(y_parent[parents])
}



###################################################
### code chunk number 29: PartB - Quantitative Genetics.Rnw:391-396
###################################################
plot(y_midparent, y_kid)
reg = lm(y_kid~y_midparent)
reg
abline(reg)
abline(0,h2, col="red")


###################################################
### code chunk number 30: PartB - Quantitative Genetics.Rnw:440-443
###################################################
twin.data <- read.table(paste(data_path,"twin_height_bmi.txt", sep=""), header = TRUE)
dim(twin.data)
head(twin.data)


###################################################
### code chunk number 31: PartB - Quantitative Genetics.Rnw:448-451
###################################################
table(twin.data$twin)
MZ = which(twin.data$twin == "MZ")
DZ = which(twin.data$twin == "DZ")


###################################################
### code chunk number 32: PartB - Quantitative Genetics.Rnw:455-457
###################################################
cor(twin.data$ht_t1[MZ],twin.data$ht_t2[MZ])
cor(twin.data$ht_t1[DZ],twin.data$ht_t2[DZ])


###################################################
### code chunk number 33: PartB - Quantitative Genetics.Rnw:461-463
###################################################
cor(twin.data$bmi_t1[MZ],twin.data$bmi_t2[MZ])
cor(twin.data$bmi_t1[DZ],twin.data$bmi_t2[DZ])


###################################################
### code chunk number 34: PartB - Quantitative Genetics.Rnw:473-479
###################################################
cor(twin.data$ht_t1[MZ][twin.data$sex[MZ] == 1],twin.data$ht_t2[MZ][twin.data$sex[MZ] == 1])
cor(twin.data$ht_t1[MZ][twin.data$sex[MZ] == 2],twin.data$ht_t2[MZ][twin.data$sex[MZ] == 2])


cor(twin.data$ht_t1[DZ][twin.data$sex[DZ] == 1],twin.data$ht_t2[DZ][twin.data$sex[DZ] == 1])
cor(twin.data$ht_t1[DZ][twin.data$sex[DZ] == 2],twin.data$ht_t2[DZ][twin.data$sex[DZ] == 2])


###################################################
### code chunk number 35: PartB - Quantitative Genetics.Rnw:484-493
###################################################
plot(density(c(twin.data$ht_t1[MZ][twin.data$sex[MZ] == 1],twin.data$ht_t2[MZ][twin.data$sex[MZ] == 1]), 
      adjust = 3), col = "Blue", lwd = 2, ylim = c(0, 0.06), main = "Density of Height")

lines(density(c(twin.data$ht_t1[MZ][twin.data$sex[MZ] == 2],twin.data$ht_t2[MZ][twin.data$sex[MZ] == 2]), 
      adjust = 3), col = "Pink", lwd = 2)
lines(density(c(twin.data$ht_t1[MZ], twin.data$ht_t2[MZ]), 
      adjust = 3), col = "Green", lwd = 2)

legend("topright",legend = c("Male", "Female", "Both"), col = c("Blue", "Pink", "Green"), pch=15)


###################################################
### code chunk number 36: PartB - Quantitative Genetics.Rnw:500-501
###################################################
2*(cor(twin.data$ht_t1[MZ],twin.data$ht_t2[MZ]) - cor(twin.data$ht_t1[DZ],twin.data$ht_t2[DZ]))


###################################################
### code chunk number 37: PartB - Quantitative Genetics.Rnw:508-513
###################################################
reg.MZ_ht = lm(twin.data$ht_t2[MZ] ~ twin.data$sex[MZ] + twin.data$ht_t1[MZ])
reg.DZ_ht = lm(twin.data$ht_t2[DZ] ~ twin.data$sex[DZ] + twin.data$ht_t1[DZ])

A_ht = 2*(coef(reg.MZ_ht)[3] - coef(reg.DZ_ht)[3])
A_ht


###################################################
### code chunk number 38: PartB - Quantitative Genetics.Rnw:521-523
###################################################
C_ht = coef(reg.MZ_ht)[3] - A_ht
E_ht = 1 - (A_ht + C_ht)


###################################################
### code chunk number 39: PartB - Quantitative Genetics.Rnw:534-540 (eval = FALSE)
###################################################
## cor(twin.data$bmi_t1[MZ][twin.data$sex[MZ] == 1],twin.data$bmi_t2[MZ][twin.data$sex[MZ] == 1])
## cor(twin.data$bmi_t1[MZ][twin.data$sex[MZ] == 2],twin.data$bmi_t2[MZ][twin.data$sex[MZ] == 2])
## 
## 
## cor(twin.data$bmi_t1[DZ][twin.data$sex[DZ] == 1],twin.data$bmi_t2[DZ][twin.data$sex[DZ] == 1])
## cor(twin.data$bmi_t1[DZ][twin.data$sex[DZ] == 2],twin.data$bmi_t2[DZ][twin.data$sex[DZ] == 2])


###################################################
### code chunk number 40: PartB - Quantitative Genetics.Rnw:547-548 (eval = FALSE)
###################################################
## 2*(cor(twin.data$bmi_t1[MZ],twin.data$bmi_t2[MZ]) - cor(twin.data$bmi_t1[DZ],twin.data$bmi_t2[DZ]))


###################################################
### code chunk number 41: PartB - Quantitative Genetics.Rnw:555-560 (eval = FALSE)
###################################################
## reg.MZ_bmi = lm(twin.data$bmi_t2[MZ] ~ twin.data$sex[MZ] + twin.data$bmi_t1[MZ])
## reg.DZ_bmi = lm(twin.data$bmi_t2[DZ] ~ twin.data$sex[DZ] + twin.data$bmi_t1[DZ])
## 
## A_bmi = 2*(coef(reg.MZ_bmi)[3] - coef(reg.DZ_bmi)[3])
## A_bmi


###################################################
### code chunk number 42: PartB - Quantitative Genetics.Rnw:569-571 (eval = FALSE)
###################################################
## C_bmi = coef(reg.MZ_bmi)[3] - A_bmi
## E_ht = 1 - (A_bmi + C_bmi)


###################################################
### code chunk number 43: PartB - Quantitative Genetics.Rnw:576-578
###################################################
# the end of that chapter. Command not to be run by students
# Stangle('Sweave/PartB - Quantitative Genetics.Rnw', encoding = 'utf8', output = "Additional Files For Students/Rcode/PartB - Quantitative Genetics.R")


