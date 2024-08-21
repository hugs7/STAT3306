### R code from vignette source 'Sweave/PartF - Quantitative Genetics of Disease.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: PartF - Quantitative Genetics of Disease.Rnw:19-23
###################################################
v <- 0.5                          # Proportion of cases in sample
p <- 0.2                          # Frequency of risk allele A
R <- 1.2                          # Heterozygote genotype relative risk; homozygote GRR is R^2
K <- 0.01                         # Overall risk (or probability) of disease in population


###################################################
### code chunk number 2: PartF - Quantitative Genetics of Disease.Rnw:28-31
###################################################
N <- 10000                        # Total sample size
alpha <- 5e-8                     # Level of significance that will be used in the association study
t <- qnorm(alpha/2,0,1)           # Normal distribution threshold for declaring significance


###################################################
### code chunk number 3: PartF - Quantitative Genetics of Disease.Rnw:38-40
###################################################
calculate_f0 <- function(p,R) {
    f0 <- p^2*R^2 + 2*p*(1-p)*R + (1-p)^2    # Probability of disease in those homozygote for non-risk allele aa
    return (f0)
}
f0 <- calculate_f0(p,R)

calculate_PG <- function(PG) {
    PG <- c((1-p)^2, 2*p*(1-p), p^2)  # Probability of genotypes aa, Aa, AA assumes HW equilibrium
    return (PG)
}
PG <- calculate_PG(PG)


###################################################
### code chunk number 4: PartF - Quantitative Genetics of Disease.Rnw:47-54
###################################################
calculate_p_case <- function(K) {
    PDgivG <- c(f0, f0*R, f0*R^2)     # Probability of disease given the genotype
    PDandG <- PDgivG*PG               # Probability of genotype and disease
    # PDandG <- c(f0*(1-p)^2,f0*R*2*p*(1-p),f0*R^2*p^2)
    sum(PDandG); K                    # Check equals K
    PGgivD <- PDandG/K                # Probability of genotype in people who are diseased -
                                    # compare PGgivD to PG
    pcase <- 0.5*PGgivD[2] + PGgivD[3]  # Frequency of allele A in cases
    return (pcase)
}
pcase <- calculate_p_case(K)


###################################################
### code chunk number 5: PartF - Quantitative Genetics of Disease.Rnw:61-69
###################################################
calculate_p <- function(K) {
    PNDgivG <- 1 - PDgivG               # Probability of being not diseased given the genotype
    PGandND <- PNDgivG * PG             # Probability of  genotype and not diseased 
    sum(PGandND); (1-K)               # Check equals 1-K
    PGgivND <- PGandND / (1-K)          # Probability of genotype in not diseased - compare PGgivND + PG
    pcont <- 0.5 * PGgivND[2] + PGgivND[3]# Frequency of allele A in controlsp
    #pcontchk = (p/(1-K))*(1-K*R/(1+p*(R-1)))
    # Checks
    p;K*pcase + (1-K)*pcont             # Check that population weighted average of pcase and pcontrol is p
    return (list(pcase=pcase,pcont=pcont))
}

p_res <- calculate_p(K)
p_res


###################################################
### code chunk number 6: PartF - Quantitative Genetics of Disease.Rnw:76-80
###################################################
pcase-pcont
OR <- (pcase/(1-pcase)) / (pcont/(1-pcont))
OR <- (pcase/pcont) / ((1-pcase)/(1-pcont)) # equivalent
OR


###################################################
### code chunk number 7: PartF - Quantitative Genetics of Disease.Rnw:85-92
###################################################
pbar <- v*pcase + (1-v)*pcont       # Mean allele frequency in case-control sample
NCP <- (pcase-pcont) * (pcase-pcont) * 2 * N * v * (1-v) / (pbar*(1-pbar)) 
                                  # Chi-square non-centrality parameter
pow <- pnorm(sqrt(NCP)+t)         # Power - normal distribution is sqrt of chi-square

# Check agrees with genetic power calculator: 
pow
NCP

N <- 10000
v <- 0.5
p <- 0.2
R <- 1.2
K <- 0.4
alpha <- 5e-8
Dprime <- 1


###################################################
### code chunk number 8: PartF - Quantitative Genetics of Disease.Rnw:108-121 (eval = FALSE)
###################################################
# N =10000,v = 0.5,p=0.2,R=1.2,K=0.01,alpha=5e-8
calculate_power <- function(pcont) {
    pbar <- v*pcase+(1-v)*pcont       # Mean allele frequency in case-control sample
    #NCP
    NCP <- (pcase-pcont)*(pcase-pcont)*2*N*v*(1-v)/(pbar*(1-pbar)) # Chi-square non-centrality parameter
    pow <- pnorm(sqrt(NCP)+t)         # Power - normal distribution is sqrt of chi-square
    return (list(pow=pow,NCP=NCP))
}

screened <- calculate_power(pcont)
unscreened <- calculate_power(p)

screened_power <- screened$pow
unscreened_power <- unscreened$pow
screened_NCP <- screened$NCP
unscreened_NCP <- unscreened$NCP
screened_power
unscreened_power
screened_NCP
unscreened_NCP

#K=0.01
#screened pow  = ??? ; NCP =???
#unscreened pow = ???; NCP =???
## 
## # the code above assumes controls are screened, for unscreened controls pcont=p


###################################################
### code chunk number 9: PartF - Quantitative Genetics of Disease.Rnw:130-133 (eval = FALSE)
###################################################
## K=0.15 # rerun the relevant code
## #screened pow  = ??? ; NCP =???
## #unscreened pow = ???; NCP =???


###################################################
### code chunk number 10: PartF - Quantitative Genetics of Disease.Rnw:141-176
###################################################
getpower=function(p,GRR,Ncase,Ncont,A,K){
  N=Ncase+Ncont
  v<-Ncase/N
  Thr<-qnorm(A/2,0,1)
  pcase<-p*GRR/(p*GRR+1-p)
  pcont<-(p/(1-K))*(1-K*GRR/(1+p*(GRR-1)))
  pdiff=pcase-pcont
  pbar<-pcase*v+pcont*(1-v)
  #exactly the same as GPC
  NCP<-pdiff*pdiff*2*N*v*(1-v)/(pbar*(1-pbar))
  pow=pnorm(sqrt(NCP)+Thr)
  pbar=0.5*(pcase+pcont)
  Nequiv=NCP*pbar*(1-pbar)/pdiff^2   # sample size of equal cases and controls that gives same power
  result=list(pow=pow,NCP=NCP,Nequiv=Nequiv)
  return (result)
}
res <- getpower(0.2,1.2,5000,5000,5e-8,0.01) #confirm same as programmed above
res
#-------------------------
# run this section to get a power graph
#------------------------
opar<-par(mfrow=c(1,1))
Ncase=5000
Ncont=5000
maxGRR=1.5
A<-5e-8
x<-c(1,maxGRR)
y<-c(0,1)
plot(x,y,xlim=c(1,maxGRR),ylim=c(0,1),xlab="Genotype Relative Risk",ylab="power",
     main="Power for sample size of",type="n",cex.lab=2,cex.main=1.5)
mtext(paste(Ncase,"cases and",Ncont,"controls"),3,0,cex =1.5)
curve(getpower(0.1,x,Ncase,Ncont,A,K)$pow,from=1,to=maxGRR,col=1,lty=1,lwd=2,add=TRUE)
curve(getpower(0.3,x,Ncase,Ncont,A,K)$pow,from=1,to=maxGRR,col=2,lty=1,lwd=2,add=TRUE)
curve(getpower(0.5,x,Ncase,Ncont,A,K)$pow,from=1,to=maxGRR,col=3,lty=1,lwd=2,add=TRUE)
legend(x=1.3,y=0.2,col=c(1:3),lty=1,legend=paste("RAF =",c(0.1,0.3,0.5)),box.col="white")


###################################################
### code chunk number 11: PartF - Quantitative Genetics of Disease.Rnw:206-224
###################################################
set.seed(615) # for reproducibility
m = 500                                     # number of SNPs
maf = runif(m, 0, .5)                       # random MAF for each SNP

N = 400                                     # number of individuals
x012 = t(replicate(N, rbinom(m, 2, maf)))   # n x m genotype matrix
# 
# scaling the data
scale_x012 = scale(x012)
scale_x012[which(is.na(scale_x012))] = 0
# 
h2 = 0.9  # heretability
beta = rnorm(m, 0, sqrt(h2/m))
g_parent = scale_x012 %*% beta
e_parent = rnorm(N, 0, sqrt(1-h2))
# 
# phenotype parent
y_parent = g_parent + e_parent


###################################################
### code chunk number 12: PartF - Quantitative Genetics of Disease.Rnw:229-234
###################################################
K = 0.2
t = qnorm(1-K, 0, 1)
t

D_parent = rep(0, N)
D_parent[y_parent>t] = 1 # P: phenotype
sum(D_parent==1)/length(D_parent)

###################################################
### code chunk number 13: PartF - Quantitative Genetics of Disease.Rnw:242-247 (eval = FALSE)
###################################################
## mean(y_parent[D_parent==0])
## mean(y_parent[D_parent==1])
## 
## i = dnorm(t)/K #i = z/K
## 


###################################################
### code chunk number 14: PartF - Quantitative Genetics of Disease.Rnw:259-294
###################################################

# if dad is AA and mom is AA => dad gives A, mom gives A =>kid is AA
# if dad is aa and mom is aa => dad gives a, mom gives a =>kid is aa
# if dad (or mom) is Aa, then probability 0.5 of giving a


# Given a parent individual, get one of their Alleles in the gamete
getGamete <- function(indiv) {
  indiv_gives = indiv # creates a vector of 0, or 1 . 0:a  and 1:A
    
  ind_1 = which(indiv_gives == 1) #we isolate the 1s (Aa) in the parent
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
### code chunk number 15: PartF - Quantitative Genetics of Disease.Rnw:299-324
###################################################
fullsib = NULL # where we will store the offsprings
for(i in 1:(N/2)){
  # randomly pick two parents
  parents = sample (1:N,2)
  # generate 2 offsprings from these two parents
  kid_1 = getOffspring(x012[parents[1], ], x012[parents[2], ])
  kid_2 = getOffspring(x012[parents[1], ], x012[parents[2], ])
  
  # store the new offspring with the other ones
  fullsib = rbind(fullsib, rbind(kid_1,kid_2))  
}
#
# scaling the data
scale_fullsib = scale(fullsib)
scale_fullsib[which(is.na(scale_fullsib))] = 0
# 
# phenotype y_fullsib
g_fullsib = scale_fullsib %*% beta
e_fullsib = rnorm(N, 0, sqrt(1-h2))
# 
y_fullsib = g_fullsib + e_fullsib
# 
# liability of the disease
D_fullsib = rep(0, N)
D_fullsib[y_fullsib>t] = 1 # P: phenotype


###################################################
### code chunk number 16: PartF - Quantitative Genetics of Disease.Rnw:328-352
###################################################
MZ = NULL # where we will store the offsprings
for(i in 1:(N/2)){
  # randomly pick two parents
  parents = sample (1:N,2)
  # generate twice the same offspring from these two parents
  kid = getOffspring(x012[parents[1], ], x012[parents[2], ])

  # store the new offspring with the other ones
  MZ = rbind(MZ, rbind(kid,kid))  
}
# 
# scaling the data
scale_MZ = scale(MZ)
scale_MZ[which(is.na(scale_MZ))] = 0
#  
# phenotype y_MZ
g_MZ = scale_MZ %*% beta
e_MZ = rnorm(N, 0, sqrt(1-h2))
#  
y_MZ = g_MZ + e_MZ
#  
# liability of the disease
D_MZ = rep(0, N)
D_MZ[y_MZ>t] = 1 # P: phenotype


###################################################
### code chunk number 17: PartF - Quantitative Genetics of Disease.Rnw:358-378
###################################################
# disease state of randomly selected pair of parents
parents = matrix(D_parent[sample(1:N,N, replace=TRUE)],ncol=2) 
#

# parent1 affected
affected_parent1 = which(parents[,1] == 1)
# How many parents are affected when parent1 is
a = apply(parents[affected_parent1,],1,sum)
# a==2: the second parent is affected when parent1 is
risk_affected_parent1 = sum(a==2)/length(affected_parent1)


# parent1 non-affected
non_affected_parent1 = which(parents[,1] == 0)
# how many parents are affected when parent1 is not
a = apply(parents[non_affected_parent1,],1,sum)
# a==1: the second parent is affected when parent1 is not
risk_non_affected_parent1 = sum(a==1)/length(non_affected_parent1) 

risk_affected_parent1/risk_non_affected_parent1


###################################################
### code chunk number 18: PartF - Quantitative Genetics of Disease.Rnw:384-404
###################################################
per_parent = matrix(D_fullsib,ncol=2,byrow=T) # each siblings are in row
# 

# sib1 affected
affected_sib1 = which(per_parent[,1] == 1)
# how many siblings are affected when sib1 is
a = apply(per_parent[affected_sib1,],1,sum) 

# a==2: the second sibling is affected when sib1 is
risk_affected_sib1 = sum(a==2)/length(affected_sib1)


# sib1 non-affected
non_affected_sib1 = which(per_parent[,1] == 0)
# how many siblings are affected when sib1 is not
a = apply(per_parent[non_affected_sib1,],1,sum) 
# a==1: the second sibling is affected when sib1 is not
risk_non_affected_sib1 = sum(a==1)/length(non_affected_sib1) 

risk_affected_sib1/risk_non_affected_sib1


###################################################
### code chunk number 19: PartF - Quantitative Genetics of Disease.Rnw:410-429
###################################################
per_parent = matrix(D_MZ,ncol=2,byrow=T) # each sibling are in row
#

# twin1 affected
affected_twin1 = which(per_parent[,1] == 1)
# how many twins are affected when twin1 is
a = apply(per_parent[affected_twin1,],1,sum) 
# a==2: the second twin is affected when twin1 is
risk_affected_twin1 = sum(a==2)/length(affected_twin1) 


# twin1 non-affected
non_affected_twin1 = which(per_parent[,1] == 0)
# how many twins are affected when twin1 is not
a = apply(per_parent[non_affected_twin1,],1,sum) 
# a==1: the second twin is affected when twin1 is not
risk_non_affected_twin1 = sum(a==1)/length(non_affected_twin1) 

risk_affected_twin1/risk_non_affected_twin1


###################################################
### code chunk number 20: PartF - Quantitative Genetics of Disease.Rnw:439-441
###################################################
# the end of that chapter. Command not to be run by students
Stangle('Sweave/PartF - Quantitative Genetics of Disease.Rnw', encoding = 'utf8', output = "Additional Files For Students/Rcode/PartF - Quantitative Genetics of Disease.R")


