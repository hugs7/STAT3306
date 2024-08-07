### R code from vignette source 'Week1 - Matrix Reminder.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Week1 - Matrix Reminder.Rnw:65-71
###################################################
A = matrix(c(3, 1, 2, 5), nrow = 2, ncol = 2, byrow = TRUE)
B = matrix(c(1, 2, 0, 1), nrow = 2, ncol = 2, byrow = TRUE)
C = A + B
C
D = A - B
D


###################################################
### code chunk number 2: Week1 - Matrix Reminder.Rnw:87-94
###################################################
C = matrix(c(3, 1, 2, 2, 5, 4, 1, 1, 2), nrow = 3, ncol = 3, byrow = TRUE)

a = C[1, 1]
b = C[1, 2:3]
d = C[2:3, 1] # not as expected 
d = C[2:3, 1, drop = FALSE] # to keep the same "object/structure" as C
B = C[2:3, 2:3]


###################################################
### code chunk number 3: Week1 - Matrix Reminder.Rnw:106-110
###################################################
a = matrix(c(1, 2, 3, 4), ncol = 1)
b = c(4, 5, 7, 9)
a*b
sum(a*b)


###################################################
### code chunk number 4: Week1 - Matrix Reminder.Rnw:145-157
###################################################
A = matrix(c(3, 1, 2, 5), nrow = 2, ncol = 2, byrow = TRUE)
B = matrix(c(1, 2, 0, 1), nrow = 2, ncol = 2, byrow = TRUE)

A %*% B
B %*% A

C = matrix(c(3, 1, 2, 2, 5, 4, 1, 1, 2), nrow = 3, ncol = 3, byrow = TRUE)
D = matrix(c(2, 0, 5, 0, 4, 21), nrow = 3, ncol = 2)
# C
# D
C %*% D
# D %*% C # non-conformable arguments


###################################################
### code chunk number 5: Week1 - Matrix Reminder.Rnw:174-180 (eval = FALSE)
###################################################
## C = matrix(c(3, 1, 2, 2, 5, 4, 1, 1, 2), nrow = 3, ncol = 3, byrow = TRUE)
## D = matrix(c(2, 0, 5, 0, 4, 21), nrow = 3, ncol = 2)
## C
## t(C)
## D
## t(D)


###################################################
### code chunk number 6: Week1 - Matrix Reminder.Rnw:189-191 (eval = FALSE)
###################################################
## I = diag(3)
## I


###################################################
### code chunk number 7: Week1 - Matrix Reminder.Rnw:214-223
###################################################
B = matrix(c(1, 2, 0, 1), nrow = 2, ncol = 2, byrow = TRUE)
solve(B)

B %*% solve(B)

C = matrix(c(3, 1, 2, 2, 5, 4, 1, 1, 2), nrow = 3, ncol = 3, byrow = TRUE)
solve(C)

C %*% solve(C) # comment


###################################################
### code chunk number 8: Week1 - Matrix Reminder.Rnw:267-274
###################################################
B = matrix(c(1, 2, 0, 1), nrow = 2, ncol = 2, byrow = TRUE)
diag(B)
sum(diag(B))

C = matrix(c(3, 1, 2, 2, 5, 4, 1, 1, 2), nrow = 3, ncol = 3, byrow = TRUE)
diag(C)
sum(diag(C))


###################################################
### code chunk number 9: Week1 - Matrix Reminder.Rnw:330-337
###################################################
V = matrix(c(10, -5, 10, -5, 20, 0, 10, 0, 30), nrow = 3, ncol = 3, byrow = TRUE)
c1 = matrix(c(1, -2, 5), ncol = 1)
c2 = matrix(c(0, 6, -4), ncol = 1)

t(c1) %*% V %*% c1
t(c2) %*% V %*% c2
t(c1) %*% V %*% c2


###################################################
### code chunk number 10: Week1 - Matrix Reminder.Rnw:352-361
###################################################
# simulate 1000 random numbers from a normal distribution with mean 10 and variance 2
set.seed(123) # for reproducibility
num = rnorm(1000, mean = 10, sd = sqrt(2))
mean(num) # sample mean
var(num) # sample variance


# plot the density of the normal distribution
plot(density(num), main = "Estimation of the density of X")


###################################################
### code chunk number 11: Week1 - Matrix Reminder.Rnw:366-370
###################################################
#### OR
x   <- seq(5,15,length=1000)
y   <- dnorm(x,mean=10, sd=sqrt(2))
plot(x,y, type="l", lwd=1, main ="Normal distribution with mean 10 and variance 2")


###################################################
### code chunk number 12: Week1 - Matrix Reminder.Rnw:463-465
###################################################
# the end of that chapter. Command not to be run by students
Stangle('Week1 - Matrix Reminder.Rnw', encoding = 'utf8')


