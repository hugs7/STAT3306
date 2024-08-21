
# Question 9

loci <- 15
indiv <- 200
h2 <- 0.5


af <- runif(loci)

x <- t(replicate(indiv, rbinom(loci, 2, af)))


# Check polymorphism
polymorphic = apply(x, 2, var) == 0
sum(polymorphic)

x <- scale(x)


beta <- rnorm(loci, 0, sqrt(h2/loci))
genetic <- x %*% beta


# Environmental is per person
environmental <- rnorm(indiv, 0, sqrt(1 - h2))

trait <- genetic + environmental

hist(trait, breaks = 20)