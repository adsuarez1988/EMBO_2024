thetas <- seq(0,1,0.01)

like <- c()
for (theta in thetas) like <- c(like, (theta^3)*(1-theta)^1)

plot(thetas,like, type="l")

# MLE
cat("estimate: ", thetas[which.max(like)])

tD <- 2*log(max(like)/like[6])
tD
# p-value
1-pchisq(tD,1)

# skeptical belief function
belief1 <- 1/thetas
belief1[1] <- NA
belief1 <- belief1/sum(belief1, na.rm=T)

par(mfrow=c(2,3))
plot(thetas,like,type="l")
plot(thetas,belief1,type="l")
plot(thetas,belief1*like,type="l")
cat("point estimate:", thetas[which.max(belief1*like)])


# very skeptical belief function
belief2 <- 1/thetas^3
belief2[1] <- NA
belief2 <- belief2/sum(belief2, na.rm=T)

par(mfrow=c(2,3))
plot(thetas,like,type="l")
plot(thetas,belief2,type="l")
plot(thetas,belief2*like,type="l")
cat("point estimate:", thetas[which.max(belief2*like)])



mu <- 2
tau <- 1

x <- seq(-4,10,0.01)
plot(x=x, dnorm(x=x, mean=mu, sd=tau), ylim=c(0,0.6),
     type="l", lty=1, ylab="Density", xlab=expression(theta), main="")
legend(x="topleft", legend=c(expression(pi(theta)),
                             expression(f(y~"|"~theta)), expression(p(theta~"|"~y))), lty=1:3)

plot(x=x, dnorm(x=x, mean=mu, sd=tau), ylim=c(0,0.6),
     type="l", lty=1, ylab="Density", xlab=expression(theta), main="")
legend(x="topleft", legend=c(expression(pi(theta)),
                             expression(f(y~"|"~theta)), expression(p(theta~"|"~y))), lty=1:3) # prior

# likelihood
y <- 6
sigma <- 1
points(x=x, y=dnorm(x=y, mean=x, sd=sigma), type="l", lty=2)


plot(x=x, dnorm(x=x, mean=mu, sd=tau), ylim=c(0,0.6),
     type="l", lty=1, ylab="Density", xlab=expression(theta), main="")
legend(x="topleft", legend=c(expression(pi(theta)),
                             expression(f(y~"|"~theta)), expression(p(theta~"|"~y))), lty=1:3) # prior

# likelihood
points(x=x, y=dnorm(x=y, mean=x, sd=sigma), type="l", lty=2) # likelihood

# posterior
B <- sigma^2/(sigma^2+tau^2)
postMean <- B*mu + (1-B)*y
postVar <- B*tau^2
points(x=x, y=dnorm(x=x, mean=postMean, sd=sqrt(postVar)), type="l", lty=3)



#### Increasing the sd in the prior distribution #######

mu <- 2
tau <- 2

x <- seq(-4,10,0.01)
plot(x=x, dnorm(x=x, mean=mu, sd=tau), ylim=c(0,0.6),
     type="l", lty=1, ylab="Density", xlab=expression(theta), main="")
legend(x="topleft", legend=c(expression(pi(theta)),
                             expression(f(y~"|"~theta)), expression(p(theta~"|"~y))), lty=1:3)

plot(x=x, dnorm(x=x, mean=mu, sd=tau), ylim=c(0,0.6),
     type="l", lty=1, ylab="Density", xlab=expression(theta), main="")
legend(x="topleft", legend=c(expression(pi(theta)),
                             expression(f(y~"|"~theta)), expression(p(theta~"|"~y))), lty=1:3) # prior

# likelihood
y <- 6
sigma <- 1
points(x=x, y=dnorm(x=y, mean=x, sd=sigma), type="l", lty=2)


plot(x=x, dnorm(x=x, mean=mu, sd=tau), ylim=c(0,0.6),
     type="l", lty=1, ylab="Density", xlab=expression(theta), main="")
legend(x="topleft", legend=c(expression(pi(theta)),
                             expression(f(y~"|"~theta)), expression(p(theta~"|"~y))), lty=1:3) # prior

# likelihood
points(x=x, y=dnorm(x=y, mean=x, sd=sigma), type="l", lty=2) # likelihood

# posterior
B <- sigma^2/(sigma^2+tau^2)
postMean <- B*mu + (1-B)*y
postVar <- B*tau^2
points(x=x, y=dnorm(x=x, mean=postMean, sd=sqrt(postVar)), type="l", lty=3)

 ############## Afternoon #############


install.packages("abc")
library(abc)

sims <- read.csv("mosquito-task2.csv", head=T)


# check prior distributions
x11()
par(mfrow=c(2,2))
hist(sims$N1)
hist(sims$N2)
hist(sims$T_split)
hist(sims$MigRate)

# remove simulations with NaN for some summary stats!

# find useful summary stats which correlate with T_split
cor(sims$Fst, sims$T_split)
cor(sims$dxy, sims$T_split)
cor(sims$segsites1, sims$T_split)
cor(sims$segsites2, sims$T_split)
cor(sims$pi1, sims$T_split)
cor(sims$pi2, sims$T_split)
cor(sims$tajima1, sims$T_split)
cor(sims$tajima2, sims$T_split)

obs <- read.csv("mosquito-observed.csv", head=T)

# check if simulated retained summary stats contain the observed one
quantile(sims$Fst); cat(obs$Fst)
quantile(sims$segsites1); cat(obs$segsites1)
quantile(sims$segsites2); cat(obs$segsites2)

# merge obs with retained sims to scale them
sumstats <- scale(rbind(obs[c(1,3,4)],sims[,c(5,7,8)]))

est <- abc(target=sumstats[1,], param=sims$T_split, sumstat=sumstats[-1,], tol=0.05, method="rejection")

hist(est$dist)
abline(v=max(est$dist[which(est$region)]), lty=2)

# posterior
x11()
par(mfrow=c(2,1))
hist(est$unadj.values, freq=FALSE, xlim=range(sims$T_split), col=rgb(0,0,1,1/4), main="Posterior probability", xlab="Split time")

# MAP
map <- mean(est$unadj.values)
abline(v=map, lty=2)

# confidence intervals
hpd <- quantile(x=est$unadj.values, probs=c(0.025,0.975))
abline(v=hpd, lty=3)

# prior
hist(sims$T_split, freq=FALSE, xlim=range(sims$T_split), col=rgb(1,0,0,1/4))

