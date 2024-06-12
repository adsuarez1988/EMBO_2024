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

 ############## PRACTICAL 1 #############


