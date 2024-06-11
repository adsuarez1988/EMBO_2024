N <- 50
fA <- 0.5

rbinom(1, 2*N, fA)/(2*N)
rbinom(1, 2*N, 0.53)/(2*N)


#### Different initial allele frequencies ######
# repeat 20 times each value, compare how many times the allele is fixed or extint 
Nexp <- 100
fexp <- 0.6
rbinom(1, 2*Nexp, fexp)/(2*N)
mean(0.62, 0.55, 0.59, 0.645, 0.56, 0.615)

f_0.6 <- rep(NA, 100)
f_0.6[1] <- 0.6
for (t in 1:99) f_0.6[t+1] <- rbinom(1, 2*N, f_0.6[t]) / (2*N)
 plot(x=1:100, y=f_0.6, type="l", ylim=c(0,1), lwd=2)

 
f_0.4 <- rep(NA, 100)
f_0.4[1] <- 0.4
for (t in 1:99) f_0.4[t+1] <- rbinom(1, 2*N, f_0.4[t]) / (2*N)
plot(x=1:100, y=f_0.4, type="l", ylim=c(0,1), lwd=2)


f_0.2 <- rep(NA, 100)
f_0.2[1] <- 0.2
for (t in 1:99) f_0.2[t+1] <- rbinom(1, 2*N, f_0.2[t]) / (2*N)
plot(x=1:100, y=f_0.2, type="l", ylim=c(0,1), lwd=2)


#### Changing population size #####

# initial population of 50
N=50
f_0.6_50 <- rep(NA, 50)
f_0.6_50[1] <- 0.6
for (t in 1:49) f_0.6_50[t+1] <- rbinom(1, 2*N, f_0.6_50[t]) / (2*N) 
plot(x=1:50, y=f_0.6_50, type="l", ylim=c(0,1), lwd=2)

# initial population of 200
N=200
f_0.6_200 <- rep(NA, 50)
f_0.6_200[1] <- 0.6
for (t in 1:199) f_0.6_200[t+1] <- rbinom(1, 2*N, f_0.6_200[t]) / (2*N)
plot(x=1:200, y=f_0.6_200, type="l", ylim=c(0,1), lwd=2)

