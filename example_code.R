####################################
# simulate data for model with K=1 #
####################################

# specify the number of curves
J <- 10

# specify the fixed parameters
a0 <- 0.01
a1 <- 0.8
b1off <- 0.15

# specify the random parameters for all J nuclei
alpha0j <- rnorm(J, mean=0, sd=sqrt(0.0001))
alpha1j <- rnorm(J, mean=0, sd=sqrt(0.0007))
beta1j <- rlnorm(J, meanlog = 0, sdlog = sqrt(0.01))

# create numeric vector containing the number of measurements for each nucleus
timep <- c(780,780,778,750,780,780,780,780,780,780)

# create list containing the time points at which measurements are available for all nuclei
# measurements are available every 0.15 seconds
timelist <- list()
for(i in 1:J){
  timevec <- seq(0.15,by=0.15,length.out=timep[[i]])
  timelist[[i]] <- timevec
}

# simulate measured concentrations for all nuclei
set.seed(383)
seeds <- sample(1:1000,J)

nobs <- c()
nobs_cum <- c(0)
timeall <- c()
concall <- c()

for(j in 1:J){
  set.seed(seeds[j])
  err <- rnorm(timep[j], mean=0, sd=sqrt(0.0003))
  concvec <- 1 - (a0 + alpha0j[j] + (a1 + alpha1j[j]) * exp(- b1off * beta1j[j] * timelist[[j]])) + err
  concall <- c(concall,concvec)
  nobs <- c(nobs,timep[j])
  nobs_cum <- c(nobs_cum,nobs_cum[j]+timep[j])
  timeall <- c(timeall, timelist[[j]]) 
}

# total number of time points at which measurements are available
ttotal <- sum(nobs) 

##############
# do the fit #
##############
require(frapmm)
result <- fitk1(timeall=timeall,concall=concall,nobs=nobs,ttotal=ttotal,J=J,sv=c(1,1,1,10^-4,1,10^-4,1,10^-4),size=100,burnin=0,thin=1,tune=20,nobs_cum=nobs_cum,sigma2b1=1.0,sigma2beta1=1.0,a0start=0.1,a1start=0.7,b1start=-2.0)

# similar for function fitk2 (for cases where K=2)