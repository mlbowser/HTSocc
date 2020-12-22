
## Now simulating covariates.

## Load libraries.
library(maptools)
library(jagsUI)
library(MCMCglmm)

## Load functions.
source("functions/functions.R")

## Load data.
load("../data/final_data/geodata/2020-12-16-1026_plot_circles.RData")
load("../data/final_data/geodata/2020-12-16-1026_subplot_semicircles.RData")

## Sweep net samples of terrestrial arthropods were collected in two sampling periods: from 14 to 17 June 2016 and from from 18 July to 9 August 2016
## I am going to simulate data similar to these date ranges.
s1 <- as.numeric(format(as.Date("2016-06-14"), format="%j"))
e1 <- as.numeric(format(as.Date("2016-06-17"), format="%j"))
s2 <- as.numeric(format(as.Date("2016-07-18"), format="%j"))
e2 <- as.numeric(format(as.Date("2016-08-09"), format="%j"))

sample(s1:e1, 1)

## transformation
transtime <- function(x){
 (x - (s1+e2)/2)*2/(e2-s1)
 }

t1vals <- transtime(s1:e1) 
t2vals <- transtime(s2:e2) 
  
M <- 40 ## Number of sites.
T <- 2 ## Number of time periods.
J <- 2 ## Number of visits.
y <- matrix(NA, nrow=M, ncol=J)
simreps <- 1 ## Number of simulations per parameter combination.
#psi <- c(0.25, 0.5, 0.75) ## Values of psi to simulate.
p <- 0.5 ## Values of p to simulate.
betapsi0 <- 0
betapsitime1 <- 0
est <- array(dim=c(3, 8, simreps, length(betapsitime1), length(p)))

this_psi <- 1 ## These set values are for testing.
this_p <- 1
this_sim <- 1

et <- system.time(
for (this_psi in 1:length(betapsitime1))
 {
 for (this_p in 1:length(p))
  {
  for (this_sim in 1:simreps)
   {
   ## Counter.
   cat("psi", this_psi, "p", this_p, "simulation", this_sim)
   ## Dates.
   t1 <- sample(s1:e1, M, replace=TRUE)
   t2 <- sample(s2:e2, M, replace=TRUE)
   t1t <- transtime(t1)
   t2t <- transtime(t2)
   jdate <- cbind(t1t,t2t)
   psi1 <- plogis(betapsi0[this_psi] + betapsitime1[this_psi]*t1t)
   psi2 <- plogis(betapsi0[this_psi] + betapsitime1[this_psi]*t2t)
   ## Generate presence/absence data.
   ## Realizations of occurrence.
   z <- matrix(NA, nrow=M, ncol=2)
   z[,1] <- rbinom(n=M, size=1, prob=psi1[this_psi]) 
   z[,2] <- rbinom(n=M, size=1, prob=psi2[this_psi])
   ## Generate detection/nondetection data.
   y <- array(NA, dim=c(M, T, J))   
   for (j in 1:2)
    {
    y[,1,j] <- rbinom(n=M, size=1, prob=z[,1]*p[this_p])
    } 
   for (j in 1:2)
    {
    y[,2,j] <- rbinom(n=M, size=1, prob=z[,2]*p[this_p])
    } 
   ## Insert NAs for sites sampled only once.
   y[,1,2][plotsdf$west_1==0] <- NA   
   y[,2,2][plotsdf$west_2==0] <- NA
   ## Now run model.
   str(win.data <- list(y=y, M=M, T=T, J=J, jdate=jdate))
   ## Initial values.
   zst <- apply(y, c(1,2), max, na.rm=TRUE)
   inits <- function(){list(z=zst)}
   ## Parameters monitored
   params <- c("p", "betapsi0", "betapsitime1")
   ## MCMC settings.
   ni <- 5000; nt <- 1; nb <- 1000; nc <- 1
   ## Call JAGS and summarize posteriors.
   fm2 <- jags(win.data,
    inits,
    params,
    "2020-12-21-1249_model.txt",
    n.chains=nc,
    n.thin=nt,
    n.iter=ni,
    n.burnin=nb
    )
   ## Save results.
   est[,1:7,this_sim,this_psi,this_p] <- fm2$summary[1:3,1:7]
   ## Calculate posterior modes.
   fmm <- mcmc(fm2$samples[[1]])
   est[,8,this_sim,this_psi,this_p] <- posterior.mode(fmm)[1:3]
   }
  }
 }
)
et

## That worked!

## Labelling the output to make it more clear.
dimnames(est) <- list(
 variable=c("psi", "p", "betapsi0", "betapsitime1"),
 statistic=c(colnames(fm2$summary)[1:7], "mode"),
 simulation=1:simreps,
 beta=1:length(betapsitime1), ## These are just combinations of betas.
 p=p
 )

## Now trying with some more values to see if all the iterating is working...

M <- 40 ## Number of sites.
T <- 2 ## Number of time periods.
J <- 2 ## Number of visits.
simreps <- 10 ## Number of simulations per parameter combination.
#psi <- c(0.25, 0.5, 0.75) ## Values of psi to simulate.
p <- c(0.25, 0.5, 0.75) ## Values of p to simulate.
betapsi0 <- c(-1, -1, 1, 1)
betapsitime1 <- c(-1, 1, -1, 1) ## betapsi0 and betapsitime1 must be the same length.
est <- array(dim=c(3, 8, simreps, length(betapsitime1), length(p)))

this_psi <- 1 ## These set values are for testing.
this_p <- 1
this_sim <- 1

et <- system.time(
for (this_psi in 1:length(betapsitime1))
 {
 for (this_p in 1:length(p))
  {
  for (this_sim in 1:simreps)
   {
   ## Counter.
   cat("psi", this_psi, "p", this_p, "simulation", this_sim)
   ## Dates.
   t1 <- sample(s1:e1, M, replace=TRUE)
   t2 <- sample(s2:e2, M, replace=TRUE)
   t1t <- transtime(t1)
   t2t <- transtime(t2)
   jdate <- cbind(t1t,t2t)
   psi1 <- plogis(betapsi0[this_psi] + betapsitime1[this_psi]*t1t)
   psi2 <- plogis(betapsi0[this_psi] + betapsitime1[this_psi]*t2t)
   ## Generate presence/absence data.
   ## Realizations of occurrence.
   z <- matrix(NA, nrow=M, ncol=2)
   z[,1] <- rbinom(n=M, size=1, prob=psi1[this_psi]) 
   z[,2] <- rbinom(n=M, size=1, prob=psi2[this_psi])
   ## Generate detection/nondetection data.
   y <- array(NA, dim=c(M, T, J))   
   for (j in 1:2)
    {
    y[,1,j] <- rbinom(n=M, size=1, prob=z[,1]*p[this_p])
    } 
   for (j in 1:2)
    {
    y[,2,j] <- rbinom(n=M, size=1, prob=z[,2]*p[this_p])
    } 
   ## Insert NAs for sites sampled only once.
   y[,1,2][plotsdf$west_1==0] <- NA   
   y[,2,2][plotsdf$west_2==0] <- NA
   ## Now run model.
   str(win.data <- list(y=y, M=M, T=T, J=J, jdate=jdate))
   ## Initial values.
   zst <- apply(y, c(1,2), max, na.rm=TRUE)
   inits <- function(){list(z=zst)}
   ## Parameters monitored
   params <- c("p", "betapsi0", "betapsitime1")
   ## MCMC settings.
   ni <- 5000; nt <- 1; nb <- 1000; nc <- 1
   ## Call JAGS and summarize posteriors.
   fm2 <- jags(win.data,
    inits,
    params,
    "2020-12-21-1249_model.txt",
    n.chains=nc,
    n.thin=nt,
    n.iter=ni,
    n.burnin=nb
    )
   ## Save results.
   est[,1:7,this_sim,this_psi,this_p] <- fm2$summary[1:3,1:7]
   ## Calculate posterior modes.
   fmm <- mcmc(fm2$samples[[1]])
   est[,8,this_sim,this_psi,this_p] <- posterior.mode(fmm)[1:3]
   }
  }
 }
)
et
   user  system elapsed 
  95.44   32.74  129.27 

## Labelling the output to make it more clear.
dimnames(est) <- list(
 variable=c("p", "betapsi0", "betapsitime1"),
 statistic=c(colnames(fm2$summary)[1:7], "mode"),
 simulation=1:simreps,
 beta=1:length(betapsitime1), ## These are just combinations of betas.
 p=p
 )
 
## Ok that took
tpersim <- et[1]/(simreps*length(p)*length(betapsitime1))
tpersim
user.self 
0.7953333 ## Seconds per simulation.

## So if I run 100 simreps * 3 levels of p * 9 combinations of betas, that would take
tpersim * 100*3*9
user.self 
    2147.4  ## seconds or 
tpersim * 100*3*9/60
user.self 
    35.79 ## minutes.
	
## Going for it.	
M <- 40 ## Number of sites.
T <- 2 ## Number of time periods.
J <- 2 ## Number of visits.
simreps <- 100 ## Number of simulations per parameter combination.
p <- c(0.25, 0.5, 0.75) ## Values of p to simulate.
betapsi0 <- c(-1, -1, -1, 0, 0, 0, 1, 1, 1)
betapsitime1 <- c(-1, 0, 1, -1, 0, 1, -1, 0, 1) ## betapsi0 and betapsitime1 must be the same length.
est <- array(dim=c(3, 8, simreps, length(betapsitime1), length(p)))
dimnames(est) <- list(
 variable=c("p", "betapsi0", "betapsitime1"),
 statistic=c(colnames(fm2$summary)[1:7], "mode"),
 simulation=1:simreps,
 beta=1:length(betapsitime1), ## These are just combinations of betas.
 p=p
 )

this_psi <- 1 ## These set values are for testing.
this_p <- 1
this_sim <- 1

et <- system.time(
for (this_psi in 1:length(betapsitime1))
 {
 for (this_p in 1:length(p))
  {
  for (this_sim in 1:simreps)
   {
   ## Counter.
   cat("beta", this_psi, "p", this_p, "simulation", this_sim[[1]])
   ## Dates.
   t1 <- sample(s1:e1, M, replace=TRUE)
   t2 <- sample(s2:e2, M, replace=TRUE)
   t1t <- transtime(t1)
   t2t <- transtime(t2)
   jdate <- cbind(t1t,t2t)
   psi1 <- plogis(betapsi0[this_psi] + betapsitime1[this_psi]*t1t)
   psi2 <- plogis(betapsi0[this_psi] + betapsitime1[this_psi]*t2t)
   ## Generate presence/absence data.
   ## Realizations of occurrence.
   z <- matrix(NA, nrow=M, ncol=2)
   z[,1] <- rbinom(n=M, size=1, prob=psi1[this_psi]) 
   z[,2] <- rbinom(n=M, size=1, prob=psi2[this_psi])
   ## Generate detection/nondetection data.
   y <- array(NA, dim=c(M, T, J))   
   for (j in 1:2)
    {
    y[,1,j] <- rbinom(n=M, size=1, prob=z[,1]*p[this_p])
    } 
   for (j in 1:2)
    {
    y[,2,j] <- rbinom(n=M, size=1, prob=z[,2]*p[this_p])
    } 
   ## Insert NAs for sites sampled only once.
   y[,1,2][plotsdf$west_1==0] <- NA   
   y[,2,2][plotsdf$west_2==0] <- NA
   ## Now run model.
   str(win.data <- list(y=y, M=M, T=T, J=J, jdate=jdate))
   ## Initial values.
   zst <- apply(y, c(1,2), max, na.rm=TRUE)
   inits <- function(){list(z=zst)}
   ## Parameters monitored
   params <- c("p", "betapsi0", "betapsitime1")
   ## MCMC settings.
   ni <- 5000; nt <- 1; nb <- 1000; nc <- 1
   ## Call JAGS and summarize posteriors.
   fm2 <- jags(win.data,
    inits,
    params,
    "2020-12-21-1249_model.txt",
    n.chains=nc,
    n.thin=nt,
    n.iter=ni,
    n.burnin=nb
    )
   ## Save results.
   est[,1:7,this_sim,this_psi,this_p] <- fm2$summary[1:3,1:7]
   ## Calculate posterior modes.
   fmm <- mcmc(fm2$samples[[1]])
   est[,8,this_sim,this_psi,this_p] <- posterior.mode(fmm)[1:3]
   }
  }
 }
)
et
   user  system elapsed 
2099.42  667.83 2808.19 
    
## Saving results.
save(est, file=paste0("../data/final_data/occupancy/", nowstring(), "_simulation_results.RData")) 

