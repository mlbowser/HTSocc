
## Script to start simulating data comparable to the Slikok arthropod dataset.

## Load libraries.
#library(unmarked)
#library(AHMbook)
library(maptools)
#library(rgdal)
library(raster)
library(jagsUI)
library(MCMCglmm)

## Load functions.
source("functions/functions.R")

## Load data.
load("../data/final_data/geodata/2020-12-16-1026_plot_circles.RData")
load("../data/final_data/geodata/2020-12-16-1026_subplot_semicircles.RData")

## Now trying higher numbers of simulations.
M <- 40 ## Number of sites.
J <- 2 ## Number of visits.
y <- matrix(NA, nrow=M, ncol=J)
simreps <- 1000 ## Number of simulations per parameter combination.
psi <- (1:9)/10 ## Values of psi to simulate.
p <- (1:9)/10 ## Values of p to simulate.
est <- array(dim=c(2, 8, simreps, length(psi), length(p)))

et <- system.time(
for (this_psi in 1:length(psi))
 {
 for (this_p in 1:length(p))
  {
  for (this_sim in 1:simreps)
   {
   ## Counter.
   cat("psi", this_psi, "p", this_p, "simulation", this_sim)   
   ## Generate presence/absence data.
   z <- rbinom(n=M, size=1, prob=psi[this_psi]) ## Realizations of occurrence.
   ## Generate detection/nondetection data.  
   for (j in 1:J)
    {
    y[,j] <- rbinom(n=M, size=1, prob=z*p[this_p])
    } 
   ## Insert NAs for sites sampled only once. 
   y[,2][plotsdf$west_2==0] <- NA
   ## Now run model.
   str(win.data <- list(y=y, M=M, J=J))
   ## Initial values.
   zst <- apply(y, 1, max, na.rm=TRUE)
   inits <- function(){list(z=zst)}
   ## Parameters monitored
   params <- c("psi", "p")
   ## MCMC settings.
   ni <- 5000; nt <- 1; nb <- 1000; nc <- 1
   ## Call JAGS and summarize posteriors.
   fm2 <- jags(win.data,
    inits,
    params,
    "2020-12-17-0623_model.txt",
    n.chains=nc,
    n.thin=nt,
    n.iter=ni,
    n.burnin=nb
    )
   ## Save results.
   est[,1:7,this_sim,this_psi,this_p] <- fm2$summary[1:2,1:7]
   ## Calculate posterior modes.
   fmm <- mcmc(fm2$samples[[1]])
   est[,8,this_sim,this_psi,this_p] <- posterior.mode(fmm)[1:2]
   }
  }
 }
)
et

## Labelling the output to make it more clear.
dimnames(est) <- list(
 variable=c("psi", "p"),
 statistic=c(colnames(fm2$summary)[1:7], "mode"),
 simulation=1:simreps,
 psi=psi,
 p=p
 )

## Saving results.
save(est, file=paste0("../data/final_data/occupancy/", nowstring(), "_simulation_results.RData")) 

