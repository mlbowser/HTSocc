
## Script to start simulating data comparable to the Slikok arthropod dataset.

## Load libraries.
library(unmarked)
library(AHMbook)
library(maptools)
library(rgdal)
library(raster)

## Load functions.
source("functions/functions.R")

## Load data.
load("../data/final_data/geodata/2020-12-16-1026_plot_circles.RData")
load("../data/final_data/geodata/2020-12-16-1026_subplot_semicircles.RData")

simreps <- 10
nsites <- 40
nsurveys <- 2
p <- array(dim=c(simreps, 3, 3))
estimates <- array(dim=c(2, simreps, 3, 3))

system.time( ## Time the whole thing.
for (i in 1:simreps) ## Loop i over simreps
 {
 det.prob <- runif(1, 0.01, 0.99)
 data <- simOcc(M=nsites,
  J=nsurveys,
  mean.occupancy=0.5,
  beta1=0, beta2=0,
  beta3=0,
  mean.detection=det.prob,
  time.effects=c(0,0),
  alpha1=0,
  alpha2=0,
  alpha3=0,
  sd.lp=0,
  b=0,
  show.plot=F
  )
 }
)


