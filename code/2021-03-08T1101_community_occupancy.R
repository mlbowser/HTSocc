
## Trying out the community occupancy model with bivariate species specific random effects from Kery and Royle (2016), p. 668

## As I am using it initially there are problems with not dealing with the separate time periods. This will just be a test of running the model.

## Load libraries.
library(reshape)
library(maptools)
library(jagsUI)
library(MCMCglmm)

## Load functions.
source("functions/functions.R")

## Load data.
load("../data/final_data/geodata/2020-12-16-1026_plot_circles.RData")
nlcdm <- read.csv("../data/final_data/geodata/2020-12-29-0937_plots_nlcd_classes.csv", row.names=1)
d1 <- read.csv("../data/raw_data/occurrence_data/2020-11-12-1400_occurrences.csv", stringsAsFactors=FALSE)

## Now reformat data to conform to the current model.
d1$BEGAN_DATE <- as.Date(d1$BEGAN_DATE)
d1$julday <- as.numeric(format(d1$BEGAN_DATE, "%j"))
d2 <- d1

## Reshaping to fit
d2$Rep <- NA
sle <- grepl("east half", d2$SPEC_LOCALITY)
d2$Rep[sle] <- 1
slw <- grepl("west half", d2$SPEC_LOCALITY)
d2$Rep[slw] <- 2

## Making visits 3 and 4...
sl1 <- d2$julday > 190
d2$Rep[sl1] <- d2$Rep[sl1] + 2

d2$Point <- gsub(", east half of plot", "", d2$SPEC_LOCALITY)
d2$Point <- gsub(", west half of plot", "", d2$Point)
d2$Point <- substr(d2$Point, nchar(d2$Point)-3, nchar(d2$Point))

names(d2)[which(names(d2)=="BEGAN_DATE")] <- "Date"

names(d2)[which(names(d2)=="SCIENTIFIC_NAME")] <- "Species"

data1 <- d2[,c("Point", "Date", "Species", "Rep")]
data1 <- unique(data1)

spp <- unique(data1$Species)
spp <- spp[order(spp)]
nspp <- length(spp)

M <- 40 ## Number of sites.
T <- 2 ## Number of time periods.
J <- 2 ## Number of visits.

site <- levels(as.factor(data1$Point))

## Formatting times/events. 
s1 <- as.numeric(format(as.Date("2016-06-14"), format="%j"))
e1 <- as.numeric(format(as.Date("2016-06-17"), format="%j"))
s2 <- as.numeric(format(as.Date("2016-07-18"), format="%j"))
e2 <- as.numeric(format(as.Date("2016-08-09"), format="%j"))
## transformation.
transtime <- function(x){
 (x - (s1+e2)/2)*2/(e2-s1)
 } 
event_date <- matrix(NA, nrow=M, ncol=T)
dimnames(event_date) <- list(site=site, time_period=c("early", "late"))
event_data <- unique(data1[,c("Point", "Date")])
event_data <- event_data[order(event_data$Point, event_data$Date),]
for (this_site in 1:M)
 {
 sl <- event_data$Point == site[this_site]
 event_date[this_site,1] <- as.numeric(format(event_data[sl,][1,2], "%j"))
 event_date[this_site,2] <- as.numeric(format(event_data[sl,][2,2], "%j"))
 }
event_date[,] <- apply(event_date[,], c(1,2), transtime)
jdate <- event_date

yf <- array(0, dim=c(nspp, M, T, J))
dimnames(yf) <- list(species=spp,
 site=site,
 time_period=c("early", "late"),
 subplot=c("E", "W")
 )

this_sp <- 1 # For testing. 
for (this_sp in 1:nspp)
 {
 slsp <- data1$Species == spp[this_sp]
 datasp <- data1[slsp,]
  for (this_site in 1:M)
   {
   sls <- datasp$Point == site[this_site]
   datast <- datasp[sls,]
   if (nrow(datast) > 0)
    {
	for (this_row in 1:nrow(datast))
	 {
	 if (datast$Rep[this_row]==1) {yf[this_sp,this_site,1,1] <- 1} 
	 if (datast$Rep[this_row]==2) {yf[this_sp,this_site,1,2] <- 1} 
	 if (datast$Rep[this_row]==3) {yf[this_sp,this_site,2,1] <- 1} 
	 if (datast$Rep[this_row]==4) {yf[this_sp,this_site,2,2] <- 1} 
	 }
	}
   }
 ## Insert NAs for sites sampled only once.
 yf[this_sp,,1,2][plotsdf$west_1==0] <- NA   
 yf[this_sp,,2,2][plotsdf$west_2==0] <- NA  
 }

nonforest <- nlcdm$nonforest

## Formatting data for this particular occupancy model.

J <- 4
y <- array(NA, dim=c(nspp, M, J))
dimnames(y) <- list(species=spp,
 site=site,
 sample=1:4 
 )
y[,,1:2] <- yf[,,1,] 
y[,,3:4] <- yf[,,2,]
 

win.data <- list(
 y=y,
 M=M,
 J=J,
 nspec=dim(y)[1],
 R=matrix(c(5,0,0,1), ncol=2),
 df=3
 )

## Initial values.
zst <- apply(y, c(1,2), max, na.rm=TRUE)
inits <- function()list(z=zst, Omega=matrix(c(1,0,0,1), ncol=2), eta=matrix(0, nrow=nspec, ncol=2)) 
 
## Parameters monitored.
params <- c("mu.eta", "probs", "psi", "p", "Nsite", "Nocc.fs", "Sigma", "rho")

# MCMC settings.
ni <- 2000; nt <- 10; nb <- 500; nc <- 3

## Run it.
out6 <- jags(win.data, inits, params, "2021-03-08T1039AKST_model.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=TRUE)

print(out6,3)
## There seemed to be too much correlation in the detection and occupancy values of these.

## Let's try looking at just one season (June sampling).
J=2
win.data <- list(
 y=y[,,1:2],
 M=M,
 J=J,
 nspec=dim(y)[1],
 R=matrix(c(5,0,0,1), ncol=2),
 df=3
 )

## Initial values.
zst <- apply(y, c(1,2), max, na.rm=TRUE)
inits <- function()list(z=zst, Omega=matrix(c(1,0,0,1), ncol=2), eta=matrix(0, nrow=nspec, ncol=2)) 
 
## Parameters monitored.
params <- c("mu.eta", "probs", "psi", "p", "Nsite", "Nocc.fs", "Sigma", "rho")

# MCMC settings.
ni <- 2000; nt <- 10; nb <- 500; nc <- 3

## Run it.
out6 <- jags(win.data, inits, params, "2021-03-08T1039AKST_model.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=TRUE)
## That looked much better.


