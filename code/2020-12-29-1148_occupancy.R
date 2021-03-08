
## Now I am going to try applying the latest model to real data.

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
dimnames(y) <- list(species=spp,
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
 
## Now run model.
for (this_sp in 1:1)
 {
   y <- yf[this_sp,,,]
   str(win.data <- list(y=y, M=M, T=T, J=J, jdate=jdate, nonforest=nonforest))
   ## Initial values.
   zst <- apply(y, c(1,2), max, na.rm=TRUE)
   inits <- function(){list(z=zst)}
   ## Parameters monitored
   params <- c("p", "betapsiforest0", "betapsinonforest0", "betapsitime1")
   ## MCMC settings.
   ni <- 5000; nt <- 1; nb <- 1000; nc <- 3
   ## Call JAGS and summarize posteriors.
   fm2 <- jags(win.data,
    inits,
    params,
    "2020-12-29-0953_model.txt",
    n.chains=nc,
    n.thin=nt,
    n.iter=ni,
    n.burnin=nb
    )
 }	
 
## I might monitor some derived nodes (occupancy in certain situations). For these I want the midpoints of the two time periods.
midtime1 <- transtime(mean(c(s1, e1))) 
midtime1
[1] -0.9464286

midtime2 <- transtime(mean(c(s2, e2))) 
midtime2
[1] 0.6071429

## Trying a new model.
for (this_sp in 442:442)
 {
   y <- yf[this_sp,,,]
   str(win.data <- list(y=y, M=M, T=T, J=J, jdate=jdate, nonforest=nonforest))
   ## Initial values.
   zst <- apply(y, c(1,2), max, na.rm=TRUE)
   inits <- function(){list(z=zst)}
   ## Parameters monitored
   params <- c("p",
    "betapsiforest0",
    "betapsinonforest0",
    "betapsitime1",
    "psiforestearly",
	"psinonforestearly",
	"psiforestlate",
	"psinonforestlate"
   )
   ## MCMC settings.
   ni <- 10000; nt <- 1; nb <- 2000; nc <- 4
   ## Call JAGS and summarize posteriors.
   fm2 <- jags(win.data,
    inits,
    params,
    "2020-12-30-0653_model.txt",
    n.chains=nc,
    n.thin=nt,
    n.iter=ni,
    n.burnin=nb
    )
 }	


## Now iterating through.
for (this_sp in 1:10)
 {
   outfile <- paste0("../data/final_data/occupancy/", nowstring(), "_output_sp_", substr(1000+this_sp, 2, 4), ".RData")
   y <- yf[this_sp,,,]
   str(win.data <- list(y=y, M=M, T=T, J=J, jdate=jdate, nonforest=nonforest))
   ## Initial values.
   zst <- apply(y, c(1,2), max, na.rm=TRUE)
   inits <- function(){list(z=zst)}
   ## Parameters monitored
   params <- c("p",
    "betapsiforest0",
    "betapsinonforest0",
    "betapsitime1",
    "psiforestearly",
	"psinonforestearly",
	"psiforestlate",
	"psinonforestlate"
   )
   ## MCMC settings.
   ni <- 10000; nt <- 1; nb <- 2000; nc <- 4
   ## Call JAGS and summarize posteriors.
   fm2 <- jags(win.data,
    inits,
    params,
    "2020-12-30-0653_model.txt",
    n.chains=nc,
    n.thin=nt,
    n.iter=ni,
    n.burnin=nb
    )
  save(fm2, file=outfile) 
 }	 
 
## It is obvious that this model is too complicated for most species.
## Paring this down based on earlier info about detection probability.
det1 <- read.csv("../data/final_data/occupancy/2020-11-23-1600_det1.csv", stringsAsFactors=FALSE, row.names=1)
det2 <- read.csv("../data/final_data/occupancy/2020-11-23-1600_det2.csv", stringsAsFactors=FALSE, row.names=1)

shortlist <- (1:nspp)[(det1$q5 >= 0.5) | (det2$q5 >= 0.5)]

results <- list(NULL)
for (this_sp in shortlist)
 { 
   cat(paste0("SPECIES ", this_sp, ", ", spp[this_sp])) 
   y <- yf[this_sp,,,]
   str(win.data <- list(y=y, M=M, T=T, J=J, jdate=jdate, nonforest=nonforest))
   ## Initial values.
   zst <- apply(y, c(1,2), max, na.rm=TRUE)
   inits <- function(){list(z=zst)}
   ## Parameters monitored
   params <- c("p",
    "betapsiforest0",
    "betapsinonforest0",
    "betapsitime1",
    "psiforestearly",
	"psinonforestearly",
	"psiforestlate",
	"psinonforestlate"
   )
   ## MCMC settings.
   ni <- 10000; nt <- 1; nb <- 2000; nc <- 4
   ## Call JAGS and summarize posteriors.
   fm2 <- jags(win.data,
    inits,
    params,
    "2020-12-30-0653_model.txt",
    n.chains=nc,
    n.thin=nt,
    n.iter=ni,
    n.burnin=nb
    )
  results[[this_sp]] <- fm2
 }	 
outfile <- paste0("../data/final_data/occupancy/", nowstring(), "_output.RData") 
save(results, file=outfile)

spp[shortlist[1]]
[1] "Ablabesmyia monilis"
results[[shortlist[1]]]
## Failed to converge. 

spp[shortlist[2]]
results[[shortlist[2]]]
## Convergence failure.

spp[shortlist[43]]
results[[shortlist[43]]]
plot(results[[shortlist[43]]])
## Ochlerotatus communis was the best-looking set of results yet.
d1$BEGAN_DATE[d1$SCIENTIFIC_NAME=="Ochlerotatus communis"]

spp[shortlist[56]]
results[[shortlist[56]]]
plot(results[[shortlist[56]]])
## Convergence failure for Psilidae sp. BOLD:AAF9705.

spp[shortlist[22]]
results[[shortlist[22]]]
plot(results[[shortlist[22]]])
yf[shortlist[22],,,]

spp[shortlist[30]]
results[[shortlist[30]]]
plot(results[[shortlist[30]]])
## That looked good.
yf[shortlist[30],,,]

## Checking Estrandia grandaeva, a species with low detection probability and fairly narrow CI for detectability, but which I know is common.
for (this_sp in which(spp == "Estrandia grandaeva"))
 {
   y <- yf[this_sp,,,]
   str(win.data <- list(y=y, M=M, T=T, J=J, jdate=jdate, nonforest=nonforest))
   ## Initial values.
   zst <- apply(y, c(1,2), max, na.rm=TRUE)
   inits <- function(){list(z=zst)}
   ## Parameters monitored
   params <- c("p",
    "betapsiforest0",
    "betapsinonforest0",
    "betapsitime1",
    "psiforestearly",
	"psinonforestearly",
	"psiforestlate",
	"psinonforestlate"
   )
   ## MCMC settings.
   ni <- 10000; nt <- 1; nb <- 2000; nc <- 4
   ## Call JAGS and summarize posteriors.
   fm2 <- jags(win.data,
    inits,
    params,
    "2020-12-30-0653_model.txt",
    n.chains=nc,
    n.thin=nt,
    n.iter=ni,
    n.burnin=nb
    )
 }	
## That looked pretty good.
yf[which(spp == "Estrandia grandaeva"),,,]

## Overall, though, I think that model is too complicated for most species.
