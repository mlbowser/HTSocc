
## This script reshapes the data and runs a basic multi-species, multi-season occupancy model.
## Code was borrowed from Joseph (2013).

## Load libraries.
library(reshape)
library(rjags)

## Load functions.
source("functions/functions.R")

d1 <- read.csv("../data/raw_data/occurrence_data/2020-11-12-1400_occurrences.csv", stringsAsFactors=FALSE)

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

## Now copying from Occ data script...
#How many citings for each species
data1$Species <- as.character(data1$Species)
data1$Occ <- rep(1, dim(data1)[1])
#How many citings for each species
total.count = tapply(data1$Occ, data1$Species, sum)

#Find the number of unique species
uspecies = as.character(unique(data1$Species))
uspecies <- uspecies[order(uspecies)]
#n is the number of observed species
n=length(uspecies)

#Find the number of unique sampling locations
upoints = as.character(unique(data1$Point))
#J is the number of sampled points
J=length(upoints)

#The detection/non-detection data is reshaped into a three dimensional
#array X where the first dimension, j, is the point; the second
#dimension, k, is the rep; and the last dimension, i, is the species.
junk.melt=melt(data1,id.var=c("Species", "Point", "Rep"), measure.var="Occ")
X=cast(junk.melt, Point ~ Rep ~ Species, sum)

## Which of these were sampled only once?
data1 <- data1[order(data1$Point, data1$Species, data1$Rep),]
sl2 <- data1$Rep == 2
s2 <- unique(data1[sl2,c("Point", "Rep")])
points <-  unique(data1$Point)
slna2 <- !(points %in% s2$Point)
X[slna2,2,] <- NA

sl2 <- data1$Rep == 4
s4 <- unique(data1[sl2,c("Point", "Rep")])
points <-  unique(data1$Point)
slna4 <- !(points %in% s4$Point)
X[slna4,4,] <- NA

## Reshaping this matrix.
X2 <- array(NA, dim=c(J,2,2,n))
X2[,,1,] <- X[,1:2,]
X2[,,2,] <- X[,3:4,]
X <- X2

#K is a vector of length J indicating the number of reps at each point j
#K=rep(2,J) ## Don't really need this any more.

## matrix of season values.
seas <- matrix(c(rep(1,80), rep(2,80)), nrow=40, ncol=4)

#Create the necessary arguments to run the bugs() command
#Load all the data
sp.data = list(n=n, J=J, X=X)

#Specify the parameters to be monitored
sp.params = list('u', 'v', 'd')

#Specify the initial values
    sp.inits = function() {
	zinit <- array(rbinom(n*J*2, size=1, prob=runif(n,0.1,0.9)), dim=c(J,2,n))
	for(j in 1:J)
	 {
	 for (i in 1:n)
	  {
      for (t in 1:2)
       {
	   zinit[j,t,i] <- max(X[j,,t,i], na.rm=TRUE)
	   }	  
	  }
	 }
    list(
	    u=array(runif(n*2,0.1,0.9), dim=c(n,2)), 
		v=array(runif(n*2,0.1,0.9), dim=c(n,2)),
        Z = zinit
        )
    }

#Run the model and call the results “fit”
#fit = bugs(sp.data, sp.inits, sp.params, "basicmodel.txt", debug=TRUE,
#         n.chains=3, n.iter=10000, n.burnin=5000, n.thin=5)
params <- as.character(sp.params)		 
ocmod <- jags.model(file = "model_two_season.txt", inits = sp.inits, data = sp.data, n.chains = 1)
nburn <- 1000
update(ocmod, n.iter = nburn)
out <- coda.samples(ocmod, n.iter = 20000, thin=10, variable.names = params)

#pdf(file=paste0(nowstring(), "_out.pdf"),
# width=8.5,
# height=11
# )
#plot(out)
#dev.off()

#See baseline estimates of species-specific occupancy and detection 
species.occ1 = out[[1]][,(n+1):(2*n)]
species.occ2 = out[[1]][,(2*n+1):(3*n)]
species.det1 = out[[1]][,(3*n+1):(4*n)]
species.det2 = out[[1]][,(4*n+1):(5*n)]
species.d = out[[1]][,1:n]

#Show occupancy and detection estimates 
psi1 = (species.occ1)
psi2 = (species.occ2)
p1   = (species.det1)
p2   = (species.det2)
d = species.d

occ1.matrix <- cbind(apply(psi1,2,mean),apply(psi1,2,sd),apply(psi1,2,quantile,probs=0.025),apply(psi1,2,quantile,probs=0.5),apply(psi1,2,quantile,probs=0.975))
colnames(occ1.matrix) = c("mean occupancy", "sd occupancy", "q025", "q5", "q975")
rownames(occ1.matrix) = uspecies
occ1.matrix <- as.data.frame(occ1.matrix)
occ1.matrix$diff <- occ1.matrix$q975 - occ1.matrix$q025

occ2.matrix <- cbind(apply(psi2,2,mean),apply(psi2,2,sd),apply(psi2,2,quantile,probs=0.025),apply(psi2,2,quantile,probs=0.5),apply(psi2,2,quantile,probs=0.975))
colnames(occ2.matrix) = c("mean occupancy", "sd occupancy", "q025", "q5", "q975")
rownames(occ2.matrix) = uspecies
occ2.matrix <- as.data.frame(occ2.matrix)
occ2.matrix$diff <- occ2.matrix$q975 - occ2.matrix$q025

det1.matrix <- cbind(apply(p1,2,mean),apply(p1,2,sd),apply(p1,2,quantile,probs=0.025),apply(p1,2,quantile,probs=0.5),apply(p1,2,quantile,probs=0.975))
colnames(det1.matrix) = c("mean detection", "sd detection", "q025", "q5", "q975")
rownames(det1.matrix) = uspecies
det1.matrix <- as.data.frame(det1.matrix)
det1.matrix$diff <- det1.matrix$q975 - det1.matrix$q025

det2.matrix <- cbind(apply(p2,2,mean),apply(p2,2,sd),apply(p2,2,quantile,probs=0.025),apply(p2,2,quantile,probs=0.5),apply(p2,2,quantile,probs=0.975))
colnames(det2.matrix) = c("mean detection", "sd detection", "q025", "q5", "q975")
rownames(det2.matrix) = uspecies
det2.matrix <- as.data.frame(det2.matrix)
det2.matrix$diff <- det2.matrix$q975 - det2.matrix$q025

d.matrix <- cbind(apply(d,2,mean),apply(d,2,sd),apply(d,2,quantile,probs=0.025),apply(d,2,quantile,probs=0.5),apply(d,2,quantile,probs=0.975))
colnames(d.matrix) = c("mean occ diff", "sd mean occ diff", "q025", "q5", "q975")
rownames(d.matrix) = uspecies
d.matrix <- as.data.frame(d.matrix)
d.matrix$diff <- d.matrix$q975 - d.matrix$q025

## Which species increased?
sli <- d.matrix$q025 > 0
sum(sli)

d.matrix[sli,]

## Which species decreased?
sld <- d.matrix$q975 < 0
sum(sld)

d.matrix[sld,]

## Save some of the output.
write.csv(x=occ1.matrix, file=paste0("../data/final_data/occupancy/", nowstring(), "_occ1.csv"))
write.csv(x=occ2.matrix, file=paste0("../data/final_data/occupancy/", nowstring(), "_occ2.csv"))
write.csv(x=det1.matrix, file=paste0("../data/final_data/occupancy/", nowstring(), "_det1.csv"))
write.csv(x=det2.matrix, file=paste0("../data/final_data/occupancy/", nowstring(), "_det2.csv"))
write.csv(x=d.matrix, file=paste0("../data/final_data/occupancy/", nowstring(), "_dif.csv"))

## Histograms of parameters.
image_file <- "../documents/images/basic_two_season_model_histogram_occupancy_season_1.png"
width <- 600
png(filename=image_file,
 width=width,
 height=round(width/1.618),
 pointsize=12
 )
plot(hist(occ1.matrix$q5, breaks=(0:20)/20),
 , xlim=c(0,1),
 main="Histogram of frequency of occupancy for season 1",
 xlab="Frequency of occurrence"
 )
dev.off()

image_file <- "../documents/images/basic_two_season_model_histogram_occupancy_season_2.png"
width <- 600
png(filename=image_file,
 width=width,
 height=round(width/1.618),
 pointsize=12
 )
plot(hist(occ2.matrix$q5, breaks=(0:20)/20),
 , xlim=c(0,1),
 main="Histogram of frequency of occupancy for season 2",
 xlab="Frequency of occurrence"
 )
dev.off()

image_file <- "../documents/images/basic_two_season_model_histogram_detection_season_1.png"
width <- 600
png(filename=image_file,
 width=width,
 height=round(width/1.618),
 pointsize=12
 )
plot(hist(det1.matrix$q5, breaks=(0:20)/20),
 , xlim=c(0,1),
 main="Histogram of frequency of detection probability for season 1",
 xlab="Frequency of occurrence"
 )
dev.off()

image_file <- "../documents/images/basic_two_season_model_histogram_detection_season_2.png"
width <- 600
png(filename=image_file,
 width=width,
 height=round(width/1.618),
 pointsize=12
 )
plot(hist(det2.matrix$q5, breaks=(0:20)/20),
 , xlim=c(0,1),
 main="Histogram of frequency of detection probability for season 2",
 xlab="Frequency of occurrence"
 )
dev.off()
