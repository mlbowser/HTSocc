
## Reshape data

data01 <- read.csv("../data/raw_data/occurrence_data/2020-11-12-1400_occurrences.csv", stringsAsFactors=FALSE)

## How many records were there?
nrow(data01)
#[1] 2375

## How many unique identifications were there?
length(unique(data01$SCIENTIFIC_NAME))
#[1] 975

## Are all of the species X event records unique? (It is possible that some MOTUs ended up with the same identificaions.)
nrow(unique(data01[,c("SCIENTIFIC_NAME", "SPEC_LOCALITY", "BEGAN_DATE")]))
#[1] 2375
## Yes, it appears that they are all unique.

## I think it will be most handy to use Julian dates.
data01$BEGAN_DATE <- as.Date(data01$BEGAN_DATE)
data01$julian_day <- as.numeric(format(data01$BEGAN_DATE, "%j"))

## Summarize the number of records per day.
ag01 <- aggregate(data01$GUID, by=list(data01$julian_day), length)
names(ag01) <- c("julian_day", "n_observations")

## Plot number of observations versus Julian day.
width <- 600
png(filename=paste0("../documents/images/", format(Sys.time(), format="%Y-%m-%d-%H%M"), "_observations_vs_julian_day.png"),
 width=width,
 height=round(width/1.618),
 pointsize=12
 )
plot(ag01$julian_day,
 ag01$n_observations,
 type="h",
 xlab="Julian day",
 ylab="Number of observations",
 lwd=5,
 lend=2
 )
dev.off()
## It might make more sense to plot the number of observations per unit effort.


