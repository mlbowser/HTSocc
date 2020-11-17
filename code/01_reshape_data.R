
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
aggregate(data01$GUID, by)