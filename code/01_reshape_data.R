
## Reshape data

library(knitr)

data01 <- read.csv("../data/raw_data/occurrence_data/2020-11-12-1400_occurrences.csv", stringsAsFactors=FALSE)

out_file <- "../documents/summaries/data_summaries.md"

## How many records were there?d
nrow(data01)

write("
Number of records: \\", out_file)
write(nrow(data01), out_file, append=TRUE)

## How many unique identifications were there?
length(unique(data01$SCIENTIFIC_NAME))

write("
Number of unique identifications: \\", out_file, append=TRUE)
write(length(unique(data01$SCIENTIFIC_NAME)), out_file, append=TRUE)

## Are all of the species X event records unique? (It is possible that some MOTUs ended up with the same identificaions.)
nrow(unique(data01[,c("SCIENTIFIC_NAME", "SPEC_LOCALITY", "BEGAN_DATE")]))
## Yes, it appears that they are all unique.

write("
Number of unique species X event records: \\", out_file, append=TRUE)
write(length(unique(data01$SCIENTIFIC_NAME)), out_file, append=TRUE)

## I think it will be most handy to use Julian dates.
data01$BEGAN_DATE <- as.Date(data01$BEGAN_DATE)
data01$julian_day <- as.numeric(format(data01$BEGAN_DATE, "%j"))

## Summarize the number of records per day.
ag01 <- aggregate(data01$GUID, by=list(data01$julian_day), length)
names(ag01) <- c("julian_day", "n_observations")

write("
Number of records by Julian day", out_file, append=TRUE)
write(kable(ag01), out_file, append=TRUE)

## Plot number of observations versus Julian day.
image_file <- "../documents/images/observations_vs_julian_day.png"
width <- 600
png(filename=image_file,
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

image_caption <- "Numbers of observation records by Julian day."
write(paste0("
![", image_caption, "](", gsub("../documents/", "../", image_file), ")\\
", image_caption, "
"), out_file, append=TRUE)


