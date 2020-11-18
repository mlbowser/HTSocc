
## Reshape data

library(knitr)

data01 <- read.csv("../data/raw_data/occurrence_data/2020-11-12-1400_occurrences.csv", stringsAsFactors=FALSE)

out_file <- "../documents/summaries/data_summaries.md"

write("# Data summaries

This output was written by the R script [../../code/01_reshape_data.R](../../code/01_reshape_data.R).
", out_file)

## How many records were there?
nrow(data01)

write("
Number of records: \\", out_file, append=TRUE)
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
 ylim=c(0, max(ag01$n_observations)),
 type="h",
 xlab="Julian day",
 ylab="Number of observations",
 lwd=5,
 lend=2
 )
dev.off()

image_caption <- "Numbers of observation records by Julian day."
write(paste0("
![", image_caption, "](", gsub("../documents/", "../", image_file), ")\\
", image_caption, "
"), out_file, append=TRUE)

## It might make more sense to plot the number of observations per unit effort.

data01$plot_half <- NA
data01$plot_half[grepl("east half of plot", data01$SPEC_LOCALITY)] <- "east"
data01$plot_half[grepl("west half of plot", data01$SPEC_LOCALITY)] <- "west"

data01$plot_name <- gsub(", east half of plot", "", data01$SPEC_LOCALITY)
data01$plot_name <- gsub(", west half of plot", "", data01$plot_name)
length(levels(as.factor(data01$plot_name)))

write("
Number of plots: \\", out_file, append=TRUE)
write(length(levels(as.factor(data01$plot_name))), out_file, append=TRUE)

## Summarize the number of plots surveyed per day.
plot_data <- unique(data01[,c("plot_name", "julian_day")])
ag02 <- aggregate(plot_data$plot_name, by=list(plot_data$julian_day), length)
names(ag02) <- c("julian_day", "n_plots")

write("
Number of plots surveyed per day", out_file, append=TRUE)
write(kable(ag02), out_file, append=TRUE)

## Plot number of plots surveyed versus Julian day.
image_file <- "../documents/images/plots_vs_julian_day.png"
width <- 600
png(filename=image_file,
 width=width,
 height=round(width/1.618),
 pointsize=12
 )
plot(ag02$julian_day,
 ag02$n_plots,
 ylim=c(0, max(ag02$n_plots)),
 type="h",
 xlab="Julian day",
 ylab="Number of plots surveyed",
 lwd=5,
 lend=2
 )
dev.off()

image_caption <- "Number of plots surveyed per day."
write(paste0("
![", image_caption, "](", gsub("../documents/", "../", image_file), ")\\
", image_caption, "
"), out_file, append=TRUE)

## I really should summarize the number of plot halves surveyed per day.
supblot_data <- unique(data01[,c("SPEC_LOCALITY", "julian_day")])
ag03 <- aggregate(supblot_data$SPEC_LOCALITY, by=list(supblot_data$julian_day), length)
names(ag03) <- c("julian_day", "n_subplots")

write("
Number of subplots: \\", out_file, append=TRUE)
write(length(levels(as.factor(data01$SPEC_LOCALITY))), out_file, append=TRUE)

write("
Number of unique subplots X date sampling events: \\", out_file, append=TRUE)
write(sum(ag03$n_subplots), out_file, append=TRUE)

write("
Number of subplots surveyed per day", out_file, append=TRUE)
write(kable(ag03), out_file, append=TRUE)

## Plot number of subplots surveyed versus Julian day.
image_file <- "../documents/images/subplots_vs_julian_day.png"
width <- 600
png(filename=image_file,
 width=width,
 height=round(width/1.618),
 pointsize=12
 )
plot(ag03$julian_day,
 ag03$n_subplots,
 ylim=c(0, max(ag03$n_subplots)),
 type="h",
 xlab="Julian day",
 ylab="Number of subplots surveyed",
 lwd=5,
 lend=2
 )
dev.off()

image_caption <- "Number of subplots surveyed per day."
write(paste0("
![", image_caption, "](", gsub("../documents/", "../", image_file), ")\\
", image_caption, "
"), out_file, append=TRUE)

## Now I want to look at numbers of observations per unit effort.
ag03$observations_per_supblot <- ag01$n_observations/ag03$n_subplots

write("
Mean number of observations per subplot: \\", out_file, append=TRUE)
write(round(nrow(data01)/sum(ag03$n_subplots), 2), out_file, append=TRUE)

write("
Number of subplots surveyed per day and mean number of observations per subplot", out_file, append=TRUE)
write(kable(ag03, digits=2), out_file, append=TRUE)

## Plot number of observations per unit effort over time.
image_file <- "../documents/images/observations_per_subplot_vs_julian_day.png"
width <- 600
png(filename=image_file,
 width=width,
 height=round(width/1.618),
 pointsize=12
 )
plot(ag03$julian_day,
 ag03$observations_per_supblot,
 ylim=c(0, max(ag03$observations_per_supblot)*1.1),
 type="h",
 xlab="Julian day",
 ylab="Number of observations per subplot surveyed",
 lwd=5,
 lend=2
 )
text(ag03$julian_day,
 ag03$observations_per_supblot,
 labels=ag03$n_subplots,
 pos=3,
 cex=0.7
 )
dev.off()

image_caption <- "Number of observations per unit effort over time."
write(paste0("
![", image_caption, "](", gsub("../documents/", "../", image_file), ")\\
", image_caption, "
"), out_file, append=TRUE)

## What is the range full range of number of observations per subplot sampling event?
ag04 <- aggregate(data01$GUID, by=list(data01$SPEC_LOCALITY, data01$julian_day), length)
write("
Summary of number of observations per sampling event: \\", out_file, append=TRUE)
write(kable(as.data.frame(as.matrix(summary(ag04$x))), col.names=c("value")), out_file, append=TRUE)
