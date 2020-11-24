
## Script to summarize data.

## Load libraries.
library(knitr)
library(reshape2)

## Load functions.
source("functions/functions.R")

data01 <- read.csv("../data/raw_data/occurrence_data/2020-11-12-1400_occurrences.csv", stringsAsFactors=FALSE)

out_file <- "../documents/summaries/data_summaries.md"

write("# Data summaries

This output was written by the R script [../../code/01_reshape_data.R](../../code/01_reshape_data.R).
", out_file)

write("## Sampling summaries

To be clear, 160 samples (40 terrestrial plots × 2 subplots per plot × 2 time periods) were collected in the field, but only 125 samples were submitted for High-Throughput Sequencing. All summaries below are based on data obtained from these 125 samples. 
", out_file, append=TRUE)

## Separate plot and subplot identifiers.
data01$plot_half <- NA
data01$plot_half[grepl("east half of plot", data01$SPEC_LOCALITY)] <- "east"
data01$plot_half[grepl("west half of plot", data01$SPEC_LOCALITY)] <- "west"

data01$plot_name <- gsub(", east half of plot", "", data01$SPEC_LOCALITY)
data01$plot_name <- gsub(", west half of plot", "", data01$plot_name)
data01$plot_name <- substr(data01$plot_name, nchar(data01$plot_name)-3, nchar(data01$plot_name))
length(levels(as.factor(data01$plot_name)))

write("
Number of plots: \\", out_file, append=TRUE)
write(length(levels(as.factor(data01$plot_name))), out_file, append=TRUE)

## I think it will be most handy to use Julian dates.
data01$BEGAN_DATE <- as.Date(data01$BEGAN_DATE)
data01$julian_day <- as.numeric(format(data01$BEGAN_DATE, "%j"))

## Summarize the number of plots surveyed per day.
plot_data <- unique(data01[,c("plot_name", "julian_day")])
ag02 <- aggregate(plot_data$plot_name, by=list(plot_data$julian_day), length)
names(ag02) <- c("julian_day", "n_plots")

write("
Number of plots surveyed per day.", out_file, append=TRUE)
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
Number of unique subplots × date sampling events: \\", out_file, append=TRUE)
write(sum(ag03$n_subplots), out_file, append=TRUE)

write("
Number of subplots surveyed per day.", out_file, append=TRUE)
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

## It would be good to make a table of the number of revisits at each plot.
## I will divide the dates into two time periods.
data01$time_period <- NA
data01$time_period[data01$julian_day < 190] <- 1
data01$time_period[data01$julian_day > 190] <- 2

ag05 <- unique(data01[,c("plot_name", "plot_half", "time_period")])
ag05 <- melt(ag05, id.vars=c("plot_name", "plot_half", "time_period"))
ag05 <- dcast(ag05, plot_name ~ plot_half + time_period, fun.aggregate=length)

all4 <- sum(ag05$east_1 & ag05$east_2 & ag05$west_1 & ag05$west_2)
onlye <- sum(ag05$east_1 & ag05$east_2 & !ag05$west_1 & !ag05$west_2)
wonce <- sum(ag05$east_1 & ag05$east_2 & ag05$west_1 & !ag05$west_2)

n_plots <- c(all4, onlye, wonce)
sampling_regime <- c("East and west suplots sampled in time periods 1 and 2.",
 "East subplots sampled in time periods 1 and 2; west subplots not sampled.",
 "East subplots sampled in time periods 1 and 2; west subplots sampled only in time period 1."
 )
sampling_breakdown <- as.data.frame(cbind(sampling_regime, n_plots))

write("
Breakdown of plots by sampling regime.", out_file, append=TRUE)
write(kable(sampling_breakdown), out_file, append=TRUE)

write("## Observation summaries
", out_file, append=TRUE)
 
write("
Number of observation records: \\", out_file, append=TRUE)
write(nrow(data01), out_file, append=TRUE)

## How many unique identifications were there?
length(unique(data01$SCIENTIFIC_NAME))

write("
Number of unique identifications: \\", out_file, append=TRUE)
write(length(unique(data01$SCIENTIFIC_NAME)), out_file, append=TRUE)

## Are all of the species × event records unique? (It is possible that some MOTUs ended up with the same identificaions.)
nrow(unique(data01[,c("SCIENTIFIC_NAME", "SPEC_LOCALITY", "BEGAN_DATE")]))
## Yes, it appears that they are all unique.

write("
Number of unique species × event records: \\", out_file, append=TRUE)
write(nrow(unique(data01[,c("SCIENTIFIC_NAME", "SPEC_LOCALITY", "BEGAN_DATE")])), out_file, append=TRUE)

## How many of the identifications are formally described species?
sum(!grepl(" sp\\. ", unique(data01$SCIENTIFIC_NAME)))

write("
Number of formally described species: \\", out_file, append=TRUE)
write(sum(!grepl(" sp\\. ", unique(data01$SCIENTIFIC_NAME))), out_file, append=TRUE)

## How many of these are BOLD BIN identifications?
sum(grepl("BOLD:", unique(data01$SCIENTIFIC_NAME)))

write("
Number of identifications using BOLD Barcode Index Numbers (BINs): \\", out_file, append=TRUE)
write(sum(grepl("BOLD:", unique(data01$SCIENTIFIC_NAME))), out_file, append=TRUE)

## How many of these are ASVs from this study?
sum(grepl("SlikokO", unique(data01$SCIENTIFIC_NAME)))

write("
Number of identifications using molecular operational taxonomic unit (MOTU) labels from this study: \\", out_file, append=TRUE)
write(sum(grepl("SlikokOtu", unique(data01$SCIENTIFIC_NAME))), out_file, append=TRUE)

## How many of these are other provisional names?
sl_formal <- !grepl(" sp\\. ", unique(data01$SCIENTIFIC_NAME))
sl_BIN <- grepl("BOLD:", unique(data01$SCIENTIFIC_NAME))
sl_MOTU <- grepl("SlikokOtu", unique(data01$SCIENTIFIC_NAME))
sl_other <- !sl_formal & !sl_BIN & !sl_MOTU

write("
Number of other provisional names: \\", out_file, append=TRUE)
write(sum(sl_other), out_file, append=TRUE)

## Accounting.
sum(sl_formal) + sum(sl_BIN) + sum(sl_MOTU) + sum(sl_other)

## I want to have a look at those other provisional names.
unique(data01$SCIENTIFIC_NAME)[sl_other]

## Summarize the number of records per day.
ag01 <- aggregate(data01$GUID, by=list(data01$julian_day), length)
names(ag01) <- c("julian_day", "n_observations")

write("
Number of records by Julian day.", out_file, append=TRUE)
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

## Now I want to look at numbers of observations per unit effort.
ag03$observations_per_supblot <- ag01$n_observations/ag03$n_subplots

write("
Mean number of observations per subplot: \\", out_file, append=TRUE)
write(round(nrow(data01)/sum(ag03$n_subplots), 2), out_file, append=TRUE)

write("
Number of subplots surveyed per day and mean number of observations per subplot.", out_file, append=TRUE)
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

image_caption <- "Number of observations per unit effort over time. Numbers above bars indicate numbers of subplots sampled each day."
write(paste0("
![", image_caption, "](", gsub("../documents/", "../", image_file), ")\\
", image_caption, "
"), out_file, append=TRUE)

## What is the range full range of number of observations per subplot sampling event?
ag04 <- aggregate(data01$GUID, by=list(data01$SPEC_LOCALITY, data01$julian_day), length)
write("
Summary of number of observations per sampling event.", out_file, append=TRUE)
write(kable(as.data.frame(as.matrix(summary(ag04$x))), col.names=c("value")), out_file, append=TRUE)

## Now summarizing frequncy of occurrence.
data01$presence <- 1
data02 <- melt(data01, measure.vars="presence")
freq01 <- dcast(data02, SCIENTIFIC_NAME ~ value, sum)
names(freq01) <- c("species", "frequency")
n_events <- sum(ag03$n_subplots)
freq01$frequency <- freq01$frequency/n_events

write("
Summary of overall frequency of occurrence for all species.", out_file, append=TRUE)
write(kable(as.data.frame(as.matrix(summary(freq01$frequency))), col.names=c("value"), digits=3), out_file, append=TRUE)

## Wow, how many species are represented by a single occurrence?
write("
Number of species represented by a single occurrence:\\", out_file, append=TRUE)
write(sum(freq01$frequency == 1/n_events), out_file, append=TRUE)

## Histogram of overall frequencies.
image_file <- "../documents/images/histogram_overall_frequencies.png"
width <- 600
png(filename=image_file,
 width=width,
 height=round(width/1.618),
 pointsize=12
 )
plot(hist(freq01$frequency),
 main="Histogram of overall frequency of occurrence for all species",
 xlab="Frequency of occurrence"
 )
dev.off()

image_caption <- "Histogram of overall frequencies of occurrences. This frequency was determined as the number of samples in which a species was detected divided by the total number of samples."
write(paste0("
![", image_caption, "](", gsub("../documents/", "../", image_file), ")\\
", image_caption, "
"), out_file, append=TRUE)

## Now summarizing frequncy of occurrence by plots.
freq03 <- dcast(data02, SCIENTIFIC_NAME ~ plot_name, sum)
freq03[,2:ncol(freq03)] <- apply(freq03[,2:ncol(freq03)], c(1,2), to10)
freq03$frequency <- apply(freq03[,2:ncol(freq03)], 1, mean)
 
write("
Summary of overall frequency of occurrence in terms of presence or absence at plots.", out_file, append=TRUE)
write(kable(as.data.frame(as.matrix(summary(freq03$frequency))), col.names=c("value"), digits=3), out_file, append=TRUE)

## How many species were detected at only one plot?
write("
Number of species detected at only one plot:\\", out_file, append=TRUE)
write(sum(freq03$frequency == 1/sum(n_plots)), out_file, append=TRUE)

## Histogram of overall frequencies.
image_file <- "../documents/images/histogram_frequencies_by_plot.png"
width <- 600
png(filename=image_file,
 width=width,
 height=round(width/1.618),
 pointsize=12
 )
plot(hist(freq03$frequency),
 main="Histogram of frequency of occurrence for all species by plots",
 xlab="Frequency of occurrence"
 )
dev.off()
image_caption <- "Histogram of frequncy of occurence by plot. This frequency was determined as the number of plots at which a species was detected divided by the total number of plots."
write(paste0("
![", image_caption, "](", gsub("../documents/", "../", image_file), ")\\
", image_caption, "
"), out_file, append=TRUE)


write("## Cost
", out_file, append=TRUE)

write("
HTS sequencing cost per sample:\\
$85", out_file, append=TRUE)

write("
HTS sequencing cost per observation record:\\", out_file, append=TRUE)
write(paste0("$", round(85/(nrow(data01)/sum(ag03$n_subplots)), 2)), out_file, append=TRUE)
