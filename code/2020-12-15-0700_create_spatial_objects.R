
## Script to calculate and save spatial polygons from the Slikok project.
## Much of this code is recycled from the earlier data summary work.

## Load libraries.
library(reshape2)
library(unmarked)
library(AHMbook)
library(maptools)
library(rgdal)
library(raster)
library(swfscMisc)
library(rgeos)

## Load functions.
source("functions/functions.R")

data01 <- read.csv("../data/raw_data/occurrence_data/2020-11-12-1400_occurrences.csv", stringsAsFactors=FALSE)

## Separate plot and subplot identifiers.
data01$plot_half <- NA
data01$plot_half[grepl("east half of plot", data01$SPEC_LOCALITY)] <- "east"
data01$plot_half[grepl("west half of plot", data01$SPEC_LOCALITY)] <- "west"

data01$plot_name <- gsub(", east half of plot", "", data01$SPEC_LOCALITY)
data01$plot_name <- gsub(", west half of plot", "", data01$plot_name)
data01$plot_name <- substr(data01$plot_name, nchar(data01$plot_name)-3, nchar(data01$plot_name))

## I think it will be most handy to use Julian dates.
data01$BEGAN_DATE <- as.Date(data01$BEGAN_DATE)
data01$julian_day <- as.numeric(format(data01$BEGAN_DATE, "%j"))

## Summarize the number of plots surveyed per day.
plot_data <- unique(data01[,c("plot_name", "julian_day")])
ag02 <- aggregate(plot_data$plot_name, by=list(plot_data$julian_day), length)
names(ag02) <- c("julian_day", "n_plots")

## I really should summarize the number of plot halves surveyed per day.
supblot_data <- unique(data01[,c("SPEC_LOCALITY", "julian_day")])
ag03 <- aggregate(supblot_data$SPEC_LOCALITY, by=list(supblot_data$julian_day), length)
names(ag03) <- c("julian_day", "n_subplots")

## It would be good to make a table of the number of revisits at each plot.
## I will divide the dates into two time periods.
data01$time_period <- NA
data01$time_period[data01$julian_day < 190] <- 1
data01$time_period[data01$julian_day > 190] <- 2

ag05 <- unique(data01[,c("plot_name", "plot_half", "time_period", "DEC_LAT",  "DEC_LONG")])
ag05 <- melt(ag05, id.vars=c("plot_name", "plot_half", "time_period", "DEC_LAT",  "DEC_LONG"))
ag05 <- dcast(ag05, plot_name + DEC_LAT + DEC_LONG ~ plot_half + time_period, fun.aggregate=length)

colSums(ag05[,c("west_1", "west_2")])
#west_1 west_2 
#    23     22 

## Now that we have the structure, it is time to simulate some data.
## Following the example of Kery and Royle (2016).

site_data <- ag05
names(site_data)[which(names(site_data)=="DEC_LAT")] <- "latitude"
names(site_data)[which(names(site_data)=="DEC_LONG")] <- "longitude"
coordinates(site_data) <- c("longitude", "latitude")

## Projections.
albers <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj4string(site_data) <- CRS(wgs84)
row.names(site_data) <- site_data$plot_name

## Generate circular polygons for the plots.
this_plot <- 1
circs <- list(NA)
for(this_plot in 1:40)
 {
 circ <- circle.polygon(x=site_data$longitude[this_plot],
  y=site_data$latitude[this_plot],
  radius=5.64/1000,
  sides=32,
  by.length=FALSE,
  units="km",
  poly.type="gc.earth",
  ellipsoid=datum("wgs84")
  )
 p = Polygon(circ)
 ps = Polygons(list(p),1)
 ps@ID <- site_data$plot_name[this_plot]
 circs[[this_plot]] <- ps 
 }
#ps = Polygons(circs,1)
sps = SpatialPolygons(circs)
plot(sps) 
proj4string(sps) <- CRS(wgs84)
plot(sps)
plotsdf <- SpatialPolygonsDataFrame(sps, data=site_data@data)
plotsdf$wkt <- writeWKT(plotsdf, byid = TRUE)
## Saving this object.
dump("plotsdf", file=paste0("../data/final_data/geodata/", nowstring(), "_plot_circles.R"))
save(plotsdf, file=paste0("../data/final_data/geodata/", nowstring(), "_plot_circles.RData"))
## Now exporting kml.
writeOGR(
 plotsdf, 
 dsn=paste0("../data/final_data/geodata/", nowstring(), "_plot_circles.kml"), 
 layer="plots", 
 driver="KML", 
 dataset_options=c("NameField=plot_name"), 
 overwrite_layer=TRUE
 )

## Now generate subplots.
this_plot <- 1
nsides <- 32
semis <- list(NA)
for(this_plot in 1:1)
 {
 circ <- circle.polygon(x=site_data$longitude[this_plot],
  y=site_data$latitude[this_plot],
  radius=5.64/1000,
  sides=nsides,
  by.length=FALSE,
  units="km",
  brng.limits = 90,
  poly.type="gc.earth",
  ellipsoid=datum("wgs84")
  )
 e <- Polygon(circ[c(1:(round(nsides/2 + 1)),1),]) 
 ps = Polygons(list(e),1)
 ps@ID <- site_data$plot_name[this_plot]
 semis[[this_plot]] <- ps 
 }
semisp = SpatialPolygons(semis)
plot(semisp) 
proj4string(semisp) <- CRS(wgs84)
plot(semisp)
area(semisp)
#[1] 49.72726
## Area of semicircle of radius 5.64 m
pi*(5.64)^2/2
#[1] 49.9664

That semicircle looked good. Now generating the whole set.
## Now generate subplots.
this_plot <- 1
nsides <- 32
semise <- list(NA)
semisw <- list(NA)
for(this_plot in 1:40)
 {
 circ <- circle.polygon(x=site_data$longitude[this_plot],
  y=site_data$latitude[this_plot],
  radius=5.64/1000,
  sides=nsides,
  by.length=FALSE,
  units="km",
  brng.limits = 90,
  poly.type="gc.earth",
  ellipsoid=datum("wgs84")
  )
 e <- Polygon(circ[c(1:(round(nsides/2 + 1)),1),])
 w <- Polygon(circ[c(round(nsides/2 + 1):(nsides+1),round(nsides/2 + 1)),]) 
 pse = Polygons(list(e),1)
 psw = Polygons(list(w),1)
 pse@ID <- paste0(site_data$plot_name[this_plot], "_east")
 psw@ID <- paste0(site_data$plot_name[this_plot], "_west")
 semise[[this_plot]] <- pse
 semisw[[this_plot]] <- psw 
 }
sde <- site_data@data
row.names(sde) <- paste0(site_data$plot_name, "_east")
semisep = SpatialPolygons(semise)
proj4string(semisep) <- CRS(wgs84)
semisedf <- SpatialPolygonsDataFrame(semisep, data=sde)
sdw <- site_data@data
row.names(sdw) <- paste0(site_data$plot_name, "_west")
semiswp = SpatialPolygons(semisw)
proj4string(semiswp) <- CRS(wgs84)
semiswdf <- SpatialPolygonsDataFrame(semiswp, data=sdw)
semisdf <- rbind(semisedf, semiswdf)
semisdf$plot_name <- row.names(semisdf)
semisdf$wkt <- writeWKT(semisdf, byid = TRUE) 
## Saving this object.
dump("semisdf", file=paste0("../data/final_data/geodata/", nowstring(), "_subplot_semicircles.R"))
save(semisdf, file=paste0("../data/final_data/geodata/", nowstring(), "_subplot_semicircles.RData"))
## Now exporting kml.
writeOGR(
 semisdf, 
 dsn=paste0("../data/final_data/geodata/", nowstring(), "_subplot_semicircles.kml"), 
 layer="plots", 
 driver="KML", 
 dataset_options=c("NameField=plot_name"), 
 overwrite_layer=TRUE
 )

## Now export these polygons for import into Arctos later.

## Plots.
file_name <- paste0("Slikok_project_plot_", plotsdf$plot_name, ".wkt")
write.csv(as.data.frame(cbind(plotsdf$plot_name, file_name)),
 file=paste0("../data/final_data/geodata/wkt/plots/", nowstring(), "_plot_files.csv"),
 row.names=FALSE
 )
for (this_poly in 1:nrow(plotsdf))
 {
 write(plotsdf$wkt[this_poly],
 file=paste0("../data/final_data/geodata/wkt/plots/", file_name[this_poly]))
 } 

## Subplots. 
file_name <- paste0("Slikok_project_subplot_", semisdf$plot_name, ".wkt")
write.csv(as.data.frame(cbind(semisdf$plot_name, file_name)),
 file=paste0("../data/final_data/geodata/wkt/subplots/", nowstring(), "_subplot_files.csv"),
 row.names=FALSE
 )
for (this_poly in 1:nrow(semisdf))
 {
 write(semisdf$wkt[this_poly],
 file=paste0("../data/final_data/geodata/wkt/subplots/", file_name[this_poly]))
 }  
