
to10 <- function(x)
 {
 if (x > 0) {1}
 else {0}
 }

nowstring <- function()
 {
 format(Sys.time(), format="%Y-%m-%d-%H%M")
 }

## Projections.
albers <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 