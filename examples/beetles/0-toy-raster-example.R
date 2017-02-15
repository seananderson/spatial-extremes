library(rgdal)

setwd("~/src/spatiotemporal-extremes/examples/beetles/Washington_Feb14132256_100")

ogrListLayers("Washington_Feb14132256_100.gdb")
shapes <- readOGR("Washington_Feb14132256_100.gdb", "IDS_Shapes")

head(shapes@data)
head(shapes@polygons)
head(shapes$ACRES)

plot(shapes@polygons[[1]]@Polygons[[1]]@coords)

plot(shapes)

proj4string(shapes)

# Conversation to a shapefile container
# system("ogr2ogr -f 'ESRI SHAPEFILE' test Washington_Feb14132256_100.gdb")

d <- foreign::read.dbf('test/IDS_attrib.dbf')

library(tidyverse)
d <- d %>% rename(ALLYEARS_ID = ALLYEARS_I)

shapes@data <- left_join(shapes@data,
  select(d, ALLYEARS_ID, DMG_TYPE))

# library(tmap)
# qtm(shapes, "Shape_Area") # plot the basic map

# Try a grid overlay

bb <- bbox(shapes)
cs <- c(3.28084, 3.28084)*6000  # cell size 6km x 6km (for illustration)
# 1 ft = 3.28084 m
cc <- bb[, 1] + (cs/2)  # cell offset
cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
grd

sp_grd <- SpatialGridDataFrame(grd,
  data=data.frame(id=1:prod(cd)),
  proj4string=CRS(proj4string(shapes)))

summary(sp_grd)

library("lattice")
spplot(sp_grd, "id", colorkey=FALSE,
  panel = function(...) {
    panel.gridplot(..., border="black")
    sp.polygons(shapes)
    # sp.points(poi, cex=1.5)
    # panel.text(...)
  })

library(raster)
grid <- raster(extent(shapes))
# Choose its resolution. I will use 0.5 degrees of latitude and longitude.
res(grid) <- 0.5

# Make the grid have the same coordinate reference system (CRS) as the shapefile.
proj4string(grid)<-proj4string(dryland)

# Transform this raster into a polygon and you will have a grid, but without Brazil (or your own shapefile).
gridpolygon <- rasterToPolygons(grid)

#############
# http://gis.stackexchange.com/questions/154682/grid-vs-polygon-overlay-in-r-gives-different-results-than-in-qgis

# Create raster with desired resolution or rows/columns. This is the reference raster that will be used to rasterize the species polygons.
library(raster)
bb <- bbox(shapes)
# xmin=-180; xmax=180; ymin=-90; ymax=90; bin=5
nbin <- 40
bins <- raster(extent(matrix(c(bb["x", "min"],bb["y", "min"],bb["x", "max"],bb["y", "max"]),nrow=2)),
  nrow=length(seq(bb["y", "min"], bb["y", "max"], length.out = nbin)), 
  ncol=length(seq(bb["x", "min"], bb["x", "max"], length.out = nbin)),
  crs=  proj4string(shapes))
bins[] <- 1:ncell(bins)
# Now we can create binary rasters for each species corresponding to the reference raster. We first create an empty stack and then populate it by subsetting each species and rasterizing it to the reference raster. For clarity, we assign names of each species to each raster in the stack.

spp.pa <- stack()
# for(i in unique(spp@data$genus_name) ) {
#   s <- as(spp[spp$genus_name == i,], "SpatialPolygons")
#   spp.pa <- addLayer(spp.pa, rasterize(s, bins, fun = 'count', background = 0))
# } 

rr <- rasterize(shapes, bins, getCover = TRUE)
plot(rr)
plot(shapes)

names(spp.pa) <- unique(spp@data$genus_name)
# Now we can calculate the sum (species counts at each pixel) by using the calc function on the raster stack.

spp.sum <- calc(spp.pa, fun = sum) 
plot(spp.sum)  
# If you want to get at fractional cover of a species in a given pixel you can use the getCover argument in the rasterize function. In this way you could create a conditional argument in the for loop that if a species is less than a defined percent coverage the pixel is 0. Here is a quick example that sets pixels < 30% coverage to 0 else 1.

p = 30
spp.pct <- stack()
for(i in unique(spp@data$genus_name) ) {
  s <- as(spp[spp$genus_name == i,], "SpatialPolygons")
  sfp <- rasterize(s, bins, getCover = TRUE, background = 0)
  sfp[sfp < p] <- 0
  sfp[sfp >= p] <- 1
  spp.pct <- addLayer(spp.pct, sfp)
} 
names(spp.pct) <- unique(spp@data$genus_name)

spp.pct.sum <- calc(spp.pct, fun = sum)

par(mfrow=c(2,1))
plot(spp.sum)
plot(spp.pct.sum)  


