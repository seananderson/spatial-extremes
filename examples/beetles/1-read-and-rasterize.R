# Start by downloading raw data from 
# https://foresthealth.fs.usda.gov/portal/Flex/IDS
# Original downloaded on February 13 2016
# Parameters: All years, Pacific Northwest, ArcGIS version 10.0
# Place the unzipped folder in "examples/beetles/"
# And name the folder "Pacific_Northwest"
library(rgdal)
library(tidyverse)
library(raster)
library(viridis)

gdb_folder <- "Pacific_Northwest"
id <- "mountain-pine-beetle-pnw-raster"

# 4 hours of my life to discover this is needed:
# We must be within the directory with the .grb folder 
setwd(file.path("examples/beetles", gdb_folder))

ogrListLayers(paste0(gdb_folder, ".gdb"))
ids_dat <- readOGR(paste0(gdb_folder, ".gdb"), "IDS_Shapes")

# head(ids_dat@data)
# nrow(ids_dat@data)
# plot(ids_dat[ids_dat$SURVEY_YEAR == 2004, ])

# Make a grid overlay 
# and count % cover within each cell
# Modified from:
# http://gis.stackexchange.com/questions/154682/grid-vs-polygon-overlay-in-r-gives-different-results-than-in-qgis

# Create raster with desired resolution or rows/columns
bb <- bbox(ids_dat)
nbin <- 30
bins <- raster(extent(matrix(c(bb["x", "min"], bb["y", "min"], 
      bb["x", "max"], bb["y", "max"]), nrow = 2)),
    nrow = length(seq(bb["y", "min"], bb["y", "max"], length.out = nbin)),
    ncol = length(seq(bb["x", "min"], bb["x", "max"], length.out = nbin)),
    crs =  proj4string(ids_dat)
  )
bins[] <- 1:ncell(bins)

# Now we can create rasters for each year
rr <- raster::stack()
years <- sort(unique(ids_dat@data$SURVEY_YEAR))
for(i in years) {
  message(i)
  s <- as(ids_dat[ids_dat$SURVEY_YEAR == i,], "SpatialPolygons")
  rr <- addLayer(rr, rasterize(s, bins, getCover = TRUE, background = 0))
}
names(rr) <- years

saveRDS(rr, 
  file = paste0("../", id, "-", nbin, "x", nbin, ".rds"))
rr <- readRDS(paste0("../", id, "-", nbin, "x", nbin, ".rds"))

d <- data.frame(rasterToPoints(rr))
d <- gather(d, year, cover, -x, -y)
d <- d %>% mutate(year = as.numeric(sub("X", "", year)))
d <- as_tibble(d)
message(paste0("Total number of covered cells is ", 
  nrow(filter(d, cover > 0)), "."))

g <- filter(d, cover > 0) %>% 
  ggplot(aes(x/1e5, y/1e5, fill = log10(cover))) + geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis(option = "D") +
  theme_light() +
  coord_fixed()
# g
ggsave(paste0("../", id, "-", nbin, "x", nbin, ".pdf"), width = 10, height = 9)

saveRDS(d, 
  file = paste0("../", id, "-", "dataframe", "-", nbin, "x", nbin, ".rds"))

setwd(file.path("..", "..", ".."))
