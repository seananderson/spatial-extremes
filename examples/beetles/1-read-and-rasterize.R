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
# gdb_folder <- "Northern"
id <- "mountain-pine-beetle-pnw-raster"
# id <- "mountain-pine-beetle-northern-raster"

# 4 hours of my life to discover this is needed:
# We must be within the directory with the .grb folder
setwd(file.path("examples/beetles", gdb_folder))

ogrListLayers(paste0(gdb_folder, ".gdb"))
ids_dat <- readOGR(paste0(gdb_folder, ".gdb"), "IDS_Shapes")

# head(ids_dat@data)
# nrow(ids_dat@data)
# plot(ids_dat[ids_dat$SURVEY_YEAR == 2004, ])

# Conversation to a shapefile container
# system("ogr2ogr -f 'CSV' temp Pacific_Northwest.gdb")
# system("ogr2ogr -f 'CSV' temp Northern.gdb")

a <- readr::read_csv('temp/IDS_attrib.csv')

ids_dat@data <- dplyr::left_join(ids_dat@data,
  dplyr::select(a, ALLYEARS_ID, DMG_TYPE, AGNT_NM, HOST))

# Filter to keep only mountain pine beetle with pine hosts
nrow(ids_dat@data)

# cannot get a slot ("Polygons") from an object of type "NULL"...
# so, by year:
out <- list()
years <- sort(unique(ids_dat@data$SURVEY_YEAR))
for(i in seq_along(years)) {
  message(years[i])
  this <- ids_dat[ids_dat$SURVEY_YEAR == years[i], ]
  idx <- this@data$AGNT_NM == "mountain pine beetle" &
      grepl("pine", this@data$HOST)
  out[[i]] <- this[idx, ]
}
# I'm getting an error in 2015, so we will work with 2014 and before
out[[which(years == 2015)]] <- NULL
stopifnot(length(out) == 18L)
ids_dat2 <- do.call("rbind", out)

nrow(ids_dat2@data)
# plot(ids_dat2[ids_dat2$SURVEY_YEAR==2012,])

# Make a grid overlay
# and count % cover within each cell
# Modified from:
# http://gis.stackexchange.com/questions/154682/grid-vs-polygon-overlay-in-r-gives-different-results-than-in-qgis

# Create raster with desired resolution or rows/columns
bb <- bbox(ids_dat2)
nbin <- 20
bins <- raster(extent(matrix(c(bb["x", "min"], bb["y", "min"],
      bb["x", "max"], bb["y", "max"]), nrow = 2)),
    nrow = length(seq(bb["y", "min"], bb["y", "max"], length.out = nbin)),
    ncol = length(seq(bb["x", "min"], bb["x", "max"], length.out = nbin)),
    crs =  proj4string(ids_dat2)
  )
bins[] <- 1:ncell(bins)

# Now we can create rasters for each year
rr <- raster::stack()
years <- sort(unique(ids_dat2@data$SURVEY_YEAR))
for(i in years) {
  message(i)
  s <- as(ids_dat2[ids_dat2$SURVEY_YEAR == i,], "SpatialPolygons")
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
g
ggsave(paste0("../", id, "-", nbin, "x", nbin, ".pdf"), width = 10, height = 9)

saveRDS(d,
  file = paste0("../", id, "-", "dataframe", "-", nbin, "x", nbin, ".rds"))

setwd(file.path("..", "..", ".."))
