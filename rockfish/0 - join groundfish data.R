library(dplyr)
library(lubridate)

# load triennial survey data
tri = read.csv("rockfish/data/RockfishHaulCatchIndivDataTri1977To2004_20150827Fnl.csv")
tri_haul = read.csv("rockfish/data/RockfishHaulCatchIndivDataTri1977To2004_20150827Fnl_haul.csv")
# load annual survey data
annual = read.csv("rockfish/data/RockfishHaulCatchIndivData2003To2014_20150826Fnl.csv")
annual_haul = read.csv("rockfish/data/RockfishHaulCatchIndivData2003To2014_20150826Fnl_haul.csv")
annual_haul = dplyr::rename(annual_haul, HAUL_IDENTIFIER = Haul.Identifier, DATE = Trawl.Date)

# bring in area swept + date
annual = left_join(annual, annual_haul[,c("HAUL_IDENTIFIER","Area.Swept.by.the.Net..hectares.","DATE")])

# convert lat/lon to midpoint
tri$HAUL_LATITUDE_DD = (tri[,"END_LATITUDE"] + tri[,"START_LATITUDE"])/2
tri$HAUL_LONGITUDE_DD = (tri[,"END_LONGITUDE"] + tri[,"START_LONGITUDE"])/2
tri_haul$HAUL_LATITUDE_DD = (tri_haul[,"END_LATITUDE"] + tri_haul[,"START_LATITUDE"])/2
tri_haul$HAUL_LONGITUDE_DD = (tri_haul[,"END_LONGITUDE"] + tri_haul[,"START_LONGITUDE"])/2

# sex, age, year, depth, weight, length
annual$YEAR = as.numeric(substr(annual$PROJECT_CYCLE, 7, 10))
annual_haul$YEAR = as.numeric(substr(annual_haul$PROJECT_CYCLE, 7, 10))
annual = dplyr::rename(annual, AREA_SWEPT_HA = Area.Swept.by.the.Net..hectares.)
tri = dplyr::rename(tri, AGE_YRS = AGE, SEX = SEX_DETERMINATION, YEAR = CRUISE_YR,
  DEPTH_M = BOTTOM_DEPTH)
tri_haul = dplyr::rename(tri_haul, YEAR = CRUISE.YEAR)
tri$WEIGHT_KG = tri$FISH_WEIGHT_GRAMS / 1000 # g to kg
tri$LENGTH_CM = tri$LENGTH_MM / 10 ## mm to cm

# species
tri = dplyr::rename(tri, SPECIES = COMMON_NAME, SCIENTIFIC_NAME = SPECIES_NAME)
tri$SPECIES = tolower(tri$SPECIES)
annual$SPECIES = tolower(annual$SPECIES)

# date
annual$DATE = mdy(annual$DATE)
tri$DATE = mdy(unlist(lapply(strsplit(as.character(tri$START_TIME), " "), getElement, 1)))
tri_haul$DATE = mdy(unlist(lapply(strsplit(as.character(tri_haul$START_TIME), " "), getElement, 1)))

# Create output for individual species -- total densities
annual_haul = dplyr::rename(annual_haul, darkblotched.rockfish = Sebastes.crameri,
  HAUL_LATITUDE_DD = Haul.Latitude..decimal.degrees.,
  HAUL_LONGITUDE_DD = Haul.Longitude..decimal.degrees.,
  DEPTH_M = Haul.Depth..meters.,
  AREA_SWEPT = Area.Swept.by.the.Net..hectares.)


# Create output for age - length - weight data for all species
tri = tri[,c("YEAR","DATE","AREA_SWEPT_HA","HAUL_LATITUDE_DD",
  "HAUL_LONGITUDE_DD", "LENGTH_CM", "SEX", "WEIGHT_KG", "AGE_YRS", "DEPTH_M", "SPECIES")]
annual = annual[,c("YEAR","DATE","AREA_SWEPT_HA","HAUL_LATITUDE_DD",
  "HAUL_LONGITUDE_DD", "LENGTH_CM", "SEX", "WEIGHT_KG", "AGE_YRS", "DEPTH_M", "SPECIES")]

joined = bind_rows(tri,annual)
names(joined) = tolower(names(joined))
saveRDS(joined, "data/triennial_annual_joined.rds")


