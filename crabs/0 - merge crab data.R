
# bring in crab data - restrict to >= 2003
dat = readxl::read_excel("crabs/data/NWFSC_trawlSurvey_crabs.xlsx")
dat$year = as.numeric(substr(as.character(as.data.frame(dat[,"Date-Time (YYYYmmdd)"])[,1]), 1, 4))
dat = dat[dat$year >= 2003,]

allhauls = read.csv("crabs/data/TrawlHaulChars_since2002.csv")

allhauls = dplyr::rename(allhauls, Haul_ID = operation_dim.operation_id)
names(dat)[which(names(dat) == c("Haul ID"))] = "Haul_ID"
allhauls = left_join(allhauls, dat, by = "Haul_ID")

allhauls = allhauls[,c("o2_at_gear_ml_per_l_der", "haul_longitude_dim.longitude_in_degrees",
  "haul_latitude_dim.latitude_in_degrees","date_dim.yyyymmdd","date_dim.year", "area_swept_ha_der" ,
  "Haul_ID", "salinity_at_gear_psu_der", "seafloor_depth_m_der", "temperature_at_gear_c_der",
  "temperature_at_surface_c_der", "CPUE Kg/Ha Derived")]
allhauls = dplyr::rename(allhauls, o2_mlperL = o2_at_gear_ml_per_l_der, lon = haul_longitude_dim.longitude_in_degrees,
  lat = haul_latitude_dim.latitude_in_degrees, area_swept_ha = area_swept_ha_der,
  salinity = salinity_at_gear_psu_der, bottom_depth = seafloor_depth_m_der,
  bottom_tempc = temperature_at_gear_c_der, surface_tempc = temperature_at_surface_c_der,
  year = date_dim.year)
names(allhauls)[which(names(allhauls) == c("CPUE Kg/Ha Derived"))] = "cpue_kgperha"

allhauls[is.na(allhauls$area_swept_ha),"area_swept_ha"] = 0
allhauls[is.na(allhauls$cpue_kgperha),"cpue_kgperha"] = 0
allhauls$cpue_kgperha = as.numeric(allhauls$cpue_kgperha)
allhauls$area_swept_ha = as.numeric(allhauls$area_swept_ha)

allhauls = allhauls[allhauls$area_swept_ha > 0, ]
saveRDS(allhauls, "crabs/data/joined_crab.rds")
