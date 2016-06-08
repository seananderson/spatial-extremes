library("dplyr")

if (!file.exists("generated-data/rf-haul-2003-2014.rds")) {
  haul <- readxl::read_excel("data/RockfishHaulCatchIndivData2003To2014_20150826Fnl.xlsx",
    sheet = 2, skip = 8) %>%
    as.data.frame()
  names(haul) <- tolower(names(haul))
  names(haul) <- sub(" ", "_", names(haul))
  haul <- mutate(haul, yr = as.numeric(sub("Cycle ", "", survey_cycle))) %>%
    rename(
      haul_latitude_dd = `haul_latitude (decimal degrees)`,
      haul_longitude_dd = `haul_longitude (decimal degrees)`,
      haul_depth = `haul_depth (meters)`,
      trawl_duration = `trawl_duration (hours)`,
      area_swept = `area_swept by the net (hectares)`)
  saveRDS(haul, file = "generated-data/rf-haul-2003-2014.rds")
}

if (!file.exists("generated-data/rf-ind-2003-2014.rds")) {
  ind <- readxl::read_excel("data/RockfishHaulCatchIndivData2003To2014_20150826Fnl.xlsx",
    sheet = 3, skip = 7) %>%
    as.data.frame()
  names(ind) <- tolower(names(ind))
  ind <-  mutate(ind, yr = as.numeric(sub("Cycle ", "", project_cycle)))
  ind <- left_join(ind,
    select(haul, yr, haul_identifier, trawl_date, haul_latitude_dd,
      haul_longitude_dd, trawl_duration, area_swept))
  saveRDS(ind, file = "generated-data/rf-ind-2003-2014.rds")
}

haul <- readRDS("generated-data/rf-haul-2003-2014.rds")
ind <- readRDS("generated-data/rf-ind-2003-2014.rds")
