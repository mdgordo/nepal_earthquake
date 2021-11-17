### data and documentation: https://microdata.worldbank.org/index.php/catalog/3705/data-dictionary
### report: https://elibrary.worldbank.org/doi/pdf/10.1596/33365 

library(raster)
library(sf)
library(tidyverse)
library(readstata13)

setwd("data")

### quake district designations: https://www.preventionweb.net/files/63268_18moderatelyaffecteddistrictsovervi.pdf
df.districts <- data.frame("district_name" = c("Gorkha", "Dhading", "Rasuwa", "Nuwakot", "Sindhupalchok", "Dolakha", "Ramechhap",
                                               "Kathmandu", "Bhaktapur", "Lalitpur", "Kabhrepalanchok", "Okhaldhunga", "Sindhuli", "Makwanpur",
                                               "Lamjung", "Tanahun", "Chitwan", "Solukhumbu", "Khotang", 
                                               "Kaski", "Parbat", "Syangja", "Palpa", "Gulmi", "Baglung",
                                               "Myagdi", "Arghakhanchi", "Nawalparasi", "Bhojpur", "Dhankuta", "Sankhuwasabha"),
                           "district" = c(36, 30, 29, 28, 23, 22, 21, 27, 26, 25, 24, 12, 20, 31, 37, 38, 35, 11, 13,
                                          40, 44, 39, 47, 46, 45, 43, 51, 48, 10, 7, 9), 
                           "dist_14" = c(rep(1, times = 14), rep(0, times = 17)),
                           "designation" = c(rep("severe", 7), rep("crisis", 7), rep("heavy", 5), rep("hit", 6), rep("slight", 6)))

df.districts$designation <- factor(df.districts$designation, levels = c("severe", "crisis", "heavy", "hit", "slight", "none"))

df.wards <- read_csv("NPL_2016-2018_HRVS_v01_M_STATA12/HRVS_gps.csv")
wards.sf <- st_as_sf(df.wards, coords = c("long", "lat"), crs = 4326)

epicenter <- st_sfc(st_point(x = c(84.731, 28.23)), crs = 4326)

df.shp <- st_read("shapefiles/NPL_adm", "NPL_adm3")

df.shp <- mutate(df.shp, district_name = if_else(NAME_3=="Chitawan", "Chitwan",
                                                 if_else(NAME_3=="Tanahu", "Tanahun",
                                                         if_else(NAME_3=="Kavrepalanchok","Kabhrepalanchok",as.character(NAME_3)))))

df.shp <- merge(df.shp, df.districts, by = "district_name", all = TRUE)
df.shp <- mutate(df.shp, dist_14 = if_else(is.na(dist_14),0,dist_14),
                 designation = if_else(is.na(designation),"none",as.character(designation)))
severe <- st_union(filter(df.shp, designation == "severe"))
crisis <- st_union(filter(df.shp, designation == "crisis"))
heavy <- st_union(filter(df.shp, designation == "heavy"))
dist_14 <- st_cast(st_union(severe, crisis), "LINESTRING")
placebo_border1 <- st_cast(st_union(filter(df.shp, district_name %in% c("Rasuwa", "Nuwakot", "Sindhupalchok"))), "LINESTRING")
placebo_border2 <- st_cast(st_union(st_union(severe, crisis), heavy), "LINESTRING")

### distance to epicenter
df.wards$distance_epicenter <- as.numeric(st_distance(wards.sf, epicenter))/1000

### creating border segments and finding closest for each ward
seg1 <- st_union(st_intersection(st_buffer(dist_14, 5), st_cast(st_union(df.shp$geometry[df.shp$district_name %in% c("Chitwan", "Manang", "Lamjung", "Tanahun")]), "MULTILINESTRING")))
seg2 <- st_union(st_intersection(st_buffer(dist_14, 5), st_cast(st_union(df.shp$geometry[df.shp$district_name %in% c("Makwanpur", "Sindhuli")]), "MULTILINESTRING")))
seg1 <- st_difference(seg1, st_buffer(seg2, 900))
seg3 <- st_union(st_intersection(st_buffer(dist_14, 5), st_cast(st_union(df.shp$geometry[df.shp$district_name %in% c("Solukhumbu", "Khotang")]), "MULTILINESTRING")))
border14_segments <- st_sf(segment = c("seg1", "seg2", "seg3"), geometry = c(seg1, seg2, seg3))

placebo1a <- st_union(st_intersection(st_buffer(placebo_border1, 5), st_cast(st_union(df.shp$geometry[df.shp$district_name=="Dhading"]), "MULTILINESTRING")))
placebo1b <- st_union(st_intersection(st_buffer(placebo_border1, 5), st_cast(st_union(df.shp$geometry[df.shp$district_name %in% c("Dolakha", "Kabhrepalanchok")]), "MULTILINESTRING")))
placebo1_segments <- st_sf(segment = c("placebo1a", "placebo1b"), geometry = c(placebo1a, placebo1b))

placebo2a <- st_union(st_intersection(st_buffer(placebo_border2, 5), st_cast(st_union(df.shp$geometry[df.shp$district_name %in% c("Lamjung", "Tanahun", "Chitwan")]), "MULTILINESTRING")))
placebo2b <- st_union(st_intersection(st_buffer(placebo_border2, 5), st_cast(st_union(df.shp$geometry[df.shp$district_name %in% c("Solukhumbu", "Khotang")]), "MULTILINESTRING")))
placebo2_segments <- st_sf(segment = c("placebo2a", "placebo2b"), geometry = c(placebo2a, placebo2b))

df.wards$border14_segment <- st_nearest_feature(wards.sf, border14_segments)
df.wards$border_segment13 <- st_nearest_feature(wards.sf, border14_segments[border14_segments$segment!="seg2",])
df.wards$placebo1_segment <- st_nearest_feature(wards.sf, placebo1_segments)
df.wards$placebo2_segment <- st_nearest_feature(wards.sf, placebo2_segments)
### distance to border 14 segments doesn't include distance to border with China
df.wards$dist_2_14 <- as.numeric(st_distance(wards.sf, st_union(border14_segments)))/1000
df.wards$dist_2_seg1 <- as.numeric(st_distance(wards.sf, seg1))/1000
df.wards$dist_2_seg13 <- as.numeric(st_distance(wards.sf, st_union(border14_segments[border14_segments$segment!="seg2",])))/1000
df.wards$dist_2_placebo1 <- as.numeric(st_distance(wards.sf, st_union(placebo1_segments)))/1000
df.wards$dist_2_placebo2 <- as.numeric(st_distance(wards.sf, st_union(placebo2_segments)))/1000

### Elevation data
raster_sf <- function(r){
  r <- projectRaster(r, crs = crs(wards.sf))
  r <- as.data.frame(r, xy=TRUE)
}

### These rasters have longitude as x and lat as y, hence the flip below
npl.alt <- getData("alt", country = "NPL", mask = TRUE, path = paste(getwd(), "/elevation", sep = ""))
npl.ter <- raster_sf(terrain(npl.alt, opt = c("slope", "aspect"), unit = "degrees"))
npl.alt <- raster_sf(npl.alt)

### USGS Shake Raster
npl.shk <- raster("USGS/raster/pga_mean.flt")
crs(npl.shk) <- crs(wards.sf)
npl.shk <- raster_sf(crop(npl.shk, y = df.shp))

ggplot() + geom_sf(data = df.shp) + geom_sf(data = border14_segments, aes(color = segment)) +
  geom_sf(data = wards.sf) +
  geom_raster(data = npl.alt, aes(x = y, y = x, fill = NPL_msk_alt), alpha = .5) +
  scale_fill_viridis_c() + theme_light()

npl.shk <- st_as_sf(npl.shk, coords = c("x", "y")) %>% filter(!is.na(pga_mean))
st_crs(npl.shk) <- st_crs(wards.sf)
npl.alt <- st_as_sf(npl.alt, coords = c("y", "x")) %>% filter(!is.na(NPL_msk_alt))
st_crs(npl.alt) <- st_crs(wards.sf)
npl.ter <- st_as_sf(npl.ter, coords = c("y", "x")) %>% filter(!is.na(slope))
st_crs(npl.ter) <- st_crs(wards.sf)

df.wards$shake_pga <- npl.shk$pga_mean[st_nearest_feature(wards.sf, npl.shk)]
df.wards$elevation <- npl.alt$NPL_msk_alt[st_nearest_feature(wards.sf, npl.alt)]
df.wards$slope <- npl.ter$slope[st_nearest_feature(wards.sf, npl.ter)]
df.wards$aspect <- npl.ter$aspect[st_nearest_feature(wards.sf, npl.ter)]

### Merge

df.wards <- merge(df.wards, df.districts, by.x = c("district", "distname"), by.y = c("district", "district_name"), all.x = TRUE)
df.wards <- mutate(df.wards, dist_14 = if_else(is.na(dist_14),0,dist_14),
                designation = if_else(is.na(designation),"none",as.character(designation))) %>%
            mutate(dist_2_14 = if_else(dist_14==1, dist_2_14, -1*dist_2_14),
                   dist_2_seg1 = if_else(dist_14==1, dist_2_seg1, -1*dist_2_seg1),
                   dist_2_seg13 = if_else(dist_14==1, dist_2_seg13, -1*dist_2_seg13),
                   dist_2_placebo1 = if_else(distname %in% c("Rasuwa", "Nuwakot", "Sindhupalchok"), dist_2_placebo1, -1*dist_2_placebo1),
                   dist_2_placebo2 = if_else(designation %in% c("severe", "crisis", "heavy"), dist_2_14, -1*dist_2_14))

rm(border14_segments, crisis, heavy, dist_14, epicenter, seg1, seg2, seg3, severe, placebo1_segments, placebo2_segments, placebo1a, placebo2a, placebo1b, placebo2b, placebo_border1, placebo_border2)

##### construct weights - email from Thomas Walker ######
df.wards <- df.wards %>% mutate(strata = if_else(district %in% c(69,70,74), 21,
                                           if_else(district %in% c(53,54,59:61), 22,
                                                   if_else(district %in% c(36:39,43,45:47), 23,
                                                           if_else(district %in% c(20,24,28,30,31), 24,
                                                                   if_else(district %in% c(3,7,10,12:14), 25,
                                                                           if_else(district==71, 31,
                                                                                   if_else(district %in% c(56,57), 32,
                                                                                           if_else(district %in% c(48,49), 33,
                                                                                                   if_else(district %in% c(17:19,33,34), 34,
                                                                                                           if_else(district %in% c(4:6,15), 35, 11)))))))))),
                          wt_hh = if_else(district %in% c(69,70,74), 547.5259259,
                                          if_else(district %in% c(53,54,59:61), 618.2705882,
                                                  if_else(district %in% c(36:39,43,45:47), 754.5958333,
                                                          if_else(district %in% c(20,24,28,30,31),487.4611111,
                                                                  if_else(district %in% c(3,7,10,12:14), 544.4632479,
                                                                          if_else(district==71, 634.0851852,
                                                                                  if_else(district %in% c(56,57), 666.9111111,
                                                                                          if_else(district %in% c(48,49), 631.6152381,
                                                                                                  if_else(district %in% c(17:19,33,34), 990.5319444,
                                                                                                          if_else(district %in% c(4:6,15), 898.6458333, 585.5533333)))))))))))

### Process survey data
setwd("NPL_2016-2018_HRVS_v01_M_STATA12")
folders <- c("Wave 1 - Household", "Wave 2 - Household", "Wave 3 - Household")

### function for finding respondents that answered "don't know"
finddks <- function(data, vars){
  df <- filter_at(data, vars(vars), all_vars(!is.na(.)))
  df <- filter_at(df, vars(vars), any_vars(. %in% c(998,999)))
  return(df$hhid)
}

timeconv <- function(x){
  t <- as.numeric(as.numeric(unlist(strsplit(x, "-")))%*%c(24, 1, 1/60))
  return(t)
}

for (i in c(1:3)) {
      setwd(folders[i])

      df.hh <- read.dta13("Section_0.dta")
      
      if (folders[i]=="Wave 1 - Household") {
        df.hh$split_hh <- 0
        df.hh$split_from <- NA
        df.hh$s00q00b <- NA
      }
      
      if (folders[i]=="Wave 2 - Household") {
        df.hh <- rename(df.hh, s00q01 = district_code,
                        s00q02a = vdc_code,
                        s00q02b = vdc,
                        s00q03 = ward) %>%
          mutate(vdc = s00q02b)
      }
      
      if (folders[i]=="Wave 3 - Household") {
        df.hh$district77 <- NULL
        df.hh$s00q01a <- NULL
        df.hh$s00q03a <- NULL
        df.hh$s00q03b <- NULL
        df.hh$s00q03c <- NULL
      }
      
      df.hh <- df.hh %>%
        select(!c(s00q00b, s00q01, s00q02b, psu)) %>%
        rename(vdc_code = s00q02a,
               ward = s00q03,
               interviewer = s00q07,
               hh_head = s00q09,
               ethnicity = s00q15,
               language = s00q16,
               religion = s00q17) %>%
        ## caste sources: https://en.wikipedia.org/wiki/Ethnic_groups_in_Nepal and https://dhsprogram.com/pubs/pdf/FA58/FA58.pdf 
        mutate(high_caste = if_else(ethnicity %in% c(12,13,14,19,48,71,81,84,95), 1, 0),
               caste_recode = if_else(ethnicity %in% c(12,13,14,19,48,81,84,95), "Brahman/Chhetri", 
                                      if_else(ethnicity==71, "Newar", 
                                              if_else(ethnicity %in% c(3,9,16,20,22,28,30,32,33,38,46,50,69,85,92,102), "Dalit", "Other"))),
               region = if_else(district %in% c(75,68,67,66,65,64,63,62,42,41,29,23,22,11,9,1),"Mountain",
                                if_else(district %in% c(72,71,58,59,57,56,50,49,48,35,34,31,33,32,20,19,18,17,16,15,14,6,5,4),"Terai","Hill")),
               split_from = if_else(is.na(split_from), hhid, as.integer(split_from)))
      
      ### roster
      df.hrvs_1 <- read.dta13("Section_1.dta")
      df.roster <- df.hrvs_1 %>%
        mutate(s01q03 = if_else(is.na(s01q03), s01q03a/12, s01q03)) %>%
        mutate(under5 = if_else(s01q03 < 6, 1, 0),
               under18 = if_else(s01q03 < 19, 1, 0),
               over65 = if_else(s01q03 > 64, 1, 0),
               left_fam = if_else(s01q00a==2 & !is.na(s01q00a), 1, 0),
               born_away = if_else(s01q05a==2, 1, 0),
               born_away_head = if_else(s01q05a==2 & s01q01==1, 1, 0),
               femalehh = if_else(s01q02!=2 & s01q01==1, 1, 0),
               currentlyliving = if_else(s01q06a==1, 1, 0),
               extended_fam = if_else(s01q01 %in% c(1, 2, 3), 0, 1),
               age_hh = if_else(s01q01==1,s01q03,0)) %>%
        filter(left_fam==0) %>%
        group_by(hhid) %>%
        summarize(hhmembers = n(),
                  femalehh = sum(femalehh, na.rm = TRUE),
                  age_hh = sum(age_hh, na.rm = TRUE),
                  under5 = sum(under5, na.rm = TRUE),
                  under18 = sum(under18, na.rm = TRUE),
                  over65 = sum(over65, na.rm = TRUE),
                  currentlyliving = sum(currentlyliving, na.rm = TRUE),
                  extended_members = sum(extended_fam, na.rm = TRUE),
                  born_elsewhere = sum(born_away, na.rm = TRUE),
                  born_away_head = sum(born_away_head, na.rm = TRUE))
      df.hh <- merge(df.hh, df.roster, by = "hhid", all = TRUE)
      rm(df.roster)
      
      ### education
      df.hrvs_2 <- read.dta13("Section_2.dta")
      df.education <- df.hrvs_2 %>% 
              mutate(highest_ed = if_else(is.na(s02q02) | s02q02 %in% c("16","17"), 0, s02q02),
                     school_costs = rowSums(.[,c("s02q09a","s02q09b","s02q09c","s02q09d","s02q09e","s02q09f")], na.rm = TRUE),
                     attending = if_else(s02q01==3, 1, 0)) %>%
              group_by(hhid) %>%
              summarize(highest_ed = max(highest_ed),
                        school_costs = sum(school_costs, na.rm = TRUE),
                        scholarships = sum(s02q08, na.rm = TRUE),
                        children_attending = sum(attending, na.rm = TRUE)) %>%
              mutate(class5 = if_else(highest_ed<6, 1, 0),
                     class10 = if_else(highest_ed<11, 1, 0))
      df.hh <- merge(df.hh, df.education, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, school_costs = if_else(is.na(school_costs),0,school_costs),
                      scholarships = if_else(is.na(scholarships),0,scholarships),
                      children_attending = if_else(is.na(children_attending),0,children_attending))
      df.hh$school_costs[df.hh$hhid %in% finddks(df.hrvs_2, c("s02q09a","s02q09b","s02q09c","s02q09d","s02q09e","s02q09f"))] <- NA
      df.hh$scholarships[df.hh$hhid %in% finddks(df.hrvs_2, c("s02q08"))] <- NA
      rm(df.education, df.hrvs_2)
      
      ### health      
      df.hrvs_3 <- read.dta13("Section_3.dta")
      df.health <- df.hrvs_3 %>% 
              mutate(health_costs = rowSums(.[,c("s03q06a","s03q06b","s03q06c","s03q06d","s03q06e")], na.rm = TRUE),
                     days_ill = rowSums(.[c("s03q03_1", "s03q03_2", "s03q03_3", "s03q03_4", "s03q03_5", "s03q03_6", 
                                            "s03q03_7", "s03q03_8", "s03q03_9", "s03q03_10", "s03q03_11", "s03q03_12")], na.rm = TRUE)) %>%
              group_by(hhid) %>%
              summarize(health_costs = sum(health_costs, na.rm = TRUE),
                        days_ill = sum(days_ill, na.rm = TRUE))
      df.hh <- merge(df.hh, df.health, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, health_costs = if_else(is.na(health_costs),0,health_costs),
                      days_ill = if_else(is.na(days_ill),0,days_ill))
      df.hh$health_costs[df.hh$hhid %in% finddks(df.hrvs_3, c("s03q06a","s03q06b","s03q06c","s03q06d","s03q06e"))] <- NA
      df.hh$health_costs[df.hh$hhid %in% finddks(df.hrvs_3, c("s03q03_1", "s03q03_2", "s03q03_3", "s03q03_4", "s03q03_5", "s03q03_6", 
                                                              "s03q03_7", "s03q03_8", "s03q03_9", "s03q03_10", "s03q03_11", "s03q03_12"))] <- NA
      rm(df.health, df.hrvs_3)
      
      ### housing
      df.hrvs_4 <- read.dta13("Section_4.dta")
      df.hrvs_4$time_to_market <- sapply(df.hrvs_4$s04q28c, timeconv)
      df.hrvs_4$time_to_bank <- sapply(df.hrvs_4$s04q29c, timeconv)
      df.hrvs_4$time_to_school <- sapply(df.hrvs_4$s04q32c, timeconv)
      df.hrvs_4$time_to_health <- sapply(df.hrvs_4$s04q34c, timeconv)
      df.housing <- df.hrvs_4 %>% 
              mutate(utilities = rowSums(.[,c("s04q21", "s04q23")], na.rm = TRUE) + 12*s04q24b,
                     steep_slope = if_else(s04q01==3, 1, 0),
                     unimproved_roof = if_else(s04q17 %in% c(1, 2, 8), 1, 0),
                     unimproved_walls = if_else(s04q15 %in% c(2, 4, 6, 7), 1, 0)) %>%
              group_by(hhid) %>%
              summarize(always_lived_house = sum(s04q09a==1, na.rm = TRUE),
                        always_lived_vdc = sum(s04q09b!=2, na.rm = TRUE),
                        always_lived_dist = sum(s04q09a!=2, na.rm = TRUE),
                        years_ago_built = sum(s04q18, na.rm = TRUE),
                        rent_earned = sum(s04q04b, na.rm = TRUE),
                        rent_paid = sum(s04q06, na.rm = TRUE),
                        home_value = sum(s04q19, na.rm = TRUE),
                        utilities_paid = sum(utilities, na.rm = TRUE),
                        steep_slope = sum(steep_slope, na.rm = TRUE), 
                        walls = sum(s04q15, na.rm = TRUE),
                        temp_earth_housing = if_else(unimproved_roof==1 | unimproved_walls==1, 1, 0), 
                        foundation = sum(s04q16, na.rm = TRUE),
                        roof = sum(s04q17, na.rm = TRUE),
                        fuel = sum(s04q25, na.rm = TRUE),
                        stove = sum(s04q26, na.rm = TRUE),
                        time_to_market = sum(time_to_market, na.rm = TRUE),
                        time_to_bank = sum(time_to_bank, na.rm = TRUE),
                        time_to_school = sum(time_to_school, na.rm = TRUE),
                        time_to_health = sum(time_to_health, na.rm = TRUE))                    
      df.hh <- merge(df.hh, df.housing, by = "hhid", all = TRUE)
      df.hh$rent_earned[df.hh$hhid %in% finddks(df.hrvs_4, c("s04q04b"))] <- NA
      df.hh$rent_paid[df.hh$hhid %in% finddks(df.hrvs_4, c("s04q06"))] <- NA
      df.hh$home_value[df.hh$hhid %in% finddks(df.hrvs_4, c("s04q19"))] <- NA
      df.hh$utilities_paid[df.hh$hhid %in% finddks(df.hrvs_4, c("s04q21", "s04q23", "s04q24b"))] <- NA
      rm(df.housing, df.hrvs_4)
      
      ### food production
      df.hrvs_5a <- read.dta13("Section_5a.dta")
      
      df.food <- df.hrvs_5a %>% 
              mutate(kgunits = ifelse(s05q05==2, s05q04/1000, ifelse(s05q05 %in% c(1,4), s05q04, NA)),
                     chicken_price = ifelse(foodid=="Chicken", s05q06/kgunits, NA),
                     rice_price = ifelse(foodid=="Rice", s05q06/kgunits, NA),
                     lentil_price = ifelse(foodid=="Lentil (Black gram)", s05q06/kgunits, NA),
                     mutton_price = ifelse(foodid=="Mutton", s05q06/kgunits, NA),
                     sugar_price = ifelse(foodid=="Sugar", s05q06/kgunits, NA),
                     potato_price = ifelse(foodid=="Potatoes", s05q06/kgunits, NA)) %>%
              group_by(hhid) %>%
              summarize(food_home_production = sum(s05q03, na.rm = TRUE),
                        food_market = sum(s05q06, na.rm = TRUE),
                        food_inkind = sum(s05q09, na.rm = TRUE),
                        chicken_price = mean(chicken_price, na.rm = TRUE),
                        rice_price = mean(rice_price, na.rm = TRUE),
                        lentil_price = mean(lentil_price, na.rm = TRUE),
                        mutton_price = mean(mutton_price, na.rm = TRUE),
                        sugar_price = mean(sugar_price, na.rm = TRUE),
                        potato_price = mean(potato_price, na.rm = TRUE))
      ### Trim miscoded outliers
      df.food <- mutate(df.food, lentil_price = if_else(lentil_price>=10000, as.numeric(NA), lentil_price),
                        chicken_price = if_else(chicken_price>=10000, as.numeric(NA), chicken_price),
                        mutton_price = if_else(mutton_price>=10000, as.numeric(NA), mutton_price),
                        sugar_price = if_else(sugar_price>=10000, as.numeric(NA), sugar_price))
      df.hh <- merge(df.hh, df.food, by = "hhid", all = TRUE)
      df.hh$food_home_production[df.hh$hhid %in% finddks(df.hrvs_5a, c("s05q03"))] <- NA
      df.hh$food_market[df.hh$hhid %in% finddks(df.hrvs_5a, c("s05q06"))] <- NA
      df.hh$food_inkind[df.hh$hhid %in% finddks(df.hrvs_5a, c("s05q09"))] <- NA
      rm(df.food, df.hrvs_5a)
      
      ### For wave 1 s05q10 1 = yes, 2 = no; for all waves s05q16 1 = never - confirmed in email from Jui Shrestha, documentation is wrong
      df.hrvs_5b <- read.dta13("Section_5b_corrected.dta")
      df.hunger <- df.hrvs_5b %>%
        mutate(felt_hunger = if_else(s05q10==1, 1, 0),
               skipped_meal = if_else(s05q16!=1, 1, 0)) %>%
        group_by(hhid) %>%
        summarize(felt_hunger = sum(felt_hunger, na.rm = TRUE),
                  skipped_meal = sum(skipped_meal, na.rm = TRUE))
      df.hh <- merge(df.hh, df.hunger, by = "hhid", all = TRUE)
      rm(df.hunger, df.hrvs_5b)
      
      ### durables
      df.hrvs_6a <- read.dta13("Section_6a.dta")
      df.nonfood <- df.hrvs_6a %>%
              mutate(energy = ifelse(nonfoodid %in% c("Wood (bundle wood, logwood, sawdust)", "Kerosene oil", "Coal, charcoal", 
                                                      "Cylinder gas (LPG)", "Matches, candles, lighters., lanterns, etc"), s06q01b, 0),
                     transportation = ifelse(nonfoodid %in% c("Public transportation (buses, taxis, rickshaws, train tickets", 
                                                              "Petrol, diesel, motor oil (for personal vehicle only)"), s06q01b, 0),
                     clothing_cleaning_home = ifelse(nonfoodid %in% c("Ready-made clothing and apparel", "Personal care items like shampoo, cosmetics, soap",
                                                                 "Shoes, slippers, sandals, etc.", "Dry cleaning and washing expenses",
                                                                 "Personal services (haircuts, shaving, shoeshine)", "Light bulbs, shades, batteries, etc",
                                                                 "Household cleaning articles (soap, bleach, washing powder)", 
                                                                 "Wages paid to watchman, servant, gardener, driver, etc"), s06q01b, 0),
                     entertainment_other = ifelse(nonfoodid %in% c("Other frequent expenses not mentioned", "Newspapers, books, stationery supplies(except educational expenses)",
                                                              "Pocket money to children", "Entertainment (cinema, CD/cassette rentals, etc.)"), s06q01b, 0)) %>%
              group_by(hhid) %>%
              summarize(energy = sum(energy, na.rm = TRUE),
                        transportation = sum(transportation, na.rm = TRUE),
                        clothing_cleaning_home = sum(clothing_cleaning_home, na.rm = TRUE),
                        entertainment_other = sum(entertainment_other, na.rm = TRUE))
      df.hh <- merge(df.hh, df.nonfood, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, energy = if_else(is.na(energy),0,energy),
                      transportation = if_else(is.na(transportation),0,transportation),
                      clothing_cleaning_home = if_else(is.na(clothing_cleaning_home),0,clothing_cleaning_home),
                      entertainment_other = if_else(is.na(entertainment_other),0,entertainment_other))
      df.hh$energy[df.hh$hhid %in% finddks(df.hrvs_6a, c("s06q01b"))] <- NA
      df.hh$transportation[df.hh$hhid %in% finddks(df.hrvs_6a, c("s06q01b"))] <- NA
      df.hh$clothing_cleaning_home[df.hh$hhid %in% finddks(df.hrvs_6a, c("s06q01b"))] <- NA
      df.hh$entertainment_other[df.hh$hhid %in% finddks(df.hrvs_6a, c("s06q01b"))] <- NA
      rm(df.nonfood, df.hrvs_6a)
      
      df.hrvs_6b <- read.dta13("Section_6b.dta")
      df.nonfood2 <- df.hrvs_6b %>%
              mutate(taxes = ifelse(nonfoodid=="Income taxes, land taxes, housing and property taxes", s06q02, 0),
                     ceremonies = ifelse(nonfoodid %in% c("Expenditure on religious ceremonies", "Marriages, births, and other ceremonies", 
                                                          "Funeral and death related expenses"), s06q02, 0),
                     durables = ifelse(nonfoodid %in% c("Expenditure on religious ceremonies", "Marriages, births, and other ceremonies", 
                                                        "Funeral and death related expenses", "Income taxes, land taxes, housing and property taxes", 
                                                        "Postal expenses, telegrams, fax, telephone", "Legal expenses and insurance (life, car, etc", 
                                                        "Excursion, holiday, (including travel and lodging)"), 0, s06q02),
                     other_expenses = ifelse(nonfoodid %in% c("Postal expenses, telegrams, fax, telephone", "Legal expenses and insurance (life, car, etc", 
                                                             "Excursion, holiday, (including travel and lodging)"), s06q02, 0)) %>%
              group_by(hhid) %>%
              summarize(durables_consumption = sum(durables, na.rm = TRUE),
                        ceremonial_expenses = sum(ceremonies, na.rm = TRUE),
                        taxes = sum(taxes, na.rm = TRUE),
                        other_expenses = sum(other_expenses, na.rm = TRUE))
      df.hh <- merge(df.hh, df.nonfood2, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, durables_consumption = if_else(is.na(durables_consumption),0,durables_consumption),
                      ceremonial_expenses = if_else(is.na(ceremonial_expenses),0,ceremonial_expenses),
                      taxes = if_else(is.na(taxes),0,taxes),
                      other_expenses = if_else(is.na(other_expenses),0,other_expenses))
      df.hh$durables_consumption[df.hh$hhid %in% finddks(df.hrvs_6b, c("s06q02"))] <- NA
      df.hh$ceremonial_expenses[df.hh$hhid %in% finddks(df.hrvs_6b, c("s06q02"))] <- NA
      df.hh$taxes[df.hh$hhid %in% finddks(df.hrvs_6b, c("s06q02"))] <- NA
      df.hh$other_expenses[df.hh$hhid %in% finddks(df.hrvs_6b, c("s06q02"))] <- NA
      rm(df.nonfood2, df.hrvs_6b)
      
      df.hrvs_6c <- read.dta13("Section_6c.dta")
      df.durables <- df.hrvs_6c %>%
              group_by(hhid) %>%
              summarize(durables_stock = sum(s06q03b, na.rm = TRUE))
      df.hh <- merge(df.hh, df.durables, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, durables_stock = if_else(is.na(durables_stock),0,durables_stock))
      df.hh$durables_stock[df.hh$hhid %in% finddks(df.hrvs_6c, c("s06q03b"))] <- NA
      rm(df.durables, df.hrvs_6c)
      
      df.hrvs_6d <- read.dta13("Section_6d.dta")
      df.baskets <- df.hrvs_6d %>%
              group_by(hhid) %>%
              summarize(craft_consumption = sum(s06q04c, na.rm = TRUE))
      df.hh <- merge(df.hh, df.baskets, by = "hhid")
      df.hh$craft_consumption[df.hh$hhid %in% finddks(df.hrvs_6d, c("s06q04c"))] <- NA
      rm(df.baskets, df.hrvs_6d)
      
      ### wage income and self employment
      df.hrvs_7 <- read.dta13("Section_7.dta")
      df.hrvs_8 <- read.dta13("Section_8.dta")
      ### no match to age data for ~ 600 individuals out of 14,000
      df.child_labor = merge(df.hrvs_7, df.hrvs_1[,c("district", "vdc", "psu", "hhid", "member_id", "s01q03")], by = c("district", "vdc", "psu", "hhid", "member_id"), all.x = TRUE) %>%
              mutate(child_labor = if_else(s01q03<18, 1, 0),
                     child_wage_labor = if_else(s01q03<18 & s07q06 %in% c(1,2), 1, 0),
                     child_days_worked = if_else(child_labor==1, rowSums(.[,c("s07q04_1", "s07q04_2", "s07q04_3","s07q04_4","s07q04_5","s07q04_6","s07q04_7","s07q04_8","s07q04_9",
                                                                              "s07q04_10","s07q04_11","s07q04_12")], na.rm = TRUE), 0),
                     child_wage_days = if_else(child_wage_labor==1, rowSums(.[,c("s07q04_1", "s07q04_2", "s07q04_3","s07q04_4","s07q04_5","s07q04_6","s07q04_7","s07q04_8","s07q04_9",
                                                                            "s07q04_10","s07q04_11","s07q04_12")], na.rm = TRUE), 0)) %>%
              group_by(hhid) %>%
              summarize(child_labor = sum(child_labor, na.rm = TRUE),
                        child_wage_labor = sum(child_wage_labor, na.rm = TRUE),
                        child_days_worked = sum(child_days_worked, na.rm = TRUE),
                        child_wage_days = sum(child_wage_days, na.rm = TRUE)) %>%
              mutate(child_labor = if_else(child_labor>0,1,0),
                     child_wage_labor = if_else(child_wage_labor>0,1,0))
      df.hh <- merge(df.hh, df.child_labor, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, child_labor = if_else(is.na(child_labor),0,child_labor),
                      child_wage_labor = if_else(is.na(child_wage_labor),0,child_wage_labor),
                      child_days_worked = if_else(is.na(child_days_worked),0,child_days_worked),
                      child_wage_days = if_else(is.na(child_wage_days),0,child_wage_days))
      
      df.wages <- filter(df.hrvs_7, s07q06 %in% c(1,2)) %>% select(!joblist)
      ### the flow is tricky here - about 200 jobs no match b/t section 7 and 8 - get annual wage of 0 but positive days worked
      df.wages <- merge(df.hrvs_8, df.wages, by.x = c("district", "vdc", "psu", "hhid", "member_id", "wagejobid"), 
                         by.y = c("district", "vdc", "psu", "hhid", "member_id", "jobid"), all = TRUE) %>%
              mutate(daysworked = rowSums(.[,c("s07q04_1", "s07q04_2", "s07q04_3","s07q04_4","s07q04_5","s07q04_6","s07q04_7","s07q04_8","s07q04_9",
                                               "s07q04_10","s07q04_11","s07q04_12")], na.rm = TRUE),
                     monthsworked = rowSums(.[,c("s07q03_1", "s07q03_2", "s07q03_3","s07q03_4","s07q03_5","s07q03_6","s07q03_7","s07q03_8","s07q03_9",
                                               "s07q03_10","s07q03_11","s07q03_12")], na.rm = TRUE),
                     annual_wage = if_else(s08q05==1, s08q06*daysworked,
                                           if_else(s08q05==3, s08q13,
                                                   if_else(s08q05==2, 
                                                           if_else(s08q09==1, s08q10, monthsworked*(s08q12a + s08q12b + s08q12c + s08q12d + s08q12e)), as.numeric(NA))))) %>%
              group_by(hhid) %>%
              summarize(wage_labor_days = sum(daysworked, na.rm = TRUE),
                        annual_wages = sum(annual_wage, na.rm = TRUE),
                        wage_workers = as.numeric(n()))
      df.hh <- merge(df.hh, df.wages, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, wage_labor_days = if_else(is.na(wage_labor_days),0,wage_labor_days),
                      annual_wages = if_else(is.na(annual_wages),0,annual_wages),
                      wage_workers = if_else(is.na(wage_workers),0,wage_workers))
      df.hh$wage_labor_days[df.hh$hhid %in% finddks(df.hrvs_7, c("s07q04_1", "s07q04_2", "s07q04_3","s07q04_4","s07q04_5","s07q04_6","s07q04_7","s07q04_8","s07q04_9",
                                                                 "s07q04_10","s07q04_11","s07q04_12"))] <- NA
      df.hh$annual_wages[df.hh$hhid %in% finddks(df.hrvs_8, c("s08q06", "s08q13", "s08q10", "s08q12a", "s08q12b", "s08q12c", "s08q12d", "s08q12e"))] <- NA
      
      df.self <- filter(df.hrvs_7, s07q06 %in% c(3,4)) %>%
              mutate(daysworked = rowSums(.[,c("s07q04_1", "s07q04_2", "s07q04_3","s07q04_4","s07q04_5","s07q04_6","s07q04_7","s07q04_8","s07q04_9",
                                               "s07q04_10","s07q04_11","s07q04_12")], na.rm = TRUE)) %>%
              group_by(hhid) %>%
              summarize(self_emp_days = sum(daysworked, na.rm = TRUE),
                        self_emp_workers = as.numeric(n()))
      df.hh <- merge(df.hh, df.self, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, self_emp_days = if_else(is.na(self_emp_days),0,self_emp_days),
                      self_emp_workers = if_else(is.na(self_emp_workers),0,self_emp_workers))
      df.hh$self_emp_days[df.hh$hhid %in% finddks(df.hrvs_7, c("s07q04_1", "s07q04_2", "s07q04_3","s07q04_4","s07q04_5","s07q04_6","s07q04_7","s07q04_8","s07q04_9",
                                                                 "s07q04_10","s07q04_11","s07q04_12"))] <- NA
      rm(df.wages, df.hrvs_7, df.hrvs_8, df.self, df.hrvs_1)
      
      ### Farming
      df.hrvs_9a1 <- read.dta13("Section_9a1.dta")
      df.agrents <- df.hrvs_9a1 %>% 
              mutate(landrents = rowSums(.[,c("s09q08a","s09q08b")], na.rm = TRUE) + rowSums(.[,c("s09q12a","s09q12b")], na.rm = TRUE)) %>%
              group_by(hhid) %>%
              summarize(landvalue = sum(s09q06, na.rm = TRUE),
                        landrents = sum(landrents, na.rm = TRUE),
                        n_plots = as.numeric(n()),
                        plot_area = sum(area_sqm, na.rm = TRUE))
      df.hh <- merge(df.hh, df.agrents, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, landowner = if_else(is.na(landvalue), 0, 1), 
                      landvalue = if_else(is.na(landvalue),0,landvalue),
                      landrents = if_else(is.na(landrents),0,landrents),
                      n_plots = if_else(is.na(n_plots),0,n_plots),
                      plot_area = if_else(is.na(plot_area),0,plot_area))
      df.hh$landvalue[df.hh$hhid %in% finddks(df.hrvs_9a1, c("s09q06"))] <- NA
      df.hh$landrents[df.hh$hhid %in% finddks(df.hrvs_9a1, c("s09q08a","s09q08b", "s09q12a","s09q12b"))] <- NA
      df.hh$plot_area[df.hh$hhid %in% finddks(df.hrvs_9a1, c("area_sqm"))] <- NA
      rm(df.agrents, df.hrvs_9a1)
      
      df.hrvs_9a2 <- read.dta13("Section_9a2.dta")
      df.agrentspaid <- df.hrvs_9a2 %>%
              group_by(hhid) %>%
              summarize(landrent_paid_cash = sum(s09q18, na.rm = TRUE),
                        landrent_paid_inkind = sum(s09q20, na.rm = TRUE))
      df.hh <- merge(df.hh, df.agrentspaid, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, landrent_paid_cash = if_else(is.na(landrent_paid_cash),0,landrent_paid_cash),
                      landrent_paid_inkind = if_else(is.na(landrent_paid_inkind),0,landrent_paid_inkind))
      df.hh$landrent_paid_cash[df.hh$hhid %in% finddks(df.hrvs_9a2, c("s09q18"))] <- NA
      df.hh$landrent_paid_inkind[df.hh$hhid %in% finddks(df.hrvs_9a2, c("s09q20"))] <- NA
      rm(df.agrentspaid, df.hrvs_9a2)
      
      df.hrvs_9a3 <- read.dta13("Section_9a3.dta")
      df.landassets <- df.hrvs_9a3 %>% 
              mutate(freetransfers = if_else(is.na(s09q31b),0,s09q31b) - if_else(is.na(s09q28b),0,s09q28b)) %>%
              group_by(hhid) %>% 
              summarize(land_sales = sum(s09q28a, na.rm = TRUE),
                        land_sales_area = sum(area_sqm_sold, na.rm = TRUE),
                        land_purchased = sum(s09q31a, na.rm = TRUE),
                        land_purchased_area = sum(area_sqm_bought, na.rm = TRUE),
                        land_net_free_transfers = sum(freetransfers, na.rm = TRUE))
      df.hh <- merge(df.hh, df.landassets, by = "hhid", all = TRUE)
      df.hh$land_sales[df.hh$hhid %in% finddks(df.hrvs_9a3, c("s09q28a"))] <- NA
      df.hh$land_sales_area[df.hh$hhid %in% finddks(df.hrvs_9a3, c("area_sqm_sold"))] <- NA
      df.hh$land_purchased[df.hh$hhid %in% finddks(df.hrvs_9a3, c("s09q31a"))] <- NA
      df.hh$land_purchased_area[df.hh$hhid %in% finddks(df.hrvs_9a3, c("area_sqm_bought"))] <- NA
      df.hh$land_net_free_transfers[df.hh$hhid %in% finddks(df.hrvs_9a3, c("s09q31b", "s09q28b"))] <- NA
      rm(df.landassets, df.hrvs_9a3)
      
      df.hrvs_9b1 <- read.dta13("Section_9b1.dta")
      df.wetseason_sales <- df.hrvs_9b1 %>% 
              mutate(sales = s09q41e*s09q42) %>%
              group_by(hhid) %>%
              summarize(wet_ag_sales = sum(sales, na.rm = TRUE))
      df.hh <- merge(df.hh, df.wetseason_sales, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, wet_ag_sales = if_else(is.na(wet_ag_sales),0,wet_ag_sales))
      df.hh$wet_ag_sales[df.hh$hhid %in% finddks(df.hrvs_9b1, c("s09q41e", "s09q42"))] <- NA
      rm(df.wetseason_sales, df.hrvs_9b1)
      
      df.hrvs_9b2 <- read.dta13("Section_9b2.dta")
      df.dryseason_sales <- df.hrvs_9b2 %>% 
              mutate(sales = s09q50e*s09q51) %>%
              group_by(hhid) %>%
              summarize(dry_ag_sales = sum(sales, na.rm = TRUE))
      df.hh <- merge(df.hh, df.dryseason_sales, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, dry_ag_sales = if_else(is.na(dry_ag_sales),0,dry_ag_sales))
      df.hh$dry_ag_sales[df.hh$hhid %in% finddks(df.hrvs_9b2, c("s09q50e", "s09q51"))] <- NA
      rm(df.dryseason_sales, df.hrvs_9b2)
      
      df.hrvs_9c <- read.dta13("Section_9c.dta")
      df.aginputcosts <- df.hrvs_9c %>% 
              mutate(agcosts = rowSums(.[,c("s09q52b","s09q52d","s09q52f","s09q52h","s09q52j",
                                            "s09q53b","s09q53d","s09q53f","s09q53h","s09q53j",
                                            "s09q55b","s09q55d","s09q55f","s09q55h","s09q55j",
                                            "s09q55l","s09q55n","s09q55p","s09q55r")], na.rm = TRUE),
                     equip_rental_income = rowSums(.[,c("s09q54b","s09q54d","s09q54f")], na.rm = TRUE)) %>%
        group_by(hhid) %>%
              summarize(ag_costs = sum(agcosts, na.rm = TRUE),
                        equip_rental_income = sum(equip_rental_income, na.rm = TRUE))
      df.hh <- merge(df.hh, df.aginputcosts, by = "hhid", all = TRUE)
      df.hh$ag_costs[df.hh$hhid %in% finddks(df.hrvs_9c, c("s09q52b","s09q52d","s09q52f","s09q52h","s09q52j",
                                                            "s09q53b","s09q53d","s09q53f","s09q53h","s09q53j",
                                                            "s09q55b","s09q55d","s09q55f","s09q55h","s09q55j",
                                                            "s09q55l","s09q55n","s09q55p","s09q55r"))] <- NA
      df.hh$equip_rental_income[df.hh$hhid %in% finddks(df.hrvs_9c, c("s09q54b","s09q54d","s09q54f"))] <- NA
      rm(df.aginputcosts, df.hrvs_9c)
      
      df.hrvs_9d <- read.dta13("Section_9d.dta")
      df.livestock <- df.hrvs_9d %>%
              group_by(hhid) %>%
              summarize(livestock_value = sum(s09q57b, na.rm = TRUE),
                        livestock_sales = sum(s09q60b, na.rm = TRUE),
                        livestock_purchases = sum(s09q61b, na.rm = TRUE))
      df.hh <- merge(df.hh, df.livestock, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, livestock_value = if_else(is.na(livestock_value),0,livestock_value),
                      livestock_sales = if_else(is.na(livestock_sales),0,livestock_sales),
                      livestock_purchases = if_else(is.na(livestock_purchases),0,livestock_purchases))
      df.hh$livestock_value[df.hh$hhid %in% finddks(df.hrvs_9d, c("s09q57b"))] <- NA
      df.hh$livestock_sales[df.hh$hhid %in% finddks(df.hrvs_9d, c("s09q60b"))] <- NA
      df.hh$livestock_purchases[df.hh$hhid %in% finddks(df.hrvs_9d, c("s09q61b"))] <- NA
      rm(df.livestock, df.hrvs_9d)
      
      df.hrvs_9e <- read.dta13("Section_9e.dta")
      df.livestock_income <- df.hrvs_9e %>% 
        mutate(livestock_inc = rowSums(.[,c("s09q62a","s09q62b","s09q62c","s09q62d","s09q62e",
                                            "s09q62f","s09q62g","s09q62h","s09q62i")], na.rm = TRUE),
               livestock_costs = rowSums(.[,c("s09q63a","s09q63b","s09q63c","s09q63d")], na.rm = TRUE)) %>%
      group_by(hhid) %>%
              summarize(livestock_income = sum(livestock_inc, na.rm = TRUE),
                        livestock_costs = sum(livestock_costs, na.rm = TRUE))
      df.hh <- merge(df.hh, df.livestock_income, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, livestock_income = if_else(is.na(livestock_income),0,livestock_income),
                      livestock_costs = if_else(is.na(livestock_costs),0,livestock_costs))
      df.hh$livestock_income[df.hh$hhid %in% finddks(df.hrvs_9e, c("s09q62a","s09q62b","s09q62c","s09q62d","s09q62e",
                                                                   "s09q62f","s09q62g","s09q62h","s09q62i"))] <- NA
      df.hh$livestock_costs[df.hh$hhid %in% finddks(df.hrvs_9e, c("s09q63a","s09q63b","s09q63c","s09q63d"))] <- NA
      rm(df.livestock_income, df.hrvs_9e)
      
      df.hrvs_9f <- read.dta13("Section_9f.dta")
      df.equipment <- df.hrvs_9f %>% 
              group_by(hhid) %>%
              summarize(equip_stock = sum(s09q66, na.rm = TRUE),
                        equip_sales = sum(s09q68, na.rm = TRUE),
                        equip_purchases = sum(s09q70, na.rm = TRUE))
      df.hh <- merge(df.hh, df.equipment, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, equip_stock = if_else(is.na(equip_stock),0,equip_stock),
                      equip_sales = if_else(is.na(equip_sales),0,equip_sales),
                      equip_purchases = if_else(is.na(equip_purchases),0,equip_purchases))
      df.hh$equip_stock[df.hh$hhid %in% finddks(df.hrvs_9f, c("s09q66"))] <- NA
      df.hh$equip_sales[df.hh$hhid %in% finddks(df.hrvs_9f, c("s09q68"))] <- NA
      df.hh$equip_purchases[df.hh$hhid %in% finddks(df.hrvs_9f, c("s09q70"))] <- NA
      rm(df.equipment, df.hrvs_9f)
      
      df.hrvs_10 <- read.dta13("Section_10.dta")
      df.business_income <- mutate(df.hrvs_10, business_expenses = s10q05 + s10q06 + s10q07 + s10q08) %>%
              group_by(hhid) %>%
              summarize(business_revenues = sum(s10q04, na.rm = TRUE),
                        business_expenses = sum(business_expenses, na.rm = TRUE),
                        business_investment = sum(s10q09, na.rm = TRUE),
                        business_asset_sales = sum(s10q10, na.rm = TRUE))
      df.hh <- merge(df.hh, df.business_income, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, business_revenues = if_else(is.na(business_revenues),0,business_revenues),
                      business_expenses = if_else(is.na(business_expenses),0,business_expenses),
                      business_investment = if_else(is.na(business_investment),0,business_investment),
                      business_asset_sales = if_else(is.na(business_asset_sales),0,business_asset_sales))
      df.hh$business_revenues[df.hh$hhid %in% finddks(df.hrvs_10, c("s10q04"))] <- NA
      df.hh$business_expenses[df.hh$hhid %in% finddks(df.hrvs_10, c("s10q05", "s10q06", "s10q07", "s10q08"))] <- NA
      df.hh$business_investment[df.hh$hhid %in% finddks(df.hrvs_10, c("s10q09"))] <- NA
      df.hh$business_asset_sales[df.hh$hhid %in% finddks(df.hrvs_10, c("s10q10"))] <- NA
      rm(df.business_income, df.hrvs_10)
      
      df.hrvs_11 <- read.dta13("Section_11.dta")
      df.remittance_income <- mutate(df.hrvs_11, migration_costs = if_else(s11q03 < 13, s11q08b, 0),
                                     pastyear = if_else(s11q03 < 13, 1, 0),
                                     overseas = if_else(s11q02==2,1,0),
                                     not_roster = if_else(s11q01a==2,1,0)) %>%
              group_by(hhid) %>%
              summarize(n_migrants = as.numeric(n()),
                        tot_hhmembers = sum(not_roster, na.rm = TRUE),
                        migrants_overseas = sum(overseas, na.rm = TRUE),
                        migrants_past_year = sum(pastyear, na.rm = TRUE),
                        remittance_income = sum(s11q07c, na.rm = TRUE),
                        migration_costs = sum(migration_costs, na.rm = TRUE))
      df.hh <- merge(df.hh, df.remittance_income, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, n_migrants = if_else(is.na(n_migrants),0,n_migrants),
                      tot_hhmembers = if_else(is.na(tot_hhmembers),0,tot_hhmembers) + hhmembers,
                      migrants_overseas = if_else(is.na(migrants_overseas),0,migrants_overseas),
                      migrants_past_year = if_else(is.na(migrants_past_year),0,migrants_past_year),
                      remittance_income = if_else(is.na(remittance_income),0,remittance_income),
                      migration_costs = if_else(is.na(migration_costs),0,migration_costs),
                      prev_migrants = n_migrants - migrants_past_year)
      df.hh$remittance_income[df.hh$hhid %in% finddks(df.hrvs_11, c("s11q07c"))] <- NA
      df.hh$migration_costs[df.hh$hhid %in% finddks(df.hrvs_11, c("s11q08b"))] <- NA
      rm(df.remittance_income, df.hrvs_11)
      
      df.hrvs_12a <- read.dta13("Section_12a.dta")
      df.credit1 <- mutate(df.hrvs_12a, loans_taken = if_else(s12q05 < 13, s12q06, 0)) %>%
              group_by(hhid) %>%
              summarize(loans_taken_past_year = sum(loans_taken, na.rm = TRUE),
                        loans_total = sum(s12q06, na.rm = TRUE),
                        loan_payments = sum(s12q11, na.rm = TRUE),
                        avg_interest = weighted.mean(s12q09, s12q06, na.rm = TRUE),
                        avg_interest_past_year = weighted.mean(s12q09, loans_taken, na.rm = TRUE))                  ### payments are not necessarily last 12 months
      df.hh <- merge(df.hh, df.credit1, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, loans_taken_past_year = if_else(is.na(loans_taken_past_year),0,loans_taken_past_year),
                      loans_total = if_else(is.na(loans_total),0,loans_total),
                      loan_payments = if_else(is.na(loan_payments),0,loan_payments),
                      prev_loans_taken = loans_total - loans_taken_past_year)
      df.hh$loans_taken_past_year[df.hh$hhid %in% finddks(df.hrvs_12a, c("s12q06"))] <- NA
      df.hh$loans_total[df.hh$hhid %in% finddks(df.hrvs_12a, c("s12q06"))] <- NA
      df.hh$loan_payments[df.hh$hhid %in% finddks(df.hrvs_12a, c("s12q11"))] <- NA
      df.hh$avg_interest[df.hh$hhid %in% finddks(df.hrvs_12a, c("s12q09"))] <- NA
      df.hh$avg_interest_past_year[df.hh$hhid %in% finddks(df.hrvs_12a, c("s12q09"))] <- NA
      rm(df.credit1, df.hrvs_12a)
      
      df.hrvs_12b <- read.dta13("Section_12b.dta")
      df.credit2 <- mutate(df.hrvs_12b, loans_made = if_else(s12q16 < 13, s12q17, 0)) %>%
              group_by(hhid) %>%
              summarize(loans_made_past_year = sum(loans_made, na.rm = TRUE),
                        loans_made_total = sum(s12q17, na.rm = TRUE),
                        loan_payments_received = sum(s12q21, na.rm = TRUE),
                        avg_interest_charged = weighted.mean(s12q19, s12q17, na.rm = TRUE),
                        avg_interest_charged_past_year = weighted.mean(s12q19, loans_made, na.rm = TRUE))        
      df.hh <- merge(df.hh, df.credit2, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, loans_made_past_year = if_else(is.na(loans_made_past_year),0,loans_made_past_year),
                      loans_made_total = if_else(is.na(loans_made_total),0,loans_made_total),
                      loan_payments_received = if_else(is.na(loan_payments_received),0,loan_payments_received),
                      prev_loans_made = loans_made_total - loans_made_past_year)
      df.hh$loans_made_past_year[df.hh$hhid %in% finddks(df.hrvs_12b, c("s12q17"))] <- NA
      df.hh$loans_made_total[df.hh$hhid %in% finddks(df.hrvs_12b, c("s12q17"))] <- NA
      df.hh$loan_payments_received[df.hh$hhid %in% finddks(df.hrvs_12b, c("s12q21"))] <- NA
      df.hh$avg_interest_charged[df.hh$hhid %in% finddks(df.hrvs_12b, c("s12q19"))] <- NA
      df.hh$avg_interest_charged_past_year[df.hh$hhid %in% finddks(df.hrvs_12b, c("s12q19"))] <- NA
      rm(df.credit2, df.hrvs_12b)
      
      df.hrvs_12c <- read.dta13("Section_12c.dta")
      df.assets <- df.hrvs_12c %>% 
              mutate(fin_assest = if_else(assetid %in% c("Other Cash in Hand", "Bank (Current/Savings Account)", "Fixed Deposit",
                                                         "Stocks, Shares, Treasury Bills", "Employee Provident Fund/Citizen Investment Fund"), s12q23, 0),
                     savings_group = if_else(assetid=="Saving Cooperative/Savings group", s12q23, 0),
                     insurance_assets = if_else(assetid %in% c("Private health insurance", "Livestock and other agriculture related insurance", "Government health insurance",
                                                               "Life Insurance"), s12q23, 0),
                     cap_gains = if_else(s12q24 %in% c(998,999), as.numeric(NA), s12q24)) %>%
              group_by(hhid) %>%
              summarize(financial_assets = sum(fin_assest, na.rm = TRUE),
                        savings_group = sum(savings_group, na.rm = TRUE),
                        insurance_assets = sum(insurance_assets, na.rm = TRUE),
                        cap_gains = sum(cap_gains, na.rm = TRUE))
      df.hh <- merge(df.hh, df.assets, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, financial_assets = if_else(is.na(financial_assets),0,financial_assets),
                      savings_group = if_else(is.na(savings_group),0,savings_group),
                      insurance_assets = if_else(is.na(insurance_assets),0,insurance_assets),
                      cap_gains = if_else(is.na(cap_gains),0,cap_gains))
      df.hh$financial_assets[df.hh$hhid %in% finddks(df.hrvs_12c, c("s12q23"))] <- NA
      df.hh$savings_group[df.hh$hhid %in% finddks(df.hrvs_12c, c("s12q23"))] <- NA
      df.hh$insurance_assets[df.hh$hhid %in% finddks(df.hrvs_12c, c("s12q23"))] <- NA
      df.hh$cap_gains[df.hh$hhid %in% finddks(df.hrvs_12c, c("s12q24"))] <- NA
      rm(df.assets, df.hrvs_12c)
      
      df.hrvs_12d <- read.dta13("Section_12d.dta")
      df.pension <- df.hrvs_12d %>%
              group_by(hhid) %>%
              summarize(pension = sum(s12q27, na.rm = TRUE))
      df.hh <- merge(df.hh, df.pension, by = "hhid", all = TRUE)
      df.hh$pension[df.hh$hhid %in% finddks(df.hrvs_12d, c("s12q27"))] <- NA    
      rm(df.pension, df.hrvs_12d)
      
      df.hrvs_13a <- read.dta13("Section_13a.dta")
      df.gifts1 <- df.hrvs_13a %>% 
              group_by(hhid) %>%
              summarize(gifts_given_cash = sum(s13q07a, na.rm = TRUE),
                        gifts_given_inkind = sum(s13q07b, na.rm = TRUE))
      df.hh <- merge(df.hh, df.gifts1, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, gifts_given_cash = if_else(is.na(gifts_given_cash),0,gifts_given_cash),
                      gifts_given_inkind = if_else(is.na(gifts_given_inkind),0,gifts_given_inkind))
      df.hh$gifts_given_cash[df.hh$hhid %in% finddks(df.hrvs_13a, c("s13q07a"))] <- NA 
      df.hh$gifts_given_inkind[df.hh$hhid %in% finddks(df.hrvs_13a, c("s13q07b"))] <- NA   
      rm(df.gifts1, df.hrvs_13a)
      
      df.hrvs_13b <- read.dta13("Section_13b.dta")
      df.gifts2 <- df.hrvs_13b %>% 
              group_by(hhid) %>%
              summarize(gifts_received_cash = sum(s13q16a, na.rm = TRUE),
                        gifts_received_inkind = sum(s13q16b, na.rm = TRUE))
      df.hh <- merge(df.hh, df.gifts2, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, gifts_received_cash = if_else(is.na(gifts_received_cash),0,gifts_received_cash),
                      gifts_received_inkind = if_else(is.na(gifts_received_inkind),0,gifts_received_inkind))
      df.hh$gifts_received_cash[df.hh$hhid %in% finddks(df.hrvs_13b, c("s13q16a"))] <- NA 
      df.hh$gifts_received_inkind[df.hh$hhid %in% finddks(df.hrvs_13b, c("s13q16b"))] <- NA   
      rm(df.gifts2, df.hrvs_13b)
      
      df.hrvs_13c <- read.dta13("Section_13c.dta")
      df.NGOs <- df.hrvs_13c %>% 
              group_by(hhid) %>%
              summarize(NGO_cash = sum(s13q19a, na.rm = TRUE),
                        NGO_inkind = sum(s13q19c, na.rm = TRUE))
      df.hh <- merge(df.hh, df.NGOs, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, NGO_cash = if_else(is.na(NGO_cash),0,NGO_cash),
                      NGO_inkind = if_else(is.na(NGO_inkind),0,NGO_inkind))
      df.hh$NGO_cash[df.hh$hhid %in% finddks(df.hrvs_13c, c("s13q19a"))] <- NA 
      df.hh$NGO_inkind[df.hh$hhid %in% finddks(df.hrvs_13c, c("s13q19c"))] <- NA   
      rm(df.NGOs, df.hrvs_13c)
      
      df.hrvs_13d <- read.dta13("Section_13d.dta")
      df.donations <- df.hrvs_13d %>% 
              group_by(hhid) %>%
              summarize(donations_made = sum(s13q23, na.rm = TRUE))
      df.hh <- merge(df.hh, df.donations, by = "hhid", all = TRUE)
      df.hh$donations_made[df.hh$hhid %in% finddks(df.hrvs_13d, c("s13q23"))] <- NA   
      rm(df.donations, df.hrvs_13d)
      
      df.hrvs_14a <- read.dta13("Section_14a.dta")
      df.public <- df.hrvs_14a %>% mutate(quake_aid = if_else(pubcashid == "Earthquake Relief from government", s14q04b, 0),
                                          quake_aid_NGO = if_else(pubcashid == "Earthquake Relief from non-governmental sources", s14q04b, 0)) %>%
              group_by(hhid) %>%
              summarize(quake_aid = sum(quake_aid, na.rm = TRUE),
                        quake_aid_NGO = sum(quake_aid_NGO, na.rm = TRUE),
                        total_public_asst = sum(s14q04b, na.rm = TRUE),
                        non_quake_aid = total_public_asst - (quake_aid + quake_aid_NGO))
      df.hh <- merge(df.hh, df.public, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, quake_aid = if_else(is.na(quake_aid),0,quake_aid),
                      quake_aid_NGO = if_else(is.na(quake_aid_NGO),0,quake_aid_NGO),
                      total_public_asst = if_else(is.na(total_public_asst),0,total_public_asst),
                      non_quake_aid = if_else(is.na(non_quake_aid),0,non_quake_aid))
      df.hh$quake_aid[df.hh$hhid %in% finddks(df.hrvs_14a, c("s14q04b"))] <- NA   
      df.hh$quake_aid_NGO[df.hh$hhid %in% finddks(df.hrvs_14a, c("s14q04b"))] <- NA   
      df.hh$total_public_asst[df.hh$hhid %in% finddks(df.hrvs_14a, c("s14q04b"))] <- NA
      df.hh$non_quake_aid[df.hh$hhid %in% finddks(df.hrvs_14a, c("s14q04b"))] <- NA      
      rm(df.public, df.hrvs_14a)
      
      df.hrvs_14b <- read.dta13("Section_14b.dta")
      if (folders[i]=="Wave 3 - Household") {
            df.hrvs_14b <- rename(df.hrvs_14b,  s14q13b_q = s14q13b_i)
      }
      df.inkind <- df.hrvs_14b %>%
              group_by(hhid) %>%
              summarize(public_asst_inkind = sum(s14q13b_q, na.rm = TRUE)) 
      df.hh <- merge(df.hh, df.inkind, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, public_asst_inkind = if_else(is.na(public_asst_inkind),0,as.numeric(public_asst_inkind)))
      df.hh$public_asst_inkind[df.hh$hhid %in% finddks(df.hrvs_14b, c("s14q13b_q"))] <- NA   
      rm(df.inkind, df.hrvs_14b)
      
      df.hrvs_14c <- read.dta13("Section_14c.dta")
      df.workaid <- df.hrvs_14c %>% 
              mutate(wages = s14q19*s14q22a) %>%
              group_by(hhid) %>%
              summarize(work_aid_days = sum(s14q19, na.rm = TRUE),
                        work_aid_wages = sum(wages, na.rm = TRUE))
      df.hh <- merge(df.hh, df.workaid, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, work_aid_days = if_else(is.na(work_aid_days),0,as.numeric(work_aid_days)),
                      work_aid_wages = if_else(is.na(work_aid_wages),0,as.numeric(work_aid_wages)))
      df.hh$work_aid_days[df.hh$hhid %in% finddks(df.hrvs_14c, c("s14q19"))] <- NA   
      df.hh$work_aid_wages[df.hh$hhid %in% finddks(df.hrvs_14c, c("s14q22a", "s14q19"))] <- NA   
      rm(df.workaid, df.hrvs_14c)
      
      ### shock module
      df.hrvs_15a <- read.dta13("Section_15a.dta")
      df.quakebalance <- df.hrvs_15a %>% 
              mutate(quake_losses = if_else(shockid %in% c("Earthquake", "Landslide"),s15q03,0),
                     quake = if_else(shockid %in% c("Earthquake", "Landslide"),1,0),
                     riot_losses = if_else(shockid=="Riots/Blockage",s15q03,0),
                     riot = if_else(shockid=="Riots/Blockage",1,0),
                     price_shock_losses = if_else(shockid %in% c("Unexpected HIgher Prices", "Fuel Shortage"),s15q03,0),
                     price_shock = if_else(shockid %in% c("Unexpected HIgher Prices", "Fuel Shortage"),1,0),
                     livestock_farm_losses = if_else(shockid %in% c("Pests and Plant Diseases", "Lifestock Loss", "Post Harvest Loss"),s15q03,0),
                     livestock_farm_shock = if_else(shockid %in% c("Pests and Plant Diseases", "Lifestock Loss", "Post Harvest Loss"),1,0),
                     illness_injury_losses = if_else(shockid %in% c("Break up of family", "Disease or injury of family member", "Death of Family Member"),s15q03,0),
                     illness_injury_shock = if_else(shockid %in% c("Break up of family", "Disease or injury of family member", "Death of Family Member"),1,0),
                     job_default_losses = if_else(shockid %in% c("Loss of a regular job of a household member", "Failure or bankruptcy",
                                                                    "Loss of contract or default by creditor", "Withdrawal of government assistance"),s15q03,0),
                     job_default_shock = if_else(shockid %in% c("Loss of a regular job of a household member", "Failure or bankruptcy",
                                                                   "Loss of contract or default by creditor", "Withdrawal of government assistance"),1,0),,
                     violence_losses = if_else(shockid %in% c("Forced Displacement", "Theft"),s15q03,0),
                     violence_shock = if_else(shockid %in% c("Forced Displacement", "Theft"),1,0),
                     other_nat_disaster_losses = if_else(shockid %in% c("Hail/Lightening", "Flood", "Fire", "Drought"),s15q03,0),
                     other_nat_disaster = if_else(shockid %in% c("Hail/Lightening", "Flood", "Fire", "Drought"),1,0),
                     dissaving = if_else(s15q05a==1, 1, 0),
                     cut_food = if_else(s15q07a==1, 1, 0), 
                     cut_nonfood = if_else(s15q08a==1, 1, 0),  
                     school_interrupt = if_else(s15q09a==1, 1, 0)) %>%
              group_by(hhid) %>%
              summarize(quake_losses = sum(quake_losses, na.rm = TRUE),
                        quake = sum(quake, na.rm = TRUE),
                        riot_losses = sum(riot_losses, na.rm = TRUE),
                        riot = sum(riot, na.rm = TRUE),
                        price_shock_losses = sum(price_shock_losses, na.rm = TRUE),
                        price_shock = sum(price_shock, na.rm = TRUE),
                        livestock_farm_losses = sum(livestock_farm_losses, na.rm = TRUE),
                        livestock_farm_shock = sum(livestock_farm_shock, na.rm = TRUE),
                        illness_injury_losses = sum(illness_injury_losses, na.rm = TRUE),
                        illness_injury_shock = sum(illness_injury_shock, na.rm = TRUE),
                        job_default_losses = sum(job_default_losses, na.rm = TRUE),
                        job_default_shock = sum(job_default_shock, na.rm = TRUE),
                        violence_losses = sum(violence_losses, na.rm = TRUE),
                        violence_shock = sum(violence_shock, na.rm = TRUE),
                        other_nat_disaster = sum(other_nat_disaster, na.rm = TRUE),
                        other_nat_disaster_losses = sum(other_nat_disaster_losses, na.rm = TRUE),
                        dissaving = sum(dissaving, na.rm = TRUE),
                        cut_food = sum(cut_food, na.rm = TRUE),
                        cut_nonfood = sum(cut_nonfood, na.rm = TRUE),
                        school_interrupt = sum(school_interrupt, na.rm = TRUE))
      df.hh <- merge(df.hh, df.quakebalance, by = "hhid", all = TRUE)
      df.hh <- mutate(df.hh, quake_losses = if_else(is.na(quake_losses),0,quake_losses),
                      quake = if_else(is.na(quake) | quake==0,0,1),
                      riot_losses = if_else(is.na(riot_losses),0,riot_losses),
                      riot = if_else(is.na(riot) | riot==0,0,1),
                      price_shock_losses = if_else(is.na(price_shock_losses),0,price_shock_losses),
                      price_shock = if_else(is.na(price_shock) | price_shock==0,0,1),
                      livestock_farm_losses = if_else(is.na(livestock_farm_losses),0,livestock_farm_losses),
                      livestock_farm_shock = if_else(is.na(livestock_farm_shock) | livestock_farm_shock==0,0,1),
                      illness_injury_losses = if_else(is.na(illness_injury_losses),0,illness_injury_losses),
                      illness_injury_shock = if_else(is.na(illness_injury_shock) | illness_injury_shock==0,0,1),
                      job_default_losses = if_else(is.na(job_default_losses),0,job_default_losses),
                      job_default_shock = if_else(is.na(job_default_shock) | job_default_shock==0,0,1),
                      violence_shock = if_else(is.na(violence_shock) | violence_shock==0,0,1),
                      violence_losses = if_else(is.na(violence_losses),0,violence_losses),
                      other_nat_disaster = if_else(is.na(other_nat_disaster) | other_nat_disaster==0,0,1),
                      other_nat_disaster_losses = if_else(is.na(other_nat_disaster_losses),0,other_nat_disaster_losses),
                      dissaving = if_else(is.na(dissaving) | dissaving==0,0,1),
                      cut_food = if_else(is.na(cut_food) | cut_food==0,0,1),
                      cut_nonfood = if_else(is.na(cut_nonfood) | cut_nonfood==0,0,1),
                      school_interrupt = if_else(is.na(school_interrupt) | school_interrupt==0,0,1))
      rm(df.quakebalance, df.hrvs_15a)
      
      if (i==1) {
        df.hrvs_15b <- read.dta13("Section_15b.dta")
        df.emotion <- df.hrvs_15b %>%
          mutate(emotion_good = if_else(s15q15c %in% c("Good", "Very Good"), 1, 0)) %>%
          group_by(hhid) %>%
          summarize(emotion_good = sum(emotion_good, na.rm = TRUE))
        df.hh <- merge(df.hh, df.emotion, by = "hhid", all = TRUE)
        rm(df.emotion, df.hrvs_15b)
      } else {
        df.hh$emotion_good <- NA
      }
      
      ######## merge survey and ward dataframes ############# - land rents paid
      df.hh <- mutate(df.hh, income = annual_wages + landrents + rent_earned + wet_ag_sales + dry_ag_sales + food_home_production*52 + livestock_income + business_revenues + 
                          equip_rental_income + if_else(is.na(cap_gains),0,cap_gains) - business_expenses - livestock_costs - ag_costs - landrent_paid_cash - landrent_paid_inkind,
                      income_gross = annual_wages + landrents + rent_earned + wet_ag_sales + dry_ag_sales + food_home_production*52 + livestock_income + business_revenues + 
                        equip_rental_income + if_else(is.na(cap_gains),0,cap_gains),
                      inputs = business_expenses + livestock_costs + ag_costs + landrent_paid_cash + landrent_paid_inkind,
                      income_pc = income/currentlyliving,
                      ag_livestock_income = wet_ag_sales + dry_ag_sales + food_home_production*52 + livestock_income - livestock_costs - ag_costs - landrent_paid_cash - landrent_paid_inkind,
                      business_income = business_revenues - business_expenses,
                      consumption = food_home_production*52 + food_market*52 + food_inkind*52 + durables_consumption + energy + utilities_paid + rent_paid + transportation + 
                          clothing_cleaning_home + entertainment_other + other_expenses + craft_consumption,
                      consumption_pc = consumption/currentlyliving,
                      non_durables = food_home_production*52 + food_market*52 + food_inkind*52 + energy + utilities_paid + rent_paid + transportation + 
                          clothing_cleaning_home + entertainment_other + other_expenses + craft_consumption,
                      food_consumption = food_home_production*52 + food_market*52 + food_inkind*52,
                      food_consumption_pc = food_consumption/currentlyliving,
                      energy_utilities = energy + utilities_paid,
                      other_consumption = clothing_cleaning_home + entertainment_other + other_expenses + craft_consumption + rent_paid,
                      asset_stock = home_value + durables_stock + landvalue + equip_stock + livestock_value + financial_assets + savings_group + insurance_assets,
                      productive_assets = landvalue + equip_stock + livestock_value,
                      investments = land_purchased + business_investment + livestock_purchases + equip_purchases,
                      asset_sales = land_sales + business_asset_sales + livestock_sales + equip_sales,
                      net_asset_investment = investments - asset_sales,
                      net_loans = loans_taken_past_year - loans_made_past_year,
                      net_loan_payments = loan_payments - loan_payments_received,
                      prev_loans = prev_loans_taken - prev_loans_made,
                      pub_transfers = pension + non_quake_aid + public_asst_inkind + work_aid_wages,
                      inf_transfers = gifts_received_cash + gifts_received_inkind,
                      NGO_transfers = NGO_cash + NGO_inkind,
                      total_income = income_gross + pub_transfers + inf_transfers + NGO_transfers + remittance_income + quake_aid,
                      total_income_pc = total_income/currentlyliving,
                      market_shocks = riot_losses + price_shock_losses,
                      natural_shocks = livestock_farm_losses + other_nat_disaster_losses,
                      idiosyncratic_shocks = illness_injury_losses + violence_losses + job_default_losses,
                      market_shocks_bin = if_else(riot + price_shock>0,1,0),
                      natural_shocks_bin = if_else(livestock_farm_shock + other_nat_disaster>0,1,0),
                      idiosyncratic_shocks_bin = if_else(illness_injury_shock + violence_shock + job_default_shock>0,1,0),
                      total_shocks = market_shocks + natural_shocks,
                      ag_1season = if_else(wet_ag_sales + dry_ag_sales>0 & (wet_ag_sales==0 | dry_ag_sales==0) & livestock_income==0 & annual_wages==0 & business_revenues==0, 1, 0),
                      ag_2season = if_else(wet_ag_sales>0 & dry_ag_sales>0 & livestock_income==0 & annual_wages==0 & business_revenues==0, 1, 0),
                      livestock_only = if_else(wet_ag_sales==0 & dry_ag_sales==0 & livestock_income>0 & annual_wages==0 & business_revenues==0, 1, 0),
                      ag_livestock = if_else(wet_ag_sales + dry_ag_sales>0 & livestock_income>0 & annual_wages==0 & business_revenues==0, 1, 0),
                      ag_ls_wages = if_else(wet_ag_sales + dry_ag_sales + livestock_income>0 & annual_wages>0 & business_revenues==0, 1, 0),
                      ag_ls_business = if_else(wet_ag_sales + dry_ag_sales + livestock_income>0 & annual_wages==0 & business_revenues>0, 1, 0),
                      wage_only = if_else(wet_ag_sales + dry_ag_sales + livestock_income==0 & annual_wages>0 & business_revenues==0, 1, 0),
                      business_only = if_else(wet_ag_sales + dry_ag_sales + livestock_income + annual_wages==0 & business_revenues>0, 1, 0),
                      wage_bus = if_else(wet_ag_sales + dry_ag_sales + livestock_income==0 & annual_wages>0 & business_revenues>0, 1, 0),
                      ag_wage_bus = if_else(wet_ag_sales + dry_ag_sales + livestock_income>0 & annual_wages>0 & business_revenues>0, 1, 0), 
                      no_inc = if_else(wet_ag_sales + dry_ag_sales + livestock_income + annual_wages + business_revenues==0, 1, 0),
                      ag_ls_all = if_else(wet_ag_sales + dry_ag_sales + livestock_income>0 & annual_wages==0 & business_revenues==0, 1, 0),
                      diversified = if_else(ag_ls_wages==1 | ag_ls_business==1 | wage_bus==1 | ag_wage_bus==1, 1, 0))

      df.hh <- merge(df.hh, df.wards, by = c("district", "vdc", "ward")) 
      
      path <- paste(dirname(dirname(getwd())), "/processed/", folders[i], ".csv", sep = "")
      write_csv(df.hh, path = path)
      rm(df.hh, path)
      setwd("..")
}

setwd("..")
setwd("processed")
df.wave1 <- read_csv("Wave 1 - Household.csv")
df.wave2 <- read_csv("Wave 2 - Household.csv")
df.wave3 <- read_csv("Wave 3 - Household.csv")

df.wave1$wave <- 1
df.wave2$wave <- 2
df.wave3$wave <- 3

df.hh <- rbind(df.wave1, df.wave2, df.wave3)

### Create variables for each round of aid, cumulative aid, and district level aid, ward level avg losses, categorical distance buckets, and loss dummy

df.hh <- df.hh %>% mutate(quake_aid_bin = if_else(quake_aid > 0, 1, 0),
                          emergency_aid = if_else(quake_aid_bin==1 & wave==1, 1, 0),
                          reconstruction_aid_bin = if_else(quake_aid_bin==1 & wave!=1, 1, 0),
                          reconstruction_aid = if_else(quake_aid_bin==1 & wave!=1, quake_aid, 0),
                          gorkha_loss = if_else(quake_losses>0 & wave==1, 1, 0),
                          gorkha_loss_amt = if_else(wave==1, quake_losses, 0)) %>%
  group_by(hhid) %>%
  mutate(quake_aid_lag = if_else(is.na(lag(quake_aid_bin, order_by = wave)), 0, lag(quake_aid_bin, order_by = wave)),
         lag_remitt = lag(remittance_income, order_by = wave),
         lag_loan = lag(loans_taken_past_year, order_by = wave))

df.aid <- df.hh %>% arrange(wave) %>%
  group_by(hhid) %>% 
  summarize(aid_total = sum(quake_aid),
          var_cons = var(log(consumption)),
          avg_cons = mean(log(consumption)),
          avg_land = mean(landvalue),
          var_cons_pc = var(log(consumption_pc)),
          avg_cons_pc = mean(log(consumption_pc)),
          gorkha_loss_ever = sum(gorkha_loss))

df.ward_losses <- df.hh %>% 
  group_by(district, vdc, ward) %>%
  summarize(ward_avg_losses = mean(gorkha_loss_amt, na.rm=TRUE),
            ward_frac_losses = mean(gorkha_loss, na.rm = TRUE))

df.hh <- merge(df.hh, df.ward_losses, by = c("district", "vdc", "ward"))
df.hh <- merge(df.hh, df.aid, by = "hhid") %>%
  mutate(received_aid = if_else(aid_total>0, 1, 0),
         caste_yr = paste(ethnicity, wave, sep = "_"),
         vdc_yr = paste(vdc, wave, sep = "_"),
         cons_shock = consumption - avg_cons) %>%
  arrange(wave) %>% group_by(hhid) %>%
  mutate(aid_cumulative = cumsum(quake_aid),
         recon_aid_cum = cumsum(reconstruction_aid),
         aid_cumulative_bin = as.numeric(aid_cumulative>0),
         recon_aid_cum_bin = as.numeric(recon_aid_cum>0),
         aid_cumulative0000s = aid_cumulative/10000) %>% ungroup()

df.hh$consumption_qtle <- cut(df.hh$avg_cons, quantile(df.hh$avg_cons, seq(0, 1, by = .2), na.rm = TRUE), labels = FALSE, include.lowest = TRUE)
df.hh$land_qtle <- cut(df.hh$avg_land, quantile(df.hh$avg_land, c(0,.5,1), na.rm = TRUE), labels = FALSE, include.lowest = TRUE)

df.hh$caste_recode <- factor(df.hh$caste_recode, levels = c("Other", "Brahman/Chhetri", "Dalit", "Newar", NA))

write_csv(df.hh, "full_panel.csv")

