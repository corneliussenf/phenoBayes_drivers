
library(tidyverse)
library(RCurl)
library(raster)
library(ggmap)
library(lubridate)

#### Get extent of study region ####

loc <- shapefile("paper/data/site/BFNP_simplified.shp")
loc_latlong <- spTransform(loc, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

extent <- extent(loc_latlong)

extent <- data.frame(xmin = extent[1] - 0.1,
                     xmax = extent[2] + 0.1,
                     ymin = extent[3] - 0.1,
                     ymax =  extent[4] + 0.1)

ggplot() +
  geom_polygon(data = fortify(loc_latlong), aes(x = long, y = lat, group = group)) +
  geom_rect(data = extent, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.2)

#### Select DWD stations ####

getDateFromDWD <- function(x){
  d <- as.Date(paste(substring(x, 1, 4),
                     substring(x, 5, 6),
                     substring(x, 7, 8), sep = "-"),
               format = "%Y-%m-%d")
  return(d)
}

stations <- read.csv("paper/data/dwd_climate/KL_Tageswerte_Beschreibung_Stationen.csv", sep = ";")
stations <- stations[-nrow(stations),]
stations$date_start <- getDateFromDWD(stations$von_datu)
stations$date_end <- getDateFromDWD(stations$m.bis_datum)

loc_check <- (findInterval(stations$geoBreite, extent[,c("ymin", "ymax")]) == 1 &
                findInterval(stations$geoLaenge, extent[,c("xmin", "xmax")]) == 1)
time_check <- year(stations$date_start) <= 1985 & year(stations$date_end) >= 2016

stations_included <- stations[loc_check & time_check, ]

#map <- get_map(location = c(lon = mean(as.double(extent[1,1:2])), lat = mean(as.double(extent[1,3:4]))), 
#               zoom = 7, maptype = "satellite", source = "google")

#ggmap(map) +
ggplot() + 
  geom_rect(data = extent, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_polygon(data = fortify(loc_latlong), aes(x = long, y = lat, group = group), fill = NA, col = "black") +
  geom_point(data = stations_included, aes(x = geoLaenge, y = geoBreite), shape = 2, col = "red")

#### Download and unzip DWD data ####

url <- "ftp://ftp-cdc.dwd.de/pub/CDC/observations_germany/climate/daily/kl/historical/"

filenames <- getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
files <- grep("*.zip$", strsplit(filenames, "\n", fixed = FALSE)[[1]], value = TRUE)
files_include <- files[substring(files, 15, 19) %in% str_pad(stations_included$Stations_id, 5, pad = "0")]

for(i in files_include){
  print(paste0("Downloading ", i))
  destfile <- paste0("paper/data/dwd_climate/", i)
  download.file(paste0(url, i), destfile = destfile)
  print(paste0("Un-zipping ", i))
  unzip(destfile, exdir = gsub(".zip", "", destfile))
  file.remove(destfile)
}

#### Load climate data ####

files <- list.files("data/dwd_climate/", pattern = "produkt_klima_tag*", recursive = TRUE, full.names = TRUE)

dat <- lapply(files, function(x){
  d <- read.table(x, sep = ";", header = TRUE) 
  d[d == -999] <- NA
  d$date <- as.Date(with(d, paste0(substring(MESS_DATUM, 1, 4), "-", 
                                   substring(MESS_DATUM, 5, 6), "-",
                                   substring(MESS_DATUM, 7, 8))), format = "%Y-%m-%d")
  return(d)
})

dat <- do.call("rbind", dat)

dat$location_id <- as.integer(as.character(dat$STATIONS_ID))
dat$year <- year(as.Date(dat$date))
dat$month <- month(as.Date(dat$date))
dat$doy <- yday(as.Date(dat$date))

# dat %>%
#   group_by(STATIONS_ID) %>%
#   summarize(start = min(date), 
#             end = max(date))

### Statistics

# Pre-Season mean temperature

daily_temp_prec <- dat %>%
  filter(year %in% 1985:2015) %>%
  group_by(year, month, doy) %>%
  summarize(tmean = mean(TMK, na.rm = TRUE),
            tmax = mean(TXK, na.rm = TRUE),
            tmin = mean(TNK, na.rm = TRUE),
            prec = mean(RSK, na.rm = TRUE))

climate_vars <- daily_temp_prec %>%
  filter(month %in% c(4:5)) %>%
  group_by(year) %>%
  summarize(preseason_temp = mean(tmean)) %>%
  mutate(preseason_temp_anomaly = (preseason_temp - mean(preseason_temp)) / sd(preseason_temp))

spring_chilling <- lapply(unique(daily_temp_prec$year), function(y){
  t1 <- filter(daily_temp_prec, year == y & month >= 11)$tmean 
  t2 <- filter(daily_temp_prec, year == y + 1 & month <= 5)$tmean
  cd <- sum(t1 > 0 & t1 <= 5) + sum(t2 > 0 & t2 <= 5)
})
climate_vars$chilling <- unlist(spring_chilling)
climate_vars$chilling_anomaly <- (climate_vars$chilling - mean(climate_vars$chilling)) / sd(climate_vars$chilling)

spring_prec <- daily_temp_prec %>%
  filter(month %in% c(4:5)) %>%
  group_by(year) %>%
  summarize(preseason_prec = mean(prec))

climate_vars$preseason_prec <- spring_prec$preseason_prec
climate_vars$preseason_prec_anomaly <- (climate_vars$preseason_prec - mean(climate_vars$preseason_prec)) / sd(climate_vars$preseason_prec)

save(climate_vars, file = "paper/data/dwd_climate/climate_dwd_spring.RData")

# ggplot(cliamte_vars, aes(x = year, y = preseason_prec_anomaly)) +
#   geom_bar(stat = "identity")


