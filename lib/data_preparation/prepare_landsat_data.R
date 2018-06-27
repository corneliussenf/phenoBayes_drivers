
setwd("Q:/PRIME/phenology_bwald/model_development/processmodel/")

library(dplyr)
library(lubridate)
library(foreign)
library(rgdal)
library(raster)
library(maptools)

# Sample size

n <- 250

# Sampling -------------------------------------------------------------

habitat <- shapefile("../basemodel/data/bwald/Habitatkarte_utm33.shp")
np <- shapefile("../basemodel/data/bwald/BFNP.shp")
np <- spTransform(np, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# bl_old <- subset(habitat, habitat@data$KLS_DEU == "Laubbestand - alt" & habitat@data$Shape_area > (5*30)^2)
con_old <- subset(habitat, habitat@data$KLS_DEU == "Nadelbestand - alt" & habitat@data$Shape_area > (5*30)^2)
# bl_old_latlong <- spTransform(bl_old, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
con_old_latlong <- spTransform(con_old, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# sample_bl_old <- spsample(bl_old, n, type = "random")
# sample_bl_old_latlong <- spTransform(sample_bl_old, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

sample_con_old <- spsample(con_old, n, type = "random")
sample_con_old_latlong <- spTransform(sample_con_old, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Extract Landsat data ----------------------------------------------------

extract_bands <- function(band, loc, image_path){
  
  lt <- list.files(path = image_path, pattern = glob2rx(paste0("LT*sr_band", band, "_*.vrt")), recursive = TRUE, full.names = TRUE)
  le <- list.files(path = image_path, pattern = glob2rx(paste0("LE*sr_band", band, "_*.vrt")), recursive = TRUE, full.names = TRUE)
  if(band < 7){
    lc <- list.files(path = image_path, pattern = glob2rx(paste0("LC*sr_band", band + 1, "_*.vrt")), recursive = TRUE, full.names = TRUE)
  }else{
    lc <- list.files(path = image_path, pattern = glob2rx(paste0("LC*sr_band", band, "_*.vrt")), recursive = TRUE, full.names = TRUE)
  }
  l <- c(lt, le, lc)
  
  s <- lapply(l, function(x)tryCatch(raster(x), error = function(e)return(NA)))
  s_correct <- unlist(lapply(s, function(x)class(x)=="RasterLayer"))
  s <- s[s_correct]
  
  out <- extract(stack(s), loc, df = TRUE)
  out$x <- coordinates(loc)[, 1]
  out$y <- coordinates(loc)[, 2]
  
  m <- reshape2::melt(out, id.vars = c("ID", "x", "y"))
  m$sceneid <- substring(m$variable, 1, 21)
  m <- m[,-which(names(m) == "variable")]
  
  return(m)
}

samples <- lapply(c(1, 2, 3, 4, 5, 7), extract_bands, loc = sample_con_old, 
                  image_path = "Q:/Geomultisens/workspace/bwald/landsat/tiles/33/33UUQ/")

dat <- plyr::join_all(samples, by =  c("ID", "x", "y", "sceneid"))

names(dat) <- c("pixel", "x", "y", "band1", "sceneid", "band2", "band3", "band4", "band5", "band7")
dat <- dat[,c("pixel", "sceneid", "x", "y", paste0("band", 1:5), "band7")]

dat$year <- as.integer(substring(dat$sceneid, 10, 13))
dat <- subset(dat, dat$year < 2016 & dat$year > 1984)
dat$doy <- as.integer(substring(dat$sceneid, 14, 16))
dat$date <- as.Date(paste(dat$year, dat$doy, sep = "-"), format = "%Y-%j")
dat$year_index <- dat$year - (min(dat$year) - 1)
dat$station <- substring(dat$sceneid, 17, 19)
dat$archive <- ifelse(dat$station == "ESA", "ESA", "USGS")
dat$sensor <- factor(substring(dat$sceneid, 1, 3))

# Remove clouds and filter for remaining artefacts

extract_clouds <- function(loc, image_path){
  
  l <- list.files(path = image_path, pattern = glob2rx("*cfmask_3*.vrt"), recursive = TRUE, full.names = TRUE)
  
  s <- lapply(l, function(x)tryCatch(raster(x), error = function(e)return(NA)))
  s_correct <- unlist(lapply(s, function(x)class(x)=="RasterLayer"))
  s <- s[s_correct]
  
  out <- extract(stack(s), loc, df = TRUE)
  out$x <- coordinates(loc)[, 1]
  out$y <- coordinates(loc)[, 2]
  
  m <- reshape2::melt(out, id.vars = c("ID", "x", "y"))
  m$sceneid <- substring(m$variable, 1, 21)
  m <- m[,-which(names(m) == "variable")]
  m <- m[,-which(names(m) == "ID")]
  
  return(m)
}

clouds <- extract_clouds(sample_con_old, "Q:/Geomultisens/workspace/bwald/landsat/tiles/33/33UUQ/")

clouds$cloud <- ifelse(clouds$value == 0, 0, 1)
clouds$cloud_value <- clouds$value
clouds <- clouds[,-which(names(clouds) == "value")]

dat <- full_join(dat, clouds, by = c("x", "y", "sceneid"))

dat <- subset(dat, dat$cloud == 0)
dat <- subset(dat, !(dat$band1 < 0 | dat$band3 < 0 | dat$band4 < 0))
dat <- na.omit(dat)

# Calculate spectral indices

dat <- mutate(dat, ndvi = (band4 - band3) / (band4 + band3))
dat <- mutate(dat, evi = (band4/10000 - band3/10000) / (band4/10000 + 6 * band3/10000 - 7.5 * band1/10000 + 1))

# Save data

save(dat, file = "data/landsat_coniferous.RData")

# Save locations

loc <- dat %>% group_by(pixel) %>% summarize(x = mean(x), y = mean(y))

coordinates(loc) <- ~ x+y

loc@data <- as.data.frame(loc@data)

writeOGR(obj = loc, dsn = "data", 
         layer = "landsat_locations_coniferous", 
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)

