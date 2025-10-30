rm(list = ls())

#Setup ----
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)
library(geodata) #geodata::travel_time, geodata::crop_spam

dir_path = "/maps/epr26/captain_brazil/"
#proj_path = "/maps/epr26/captain_brazil/af_10km/"
proj_path = "/maps/epr26/captain_brazil/ideal_250m/"
bioclim = rast(paste0(proj_path, "rasters/bioclim.tif"))
aoi = vect(paste0(proj_path, "aoi.geojson"))
aoi_proj = vect(paste0(proj_path, "aoi_proj.geojson"))
proj_res = 250 #in meters

#Create or read AOI mask
if(!file.exists(paste0(proj_path, "rasters/aoi_mask.tif"))) {
  aoi_mask = bioclim[[1]] %>%
    mask(aoi_proj) %>%
    classify(rcl = matrix(c(-Inf, Inf, 1, NA, NA, 0), ncol = 3, byrow = T), right = NA, others = 0) #turn non Na values to 1
  writeRaster(aoi_mask, paste0(proj_path, "rasters/aoi_mask.tif"), overwrite = T)
} else {
  aoi_mask = rast(paste0(proj_path, "rasters/aoi_mask.tif"))
}

#Load data of travel time to nearest city
accessibility_orig = geodata::travel_time(to = "city", size = 9, up = T, path = dir_path)
accessibility = accessibility_orig %>%
  crop(ext(aoi)) %>%
  project(crs(aoi_proj), res = proj_res) #reproject to defined resolution
writeRaster(accessibility, paste0(proj_path, "rasters/accessibility.tif"), overwrite = T)

#Load data of production per ha for all crops using all technologies
#val_prod = rast(paste0(path, "rasters/spam_2010_val_prod_per_area_proj.tif"))
if(!file.exists(paste0(dir_path, "crop_spam/spam/spam2010V2r0_global_V_agg_VP_CR_AR_A.tif"))) {
  geodata::crop_spam(crop = "maize", var = "val_prod", path = dir_path)
}
val_prod = rast(paste0(dir_path, "crop_spam/spam/spam2010V2r0_global_V_agg_VP_CR_AR_A.tif")) %>%
  crop(ext(buffer(aoi, 20000))) %>% #20-km buffer to avoid edge effects
  focal(w = 3, fun = function(x) {if (is.na(x[5])) mean(x, na.rm = T) else x[5]}) %>% #interpolate NAs using 3x3 moving window
  project(crs(aoi_proj), res = proj_res) %>% #reproject to defined resolution
  crop(ext(aoi_proj)) #crop again
writeRaster(val_prod, paste0(proj_path, "rasters/val_prod.tif"), overwrite = T)

#Retrieve or read elevation data
elevation_orig = geodata::elevation_global(res = 0.5, path = paste0(dir_path, "elevation"), mask = T)
elevation = elevation_orig %>%
  crop(ext(aoi)) %>%
  project(crs(aoi_proj), res = proj_res) #reproject to defined resolution
writeRaster(elevation, paste0(proj_path, "rasters/elevation.tif"), overwrite = T)

#scale values to 0-1 range
accessibility_range = range(extract(accessibility, aoi_proj)[, 2], na.rm = T)
accessibility_scaled = (accessibility - accessibility_range[1]) / (accessibility_range[2] - accessibility_range[1])
writeRaster(accessibility_scaled, paste0(path, "rasters/accessibility_scaled.tif"), overwrite = T)

val_prod_range = range(extract(val_prod, aoi_proj)[, 2], na.rm = T)
val_prod_scaled = (val_prod - val_prod_range[1]) / (val_prod_range[2] - val_prod_range[1])
writeRaster(val_prod_scaled, paste0(path, "rasters/val_prod_scaled.tif"), overwrite = T)

elevation_range = range(extract(elevation, aoi_proj)[, 2], na.rm = T)
elevation_scaled = (elevation - elevation_range[1]) / (elevation_range[2] - elevation_range[1])
writeRaster(elevation_scaled, paste0(path, "rasters/elevation_scaled.tif"), overwrite = T)

cost = (val_prod_scaled + accessibility_scaled + elevation_scaled) / 3
cost_na = is.na(cost) & !(aoi_mask == 0) #find cells that are NA *inside* the vector
cost_out = (aoi_mask == 0)
cost_rcl = cost
cost_rcl[cost_out] = NA #replace those with 1 in the original raster
cost_rcl[cost_na] = 1 #replace those with 1 in the original raster
writeRaster(cost_rcl, paste0(path, "rasters/cost_reclassified.tif"), overwrite = T)
