rm(list = ls())

#Setup ----
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)
library(geodata) #geodata::travel_time, geodata::crop_spam

path = "/maps/epr26/sdm_captain_out/"
bioclim = rast(paste0(path, "rasters/bioclim_reduced.tif"))
aoi_proj = vect(paste0(path, "atlantic_forest_global_200.geojson")) %>%
  project("EPSG:3857")
af_mask = rast(paste0(path, "rasters/af_mask.tif"))

#Load data of travel time to nearest city
accessibility = rast(paste0(path, "rasters/travel_time_to_cities_u9_proj.tif"))
#Created by:
# accessibility = geodata::travel_time(to = "city", size = 9, up = T,
#                                      path = paste0(path, "travel_time/"))
#gdalwarp -te -58 -34 -34 -3 -r bilinear /maps/epr26/sdm_captain_out/travel_time/travel/travel_time_to_cities_u9.tif \
#   /maps/epr26/sdm_captain_out/rasters/travel_time_to_cities_u9_cropped.tif
#gdalwarp -t_srs EPSG:3857 /maps/epr26/sdm_captain_out/rasters/travel_time_to_cities_u9_cropped.tif \
#   /maps/epr26/sdm_captain_out/rasters/travel_time_to_cities_u9_proj.tif


#Load data of production per ha for all crops using all technologies
val_prod = rast(paste0(path, "rasters/spam_2010_val_prod_per_area_proj.tif"))
#Created by:
# geodata::crop_spam(crop = "maize", var = "val_prod",
#                    path = "/maps/epr26/sdm_captain_out/crop_spam")
# gdalwarp -te -58 -34 -34 -3 -r bilinear /maps/epr26/sdm_captain_out/crop_spam/spam/spam2010V2r0_global_V_agg_VP_CR_AR_A.tif \
#   /maps/epr26/sdm_captain_out/rasters/spam_2010_val_prod_per_area_cropped.tif
# gdalwarp -t_srs EPSG:3857 /maps/epr26/sdm_captain_out/rasters/spam_2010_val_prod_per_area_cropped.tif \
#   /maps/epr26/sdm_captain_out/rasters/spam_2010_val_prod_per_area_proj.tif

#Retrieve or read elevation data
elevation = elevation_global(res = 0.5, path = paste0(path, "elevation"), mask = T) %>%
    project("EPSG:3857")

#resample to 10-km resolution and scale values to 0-1 range
accessibility_resamp = resample(accessibility, bioclim[[1]], method = "bilinear",
                                filename = paste0(path, "rasters/travel_time_to_cities_u9_proj_resamp.tif"), overwrite = T)
accessibility_range = range(extract(accessibility_resamp, aoi_proj)[, 2], na.rm = T)
accessibility_scaled = (accessibility_resamp - accessibility_range[1]) / (accessibility_range[2] - accessibility_range[1])
writeRaster(accessibility_scaled, paste0(path, "rasters/accessibility_scaled.tif"), overwrite = T)

val_prod_resamp = resample(val_prod, bioclim[[1]], method = "bilinear",
                           filename = paste0(path, "rasters/spam_2010_val_prod_per_area_proj_resamp.tif"), overwrite = T)
val_prod_range = range(extract(val_prod_resamp, aoi_proj)[, 2], na.rm = T)
val_prod_scaled = (val_prod_resamp - val_prod_range[1]) / (val_prod_range[2] - val_prod_range[1])
writeRaster(val_prod_scaled, paste0(path, "rasters/val_prod_scaled.tif"), overwrite = T)

elevation_resamp = resample(elevation, bioclim[[1]], method = "bilinear",
                            filename = paste0(path, "rasters/elevation_resamp.tif"), overwrite = T)
elevation_range = range(extract(elevation_resamp, aoi_proj)[, 2], na.rm = T)
elevation_scaled = (elevation_resamp - elevation_range[1]) / (elevation_range[2] - elevation_range[1])
writeRaster(elevation_scaled, paste0(path, "rasters/elevation_scaled.tif"), overwrite = T)

cost = (val_prod_scaled + accessibility_scaled + elevation_scaled) / 3
cost_na = is.na(cost) & !(af_mask == 0) #find cells that are NA *inside* the vector
cost_out = (af_mask == 0)
cost_rcl = cost
cost_rcl[cost_out] = NA #replace those with 1 in the original raster
cost_rcl[cost_na] = 1 #replace those with 1 in the original raster
writeRaster(cost_rcl, paste0(path, "rasters/cost_reclassified.tif"), overwrite = T)
