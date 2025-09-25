rm(list = ls())

#Setup ----
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)
library(future) #parallelise lapply() : future_lapply()
library(future.apply) #parallelise lapply(): future_lapply()

#Set parallelise plan
# parallel::detectCores(logical = F) #256
plan(multicore, workers = 16)
Sys.setenv(OMP_NUM_THREADS = 4, OPENBLAS_NUM_THREADS = 4)
#maximum number of threads used by OpenMP and OpenBLAS, usually under the hood for R

path = "/maps/epr26/sdm_captain_out/"
log_file = paste0(path, "log_disturbance_score.txt")
aoi_proj = vect(paste0(path, "atlantic_forest_global_200.geojson")) %>%
  project("EPSG:3857")
bioclim = rast(paste0(path, "rasters/bioclim_reduced.tif"))
samp_size_df = read.csv(paste0(path, "species_sample_size.csv"), header = T)
n_sp = nrow(samp_size_df)

#Create or read raster mask for the Atlantic Forest ecoregion
if(!file.exists(paste0(path, "rasters/af_mask.tif"))) {
  af_mask = bioclim[[1]] %>%
    mask(aoi_proj) %>%
    classify(rcl = matrix(c(-Inf, Inf, 1, NA, NA, 0), ncol = 3, byrow = T), right = NA, others = 0) #turn non Na values to 1
  writeRaster(af_mask, paste0(path, "rasters/af_mask.tif"), overwrite = T)
} else {
  af_mask = rast(paste0(path, "rasters/af_mask.tif"))
}

#Create or read 10-km grids over which to calculate percentage cover of each land cover class
if(!file.exists(paste0(path, "grid_vect.geojson"))) {
  grid = rast(bioclim[[1]])
  names(grid) = "grid_id"
  grid_vect = as.polygons(grid)
  values(grid_vect) = seq_len(ncell(grid)) #assign unique IDs
  writeVector(grid_vect, paste0(path, "grid_vect.geojson"), overwrite = T)
} else {
  grid_vect = vect(paste0(path, "grid_vect.geojson"))
}


#Set up function to calculate percentage covers in each grid ----
PercCover = function(i) {
  a_grid = Sys.time()
  grid_i = vect(paste0(path, "grid_vect.geojson"))[i]

  #Load cropped and reprojected land cover raster
  lc_i = rast(paste0(path, "rasters/brazil_coverage_2024_af_proj.tif")) %>%
    crop(grid_i, snap = "out")
  #Created by:
  # gdalwarp -te -58 -34 -34 -3 -r near /maps/epr26/sdm_captain_out/rasters/brazil_coverage_2024.tif \
  #   /maps/epr26/sdm_captain_out/rasters/brazil_coverage_2024_af.tif
  # gdalwarp -t_srs EPSG:3857 /maps/epr26/sdm_captain_out/rasters/brazil_coverage_2024_af.tif \
  #   /maps/epr26/sdm_captain_out/rasters/brazil_coverage_2024_af_proj.tif

  #Load secondary vegetation age raster
  sf_i = rast(paste0(path, "rasters/secondary_vegetation_age_2023_proj.tif")) %>%
    crop(grid_i, snap = "out")
  #Created by:
  # gdalwarp -te -58 -34 -34 -3 -r bilinear /maps/epr26/sdm_captain_out/rasters/mapbiomas_collection90_secondary_vegetation_age_v1-secondary_vegetation_age_2023.tif \
  #   /maps/epr26/sdm_captain_out/rasters/secondary_vegetation_age_2023_cropped.tif
  # gdalwarp -t_srs EPSG:3857 /maps/epr26/sdm_captain_out/rasters/secondary_vegetation_age_2023_cropped.tif \
  #   /maps/epr26/sdm_captain_out/rasters/secondary_vegetation_age_2023_proj.tif

  sf_mask = (sf_i != 0) #create a mask for when secondary forest age > 0
  lc_i[sf_mask] = 76 #assign new code for secondary forest

  #get pixel count and proportions
  tb = table(values(lc_i))
  prop = as.numeric(tb / sum(tb) * 100)

  #fill in zeros for missing classes and standardise output
  missing_classes = setdiff(0:76, as.numeric(names(tb))) #0-75 are land cover classes
  if(length(missing_classes) > 0) {
    prop = c(prop, rep(0, length(missing_classes)))
    names(prop) = paste0("c", str_pad(c(names(tb), as.character(missing_classes)), 2, side = "left", pad = "0"))
  }
  prop = prop[order(names(prop))] %>%
    t() %>%
    as.data.frame() %>%
    mutate(grid = i)
  write.csv(prop, paste0(path, "props/land_cover_percentage_grid_", i, ".csv"), row.names = F)
  b_grid = Sys.time()
  cat("Completed grid", i, ":", difftime(b_grid, a_grid, units = "secs"), "seconds\n")
  return(prop)
}

#Run calculations ----
PercCover_safe = function(i) {
  tryCatch({
    PercCover(i)
    cat("Success for grid ", i, "\n", file = log_file, append = T)
  }, error = function(e) {
    cat("Error in grid ", i, ":", e$message, "\n", file = log_file, append = T)
  })
}

future_lapply(seq_along(grid_vect), PercCover_safe,
              future.packages = c("terra"),
              future.seed = T)


#Collate results ----
prop_list = vector("list", nrow(grid_vect))
for(i in seq_along(grid_vect)) {
  prop_list[[i]] = read.csv(paste0(path, "props/land_cover_percentage_grid_", i, ".csv"), header = T)
  if (!identical(colnames(prop_list[[i]]), c(paste0("c", str_pad(0:76, 2, side = "left", pad = "0")), "grid"))) {
    cat("Column names of grid", i, "is wrong, it is:", colnames(prop_list[[i]]), "\n")
  }
}
prop_df = bind_rows(prop_list)
write.csv(prop_df, paste0(path, "land_cover_percentage.csv"), row.names = F)

#Reclassify and calculate disturbance score ----
prop_df = read.csv(paste0(path, "land_cover_percentage.csv"), header = T)

prop_rcl = prop_df %>%
  mutate(forest = c01 + c03 + c05 + c06 + c49,
         grassland = c04 + c12 + c29 + c50,
         water = c11 + c32 + c26 + c33,
         beach = c23,
         forest_sec = c76,
         plantation = c09,
         cultivation = c14 + c15 + c18 + c19 + c39 + c20 + c40 + c62 + c41 + c36 + c46 + c47 + c35 + c48 + c21 + c31,
         urban = c24 + c30 + c75 + c25) %>%
  mutate(other = 100 - forest - grassland - water - beach - forest_sec - plantation - cultivation - urban) %>%
  dplyr::select(grid, forest, grassland, water, beach, forest_sec, plantation, cultivation, urban, other) %>%
  mutate(score = (forest_sec * 0.3 + plantation * 0.5 + cultivation * 0.7 + urban * 1 + other * 0.5) / 100)
write.csv(prop_rcl, paste0(path, "land_cover_percentage_reclassified.csv"), row.names = F)

grid_vect$score = prop_rcl$score
disturbance = rasterize(grid_vect, bioclim[[1]], field = "score")
varnames(disturbance) = "disturbance"
writeRaster(disturbance, paste0(path, "rasters/disturbance_score.tif"), overwrite = T)
