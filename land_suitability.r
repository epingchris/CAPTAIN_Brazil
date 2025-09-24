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
log_file = paste0(path, "log_land_cover.txt")
aoi = vect(paste0(path, "atlantic_forest_global_200.geojson")) %>%
  project("EPSG:4326")
aoi_proj = project(aoi, "EPSG:3857")
bioclim = rast(paste0(path, "rasters/bioclim_reduced.tif"))

#Create raster mask for the Atlantic Forest ecoregion
if(!file.exists(paste0(path, "rasters/af_mask.tif"))) {
  af_mask = bioclim[[1]] %>%
    mask(aoi_proj) %>%
    classify(rcl = matrix(c(-Inf, Inf, 1, NA, NA, 0), ncol = 3, byrow = T), right = NA, others = 0) #turn non Na values to 1
  writeRaster(af_mask, paste0(path, "rasters/af_mask.tif"), overwrite = T)
} else {
  af_mask = rast(paste0(path, "rasters/af_mask.tif"))
}

#Create or load cropped and reprojected land cover raster
if(!file.exists(paste0(path, "rasters/brazil_coverage_2024_af_proj.tif"))) {
  lc_af = rast(paste0(path, "rasters/brazil_coverage_2024_af.tif"))
  lc_proj = project(lc_af, "EPSG:3857", method = "near",
                    filename = paste0(path, "rasters/brazil_coverage_2024_af_proj.tif"), overwrite = T)
} else {
  lc_proj = rast(paste0(path, "rasters/brazil_coverage_2024_af_proj.tif"))
}

#Create 10-km grids over which to calculate percentage cover of each land cover class
if(!file.exists(paste0(path, "grid_vect.geojson"))) {
  grid = rast(bioclim[[1]])
  names(grid) = "grid_id"
  grid_vect = as.polygons(grid)
  values(grid_vect) = seq_len(ncell(grid)) #assign unique IDs
  writeVector(grid_vect, paste0(path, "grid_vect.geojson"), overwrite = T)
} else {
  grid_vect = vect(paste0(path, "grid_vect.geojson"))
}

#Set up function to calculate percentage covers in each grid
PercCover = function(i) {
  a_grid = Sys.time()
  lc_proj = rast(paste0(path, "rasters/brazil_coverage_2024_af_proj.tif"))
  grid_i = vect(paste0(path, "grid_vect.geojson"))[i]
  lc_extract = extract(lc_proj, grid_i)
  tb = table(lc_extract$brazil_coverage_2024)
  prop = as.numeric(tb / sum(tb) * 100)
  names(prop) = paste0("class_", names(tb))

  #fill in zeros for missing classes
  missing_classes = setdiff(0:75, as.numeric(names(tb))) #0-75 are land cover classes
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
}

#Run models ----
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

prop_list = vector("list", length(grid_vect))
for(i in seq_along(grid_vect)) {
  prop_i = read.csv(paste0(path, "props/land_cover_percentage_grid_", i, ".csv"), header = T)
  colnames_to_fix = colnames(prop_i)[1:76]
  colnames_fixed = paste0("c", str_pad(gsub("class_", "", colnames_to_fix, fixed = T), 2, side = "left", pad = "0"))
  colnames(prop_i) = c(colnames_fixed, "grid")
  prop_i = prop_i[order(names(prop_i))]
  prop_list[[i]] = prop_i
  if(ncol(prop_i) != 77) {
    cat("Warning: grid", i, "has", ncol(prop_i), "columns instead of 77\n", file = log_file, append = T)
  } 
  cat("Read grid", i, "\n")
}
prop_df = bind_rows(prop_list)
write.csv(prop_df, paste0(path, "land_cover_percentage.csv"), row.names = F)

prop_df_new = prop_df %>%
  mutate(forest = c01 + c03 + c06,
         agriculture = c14 + c18 + c19 + c39 + c20 + c40 + c62 + c41 + c36 + c46 + c47 + c35 + c48,
         pasture = c15,
         plantation = c09,
         grassland = c04 + c12,
         wetland = c05 + c49 + c11 + c32 + c50) %>%
  mutate(other = 100 - forest - agriculture - pasture - plantation - grassland - wetland) %>%
  dplyr::select(grid, forest, agriculture, pasture, plantation, grassland, wetland, other) %>%
  mutate(score = (agriculture * 1 + pasture * 1 + plantation * 0.5) / 100)
write.csv(prop_df_new, paste0(path, "land_cover_percentage_reclassified.csv"), row.names = F)

grid_vect$score = prop_df_new$score
land_suitability = rasterize(grid_vect, bioclim[[1]], field = "score")
varnames(land_suitability) = "land_cover"
writeRaster(land_suitability, paste0(path, "rasters/land_suitability_score.tif"), overwrite = T)
