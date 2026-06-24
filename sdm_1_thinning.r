rm(list = ls())

#Setup ----
library(renv)
library(dplyr)
library(terra)
library(tidyterra)
library(geodata) #world, worldclim_global
library(fuzzyjoin)
library(ENMTools) #raster.cor.plot, raster.cor.matrix, trimdupes.by.raster
library(flexsdm)
library(sf)

on_cluster = nchar(Sys.getenv("SCRATCH_PATH")) > 0

if (on_cluster) {
  cat("=== Running on CLUSTER ===\n")
  
  in_path = paste0(Sys.getenv("PROJECT_PATH"), "/data_in/")
  out_path = paste0(Sys.getenv("PROJECT_PATH"), "/CAPTAIN_Brazil_out/")
  scratch_path = Sys.getenv("SCRATCH_PATH")
  
  cat("Cluster environment detected:\n")
  cat("  PROJECT_PATH:", Sys.getenv("PROJECT_PATH"), "\n")
  cat("  SCRATCH_PATH:", Sys.getenv("SCRATCH_PATH"), "\n")
  
} else {
  cat("=== Running LOCALLY ===\n")
  cat("  Working directory:", getwd(), "\n")

  in_path = "../DATA/CAPTAIN_Brazil_in/"
  out_path = "../DATA/CAPTAIN_Brazil_out/"
  scratch_path = tempdir()
}

res_val = 30 #30-m resolution
crs_val = "EPSG:3857"
aoi_filepath = paste0(in_path, "atlantic_forest_global_200.geojson")
env_filepath = paste0(in_path, "env_final.tif")
sp_filepath = paste0(in_path, "200k/threatened_occurrences_for_sdm_200k.csv")

#Read AOI shapefile
aoi = vect(aoi_filepath) |> project(crs_val)
#aoi_bbox = as.polygons(ext(aoi), crs = crs(aoi))
#writeVector(aoi_bbox, paste0(save_path, "aoi_bbox.geojson"), overwrite = T)

#Read environmental data
env = terra::rast(env_filepath)

#Read species occurrence data
sp_occ = read.csv(sp_filepath) |>
  mutate(x = longitude, y = latitude, index = row_number())
sp_occ_vect = sp_occ |>
  vect(geom = c("x", "y"), crs = "EPSG:4326") |>
  project(crs_val)
tax_df = as.data.frame(table(sp_occ$binomial)) |>
  rename(tax = Var1, count = Freq)
n_sp = nrow(tax_df)


#Perform thinning and examine sample size
sp_occ_list = vector("list", n_sp)
samp_size_df = data.frame(index = numeric(), sp_name = character(),
                          original = numeric(), thinned = numeric(), thin_perc = numeric(),
                          data_used = character(), n_used = numeric(), flag = character())

for(i in seq_len(n_sp)) {
  a = Sys.time()
  sp_name = as.character(tax_df$binomial[i])
  sp_occ_sel = sp_occ_bbox[sp_occ_bbox$tax == sp_name, ]
  n_orig = nrow(sp_occ_sel)
  
  #optional: geographical distributions of occurrence data and features that may cause spatial biases
  #can be explored using visualization tools in the ‘sampbias' package
  
  #perform spatial-grid thinning for abundant species
  n_thin = NA
  if(n_orig >= 30) {
    thin_method = "trimdupes" #other option: "occfilt" 
    if(thin_method == "trimdupes") {
      sp_occ_thin = ENMTools::trimdupes.by.raster(sp_occ_sel, bioclim) #removes duplicates based on raster cells
      n_thin = nrow(sp_occ_thin)
    } else {
      sp_occ_thin = flexsdm::occfilt_geo(data = crds(sp_occ_sel) %>% as.data.frame(),
                                         x = "x", y = "y",
                                         env_layer = bioclim,
                                         method = c("cellsize", 1),
                                         prj = crs(bioclim)) %>% #
        vect(geom = c("x", "y"), crs = crs(bioclim))
      n_thin = nrow(sp_occ_thin)
    }
  }

  #add attributes back in
  if(n_orig < 30 | n_thin < 30) { #rare species, do not thin
    use = "original"
    n_used = n_orig
    sp_occ_used = sp_occ_sel
  } else {
    use = "thinned"
    n_used = n_thin
    sp_occ_used = sp_occ_sel[geom(sp_occ_sel) %in% geom(sp_occ_thin)]
  }
  
  #flag data point abundance
  samp_size_flag = ifelse(n_used >= 30, "abundant", ifelse(n_used >= 15, "sparse", "insufficient"))
  samp_size_df[i, ] = data.frame(index = i, sp_name = sp_name,
                                 original = n_orig, thinned = n_thin, thin_perc = round((n_thin / n_orig) * 100, 1),
                                 data_used = use, n_used = n_used, flag = samp_size_flag)
  sp_occ_list[[i]] = sp_occ_used$index

  b = Sys.time()
  cat(i, "-", sp_name, ":", b - a, "s\n")
}

write.csv(samp_size_df, paste0(in_path, "species_sample_size.csv"), row.names = F)
saveRDS(sp_occ_list, paste0(in_path, "species_occurrence_thinned.rds"))