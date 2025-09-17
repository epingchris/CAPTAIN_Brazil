rm(list = ls())

#Setup ----
library(parallel)
library(future) #parallelise lapply() : future_lapply()
library(future.apply) #parallelise lapply(): future_lapply()
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)
library(fuzzyjoin)
library(ENMTools) #raster.cor.plot, raster.cor.matrix, trimdupes.by.raster
library(flexsdm)
library(concaveman) #concaveman
library(sf)

save_path = "/maps/epr26/sdm_captain_out/"

#Set optional user-selected project(s) to run
args = commandArgs(trailingOnly = T)
samp_size_df = read.csv(paste0(save_path, "species_sample_size.csv"), header = T)
n_sp = nrow(samp_size_df)

if (length(args) == 0) {
  #all projects in project_var.csv
  n0 = 1
  n1 = n_sp
} else if (length(args) == 1) {
  n0 = 1
  n1 = as.numeric(args)
} else if (length(args) == 2) {
  n0 = as.numeric(args[1])
  n1 = as.numeric(args[2])
} else {
  stop("Error: incorrect number of arguments")
}
cat("Running species from", n0, "to", n1, "\n")

#Set parallelise plan
# parallel::detectCores(logical = F) #256
# plan(multisession, workers = 20)

#Get data
aoi_proj = vect(paste0(save_path, "atlantic_forest_global_200.geojson")) %>% project("EPSG:3857")
aoi_bbox = vect(paste0(save_path, "atlantic_forest_global_200_bbox.geojson"))
land = vect(paste0(save_path, "aoi_land.geojson"))
sp_occ_bbox = vect(paste0(save_path, "SpeciesOccurrenceData_bbox.geojson"))
bioclim = rast(paste0(save_path, "bioclim_reduced.tif"))
sp_occ_list = readRDS(paste0(save_path, "species_occurrence_thinned.rds"))
keep_vars = c(1, 2, 7, 12, 15, 18, 19) #indices of selected bioclim variables after collinearity test

#calculate range of bioclim variables across the AOI
range_bioclim = apply(values(bioclim), 2, function(x) range(x, na.rm = T)) %>%
  t() %>%
  as.data.frame() %>%
  rename(min = V1, max = V2) %>%
  mutate(range = max - min)

#Run models ----
a0 = Sys.time()
#future_lapply(seq(n0, n1, 1), function(i) {
for(i in seq(n0, n1, 1)) { #change to n_sp for full run
  a = Sys.time()

  samp_size_info = samp_size_df[i, ]
  sp_name = samp_size_info$sp_name
  samp_size = samp_size_info$n_used
  
  if(samp_size < 2) { #skip species with less than 2 occurrence points
    cat("Skipping species", i, "with less than 2 occurrence points\n")
    next
  } else {
    cat("Running species", i, "with", samp_size, "occurrence points:", sp_name, "\n")
  }

  i_pad = str_pad(i, 4, side = "left", pad = "0") #haha
  sp_occ_sel = sp_occ_bbox[sp_occ_bbox$tax == sp_name, ]
  sp_occ_used = sp_occ_bbox[sp_occ_bbox$index %in% sp_occ_list[[i]], ]
  #set training extent (other option: flexsdm::calib_area(), more simplistic)

  if(file.exists(paste0(save_path, "ashape/ashape_", i, ".geojson"))) {
    ashape = vect(paste0(save_path, "ashape/ashape_", i, ".geojson"))
    cat("Alpha shape loaded from file for species", i, "\n")  
  } else {
    cat("Creating new alpha shape for species", i, "\n")
    ashape = st_as_sf(sp_occ_used) %>%
      concaveman::concaveman(concavity = 3) %>%
      st_cast("POLYGON") %>%
      vect()
  
    #buffer with 2 x median inter-point distance, then crop to land
    dist = terra::distance(sp_occ_used, unit = "m", method = "geo") %>% as.matrix()
    diag(dist) = NA
    dist[upper.tri(dist)] = NA
    interdist = median(as.vector(dist), na.rm = T)
    ashape_buff = terra::buffer(ashape, interdist * 2) %>% #width unit is in meter!
      intersect(land)
    cat("New alpha shape created for species", i, "\n")
    writeVector(ashape_buff, paste0(save_path, "ashape/ashape_buff_", i, ".geojson"), overwrite = T)
  }
  #create alpha shape with alpha = 3
  
  #sample background points (other option: flexsdm::sample_background)
  #check if background points already sampled: g file exists
  if (file.exists(paste0(save_path, "bg/bg_", i, ".geojson"))) {
    bg = vect(paste0(save_path, "bg/bg_", i, ".geojson"))
    cat("Background points loaded from file for species", i, "\n")
  } else {
    bg = spatSample(ashape_buff, 10000, "random")
    cat("Sampling new background points for species", i, "\n")
    writeVector(bg, paste0(save_path, "bg/bg_", i, ".geojson"), overwrite = T)
  }
  
  #visualize to verify
  plot_sample = ggplot() +
    geom_spatraster(data = bioclim[[1]]) +
    geom_spatvector(data = aoi_proj, color = "red", alpha = 0, linewidth = 0.5) +
    geom_spatvector(data = bg, color = "blue", size = 0.05) +
    geom_spatvector(data = sp_occ_sel, color = "orange", size = 1) +
    geom_spatvector(data = sp_occ_used, color = "green", size = 0.5) +
    geom_spatvector(data = ashape, color = "white", alpha = 0) +
    geom_spatvector(data = ashape_buff, color = "red", alpha = 0) +
    scale_fill_continuous(type = "viridis") + #scale_fill_grass_c somehow doesn't work
    theme_bw()
  ggsave(filename = paste0(save_path, "/plots/plot_sample_", i_pad, ".png"),
          plot = plot_sample, width = 6, height = 8, units = "in", dpi = 300)

  #construct SDM model input data:
  #extract environmental variables at background points and species occurrence points
  bg_var = extract(bioclim, bg, ID = F, xy = T) %>%
    filter(complete.cases(.)) %>%
    mutate(pb = 0)
  sp_occ_var = extract(bioclim, sp_occ_used, ID = F, xy = T) %>%
    filter(complete.cases(.)) %>%
    mutate(pb = 1)
  sdm_data_xy = rbind(sp_occ_var, bg_var)
  
  #randomly partition data into training-testing datasets for cross-validation
  sdm_part = flexsdm::part_random(
      data = sdm_data_xy,
      pr_ab = "pb",
      method = c(method= "kfold", folds = ifelse(samp_size < 30, "2", "4")))
  sdm_part_vect = sdm_part %>%
      mutate(.part = as.factor(.part)) %>%
      vect(geom = c("x", "y"), crs = crs(bioclim))
  plot_part = ggplot() +
    geom_spatvector(data = land, alpha = 0) +
    geom_spatvector(data = filter(sdm_part_vect, pb == 0), aes(col = .part), cex = 0.5) +
    geom_spatvector(data = filter(sdm_part_vect, pb == 1), aes(col = .part), cex = 1) +
    scale_color_manual(values = c("red", "green", "blue", "purple")) +
    theme_bw()
  ggsave(filename = paste0(save_path, "/plots/plot_part_", i_pad, ".png"),
          plot = plot_part, width = 6, height = 8, units = "in", dpi = 300)
  sdm_data = sdm_data_xy %>%
    mutate(part = sdm_part$.part) %>%
    dplyr::select(!c(x, y))
  
  #Calculate range coverage of background points
  range_bg = apply(bg_var[, 1:7], 2, range) %>%
    t() %>%
    as.data.frame() %>%
    rename(min = V1, max = V2) %>%
    mutate(range = max - min)

  range_coverage = range_bg$range / range_bioclim$range %>%
    t() %>%
    as.data.frame()
  rownames(range_coverage) = i
  colnames(range_coverage) = paste0("wc2.1_5m_bio_", keep_vars)

  if(samp_size >= 15) { #3781 species with more than 7 occurrence points per partition
    #GAM model
    m_gam = flexsdm::fit_gam(
      data = sdm_data,
      response = "pb",
      predictors = paste0("wc2.1_5m_bio_", keep_vars),
      partition = "part",
      thr = "max_sens_spec")
  } else if (samp_size > 1) {
    #GAM model: ensemble of small models approach
    m_gam = flexsdm::esm_gam(
      data = sdm_data,
      response = "pb",
      predictors = paste0("wc2.1_5m_bio_", keep_vars),
      partition = "part",
      thr = "max_sens_spec")
  } else {
    m_gam = NULL
  }

  #generate predictions
  if(!is.null(m_gam)) {
    pred_gam = flexsdm::sdm_predict(
      models = m_gam,
      pred = bioclim,
      predict_area = land
    )
  } else {
    pred_gam = list(NULL)
  }

  if (!is.null(pred_gam[[1]])) {
    writeRaster(pred_gam[[1]], paste0(save_path, "/rasters/pred_gam_", i, ".tif"), overwrite = T)
    plot_gam = ggplot() +
      geom_spatraster(data = pred_gam[[1]]) +
      geom_spatvector(data = sp_occ_used, cex = 0.5, col = "red") +
      scale_fill_continuous(limits = c(0, 1), type = "viridis") +
      labs(title = "GAM", fill = "Probability") +
      theme_bw()
    ggsave(filename = paste0(save_path, "/plots/plot_m_gam_", i_pad, ".png"),
            plot = plot_gam, width = 6, height = 8, units = "in", dpi = 300)
  } else {
    cat("GAM model failed for species", sp_name, "\n")
  }

  b = Sys.time()
  cat("Model completed for species", i, ", run time:", as.numeric(b - a, units = "secs"), "secs\n")

  sdm_out = list(samp_size_info = samp_size_info,
                 sdm_data_xy = sdm_data_xy,
                 sdm_data = sdm_data,
                 range_coverage = range_coverage,
                 model = m_gam)
  saveRDS(sdm_out, paste0(save_path, "/models/sdm_model_outputs_", i, ".rds"))

}

b0 = Sys.time()
cat("Total time:", as.numeric(b0 - a0, units = "secs"), "secs\n")
