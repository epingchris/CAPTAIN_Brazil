rm(list = ls())

#Execute with: Rscript sdm_2_model.r [n0 n1] 2>&1 | tee output_sdm_model_n0_n1.txt

#Setup ----
#library(parallel)
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

dir_path = "/maps/epr26/captain_brazil/"
proj_folder = "ideal_250m" #af_10km, ideal_250m
proj_path = paste0(dir_path, proj_folder, "/")
if (!dir.exists(proj_path)) dir.create(proj_path)
if (!dir.exists(paste0(proj_path, "models/"))) dir.create(paste0(proj_path, "models/"))
if (!dir.exists(paste0(proj_path, "plots_glm/"))) dir.create(paste0(proj_path, "plots_glm/"))
if (!dir.exists(paste0(proj_path, "preds_glm/"))) dir.create(paste0(proj_path, "preds_glm/"))
if (!dir.exists(paste0(proj_path, "plots/"))) dir.create(paste0(proj_path, "plots/"))
if (!dir.exists(paste0(proj_path, "ashape/"))) dir.create(paste0(proj_path, "ashape/"))
if (!dir.exists(paste0(proj_path, "bg/"))) dir.create(paste0(proj_path, "bg/"))


#Set optional user-selected project(s) to run
args = commandArgs(trailingOnly = T)
sp_info = read.csv(paste0(proj_path, "species_info.csv"), header = T)
n_sp = nrow(sp_info)

if (length(args) == 0) {
  #all projects in species_info.csv
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

log_file = paste0(dir_path, "logs/output_sdm_model_", n0, "_", n1, ".txt")

#Set parallelise plan
# parallel::detectCores(logical = F) #256
plan(multicore, workers = 20)

Sys.setenv(OMP_NUM_THREADS = 4, OPENBLAS_NUM_THREADS = 4)
#maximum number of threads used by OpenMP and OpenBLAS, usually under the hood for R


#Define model ----
run_sdm = function(i) {
  a = Sys.time()

  #Get data
  dir_path = "/maps/epr26/captain_brazil/"
  proj_path = "/maps/epr26/captain_brazil/ideal_250m/"

  aoi_proj = vect(paste0(dir_path, "aoi_proj.geojson"))
  land = vect(paste0(dir_path, "aoi_land.geojson"))
  sp_occ_bbox = vect(paste0(dir_path, "spocc_bbox.geojson"))

  bioclim = rast(paste0(proj_path, "rasters/bioclim.tif"))
  keep_vars = c(1, 2, 7, 12, 15, 18, 19) #indices of selected bioclim variables after collinearity test
  names(bioclim) = paste0("wc_", keep_vars)
  sp_info = read.csv(paste0(proj_path, "species_info.csv"), header = T)
  sp_occ_list = readRDS(paste0(proj_path, "species_occurrence_thinned.rds"))
  land_sub = vect(paste0(proj_path, "aoi_land.geojson"))

  sp_info_i = sp_info[i, ]
  sp_name = sp_info_i$tax
  samp_size = sp_info_i$n_used
  sp_ind = sp_info_i$sp_ind

  i_pad = str_pad(sp_ind, 4, side = "left", pad = "0") #haha
  sp_occ_sel = sp_occ_bbox[sp_occ_bbox$tax == sp_name, ]
  sp_occ_used = sp_occ_bbox[sp_occ_bbox$index %in% sp_occ_list[[i]], ]

  if(samp_size <= 2) { #skip species with less than or equal to 2 occurrence points
    cat("Skipping species", sp_ind, "with less than or equal to 2 occurrence points\n", file = log_file, append = T)
    sdm_out = list(sp_info = sp_info_i,
                   sdm_data_xy = NULL,
                   sdm_data = NULL,
                  #  range_coverage = NULL,
                   model = NULL)
    saveRDS(sdm_out, paste0(proj_path, "/models/sdm_model_outputs_", sp_ind, ".rds"))
    return(NULL)
  } else {
    if (file.exists(paste0(proj_path, "preds_glm/pred_", sp_ind, ".tif")) &
        file.exists(paste0(proj_path, "/models/sdm_model_outputs_", sp_ind, ".rds"))) {
        cat("Skipping species", sp_ind, "as model output already exists:", sp_name, "\n", file = log_file, append = T)
        return(NULL)
    } else {
      cat("Running species", sp_ind, "with", samp_size, "occurrence points:", sp_name, "\n", file = log_file, append = T)
    }
  }

  #calculate range of bioclim variables across the AOI
  # range_bioclim = apply(values(bioclim), 2, function(x) range(x, na.rm = T)) %>%
  #   t() %>%
  #   as.data.frame() %>%
  #   rename(min = V1, max = V2) %>%
  #   mutate(range = max - min)

  #set training extent (other option: flexsdm::calib_area(), more simplistic)
  if(file.exists(paste0(proj_path, "ashape/ashape_", sp_ind, ".geojson")) &&
     file.exists(paste0(proj_path, "ashape/ashape_buff_", sp_ind, ".geojson"))) {
    ashape = vect(paste0(proj_path, "ashape/ashape_", sp_ind, ".geojson"))
    ashape_buff = vect(paste0(proj_path, "ashape/ashape_buff_", sp_ind, ".geojson"))
    cat("Alpha shape loaded from file for species", sp_ind, "\n", file = log_file, append = T)
  } else {
    cat("Creating new alpha shape for species", sp_ind, "\n", file = log_file, append = T)
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
    cat("New alpha shape created for species", sp_ind, "\n", file = log_file, append = T)
    writeVector(ashape, paste0(proj_path, "ashape/ashape_", sp_ind, ".geojson"), overwrite = T)
    writeVector(ashape_buff, paste0(proj_path, "ashape/ashape_buff_", sp_ind, ".geojson"), overwrite = T)
  }
  #create alpha shape with alpha = 3

  #sample background points (other option: flexsdm::sample_background)
  #check if background points already sampled: g file exists
  if (file.exists(paste0(proj_path, "bg/bg_", sp_ind, ".geojson"))) {
    bg = vect(paste0(proj_path, "bg/bg_", sp_ind, ".geojson"))
    cat("Background points loaded from file for species", sp_ind, "\n", file = log_file, append = T)
  } else {
    bg = spatSample(ashape_buff, 10000, "random")
    cat("Sampling new background points for species", sp_ind, "\n", file = log_file, append = T)
    writeVector(bg, paste0(proj_path, "bg/bg_", sp_ind, ".geojson"), overwrite = T)
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
  ggsave(filename = paste0(proj_path, "/plots/plot_sample_", i_pad, ".png"),
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
  sdm_data = sdm_data_xy %>%
    mutate(part = sdm_part$.part) %>%
    dplyr::select(!c(x, y))

  #Calculate range coverage of background points
  # range_bg = apply(bg_var[, 1:7], 2, range) %>%
  #   t() %>%
  #   as.data.frame() %>%
  #   rename(min = V1, max = V2) %>%
  #   mutate(range = max - min)

  # range_coverage = range_bg$range / range_bioclim$range %>%
  #   t() %>%
  #   as.data.frame()
  # rownames(range_coverage) = i
  # colnames(range_coverage) = paste0("wc2.1_5m_bio_", keep_vars)

  if(samp_size >= 15) { #3781 species with more than 7 occurrence points per partition
    #GLM model
    m_glm = flexsdm::fit_glm(
      data = sdm_data,
      response = "pb",
      predictors = paste0("wc_", keep_vars),
      partition = "part",
      thr = "max_sens_spec")
    pred_glm = flexsdm::sdm_predict(
      models = m_glm,
      pred = bioclim,
      predict_area = land_sub
    )
  } else if (samp_size > 2) {
    #GLM model: ensemble of small models approach
    m_glm = flexsdm::esm_glm(
      data = sdm_data,
      response = "pb",
      predictors = paste0("wc2.1_5m_bio_", keep_vars),
      partition = "part",
      thr = "max_sens_spec")
    pred_glm = flexsdm::sdm_predict(
      models = m_glm$esm_model,
      pred = bioclim,
      predict_area = land_sub
    )

  } else {
    m_glm = NULL
  }

  if (!is.null(pred_glm[[1]])) {
    writeRaster(pred_glm[[1]], paste0(proj_path, "/preds_glm/pred_", sp_ind, ".tif"), overwrite = T)
    plot_pred = ggplot() +
      geom_spatraster(data = pred_glm[[1]]) +
      scale_fill_continuous(limits = c(0, 1), type = "viridis") +
      labs(title = "GLM", fill = "Probability") +
      theme_bw()
    if(proj_folder == "af_10km") {
      plot_pred = plot_pred +
        geom_spatvector(data = sp_occ_used, cex = 0.5, col = "red")
    }
    ggsave(filename = paste0(proj_path, "/plots_glm/plot_", i_pad, ".png"),
           plot = plot_pred, width = 6, height = 8, units = "in", dpi = 300)
  } else {
    cat("GLM model failed for species", sp_name, "\n", file = log_file, append = T)
  }

  b = Sys.time()
  cat("Model completed for species", sp_ind, ", run time:", difftime(b, a, units = "secs"), "secs\n", file = log_file, append = T)

  sdm_out = list(sp_info = sp_info_i,
                 sdm_data_xy = sdm_data_xy,
                 sdm_data = sdm_data,
                #  range_coverage = range_coverage,
                 model = m_glm)
  saveRDS(sdm_out, paste0(proj_path, "/models/sdm_model_outputs_", sp_ind, ".rds"))

}

#Run models ----
run_sdm_safe = function(i) {
  tryCatch({
    run_sdm(i)
    NULL  # return NULL if successful
  }, error = function(e) {
    cat("Error in species ", sp_info[1, ]$sp_ind, ": ", e$message, "\n", file = log_file, append = T)
  })
}

a0 = Sys.time()
future_lapply(seq(n0, n1, 1), run_sdm_safe,
              future.seed = T,
              future.stdout = T)
b0 = Sys.time()
cat("Total time:", difftime(b0, a0, units = "secs"), "secs\n")
