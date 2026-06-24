rm(list = ls())

# Setup ----
library(renv)
library(dplyr)
library(terra)
library(tidyterra)
library(flexsdm)
library(sf)

on_cluster = nchar(Sys.getenv("SCRATCH_PATH")) > 0

if (on_cluster) {
  cat("=== Running on CLUSTER ===\n"); flush.console()
  
  in_path = paste0(Sys.getenv("PROJECT_PATH"), "/data_in/")
  model_path = paste0(Sys.getenv("PROJECT_PATH"), "/AF_foundmod_sdm_out/")
  out_path = paste0(Sys.getenv("SCRATCH_PATH"), "/pred_endangered_60k/")
  
  cat("Cluster environment detected:\n"); flush.console()
  cat("  PROJECT_PATH:", Sys.getenv("PROJECT_PATH"), "\n"); flush.console()
  cat("  SCRATCH_PATH:", Sys.getenv("SCRATCH_PATH"), "\n"); flush.console()
  
} else {
  cat("=== Running LOCALLY ===\n"); flush.console()
  cat("  Working directory:", getwd(), "\n"); flush.console()

  in_path = "../DATA/CAPTAIN_Brazil_in/"
  out_path = "../DATA/CAPTAIN_Brazil_out/"
  scratch_path = tempdir()
}

crs_val = "EPSG:3857"
aoi_filepath = paste0(in_path, "60k_ha_corridor_priority_3857.geojson")
#aoi_filepath = paste0(in_path, "200k_ha_corridor_ideal_3857.geojson")


# Read command line arguments ----
args = commandArgs(trailingOnly = T)
i = as.numeric(args[1])
cat("Species index:", i, "\n"); flush.console()


# Read fitted model ----
if(!file.exists(paste0(model_path, "/sdm_2_fit_out_", i, ".rds"))) {
  cat("RDS object for species", i, "not found: exit script"); flush.console()
  break
} else {
  sdm_out = readRDS(paste0(model_path, "/sdm_2_fit_out_", i, ".rds"))
}

if(is.null(sdm_out$model)) {
  cat("Model empty for species", i, ": exit script"); flush.console()
  break
} else {
  m_out = sdm_out$model
}


# Read species data ----
tax_df = read.csv(paste0(in_path, "species_occurrence_count.csv"), header = T)
sp_name = tax_df[i, ]$sp_name
if(!is.null(sdm_out$sdm_data_xy)) {
  pres_data = subset(sdm_out$sdm_data_xy, pb == 1)
  pres_vect = vect(pres_data, geom = c("x", "y"), crs = crs_val) |> project("EPSG:4326")
  writeVector(pres_vect, paste0(out_path, "pres_", i, ".geojson"), overwrite = T)
  samp_size = nrow(pres_data)
} else {
  cat("No species data found for species", i, ": exit script"); flush.console()
  break
}
cat("===Species index:", i, ", name:", sp_name, ", sample size:", samp_size, "===\n"); flush.console()

if(samp_size <= 2) { #skip species with less than or equal to 2 occurrence points
  cat("Skipping species", i, "with less than 3 occurrences\n"); flush.console()
  break
}


# Read AOI ----
aoi = vect(aoi_filepath) |> project(crs_val)

# Read environmental data ----
env = rast(paste0(in_path, "env_60k.tif"))
names(env) = paste0("env_final_", 1:17)
cat("Environmental data read\n"); flush.console()


# Predict from GAM model ----
#if(samp_size >= 15) {
#  cat("Predicting using GAM model\n"); flush.console()
#  model_used = m_out$model
#} else if (samp_size > 2 & samp_size < 15) {
#  #ensemble of small models approach for species with 3-14 occurrence points
#  cat("Predicting using GAM model with ensemble of small models\n"); flush.console()
#  model_used = m_out$esm_model
#}

if(!is.null(m_out$model)) {
  cat("Predicting using GAM model\n"); flush.console()
  model_used = m_out$model
  pred_terra = terra::predict(
    object = env,
    model = model_used,
    fun = function(model, ...) mgcv::predict.gam(model, ..., type = "response"),
    na.rm = T,
    cores = 1
  )
  pred_gam = list(pred_terra)
} else if (!is.null(m_out$esm_model)) {
  #ensemble of small models approach for species with 3-14 occurrence points
  cat("Predicting using GAM model with ensemble of small models\n"); flush.console()
  model_used = m_out$esm_model
  pred_gam = flexsdm::sdm_predict(
    models = model_used,
    pred = env,
    predict_area = aoi
  )
} else {
  cat("No valid model found for species", i, ": exit script"); flush.console()
  break
}

if (!is.null(pred_gam[[1]])) {
  writeRaster(pred_gam[[1]], paste0(out_path, "pred_", i, ".tif"), overwrite = T)
} else {
  cat("Model prediction failed for species", i, ":", sp_name, "\n"); flush.console()
  break
}

#Extract performance metrics at max_sens_spec threshold
AssignCriterion = function(x) {
  switch(x,
         TSS_mean = 0.5,
         TPR_mean = 0.7,
         TNR_mean = 0.7,
         AUC_mean = 0.8,
         NA_real_)
}
eval_df = m_out$performance |>
  filter(threshold == "max_sens_spec") |>
  dplyr::select(thr_value, TPR_mean, TNR_mean, TSS_mean, AUC_mean) |>
  mutate(species = i) |>
  pivot_longer(cols = ends_with("_mean"), names_to = "metric", values_to = "value") |>
  relocate(species, .before = thr_value) |>
  mutate(criterion = sapply(metric, AssignCriterion)) |>
  mutate(exclude = value < criterion)

eval_part_df = m_out$performance_part |>
  filter(threshold == "max_sens_spec") |>
  dplyr::select(partition, thr_value, n_presences, n_absences, TPR, TNR, TSS, AUC) |>
  mutate(species = i) |>
  relocate(species, .before = thr_value)
write.csv(eval_df, paste0(out_path, "eval_", i, ".csv"), row.names = F)
write.csv(eval_part_df, paste0(out_path, "eval_part_", i, ".csv"), row.names = F)