rm(list = ls())

#Setup ----
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)

path = "/maps/epr26/captain_brazil/ideal_250m/"
bioclim = rast(paste0(path, "rasters/bioclim_reduced.tif"))
orig_file = list.files(path = paste0("/maps/epr26/captain_brazil/af_10km/preds_valid/"))
orig_index = as.numeric(gsub("\\D", "", orig_file))
sp_info = read.csv("/maps/epr26/captain_brazil/af_10km/species_info.csv", header = T)
sp_retained = read.table(paste0(path, "species_retained.txt"), header = F)$V1
sp_info_retained = subset(sp_info, index %in% sp_retained & index %in% orig_index)
write.csv(sp_info_retained, paste0(path, "sp_info_retained.csv"), row.names = F)

for(i in seq_along(orig_file)) {
  filename = orig_file[i]
  index = as.numeric(gsub("\\D", "", filename))
  if(index %in% sp_retained) {
    gam_i = rast(paste0("/maps/epr26/captain_brazil/af_10km/preds_valid/", filename))
    gam_resamp = resample(gam_i, bioclim[[1]], method = "bilinear",
                          filename = paste0(path, "resamp/", filename), overwrite = T)
    cat("Done ", i, " of ", length(orig_file), ": index ", index, "\n", sep = "")
  } else {
    cat("Skipping ", i, " of ", length(orig_file), ": index ", index, " not in region\n", sep = "")
  }
}
