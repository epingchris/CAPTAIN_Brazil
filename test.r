rm(list = ls())

#Setup ----
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)
library(geodata) #geodata::travel_time, geodata::crop_spam

path = "/maps/epr26/sdm_captain_out/"
samp_size_df = read.csv(paste0(path, "species_sample_size.csv"), header = T)
n_sp = nrow(samp_size_df)
to_run = sapply(1:n_sp, function(i) ifelse(!file.exists(paste0(path, "/models/sdm_model_outputs_", i, ".rds")), T, F))
to_run_vec = 1:n_sp
to_run_vec = to_run_vec[which(to_run)]
length(to_run_vec)

log = read.table("/maps/epr26/sdm_captain_out/logs/output_sdm_model_31_4948.txt", sep = "\n")

esm_lines = grep("Rerunning species", log$V1)
#extract the lines and get the species ID after "Rerunning species XXX"
esm_ids = sapply(esm_lines, function(i) {
  line = log$V1[i]
  match = str_extract(line, "(?<=Rerunning species )\\d+")
  if (length(match) > 0) {
    as.integer(match)
  } else {
    NA
  }
})

finished_lines = grep("Model completed", log$V1)
finished_ids = sapply(finished_lines, function(i) {
  line = log$V1[i]
  match = str_extract(line, "(?<=Model completed for species )\\d+")
  if (length(match) > 0) {
    as.integer(match)
  } else {
    NA
  }
})

esm_finished = esm_ids[which(esm_ids %in% finished_ids)]


samp_size_df = read.csv(paste0(save_path, "species_sample_size.csv"), header = T)
esm_total = subset(samp_size_df, n_used > 2 & n_used < 15)$index
esm_to_run = esm_total[which(!(esm_total %in% esm_finished))]
n_sp = nrow(samp_size_df)

for(i in 1:n_sp) {
    to_run[i] = ifelse(!file.exists(paste0(save_path, "/models/sdm_model_outputs_", i, ".rds")), T, F)
}

