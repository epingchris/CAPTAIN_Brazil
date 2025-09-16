rm(list = ls())

#Setup ----
library(httpgd) #view plots in VS Code
library(parallel)
library(future) #parallelise lapply() : future_lapply()
library(future.apply) #parallelise lapply(): future_lapply()
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)
library(geodata) #world, worldclim_global
library(fuzzyjoin)
library(predicts)
library(ENMTools) #raster.cor.plot, raster.cor.matrix, trimdupes.by.raster
library(flexsdm)
library(concaveman) #concaveman
library(sf)

hgd()
save_path = "/maps/epr26/sdm_captain_out/"

#Set optional user-selected project(s) to run
args = commandArgs(trailingOnly = T)
samp_size_df = read.csv(paste0(save_path, "species_sample_size.csv"), header = T)
n_sp = nrow(samp_size_df)

if (length(args) == 2) {
  n0 = as.numeric(args[1])
  n1 = as.numeric(args[2])
} else {
  stop("Error: incorrect number of arguments")
}
cat("Saving outputs for species from", n0, "to", n1, "\n")

#check range coverage of background points
range_coverage_df = read.csv(paste0(save_path, "range_coverage.csv"), header = T)
for(i in seq(n0, n1, 1) {
  sdm_out = readRDS(paste0(save_path, "/models/sdm_model_outputs_", i, ".rds"))
  if(!is.null(sdm_out$range_coverage)) {
    range_coverage = sdm_out$range_coverage
  } else {
    range_coverage = data.frame(matrix(NA, nrow = 1, ncol = length(keep_vars)))
    rownames(range_coverage) = i
    colnames(range_coverage) = paste0("wc2.1_5m_bio_", keep_vars)
  }
  if(range_coverage_df$bg_done[i] == 0) { #if not done yet
    range_coverage_df[i, ] = c(range_coverage, 1) #save values and mark as done
  }
})
write.csv(range_coverage_df, paste0(save_path, "range_coverage.csv"), row.names = F)

for(i in seq(n0, n1, 1)) {
  sdm_out = readRDS(paste0(save_path, "/models/sdm_model_outputs_", i, ".rds"))
  sdm_models = sdm_out$models[c("glm", "gam", "raf")]

  # Create ensemble model
  m_ens = fit_ensemble(
       models = sdm_models,
       ens_method = "meanw",
       thr_model = "max_sens_spec",
       metric = "TSS"
       )

  m_ens

  if(!is.null(sdm_out$data)) {
    samp_size = ifelse(sdm_out$samp_size$data_used == "thinning",
                       sdm_out$samp_size$thinned,
                       sdm_out$samp_size$original)
    cat("Species", i, "sample size:", samp_size, "\n")
    plot_sample = sdm_out$diagnostics[[1]]

    plot_part = sdm_out$diagnostics[[3]]
    ggsave(filename = paste0(save_path, "/plots/plot_part_", i, ".png"),
           plot = plot_part, width = 6, height = 8, units = "in", dpi = 300)

    plot_m_max = sdm_out$plots[[1]]
    ggsave(filename = paste0(save_path, "/plots/plot_m_max_", i, ".png"),
           plot = plot_m_max, width = 6, height = 8, units = "in", dpi = 300)

    plot_m_glm = sdm_out$plots[[2]]
    ggsave(filename = paste0(save_path, "/plots/plot_m_glm_", i, ".png"),
           plot = plot_m_glm, width = 6, height = 8, units = "in", dpi = 300)

    plot_m_gam = sdm_out$plots[[3]]
    ggsave(filename = paste0(save_path, "/plots/plot_m_gam_", i, ".png"),
           plot = plot_m_gam, width = 6, height = 8, units = "in", dpi = 300)

    plot_m_net = sdm_out$plots[[4]]
    ggsave(filename = paste0(save_path, "/plots/plot_m_net_", i, ".png"),
           plot = plot_m_net, width = 6, height = 8, units = "in", dpi = 300)

    plot_m_raf = sdm_out$plots[[5]]
    ggsave(filename = paste0(save_path, "/plots/plot_m_raf_", i, ".png"),
           plot = plot_m_raf, width = 6, height = 8, units = "in", dpi = 300)
    cat("Completed species", i, "\n")
  } else {
    cat("Skipped species", i, "\n")
  }
}
cat("Completed saving outputs for species from", n0, "to", n1, "\n")
