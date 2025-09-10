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


save_path = "/maps/epr26/sdm_captain_out/"

for(i in seq_along(sdm_out)) {
  sdm_out = readRDS("sdm_model_outputs_", i, ".rds")

  if(!is.null(sdm_out$sdm_data_xy)) {
    cat("Species", i, "sample size:", sdm_out[[i]]$samp_size$sample_size, "\n")
    plot_sample = sdm_out[[i]]$diagnostics$plot_sample
    ggsave(filename = paste0(save_path, "/plots/plot_sample_", i, ".png"), plot = plot_sample, width = 6, height = 8, units = "in", dpi = 300)

    plot_part = sdm_out[[i]]$diagnostics$sdm_part
    ggsave(filename = paste0(save_path, "/plots/plot_part_", i, ".png"), plot = plot_part, width = 6, height = 8, units = "in", dpi = 300)

    plot_max = sdm_out[[i]]$plots$plot_max
    ggsave(filename = paste0(save_path, "/plots/plot_m_max_", i, ".png"), plot = plot_max, width = 6, height = 8, units = "in", dpi = 300)

    plot_glm = sdm_out[[i]]$plots$plot_glm
    ggsave(filename = paste0(save_path, "/plots/plot_m_glm_", i, ".png"), plot = plot_glm, width = 6, height = 8, units = "in", dpi = 300)

    plot_gam = sdm_out[[i]]$plots$plot_gam
    ggsave(filename = paste0(save_path, "/plots/plot_m_gam_", i, ".png"), plot = plot_gam, width = 6, height = 8, units = "in", dpi = 300)

    plot_net = sdm_out[[i]]$plots$plot_net
    ggsave(filename = paste0(save_path, "/plots/plot_m_net_", i, ".png"), plot = plot_net, width = 6, height = 8, units = "in", dpi = 300)

    plot_raf = sdm_out[[i]]$plots$plot_raf
    ggsave(filename = paste0(save_path, "/plots/plot_m_raf_", i, ".png"), plot = plot_raf, width = 6, height = 8, units = "in", dpi = 300)
  }
}