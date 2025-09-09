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

#Set parallelise plan
parallel::detectCores(logical = F) #256
plan(multisession, workers = 20)

#Get data
aoi_proj = vect("atlantic_forest_global_200.geojson") %>% project("EPSG:3857")
aoi_bbox = vect("atlantic_forest_global_200_bbox.geojson")
land = vect("aoi_land.geojson")

sp_occ_bbox = vect("SpeciesOccurrenceData_bbox.geojson")
bioclim = rast("bioclim_reduced.tif")

samp_size_df = read.csv("species_sample_size.csv", header = T)
sp_occ_list = readRDS("species_occurrence_thinned.rds")

keep_vars = c(1, 2, 7, 12, 15, 18, 19) #indices of selected bioclim variables after collinearity test


#Run models ----
set.seed(1963)
sdm_out = future_lapply(seq_len(n_sp), function(i) {
  samp_size_info = samp_size_df[i, ]
  sp_name = samp_size_info$sp_name
  samp_size = samp_size_info$occfilt

  if(samp_size >= 30) {
    sp_occ_sel = sp_occ_bbox[sp_occ_bbox$tax == sp_name, ]
    sp_occ_thin = sp_occ_bbox[sp_occ_bbox$index %in% sp_occ_list[[i]], ]
    #set training extent (other option: flexsdm::calib_area(), more simplistic)

    #create alpha shape with alpha = 3
    ashape = st_as_sf(sp_occ_thin) %>%
      concaveman::concaveman(concavity = 3) %>%
      st_cast("POLYGON") %>%
      vect()
    
    #buffer with 2 x median inter-point distance, then crop to land
    dist = terra::distance(sp_occ_thin, unit = "m", method = "geo") %>% as.matrix()
    diag(dist) = NA
    dist[upper.tri(dist)] = NA
    interdist = median(as.vector(dist), na.rm = T)
    ashape_buff = terra::buffer(ashape, interdist * 2) %>% #width unit is in meter!
      intersect(land)
    
    #ecological clipping: restrict sampling to areas within a specified environmental distance
    #@@@to do@@@ (‘biomod2', ‘ENMTML', ‘flexsdm')
    
    #sample background points (other option: flexsdm::sample_background)
    bg = spatSample(ashape_buff, 1000, "random")
    
    #visualize to verify
    plot_sample = ggplot() +
      geom_spatraster(data = bioclim[[1]]) +
      geom_spatvector(data = aoi_proj, color = "red", alpha = 0, linewidth = 0.5) +
      geom_spatvector(data = bg, color = "blue", size = 0.05) +
      geom_spatvector(data = sp_occ_sel, color = "orange", size = 1) +
      geom_spatvector(data = sp_occ_thin, color = "green", size = 0.5) +
      geom_spatvector(data = ashape, color = "white", alpha = 0) +
      geom_spatvector(data = ashape_buff, color = "red", alpha = 0) +
      scale_fill_continuous(type = "viridis") #scale_fill_grass_c somehow doesn't work
    plot_sample

    #construct SDM model input data:
    #extract environmental variables at background points and species occurrence points
    bg_var = extract(bioclim, bg, ID = F, xy = T) %>%
      filter(complete.cases(.)) %>%
      mutate(pb = 0)
    sp_occ_var = extract(bioclim, sp_occ_thin, ID = F, xy = T) %>%
      filter(complete.cases(.)) %>%
      mutate(pb = 1)
    sdm_data_xy = rbind(sp_occ_var, bg_var)
    
    #partition data into training-testing using spatial blocks
    sdm_part = flexsdm::part_sblock(
        data = sdm_data_xy,
        bioclim,
        pr_ab = "pb",
        x = "x",
        y = "y",
        n_part = 5)
    sdm_part_vect = sdm_part$part %>%
      mutate(.part = as.factor(.part)) %>%
      vect(geom = c("x", "y"), crs = crs(bioclim))
    plot_part = ggplot() +
      geom_spatvector(data = land, alpha = 0) +
      geom_spatraster(data = sdm_part$grid, alpha = 0.3) +
      geom_spatvector(data = filter(sdm_part_vect, pr_ab == 0), aes(col = .part), cex = 0.5) +
      geom_spatvector(data = filter(sdm_part_vect, pr_ab == 1), aes(col = .part), cex = 1) +
      scale_fill_continuous(type = "gradient") + #scale_fill_grass_c somehow doesn't work
      scale_color_manual(values = c("red", "green", "blue", "purple",  "yellow")) +
      theme_bw()
    plot_part
    sdm_data = sdm_data_xy %>%
      mutate(part = sdm_part$part$.part) %>%
      dplyr::select(!c(x, y))
      
    #MaxEnt model
    m_max = flexsdm::fit_max(
      data = sdm_data,
      response = "pb",
      predictors = paste0("wc2.1_5m_bio_", keep_vars),
      partition = "part",
      thr = "max_sens_spec")
    
    #GLM model
    m_glm = flexsdm::fit_glm(
      data = sdm_data,
      response = "pb",
      predictors = paste0("wc2.1_5m_bio_", keep_vars),
      partition = "part",
      thr = "max_sens_spec")
    
    #GAM model
    m_gam = flexsdm::fit_gam(
      data = sdm_data,
      response = "pb",
      predictors = paste0("wc2.1_5m_bio_", keep_vars),
      partition = "part",
      thr = "max_sens_spec")

    #neural network model
    m_net = flexsdm::fit_net(
      data = sdm_data,
      response = "pb",
      predictors = paste0("wc2.1_5m_bio_", keep_vars),
      partition = "part",
      thr = "max_sens_spec")

    #random forest model
    m_raf = flexsdm::fit_raf(
      data = sdm_data,
      response = "pb",
      predictors = paste0("wc2.1_5m_bio_", keep_vars),
      partition = "part",
      thr = "max_sens_spec")
    
    
    #generate predictions
    pred_list = sdm_predict(
      models = list(m_max, m_glm, m_gam, m_net, m_raf),
      pred = bioclim,
      predict_area = land
    )

    plot_max = ggplot() +
      geom_spatraster(data = pred_list$max) +
      geom_spatvector(data = sp_occ_thin, cex = 0.5, col = "red") +
      scale_fill_continuous(limits = c(0, 1), type = "viridis") +
      labs(main = "Maxent", fill = "Probability") +
      theme_bw()
    plot_max
    
    plot_glm = ggplot() +
      geom_spatraster(data = pred_list$glm) +
      geom_spatvector(data = sp_occ_thin, cex = 0.5, col = "red") +
      scale_fill_continuous(limits = c(0, 1), type = "viridis") +
      labs(main = "GLM", fill = "Probability") +
      theme_bw()
    plot_glm
    
    plot_gam = ggplot() +
      geom_spatraster(data = pred_list$gam) +
      geom_spatvector(data = sp_occ_thin, cex = 0.5, col = "red") +
      scale_fill_continuous(limits = c(0, 1), type = "viridis") +
      labs(main = "GAM", fill = "Probability") +
      theme_bw()
    plot_gam

    plot_net = ggplot() +
      geom_spatraster(data = pred_list$net) +
      geom_spatvector(data = sp_occ_thin, cex = 0.5, col = "red") +
      scale_fill_continuous(limits = c(0, 1), type = "viridis") +
      labs(main = "GAM", fill = "Probability") +
      theme_bw()
    plot_net

    plot_raf = ggplot() +
      geom_spatraster(data = pred_list$raf) +
      geom_spatvector(data = sp_occ_thin, cex = 0.5, col = "red") +
      scale_fill_continuous(limits = c(0, 1), type = "viridis") +
      labs(main = "Random forest", fill = "Probability") +
      theme_bw()
    plot_raf
    
  } else {
    sdm_data_xy = NULL
    plot_sample = NULL
    sdm_part = NULL
    plot_part = NULL
    m_max = NULL
    m_glm = NULL
    m_gam = NULL
    m_net = NULL
    m_raf = NULL
    pred_list = NULL
    plot_max = NULL
    plot_glm = NULL
    plot_gam = NULL
    plot_net = NULL
    plot_raf = NULL
  }
  
  return(list(samp_size = samp_size_info,
              data = sdm_data_xy,
              diagnostics = list(plot_sample, sdm_part, plot_part),
              models = list(m_max, m_glm, m_gam, m_net, m_raf),
              predictions = pred_list,
              plots = list(plot_max, plot_glm, plot_gam, plot_net, plot_raf))
         )
})

saveRDS(sdm_out, "sdm_model_outputs.rds")