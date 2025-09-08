#Setup ----
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)
library(geodata) #world, worldclim_global
library(fuzzyjoin)
library(future) #parallelise lapply() : future_lapply()
library(future.apply) #parallelise lapply(): future_lapply()
library(predicts)
library(ENMTools) #raster.cor.plot, raster.cor.matrix, trimdupes.by.raster
library(flexsdm)
library(concaveman) #concaveman
library(sf)

path = "C:/Users/epr26/OneDrive - University of Cambridge/Sweden/lab_Antonelli/"

#Set parallelise plan
plan(multisession, workers = 20)

#Read AOI shapefile
aoi = vect(paste0(path, "atlantic_forest_global_200.geojson")) %>%
  project("EPSG:4326")
aoi_proj = aoi %>% project("EPSG:3857")
aoi_bbox = as.polygons(ext(aoi_proj), crs = crs(aoi_proj))
writeVector(aoi_bbox, paste0(path, "atlantic_forest_global_200_bbox.geojson"), overwrite = T)

#get world map and land boundary
worldmap = geodata::world(path = ".") %>% project("EPSG:4326") #GADM
worldmap_aoi = worldmap %>%
  project("EPSG:3857") %>%
  crop(ext(aoi_bbox))
land = aggregate(worldmap_aoi)

#bioclimatic variable names
biovars = c("Annual Mean Temperature", 
            "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
            "Isothermality (BIO2/BIO7) (×100)",
            "Temperature Seasonality (standard deviation ×100)",
            "Max Temperature of Warmest Month",
            "Min Temperature of Coldest Month",
            "Temperature Annual Range (BIO5-BIO6)",
            "Mean Temperature of Wettest Quarter",
            "Mean Temperature of Driest Quarter",
            "Mean Temperature of Warmest Quarter",
            "Mean Temperature of Coldest Quarter",
            "Annual Precipitation",
            "Precipitation of Wettest Month",
            "Precipitation of Driest Month",
            "Precipitation Seasonality (Coefficient of Variation)",
            "Precipitation of Wettest Quarter",
            "Precipitation of Driest Quarter",
            "Precipitation of Warmest Quarter",
            "Precipitation of Coldest Quarter")

#Environmental data processing ----
bioclim_paths = list.files(path = paste0(path, "wc2.1_5m_bio/"), pattern = "tif", full.names = T)
bioclim = rast(bioclim_paths) %>%
  resample(rast(extent = ext(aoi))) %>% #crop exactly to the AOI; crop() doesn't do this
  project("EPSG:3857", res = 10000) #reproject to 10-km
var_order = names(bioclim) %>% sub("wc2.1_5m_bio_", "", .) %>% as.numeric() %>% order()
bioclim = bioclim[[var_order]]

#examine collinearity
#removing redundant variables (pairwise: ‘ENMTML', ‘flexsdm', ‘modleR', ‘ntbox';
#sequential: ‘fuzzySim', ‘SDMtune', ‘usdm')
#or reducing variable dimensionality through ordination (‘ENMTML', ‘ENMTools', ‘flexsdm', ‘kuenm', ‘ntbox')
#flexsdm::correct_colinvar but there is an error
ENMTools::raster.cor.plot(bioclim) #keep 1, 2, 7, 12, 15, 18, 19
keep_vars = c(1, 2, 7, 12, 15, 18, 19)
bioclim_red = bioclim[[keep_vars]]
ENMTools::raster.cor.plot(bioclim_red)
mat_cor = ENMTools::raster.cor.matrix(bioclim_red)
diag(mat_cor) = NA
biovars[keep_vars]


#Occurrence data processing ----
sp_occ_df = readRDS(paste0(path, "SDM_CAPTAIN/SpeciesOccurrenceData.rds")) %>%
  as.data.frame() %>%
  filter(complete.cases(ddlat) & complete.cases(ddlon)) %>%
  mutate(x = ddlon, y = ddlat)
sp_occ = sp_occ_df %>%
  vect(geom = c("ddlon", "ddlat"), crs = crs(aoi))
sp_occ_proj = sp_occ %>%
  project("EPSG:3857")
writeVector(sp_occ_proj, paste0(path, "SDM_CAPTAIN/SpeciesOccurrenceData.geojson"), overwrite = T)
sp_occ_bbox = crop(sp_occ_proj, ext(aoi_bbox)) #filter species occurrence data by AOI
writeVector(sp_occ_bbox, paste0(path, "SDM_CAPTAIN/SpeciesOccurrenceData_bbox.geojson"), overwrite = T)

#visualize
ggplot() +
  geom_spatvector(data = worldmap, fill = "lightyellow") +
  geom_spatvector(data = sp_occ, color = "orange", size = 0.05) +
  coord_sf(xlim = c(-120, -5), ylim = c(-50, 40)) +
  theme_bw()
  
#examine anomalous coordinates
dim(filter(sp_occ, x > -34.793015)) #many but not all are on islands east of Brazil, 394 entries
dim(filter(sp_occ, x > -10)) #one entry, definitely wrong
dim(filter(sp_occ, x > -20 & x <= -10)) #76 entries: possibly wrong?
dim(filter(sp_occ, x > -30 & x <= -20)) #-20.5, -29.3: Ilha da Trindade; -18.x, -28~29.x: possibly wrong
dim(filter(sp_occ, x > -34.793015 & x <= -30)) #77 entries: Ilha de Fernando de Noronha

dim(filter(sp_occ, x < -85 & y <= 1)) #22 entries: Galapagos Islands
dim(filter(sp_occ, x > -75 & y > 30)) #4 entries: Bermuda Main Island
dim(filter(sp_occ, y > 38)) #3 entries: middle of the US

anomaly_coord = data.frame(x = c(-9.24255, -18.42746, -18.42299, -18.42567, -18.08748, 39.52944),
                           y = c(-8.00000, -18.42746, -29.07226, -29.08331, -28.82678, -99.15207))
anomaly = fuzzyjoin::difference_inner_join(
  sp_occ_df, anomaly_coord,
  by = c("x", "y"),
  max_dist = 1e-5
) %>%
  dplyr::select(!c("x.y", "y.y")) %>%
  rename(x = x.x, y = y.x)


#Build models ----

#examine data abundance
tax_df = as.data.frame(table(sp_occ_bbox$tax)) %>%
  rename(tax = Var1, count = Freq)
n_sp = nrow(tax_df)

samp_size_df = data.frame(original = rep(NA, n_sp), trimdupes = rep(NA, n_sp), occfilt = rep(NA, n_sp), flag = rep("", n_sp))
sp_occ_thin_list = vector("list", n_sp)

out = future_lapply(seq_len(n_sp), function(i) {
  a = Sys.time()
  sp_occ_sel = sp_occ_bbox[sp_occ_bbox$tax == tax_df$tax[i], ]
  samp_size_df$original[i] = nrow(sp_occ_sel)
  
  #geographical distributions of occurrence data and features that may cause spatial biases
  #can be explored using visualization tools in the ‘sampbias' package
  
  #spatial-grid thinning based on one point per grid cell
  sp_occ_thin_last = ENMTools::trimdupes.by.raster(sp_occ_sel, bioclim)
  sp_occ_thin = sp_occ_sel[geom(sp_occ_sel) %in% geom(sp_occ_thin_last)]
  samp_size_df$trimdupes[i] = nrow(sp_occ_thin)

  #spatial-grid thinning based on nearest-neighbor distance larger than grid cell size
  sp_occ_thin_occfilt = flexsdm::occfilt_geo(data = crds(sp_occ_sel) %>% as.data.frame(),
                                             x = "x", y = "y",
                                             env_layer = bioclim,
                                             method = c("cellsize", 1),
                                             prj = crs(bioclim)) %>%
    vect(geom = c("x", "y"), crs = crs(bioclim))
  sp_occ_thin = sp_occ_sel[geom(sp_occ_sel) %in% geom(sp_occ_thin_occfilt)]

  samp_size_df$occfilt[i] = nrow(sp_occ_thin)
  n_occfilt = nrow(sp_occ_thin)
  
  #flag data point abundance
  sp_occ_thin_list[[i]] = sp_occ_thin
  samp_size = nrow(sp_occ_thin)
  samp_size_flag = ifelse(samp_size >= 30, "abundant", ifelse(samp_size >= 15, "rare", "insufficient"))
  b = Sys.time()
  cat(i, ":", as.character(tax_df$tax[i]), ",", b - a, "\n")
}

set.seed(1963)
out = future_lapply(seq_along(tax_df$tax), function(i) {
  if(samp_size >= 30) {
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
      scale_fill_grass_c(palette = "viridis")
    plot_sample

    #construct SDM model input data:
    #extract environmental variables at background points and species occurrence points
    bg_var = extract(bioclim_red, bg, ID = F, xy = T) %>%
      filter(complete.cases(.)) %>%
      mutate(pb = 0)
    sp_occ_var = extract(bioclim_red, sp_occ_thin, ID = F, xy = T) %>%
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
      scale_fill_grass_c(palette = "grey") +
      scale_color_manual(values = c("red", "green", "blue", "purple",  "yellow")) +
      theme_bw()
    plot_part
    sdm_data = sdm_data_xy %>%
      mutate(part = sdm_part$part$.part) %>%
      dplyr::select(!c(x, y))
      
    # sdm_dirs = sdm_directory(main_dir = paste0(path, "SDM_CAPTAIN/"),
    #                          projections = NULL,
    #                          calibration_area = T,
    #                          algorithm = c("gam", "glm", "max", "raf"),
    #                          ensemble = c("mean", "meanthr"))

        
    #bioclimatic envelope model
    m_env = predicts::envelope(subset(spp_occ_var, select = -pb))
    pred_env = predict(bioclim_red, m_env)
    plot_env = ggplot() +
      geom_spatraster(data = pred_env) +
      geom_spatvector(data = sp_occ_thin, cex = 0.5, col = "red") +
      scale_fill_grass_c(limits = c(0, 1), palette = "viridis") +
      labs(main = "Bioclimatic envelope", fill = "Probability") +
      theme_bw()
    plot_env

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
    
    #random forest model
    m_raf = flexsdm::fit_raf(
      data = sdm_data,
      response = "pb",
      predictors = paste0("wc2.1_5m_bio_", keep_vars),
      partition = "part",
      thr = "max_sens_spec")
    
    
    #generate predictions
    pred_list = sdm_predict(
      models = list(m_max, m_glm, m_gam, m_raf),
      pred = bioclim,
      predict_area = land
    )

    plot_max = ggplot() +
      geom_spatraster(data = pred_list$max) +
      geom_spatvector(data = sp_occ_thin, cex = 0.5, col = "red") +
      scale_fill_grass_c(limits = c(0, 1), palette = "viridis") +
      labs(main = "Maxent", fill = "Probability") +
      theme_bw()
    plot_max
    
    plot_glm = ggplot() +
      geom_spatraster(data = pred_list$glm) +
      geom_spatvector(data = sp_occ_thin, cex = 0.5, col = "red") +
      scale_fill_grass_c(limits = c(0, 1), palette = "viridis") +
      labs(main = "GLM", fill = "Probability") +
      theme_bw()
    plot_glm
    
    plot_gam = ggplot() +
      geom_spatraster(data = pred_list$gam) +
      geom_spatvector(data = sp_occ_thin, cex = 0.5, col = "red") +
      scale_fill_grass_c(limits = c(0, 1), palette = "viridis") +
      labs(main = "GAM", fill = "Probability") +
      theme_bw()
    plot_gam

    plot_raf = ggplot() +
      geom_spatraster(data = pred_list$raf) +
      geom_spatvector(data = sp_occ_thin, cex = 0.5, col = "red") +
      scale_fill_grass_c(limits = c(0, 1), palette = "viridis") +
      labs(main = "Random forest", fill = "Probability") +
      theme_bw()
    plot_raf
    
  } else if(samp_size >= 15) {
    sdm_data_xy = NULL
    plot_sample = NULL
    sdm_part = NULL
    plot_part = NULL
    m_env = NULL
    m_max = NULL
    m_glm = NULL
    m_gam = NULL
    m_raf = NULL
    pred_env = NULL
    pred_list = NULL
    plot_env = NULL
    plot_max = NULL
    plot_glm = NULL
    plot_gam = NULL
    plot_raf = NULL
  } else {
    sdm_data_xy = NULL
    plot_sample = NULL
    sdm_part = NULL
    plot_part = NULL
    m_env = NULL
    m_max = NULL
    m_glm = NULL
    m_gam = NULL
    m_raf = NULL
    pred_env = NULL
    pred_list = NULL
    plot_env = NULL
    plot_max = NULL
    plot_glm = NULL
    plot_gam = NULL
    plot_raf = NULL
  }
  
  return(list(samp_size = samp_size_df,
              data = sdm_data_xy,
              diagnostics = list(plot_sample, sdm_part, plot_part),
              models = list(m_env, m_max, m_glm, m_gam, m_raf),
              predictions = list(pred_env, pred_list),
              plots = list(plot_env,plot_max, plot_glm, plot_gam, plot_raf))
         )
})
  