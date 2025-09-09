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

#Read AOI shapefile
aoi = vect("atlantic_forest_global_200.geojson") %>%
  project("EPSG:4326")
aoi_proj = aoi %>% project("EPSG:3857")
aoi_bbox = as.polygons(ext(aoi_proj), crs = crs(aoi_proj))
writeVector(aoi_bbox, "atlantic_forest_global_200_bbox.geojson", overwrite = T)

#get world map and land boundary
worldmap = geodata::world(path = ".") %>% project("EPSG:4326") #GADM
worldmap_aoi = worldmap %>%
  project("EPSG:3857") %>%
  crop(ext(aoi_bbox))
land = aggregate(worldmap_aoi)
writeVector(land, "aoi_land.geojson", overwrite = T)

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
bioclim_paths = list.files(path = "wc2.1_5m_bio/", pattern = "tif", full.names = T)
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
cor_plot = ENMTools::raster.cor.plot(bioclim) #keep 1, 2, 7, 12, 15, 18, 19
keep_vars = c(1, 2, 7, 12, 15, 18, 19)
bioclim_red = bioclim[[keep_vars]]
ENMTools::raster.cor.plot(bioclim_red)
mat_cor = ENMTools::raster.cor.matrix(bioclim_red)
diag(mat_cor) = NA
biovars[keep_vars]
writeRaster(bioclim_red,"bioclim_reduced.tif", overwrite = T)


#Occurrence data processing ----
sp_occ_df = readRDS("SpeciesOccurrenceData.rds") %>%
  as.data.frame() %>%
  filter(complete.cases(ddlat) & complete.cases(ddlon)) %>%
  mutate(x = ddlon, y = ddlat, index = row_number())
sp_occ = sp_occ_df %>%
  vect(geom = c("ddlon", "ddlat"), crs = crs(aoi))
sp_occ_proj = sp_occ %>%
  project("EPSG:3857")
writeVector(sp_occ_proj, "SpeciesOccurrenceData.geojson", overwrite = T)
sp_occ_bbox = crop(sp_occ_proj, ext(aoi_bbox)) #filter species occurrence data by AOI
writeVector(sp_occ_bbox, "SpeciesOccurrenceData_bbox.geojson", overwrite = T)

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


# Perform thinning and examine sample size ----
tax_df = as.data.frame(table(sp_occ_bbox$tax)) %>%
  rename(tax = Var1, count = Freq)
n_sp = nrow(tax_df)

#samp_size_df = data.frame(original = rep(NA, n_sp), trimdupes = rep(NA, n_sp), occfilt = rep(NA, n_sp), flag = rep("", n_sp))
#sp_occ_thin_list = vector("list", n_sp)

sp_occ_list = vector("list", n_sp)
samp_size_df = data.frame(index = numeric(), sp_name = character(),
                          original = numeric(), trimdupes = numeric(), occfilt = numeric(), flag = character())

for(i in seq_len(n_sp)) {
  a = Sys.time()
  sp_name = as.character(tax_df$tax[i])
  sp_occ_sel = sp_occ_bbox[sp_occ_bbox$tax == sp_name, ]
  n_orig = nrow(sp_occ_sel)
  
  #geographical distributions of occurrence data and features that may cause spatial biases
  #can be explored using visualization tools in the ‘sampbias' package
  
  #spatial-grid thinning based on one point per grid cell
  sp_occ_trimdupes = ENMTools::trimdupes.by.raster(sp_occ_sel, bioclim)
  n_trimdupes = nrow(sp_occ_trimdupes)

  #spatial-grid thinning based on nearest-neighbor distance larger than grid cell size
  sp_occ_occfilt = flexsdm::occfilt_geo(data = crds(sp_occ_sel) %>% as.data.frame(),
                                        x = "x", y = "y",
                                        env_layer = bioclim,
                                        method = c("cellsize", 1),
                                        prj = crs(bioclim)) %>%
    vect(geom = c("x", "y"), crs = crs(bioclim))
  n_occfilt = nrow(sp_occ_occfilt)
  sp_occ_occfilt_attr = sp_occ_sel[geom(sp_occ_sel) %in% geom(sp_occ_occfilt)]
  
  #flag data point abundance
  samp_size_flag = ifelse(n_occfilt >= 30, "abundant", ifelse(n_occfilt >= 15, "rare", "insufficient"))
  samp_size_df[i, ] = data.frame(index = i, sp_name = sp_name,
                                 original = n_orig, trimdupes = n_trimdupes, occfilt = n_occfilt, flag = samp_size_flag)
  sp_occ_list[[i]] = sp_occ_occfilt_attr$index

  b = Sys.time()
  cat(i, "-", sp_name, ":", b - a, "s\n")
}

samp_size_df$trimdupes_perc = round(samp_size_df$trimdupes / samp_size_df$original * 100, 1)
samp_size_df$occfilt_perc = round(samp_size_df$occfilt / samp_size_df$original * 100, 1)

write.csv(samp_size_df, "species_sample_size.csv", row.names = F)
saveRDS(sp_occ_list, "species_occurrence_thinned.rds")