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
library(ENMTools) #raster.cor.plot, raster.cor.matrix, trimdupes.by.raster
library(flexsdm)
library(sf)

hgd()
save_path = "/maps/epr26/sdm_captain_out/"

#Read AOI shapefile
aoi = vect(paste0(save_path, "atlantic_forest_global_200.geojson")) %>%
  project("EPSG:4326")
aoi_proj = aoi %>% project("EPSG:3857")
aoi_bbox = as.polygons(ext(aoi_proj), crs = crs(aoi_proj))
writeVector(aoi_bbox, paste0(save_path, "atlantic_forest_global_200_bbox.geojson"), overwrite = T)

#get world map and land boundary
worldmap = geodata::world(path = ".") %>% project("EPSG:4326") #GADM
worldmap_aoi = worldmap %>%
  project("EPSG:3857") %>%
  crop(ext(aoi_bbox))
land = aggregate(worldmap_aoi)
writeVector(land, paste0(save_path, "aoi_land.geojson"), overwrite = T)

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
#biovars[keep_vars]
writeRaster(bioclim_red, paste0(save_path, "rasters/bioclim_reduced.tif"), overwrite = T)


#Occurrence data processing ----
sp_occ_df = readRDS(paste0(save_path, "SpeciesOccurrenceData.rds")) %>%
  as.data.frame() %>%
  filter(complete.cases(ddlat) & complete.cases(ddlon)) %>%
  mutate(x = ddlon, y = ddlat, index = row_number())
sp_occ = sp_occ_df %>%
  vect(geom = c("ddlon", "ddlat"), crs = crs(aoi))
sp_occ_proj = sp_occ %>%
  project("EPSG:3857")
writeVector(sp_occ_proj, paste0(save_path, "SpeciesOccurrenceData.geojson"), overwrite = T)
sp_occ_bbox = crop(sp_occ_proj, ext(aoi_bbox)) #filter species occurrence data by AOI
writeVector(sp_occ_bbox, paste0(save_path, "SpeciesOccurrenceData_bbox.geojson"), overwrite = T)

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

sp_occ_list = vector("list", n_sp)
sp_info = data.frame(index = numeric(), sp_name = character(),
                     original = numeric(), thinned = numeric(), thin_perc = numeric(),
                     data_used = character(), n_used = numeric(), flag = character())

for(i in seq_len(n_sp)) {
  a = Sys.time()
  sp_name = as.character(tax_df$tax[i])
  sp_occ_sel = sp_occ_bbox[sp_occ_bbox$tax == sp_name, ]
  n_orig = nrow(sp_occ_sel)
  
  #optional: geographical distributions of occurrence data and features that may cause spatial biases
  #can be explored using visualization tools in the ‘sampbias' package
  
  #perform spatial-grid thinning for abundant species
  n_thin = NA
  if(n_orig >= 30) {
    thin_method = "trimdupes" #other option: "occfilt" 
    if(thin_method == "trimdupes") {
      sp_occ_thin = ENMTools::trimdupes.by.raster(sp_occ_sel, bioclim) #removes duplicates based on raster cells
      n_thin = nrow(sp_occ_thin)
    } else {
      sp_occ_thin = flexsdm::occfilt_geo(data = crds(sp_occ_sel) %>% as.data.frame(),
                                         x = "x", y = "y",
                                         env_layer = bioclim,
                                         method = c("cellsize", 1),
                                         prj = crs(bioclim)) %>% #
        vect(geom = c("x", "y"), crs = crs(bioclim))
      n_thin = nrow(sp_occ_thin)
    }
  }

  #add attributes back in
  if(n_orig < 30 | n_thin < 30) { #rare species, do not thin
    use = "original"
    n_used = n_orig
    sp_occ_used = sp_occ_sel
  } else {
    use = "thinned"
    n_used = n_thin
    sp_occ_used = sp_occ_sel[geom(sp_occ_sel) %in% geom(sp_occ_thin)]
  }
  
  #flag data point abundance
  samp_size_flag = ifelse(n_used >= 30, "abundant", ifelse(n_used >= 15, "sparse", "insufficient"))
  sp_info[i, ] = data.frame(index = i, sp_name = sp_name,
                            original = n_orig, thinned = n_thin, thin_perc = round((n_thin / n_orig) * 100, 1),
                            data_used = use, n_used = n_used, flag = samp_size_flag)
  sp_occ_list[[i]] = sp_occ_used$index

  b = Sys.time()
  cat(i, "-", sp_name, ":", b - a, "s\n")
}

write.csv(sp_info, paste0(save_path, "species_info.csv"), row.names = F)
saveRDS(sp_occ_list, paste0(save_path, "species_occurrence_thinned.rds"))


#Calculate maximum biomass per tree in a two-step process:
#1. Reverse-estimate maximum diameter from maximum height using generic model 3 in Cysneiros et al. (2020)
#https://cdnsciencepub.com/doi/full/10.1139/cjfr-2020-0060
#log(H) = 1.029 + 0.567 * log(DBH)
#DBH = exp((log(H) - 1.193) / 0.529)

#H = 50.874 * (1 - exp(-0.042 * D ^ 0.784))
#D = (log(1 - H / 50.874) / (-0.042)) ^ (1 / 0.784)
#2. Use the Chave et al. (2005) pantropical model to estimate maximum AGB (kg) per tree for each species
#https://link.springer.com/article/10.1007/s00442-005-0100-x
#AGB = exp(-29.77 + ln(WD * D^2 * H)) ~ 0.0509 * WD * D^2 * H
sp_info = read.csv(paste0(save_path, "species_info.csv"), header = T)
sp_trait = read.csv(paste0(save_path, "SpeciesInfo.csv"), header = T)

sp_info_merged = merge(sp_info,
                       sp_trait[, c("Species", "RedList_international_Category_2023", "GrowthForm", "MaximumHeight_m", "WoodSpecificGravity")],
                       by.x = "sp_name", by.y = "Species", all.x = T) %>%
  mutate(maxH_use_mean = ifelse(is.na(MaximumHeight_m), T, F),
         WSG_use_mean = ifelse(is.na(WoodSpecificGravity), T, F))

#substitute species with NA for maxH/WSG with
maxH_comm_mean = mean(sp_info_merged$MaximumHeight_m, na.rm = T)
WSG_comm_mean = mean(sp_info_merged$WoodSpecificGravity, na.rm = T)
sp_info_merged = sp_info_merged %>%
  mutate(MaximumHeight_m = ifelse(maxH_use_mean, maxH_comm_mean, MaximumHeight_m),
         WoodSpecificGravity = ifelse(WSG_use_mean, WSG_comm_mean, WoodSpecificGravity))

#estimate maximum diameter and AGB
sp_info_merged = sp_info_merged %>%
  mutate(MaximumDiameter_cm = exp((log(MaximumHeight_m) - 1.029) / 0.567)) %>% #Cysneiros et al 2020
  mutate(AGB_kg = exp(-2.977 + log(WoodSpecificGravity * (MaximumDiameter_cm ^ 2) * MaximumHeight_m))) #Chave et al 2005
write.csv(sp_info_merged, paste0(save_path, "species_info.csv"), row.names = F)
