rm(list = ls())

#Setup ----
library(httpgd) #view plots in VS Code
library(future) #parallelise lapply() : future_lapply()
library(future.apply) #parallelise lapply(): future_lapply()
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)
library(geodata) #world, worldclim_global
library(fuzzyjoin)
library(devtools)
#library(remotes)
#remotes::install_github("danlwarren/ENMTools")
library(ENMTools) #raster.cor.plot, raster.cor.matrix, trimdupes.by.raster
#remotes::install_github("sjevelazco/flexsdm")
library(flexsdm)
library(sf)

hgd()
dir_path = "/maps/epr26/captain_brazil/"

#To define: AOI path, resolution, output path
proj_aoi = "atlantic_forest_global_200.geojson"
proj_aoi_sub = "200k_ha_corridor_Ideal-area.geojson" #"atlantic_forest_global_200.geojson"
wc_res = 0.5 #10, 5, 2.5, and 0.5: 5 for af_10km and 0.5 for ideal_250m
proj_res = 250 #in meters
proj_folder = "ideal_250m" #af_10km
proj_path = paste0(dir_path, proj_folder, "/")
if (!dir.exists(proj_path)) dir.create(proj_path)
if (!dir.exists(paste0(proj_path, "rasters/"))) dir.create(paste0(proj_path, "rasters/"))
test = F #draw diagnostic plots if True

#Read AOI shapefile
aoi = vect(paste0(dir_path, proj_aoi)) %>%
  project("EPSG:4326")
aoi_proj = aoi %>% project("EPSG:3857")
aoi_bbox = as.polygons(ext(aoi_proj), crs = crs(aoi_proj))
writeVector(aoi, paste0(dir_path, "aoi.geojson"), overwrite = T)
writeVector(aoi_proj, paste0(dir_path, "aoi_proj.geojson"), overwrite = T)
writeVector(aoi_bbox, paste0(dir_path, "aoi_bbox.geojson"), overwrite = T)

#Read subsetted AOI shapefile
aoi_sub = vect(paste0(dir_path, proj_aoi_sub)) %>%
  project("EPSG:4326")
aoi_sub_proj = aoi_sub %>% project("EPSG:3857")
aoi_sub_bbox = as.polygons(ext(aoi_sub_proj), crs = crs(aoi_proj))
writeVector(aoi_sub, paste0(proj_path, "aoi.geojson"), overwrite = T)
writeVector(aoi_sub_proj, paste0(proj_path, "aoi_proj.geojson"), overwrite = T)
writeVector(aoi_sub_bbox, paste0(proj_path, "aoi_bbox.geojson"), overwrite = T)


#get world map and land boundary
worldmap = geodata::world(path = dir_path) %>% project("EPSG:4326") #GADM
if(!file.exists(paste0(dir_path, "aoi_land.geojson"))) {
  worldmap_aoi = worldmap %>%
    project(crs(aoi_proj)) %>%
    crop(ext(aoi_bbox))
  land = aggregate(worldmap_aoi)
  writeVector(land, paste0(dir_path, "aoi_land.geojson"), overwrite = T)
} else {
  land = vect(paste0(dir_path, "aoi_land.geojson"))
}

if(!file.exists(paste0(proj_path, "aoi_land.geojson"))) {
  worldmap_aoi_sub = worldmap %>%
    project(crs(aoi_sub_proj)) %>%
    crop(ext(aoi_sub_bbox))
  land_sub = aggregate(worldmap_aoi_sub)
  writeVector(land_sub, paste0(proj_path, "aoi_land.geojson"), overwrite = T)
} else {
  land_sub = vect(paste0(proj_path, "aoi_land.geojson"))
}


#bioclimatic variable names
biovars = c("Annual Mean Temperature",
            "Mean Diurnal Range",
            "Isothermality",
            "Temperature Seasonality ",
            "Max Temperature of Warmest Month",
            "Min Temperature of Coldest Month",
            "Temperature Annual Range",
            "Mean Temperature of Wettest Quarter",
            "Mean Temperature of Driest Quarter",
            "Mean Temperature of Warmest Quarter",
            "Mean Temperature of Coldest Quarter",
            "Annual Precipitation",
            "Precipitation of Wettest Month",
            "Precipitation of Driest Month",
            "Precipitation Seasonality",
            "Precipitation of Wettest Quarter",
            "Precipitation of Driest Quarter",
            "Precipitation of Warmest Quarter",
            "Precipitation of Coldest Quarter")


#Environmental data processing ----
res_text = switch(as.character(wc_res),
                  "10" = "10m",
                  "5" = "5m",
                  "2.5" = "2.5m",
                  "0.5" = "30s")
bioclim_orig = rast(paste0(dir_path, "wc_", res_text, "_sa.tif")) #res = 5 for the entire AF
var_order = names(bioclim_orig) %>% sub(".*bio_", "", .) %>% as.numeric() %>% order()
bioclim_orig = bioclim_orig[[var_order]]

#crop to AOI and reproject to defined resolution
bioclim_proj = bioclim_orig %>%
  crop(ext(aoi)) %>% #crop exactly to the bounding box of the AOI
  project(crs(aoi_proj), res = proj_res) #reproject to defined resolution
writeRaster(bioclim_proj, paste0(proj_path, "rasters/bioclim_all.tif"), overwrite = T)

bioclim_proj = rast(paste0(proj_path, "rasters/bioclim_all.tif"))


#examine collinearity
#removing redundant variables (pairwise: ‘ENMTML', ‘flexsdm', ‘modleR', ‘ntbox';
#sequential: ‘fuzzySim', ‘SDMtune', ‘usdm')
#or reducing variable dimensionality through ordination (‘ENMTML', ‘ENMTools', ‘flexsdm', ‘kuenm', ‘ntbox')
#flexsdm::correct_colinvar but there is an error
ENMTools::raster.cor.plot(bioclim_proj) #visualise: keep 1, 2, 7, 12, 15, 18, 19
keep_vars = c(1, 2, 7, 12, 15, 18, 19)
bioclim = bioclim_proj[[keep_vars]]
# mat_cor = ENMTools::raster.cor.matrix(bioclim)
# diag(mat_cor) = NA
writeRaster(bioclim, paste0(proj_path, "rasters/bioclim.tif"), overwrite = T)

if(test) {
  #plot correlation matrices
  bioclim_all_named = bioclim_all
  names(bioclim_all_named) = biovars
  cor_plot = ENMTools::raster.cor.plot(bioclim_all_named)$cor.heatmap +
    labs(x = NULL, y = NULL)
  ggsave(paste0(proj_path, "plot_bioclim_correlation_all.png"), width = 8, height = 6, units = "in", dpi = 300)

  bioclim_named = bioclim
  names(bioclim_named) = biovars[keep_vars]
  cor_plot_red = ENMTools::raster.cor.plot(bioclim_named)$cor.heatmap +
    labs(x = NULL, y = NULL)
  ggsave(paste0(proj_path, "plot_bioclim_correlation_red.png"), width = 8, height = 6, units = "in", dpi = 300)
}


#Occurrence data processing ----
if(!file.exists(paste0(dir_path, "SpeciesOccurrenceData.geojson"))) {
  sp_occ_df = readRDS(paste0(dir_path, "SpeciesOccurrenceData.rds")) %>%
    as.data.frame() %>%
    filter(complete.cases(ddlat) & complete.cases(ddlon)) %>%
    mutate(x = ddlon, y = ddlat, index = row_number())
  sp_occ = sp_occ_df %>%
    vect(geom = c("ddlon", "ddlat"), crs = crs(aoi))
  sp_occ_proj = sp_occ %>%
    project(crs(aoi_proj))
  writeVector(sp_occ_proj, paste0(dir_path, "SpeciesOccurrenceData.geojson"), overwrite = T)
} else {
  sp_occ_proj = vect(paste0(dir_path, "SpeciesOccurrenceData.geojson"))
}
if(!file.exists(paste0(dir_path, "spocc_bbox.geojson"))) {
  sp_occ_bbox = crop(sp_occ_proj, ext(aoi_bbox)) #filter species occurrence data by AOI
  writeVector(sp_occ_bbox, paste0(dir_path, "spocc_bbox.geojson"), overwrite = T)
} else {
  sp_occ_bbox = vect(paste0(dir_path, "spocc_bbox.geojson"))
}

#visualize and examine anomalous coordinates: not really needed
if(test) {
  ggplot() +
  geom_spatvector(data = worldmap, fill = "lightyellow") +
  geom_spatvector(data = sp_occ, color = "orange", size = 0.05) +
  coord_sf(xlim = c(-120, -5), ylim = c(-50, 40)) +
  theme_bw()

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
}


# Calculate maximum biomass per tree ----

#Gather information about all species
if(!file.exists(paste0(dir_path, "species_info_all.csv"))) {
  tax_df = as.data.frame(table(sp_occ_bbox$tax)) %>%
    rename(tax = Var1, count = Freq)
  n_sp = nrow(tax_df)
  sp_info_all = tax_df %>%
    mutate(sp_ind = seq_len(n_sp))
  write.csv(sp_info_all, paste0(dir_path, "species_info_all.csv"), row.names = F)
  write.table(sp_info_all$sp_ind, paste0(dir_path, "species_retained.txt"), sep = "\n", row.names = F, col.names = F)
} else {
  sp_info_all = read.csv(paste0(dir_path, "species_info_all.csv"), header = T)
  n_sp = nrow(sp_info_all)
}

#1. Reverse-estimate maximum diameter from maximum height using generic model 3 in Cysneiros et al. (2020)
#https://cdnsciencepub.com/doi/full/10.1139/cjfr-2020-0060
#log(H) = 1.029 + 0.567 * log(DBH)
#DBH = exp((log(H) - 1.193) / 0.529)

#H = 50.874 * (1 - exp(-0.042 * D ^ 0.784))
#D = (log(1 - H / 50.874) / (-0.042)) ^ (1 / 0.784)
#2. Use the Chave et al. (2014) improved pantropical model to estimate maximum AGB (kg) per tree for each species
#https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.12629
#AGB = 0.0673 * (WD * D^2 * H)^0.976

sp_trait = read.csv(paste0(dir_path, "SpeciesInfo.csv"), header = T) %>%
  dplyr::select(Species, RedList_international_Category_2023, GrowthForm, MaximumHeight_m, WoodSpecificGravity)

#calculate community mean
maxH_comm_mean = mean(sp_trait$MaximumHeight_m, na.rm = T)
WSG_comm_mean = mean(sp_trait$WoodSpecificGravity, na.rm = T)

#substitute NAs with community mean for maxH/WSG
sp_info_merged = merge(sp_info_all, sp_trait, by.x = "tax", by.y = "Species", all.x = T) %>%
  mutate(maxH_use_mean = ifelse(is.na(MaximumHeight_m), T, F),
         WSG_use_mean = ifelse(is.na(WoodSpecificGravity), T, F)) %>%
  mutate(MaximumHeight_m = ifelse(maxH_use_mean, maxH_comm_mean, MaximumHeight_m),
         WoodSpecificGravity = ifelse(WSG_use_mean, WSG_comm_mean, WoodSpecificGravity))

#estimate maximum diameter and AGB
sp_info_merged = sp_info_merged %>%
  mutate(MaximumDiameter_cm = exp((log(MaximumHeight_m) - 1.029) / 0.567)) %>% #Cysneiros et al 2020
  mutate(AGB_kg = 0.0673 * (WoodSpecificGravity * MaximumDiameter_cm ^ 2 * MaximumHeight_m) ^ 0.976) #Chave et al 2014
write.csv(sp_info_merged, paste0(dir_path, "species_info_merged.csv"), row.names = F)


#obtain species occurrence data in subsetted AOI: used to identify species for SDM
sp_occ_sub_bbox = crop(sp_occ_proj, ext(aoi_sub_bbox))
writeVector(sp_occ_sub_bbox, paste0(proj_path, "spocc_bbox.geojson"), overwrite = T)
tax_df_sub = as.data.frame(table(sp_occ_sub_bbox$tax)) %>%
  rename(tax = Var1, count = Freq)
n_sp_sub = nrow(tax_df_sub)
sp_info_sub = sp_info_all %>%
  filter(tax %in% tax_df_sub$tax)
write.table(sp_info_sub$ind, paste0(proj_path, "species_retained.txt"), sep = "\n", row.names = F, col.names = F)


#perform spatial-grid thinning for abundant species
sp_occ_list = vector("list", n_sp_sub)
sp_thinning = data.frame(sp_ind = numeric(), sp_name = character(),
                         original = numeric(), thinned = numeric(), thin_perc = numeric(),
                         data_used = character(), n_used = numeric(), flag = character())

for(i in seq_len(n_sp_sub)) {
  a = Sys.time()
  sp_ind = sp_info_sub$sp_ind[i]
  sp_name = as.character(sp_info_sub$tax[i])
  sp_occ_sel = sp_occ_bbox[sp_occ_bbox$tax == sp_name, ]
  n_orig = nrow(sp_occ_sel)

  n_thin = NA
  if(n_orig < 30) { #rare species, do not perform thinning
    use = "original"
    n_used = n_orig
    sp_occ_used = sp_occ_sel
  } else { #perform thinning
    sp_occ_thin = ENMTools::trimdupes.by.raster(sp_occ_sel, bioclim) #removes duplicates based on raster cells
    n_thin = nrow(sp_occ_thin)
    if(n_thin >= 30) { #thinned results acceptable
      use = "thinned"
      n_used = n_thin
      sp_occ_used = sp_occ_sel[geom(sp_occ_sel) %in% geom(sp_occ_thin)]
    } else { #thinned results too sparse, revert to original data
      use = "original"
      n_used = n_orig
      sp_occ_used = sp_occ_sel
    }
  }

  #flag data point abundance
  samp_size_flag = ifelse(n_used >= 30, "abundant", ifelse(n_used >= 3, "sparse", "insufficient"))
  sp_thinning[i, ] = data.frame(sp_ind = sp_ind, sp_name = sp_name,
                                original = n_orig, thinned = n_thin, thin_perc = round((n_thin / n_orig) * 100, 1),
                                data_used = use, n_used = n_used, flag = samp_size_flag)
  sp_occ_list[[i]] = sp_occ_used$index

  b = Sys.time()
  cat("Processed species ", sp_ind, " (", i, "/", n_sp_sub, "): ", round(as.numeric(b - a, units = "secs"), 2), " secs)\n", sep = "")
}

sp_info = sp_info_sub %>%
  left_join(sp_thinning, by = c("sp_ind", "tax" = "sp_name"))

write.csv(sp_info, paste0(proj_path, "species_info.csv"), row.names = F)
saveRDS(sp_occ_list, paste0(proj_path, "species_occurrence_thinned.rds"))