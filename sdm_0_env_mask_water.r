library(terra)

dir = "/cfs/klemming/projects/snic/brazilaf_sdm/data_in"

env_60k_tmp = rast(paste0(dir, "/env_60k_tmp.tif"))
env_200k_tmp = rast(paste0(dir, "/env_200k_tmp.tif"))
water_mask_60k = rast(paste0(dir, "/water_mask_60k.tif"))
water_mask_200k = rast(paste0(dir, "/water_mask_200k.tif"))

water_mask_60k_exact = resample(water_mask_60k, env_60k_tmp, method = "near")
water_mask_200k_exact = resample(water_mask_200k, env_200k_tmp, method = "near")

env_60k = mask(env_60k_tmp, water_mask_60k_exact)
env_200k = mask(env_200k_tmp, water_mask_200k_exact)

writeRaster(env_60k, paste0(dir, "/env_60k.tif"), overwrite = T,
            gdal = c("TILED=YES", "BIGTIFF=YES", "COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=4",
            "BLOCKXSIZE=512", "BLOCKYSIZE=512", "ZLEVEL=9"))
writeRaster(env_200k, paste0(dir, "/env_200k.tif"), overwrite = T,
            gdal = c("TILED=YES", "BIGTIFF=YES", "COMPRESS=DEFLATE", "PREDICTOR=3", "NUM_THREADS=4",
            "BLOCKXSIZE=512", "BLOCKYSIZE=512", "ZLEVEL=9"))