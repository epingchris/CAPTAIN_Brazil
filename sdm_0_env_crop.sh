#!/bin/bash
#SBATCH -A naiss2026-4-88
#SBATCH -N 1
#SBATCH -p shared
#SBATCH --ntasks=2
#SBATCH --mem=32G
#SBATCH -t 72:00:00
#SBATCH -J sdm_crop
#SBATCH --output=/cfs/klemming/projects/snic/brazilaf_sdm/logs/sdm_0_env_crop.log
#SBATCH --error=/cfs/klemming/projects/snic/brazilaf_sdm/logs/sdm_0_env_crop.err

# Load required modules
module load PDC/24.11
module load R/4.4.2-cpeGNU-24.11
module load gdal/3.10.0-cpeGNU-24.11
module load nano/7.2

set -euo pipefail #exit on error, undefined variable, or failed pipe

DIR="/cfs/klemming/projects/snic/brazilaf_sdm/data_in"
GTIFF_OPTS=(
  -co TILED=YES
  -co BIGTIFF=YES
  -co COMPRESS=DEFLATE
  -co PREDICTOR=3
  -co NUM_THREADS=4
  -co BLOCKXSIZE=512
  -co BLOCKYSIZE=512
  -co ZLEVEL=9
)

#crop to 200k and 60k AOI (for prediction)
gdalwarp \
  -cutline "${DIR}/200k_ha_corridor_ideal_3857.geojson" \
  -crop_to_cutline \
  -dstnodata -9999 \
  -wm 8000 \
  -overwrite \
  "${DIR}/env_land.tif" "${DIR}/env_200k_tmp.tif" \
  "${GTIFF_OPTS[@]}"

gdalwarp \
  -cutline "${DIR}/60k_ha_corridor_priority_3857.geojson" \
  -crop_to_cutline \
  -dstnodata -9999 \
  -wm 8000 \
  -overwrite \
  "${DIR}/env_land.tif" "${DIR}/env_60k_tmp.tif" \
  "${GTIFF_OPTS[@]}"