

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, glue)

g <- gc(reset = T)
rm(list = ls())


# Load data ---------------------------------------------------------------

fles <- fs::dir_ls('./rf/output/run_4/results/process')
fles <- as.character(fles)
fles <- grep('.tif$', fles, value = T)

## Raster baseline 
bsln <- terra::rast(grep('mixed_bsl', fles, value = T))
plot(bsln)

## Raster future
ftre <- terra::rast(grep('mix_ftr-modal', fles, value = T))
gcms <- terra::rast(grep('mix_ftr-gcms', fles, value = T))
