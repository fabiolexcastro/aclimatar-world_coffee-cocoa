

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
indx <- names(indx.hist) %>% str_split('_') %>% map_chr(1)
wrld <- geodata::world(resolution = 1, path = './tmpr')

## Raster data

### Current
indx.hist <- rast('./common_data/atlas_hazards/hist_indices2.tif')

### Future
indx.ftre <- rast('./common_data/atlas_hazards/ftre-gcms_indices2.tif')

## Indices
indices <- c('ndd', 'n30', 'n35', 'hsh', 'tai', 'ndw')

## Tabular data 
load('./rData/run_3_arabica-robusta/presences.rData')
pnts <- clusteredpresdata
pnts <- pnts[,c(1:3)]

pnts <- read_csv('./tbl/points/all_points_coffee-arabica.csv', show_col_types = FALSE)
pnts <- pnts[,c('Lon', 'Lat')]

### Remove duplicated by cell
cell <- terra::extract(indx.hist[[1]], pnts[,1:2], cell = T)[,3]
pnts <- mutate(pnts, cll = cell)
cell <- unique(cell)
pnts.rmv <- map(.x = cell, .f = function(i){slice(filter(pnts, cll == i), 1)})
pnts.rmv <- bind_rows(pnts.rmv)
pnts <- pnts.rmv

### ISOS
isos <- terra::extract(wrld, pnts[,1:2])
isos <- isos[,2]
pnts <- mutate(pnts, sov_a3 = isos)
isos <- table(isos) %>% as.data.frame() %>% arrange(desc(Freq)) %>% filter(Freq > 10) %>% pull(1)

# Function to classify ----------------------------------------------------
stck <- indx.ftre.avg
iso <- 'COL'
indx <- 'ndd'

clsf <- function(stck, iso, indx){
  cat('To process: ', indx, '\n')
  pnt <- filter(pnts, sov_a3 == iso)
  shp <- wrld[wrld$GID_0 == iso,]
  rst <- stck[[grep(indx, names(stck))]]
  rst <- terra::crop(rst, shp) %>% mask(shp)
  vls <- terra::extract(rst, pnt[,1:2])
  qnt <- quantile(vls[,2], c(0, 0.5, 0.75, 1), na.rm = T)
  mtr <-matrix(c(0, qnt[2], 1, qnt[2], qnt[3], 2, qnt[3], Inf, 3), byrow = T, ncol = 3)
  cls <- terra::classify(rst, mtr, include.lowest = T)
  return(cls)
}

# Future ------------------------------------------------------------------
indx.ftre.avg <- map(.x = 1:length(indices), .f = function(i){
  rst <- indx.ftre[[grep(indices[i], names(indx.ftre))]]
  rst <- mean(rst)
  names(rst) <- indices[i]
  return(rst)
})
indx.ftre.avg <- reduce(indx.ftre.avg, c)

# To classify -------------------------------------------------------------

## By each country 

### Baseline
indx.hist.cls <- map(1:length(indices), function(i){
  cat('Index: ', indices[i], '\n')
  rr <- map(.x = unique(isos), .f = function(s){
    try(
      expr = {
        cat(s, '\n')
        r <- clsf(stck = indx.hist, iso = s, indx = indices[i])          
      }
    )
  })
  rr <- rr[sapply(rr, inherits, what = "SpatRaster")]
  rr <- sprc(rr)
  rr <- mosaic(rr, fun = 'modal')
  return(rr)
})
indx.hist.cls  <- reduce(indx.hist.cls, c)

### Future
inds <- c('ndd', 'n30', 'n35', 'tai', 'ndw', 'hsh')
indx.ftre.cls <- map(1:length(inds), function(i){
  
  cat('Index: ', inds[i], '\n')
  rr <- map(.x = as.character(isos), .f = function(s){
    try(
      expr = {
        cat(s, '\n')
        r <- clsf(stck = indx.ftre.avg, iso = s, indx = inds[i])          
      }
    )
  })
  rr <- rr[sapply(rr, inherits, what = "SpatRaster")]
  rr <- sprc(rr)
  rr <- mosaic(rr, fun = 'modal')
  return(rr)
  
})
indx.ftre.cls <- reduce(indx.ftre.cls, c)
names(indx.ftre.cls) <- paste0(names(indx.ftre.cls), '_ftre')

# To write the rasters ----------------------------------------------------
plot(c(indx.hist.cls[[1]], indx.ftre.cls[[1]]))
plot(c(indx.hist.cls[[2]], indx.ftre.cls[[2]]))

terra::writeRaster(x = indx.hist.cls, filename = './common_data/stack_final/stack_indices-hist.tif', overwrite = TRUE)
terra::writeRaster(x = indx.ftre.cls, filename = './common_data/stack_final/stack_indices-ftre.tif', overwrite = TRUE)





