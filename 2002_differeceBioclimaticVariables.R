
## Fabio Castro - Llanos 
## Alliance Bioversity - CIAT 
## Difference Bioclimatic variables 

## Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, glue, dismo, geodata, gtools)
pacman::p_load(terra, fs, usdm, raster, gtools, spocc, randomForest, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth, geodata)

g <- gc(reset = T); rm(list = ls())
options(scipen = 999, warn = -1)

# Data --------------------------------------------------------------------
fles <- as.character(dir_ls('./common_data/bioclimatic_variables', regexp = '.tif$'))

## Baseline
bsln <- grep('hist', fles, value = T)
bsln.1 <- rast(grep('bios_hist', bsln, value = T))
bsln.2 <- rast(grep(paste0(c('100mm', '60mm'), collapse = '|'), fles, value = T)) %>% setNames(c('bioc_20_100', 'bioc_20_60'))
bsln.3 <- rast(grep('bios-etp', bsln, value = T))
bsln.4 <- rast(grep('balance_hist', bsln, value = T))
bsln.r <- c(bsln.1, bsln.2, bsln.3, bsln.4)

### To write the raster 
outd <- './common_data/input_bios'
terra::writeRaster(x = bsln.r, filename = glue('{outd}/bioc_hist.tif'), overwrite = TRUE)
bsln.r <- terra::rast('./common_data/input_bios/bioc_hist.tif')

## Future
ftre <- grep('all_', fles, value = T)
ftre <- map(ftre, rast)
map(ftre, ext)
ftre.r <- map(1:length(ftre), .f = function(i){
  print(i)
  ftr <- terra::resample(ftre[[i]], bsln.r, method = 'bilinear')
  return(ftr)
})
ftre.r <- reduce(ftre.r, c)
terra::writeRaster(x = ftre.r, filename = './common_data/input_bios/bioc_ftre-gcms.tif', overwrite = TRUE)
ftre.r <- terra::rast('./common_data/input_bios/bioc_ftre-gcms.tif')
gcms   <- names(ftre.r) %>% str_split('_') %>% map_chr(3) %>% unique()

# To calculate the difference ---------------------------------------------
calc.dfrn <- function(gcme){
  
  gcme <- gcms[1]
    
  ## Filtering the rasters
  cat('GCM: ', gcme, '\n')
  ftr <- ftre.r[[grep(gcme, names(ftre.r))]]
  bsl <- bsln.r
  
  ## To change the names
  names(ftr) <- c(glue('bioc_{1:19}'), glue('bioc_2060'), glue('bioc_20100'), glue('bioc_{21:33}'))
  names(bsl) <- c(glue('bioc_{1:19}'), glue('bioc_20100'), glue('bioc_2060'), glue('bioc_{21:33}'))
  
  ## To calculate the difference
  
  ### Temperature
  dfr.tas <- map(.x = 1:11, .f = function(i){
    cat('To process: ', i, '\n')
    dfr <- ftr[[i]] - bsl[[i]]
    names(dfr) <- glue('dfr_{names(dfr)}')
    return(dfr)
  })
  dfr.tas <- reduce(dfr.tas, c)
  
  ### Precipitation
  dfr.ppt <- map(.x = 12:19, .f = function(i){
    
    ## To calculate the difference
    cat('To process: ', i, '\t')
    dfr <- ftr[[i]] - bsl[[i]]
    bsl <- bsl[[i]]
    bsl[bsl == 0] <- 0.1
    dfr <- (dfr / bsl) * 100
    names(dfr) <- glue('dfr_{names(dfr)}')
    
    ##
    cat('To estimate the threshold\t')
    thr <- as.numeric(quantile(dfr[], 0.99, na.rm = T))
    dfr[dfr > thr] <- NA
    
    ##
    cat('Finish!\n')
    return(dfr)
    
  })
  dfr.ppt <- reduce(dfr.ppt, c)
  
  ### Bio 20 
  dfr.b20 <- map(.x = 20:21, .f = function(i){
    
    ## To calculate the difference
    cat('To process ', i, '\t')
    dfr <- ftr[[i]] - bsl[[i]]
   
    ##
    cat('Finish!\n')
    return(dfr)
    
  })
  dfr.b20 <- reduce(dfr.b20, c)
  
  ### Potential Evapotranspiration
  dfr.etp <- map(.x = 22:30, .f = function(i){
    
    ## To calculate the difference
    cat('To process: ', i, '\t')
    dfr <- ftr[[i]] - bsl[[i]]
    bsl <- bsl[[i]]
    bsl[bsl == 0] <- 0.1
    dfr <- (dfr / bsl) * 100
    names(dfr) <- glue('dfr_{names(dfr)}')
    
    ##
    cat('To estimate the threshold\t')
    thr <- as.numeric(quantile(dfr[], 0.99, na.rm = T))
    dfr[dfr > thr] <- NA
    
    ## Finish 
    cat('Finish!\n')
    return(dfr)
    
  })
  dfr.etp <- reduce(dfr.etp, c)

  ### Balance
  dfr.bln <- map(.x = 31:34, .f = function(i){
    
    ## To calculate the difference
    cat('To process: ', i, '\t')
    dfr <- ftr[[i]] - bsl[[i]]
    bsl <- bsl[[i]]
    bsl[bsl == 0] <- 0.1
    dfr <- (dfr / bsl) * 100
    names(dfr) <- glue('dfr_{names(dfr)}')
    
    ##
    cat('To estimate the threshold\t')
    thr <- as.numeric(quantile(dfr[], 0.99, na.rm = T))
    dfr[dfr > thr] <- NA
    
    ## Finish 
    cat('Finish!\n')
    return(dfr)
    
  })
  dfr.bln <- reduce(dfr.bln, c)
  
  ## Finish
  dfr.stk <- c(dfr.tas, dfr.ppt, dfr.b20, dfr.etp, dfr.bln)
  names(dfr.stk) <- glue('{names(dfr.stk)}_{gcme}')
  cat('Done!\n')
  return(dfr.stk)
  
}





