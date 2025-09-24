

## Fabio Castro - Llanos 
## Alliance Bioversity - CIAT 
## Postprocessing Indices

## Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, glue, dismo, geodata, gtools)
pacman::p_load(terra, fs, usdm, raster, gtools, spocc, randomForest, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth, geodata)

g <- gc(reset = T); rm(list = ls())
options(scipen = 999, warn = -1)

# Historic datasets -------------------------------------------------------

# Parameters -------------
dirs <- dir_ls('./common_data/atlas_hazards/historical')
dirs <- as.character(dirs)

# Function ---------------
calc.avrg <- function(dire, year){
  
  # dire <- dirs[5]
  # year <- 1995:2014
  
  ## To list the files 
  cat('To process: ', dire, '\n')
  fles <- as.character(dir_ls(dire))
  fles <- grep('.tif$', fles, value = T)
  indx <- basename(dire)
  
  ## Descriptive by each year
  trra <- map(.x = year, .f = function(y){
    
    ## To read as a raster
    cat('>>> Year: ', y, '\n')
    rstr <- grep(y, fles, value = T)
    rstr <- rast(rstr)
    
    ## To estimate the index average or sum
    if(indx %in% c('NDD', 'NTx30', 'NTx35')){
      cat('Index: ', indx, '\n')
      rstr <- sum(rstr)  
    } else if(indx %in% c('TAI', 'HSH')){
      cat('Index: ', indx, '\n')
      rstr <- mean(rstr)
    }
    
    ## Finish!
    cat('Done!\n')
    return(rstr)
    
  })
  
  ## Finish!\n
  trra <- reduce(trra, c)
  names(trra) <- glue('{indx}_{year}')
  cat('Finish the process!\n')
  return(trra)
  
}

# Index ------------------

## NDD -------------------
ndd.hist <- calc.avrg(dire = grep('NDD', dirs, value = T), year = 1995:2014)

## NTx30 -----------------
n30.hist <- calc.avrg(dire = grep('30', dirs, value = T), year = 1995:2014)
n35.hist <- calc.avrg(dire = grep('35', dirs, value = T), year = 1995:2014)

## HSH -------------------
hsh.hist <- calc.avrg(dire = grep('HSH', dirs, value = T), year = 1995:2014)

## TAI -------------------
tai.hist <- calc.avrg(dire = grep('TAI', dirs, value = T), year = 1995:2014)

## NDWL0 -----------------
ndw.hist <- calc.avrg(dire = grep('NDWL0', dirs, value = T), year = 1995:2014)
