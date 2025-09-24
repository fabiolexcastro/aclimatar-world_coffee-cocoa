

## Fabio Castro - Llanos 
## Alliance Bioversity - CIAT 
## August 28th / 2025

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, sf, fs, glue, tidyverse, geodata, ggspatial, RColorBrewer)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
dirs <- dir_ls('./common_data/esfg_metagrid/raw_extent') %>% as.character()
wrld <- geodata::world(resolution = 1, path = './tmpr', level = 0)

# Function ----------------------------------------------------------------
rott.conv <- function(dire){
  
  # dire <- dirs[1]
  
  ## To list the files
  cat('To process: ', basename(dire), '\n')
  fles <- dir_ls(dire, regexp = '.tif$')
  fles <- as.character(fles)
  gcme <- basename(dire)
  
  map(.x = fles, .f = function(f){
  
    ## To get the name and the variable
    cat('>>>> ', basename(f), '\n')
    name <- basename(f)
    varb <- str_split(name, '_') %>% map_chr(1)
    
    ## To read the raster
    rstr <- rast(f)
    dtes <- time(rstr)
    
    ## To convert
    if(varb == 'pr'){
      rstr <- rstr * 86400
    } else if(varb %in% c('tasmin', 'tasmax')){
      rstr <- rstr - 273.15
    } else if(varb == 'rsds'){
      rstr <- rstr * 86400 / 1000000
    }
    
    ## To rotate 
    rstr <- terra::rotate(rstr)
    rstr <- terra::crop(rstr, ext(-115, 179, -25.625, 25.625))
    
    ## To write the raster
    dout <- glue('./common_data/esfg_metagrid/raw_extent-rotate/{gcme}')
    dir_create(dout)
    terra::writeRaster(x = rstr, filename = glue('{dout}/{name}'), overwrite = TRUE)
    
    ## Finish
    rm(rstr, dout)
    gc(reset = T)
    cat('Done!\n')
      
  })
  
}

# Apply function ----------------------------------------------------------
map(dirs, rott.conv)




