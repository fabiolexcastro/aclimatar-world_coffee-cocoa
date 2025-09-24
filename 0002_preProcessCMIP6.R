

## Fabio Castro - Llanos 
## Alliance Bioversity - CIAT 
## August 25th / 2025

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, sf, fs, tidyterra, glue, tidyverse, geodata, ggspatial, RColorBrewer)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
root <- './common_data/esfg_metagrid/raw_world'
gcms <- basename(as.character(dir_ls(root, type = 'directory')))
extn <- c(-0.9375, 359.0625, -25, 25)

# Functions to use --------------------------------------------------------
gcme <- gcms[2]; prdo <- 'ssp370'
make.merge <- function(gcme, prdo){
  
  cat('To process!\n')
  dire <- paste0(root, '/', gcme)
  fles <- dir_ls(dire)
  fles <- as.character(fles)
  
  cat('To make the filter!\n')
  fles <- grep(prdo, fles, value = T)
  vars <- c('hurs','tasmin', 'tasmax', 'pr', 'rsds')
  if(prdo == 'ssp370' & gcme == 'EC-Earth3'){fles <- grep(paste0(2025:2055, collapse = '|'), fles, value = T)}
  
  cat('Map loop by each variable\n')
  msco <- map(.x = 1:length(vars), .f = function(i){
    
    cat('Variable: ', vars[i], '\n')
    rstr <- grep(vars[i], fles, value = T)
  
    trra <- map(.x = 1:length(rstr), .f = function(j){
      
      cat('Layer: ', j, '\n')
      rst <- rstr[[j]]
      rst <- rast(rst)
      
      if(prdo == 'historical'){
        
        year <- 1995:2014
        trr <- rst[[grep(paste0(year, '-', collapse = '|'), time(rst), value = F)]]
        trr <- terra::crop(trr, ext(extn))
        
      } else {
        
        year <- 2025:2055
        trr <- rst[[grep(paste0(year, '-', collapse = '|'), time(rst), value = F)]]
        trr <- terra::crop(trr, ext(extn))
        
      }
      
      return(trr)
      
    })
    
    if(length(trra) > 1){
      
      cat('Files more than 1\n')
      trra <- do.call(c, trra)
      
    } else {
      
      cat('Files just 1\n')
      trra <- trra[[1]]
    
    }
    
    dout <- glue('./common_data/esfg_metagrid/raw_extent/{gcme}')
    dir_create(dout)
    
    if(prdo == 'historical'){
      dtes <- '199550101-20141231'
    } else {
      dtes <- '20250101-20551231'
    }
    
    name <- glue('{vars[i]}_day_{gcme}_{prdo}_r1i1p1f1_{dtes}.tif')
    terra::writeRaster(x = trra, filename = glue('{dout}/{name}'), overwrite = TRUE)
    rm(trra); gc(reset = T)
    cat('File wrotten\n')
    
    
  })
  
  
  
}

# To apply the function ---------------------------------------------------
map(.x = gcms, .f = function(g){
  map(.x = c('historical', 'ssp370'), .f = function(p){
    make.merge(gcme = g, prdo = p)
  })
})

# To check the rasters cropped --------------------------------------------
dirs <- as.character(dir_ls('./common_data/esfg_metagrid/raw_extent',  type = 'directory'))
make.check <- function(dire){
  
  ## To list the files
  cat('To process: ', basename(dire), '\n')
  fles <- as.character(dir_ls(dire, regexp = '.tif$'))
  vars <- c('hurs', 'pr', 'rsds', 'tasmax', 'tasmin')
  
  ## To draw a simple map
  map(.x = vars, .f = function(v){
    
    ### Read as a raster
    cat('To process :', v, '\n')
    fls <- grep(v, fles, value = T)
    bsl <- rast(grep('historical', fls, value = T))
    ftr <- rast(grep('ssp370', fls, value = T))
    
    ### To make the plot
    png(filename = paste0('./png/plot_raster-check/', basename(dire), ' ', v, '.jpg'), units = 'in', width = 5.5, height = 3.5, res = 300)
    par(mfrow = c(2, 1))
    plot(bsl[[1]], main = paste0(basename(dire), ' ', v, ' ', time(bsl)[[1]]))
    plot(ftr[[1]], main = paste0(basename(dire), ' ', v, ' ', time(ftr)[[1]]))
    par(mfrow = c(1, 1))
    dev.off()
    
    ## Finish 
    rm(fls, bsl, ftr)
    gc(reset = T)
    cat('Done!\n')
      
  })
  
  
  
}
map(dirs, make.check)

