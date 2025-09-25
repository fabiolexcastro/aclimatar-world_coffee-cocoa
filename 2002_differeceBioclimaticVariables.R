
## Fabio Castro - Llanos 
## Alliance Bioversity - CIAT 
## Difference Bioclimatic variables 

## Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, glue, dismo, geodata, gtools)
pacman::p_load(terra, fs, usdm, raster, gtools, spocc, randomForest, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth, geodata)

g <- gc(reset = T); rm(list = ls())
options(scipen = 999, warn = -1)

# Vector data -------------------------------------------------------------
wrld <- geodata::world(resolution = 1, path = './tmpr', level = 0)
wrld <- st_as_sf(wrld)

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
dfrn <- map(gcms, calc.dfrn)
dfrn <- reduce(dfrn, c)
terra::writeRaster(x = dfrn, filename = './common_data/bioclimatic_variables/dfrn-bioc_gcms.tif', overwrite = TRUE)


# To draw the maps --------------------------------------------------------

make.map <- function(gcm){
  
  ## Filtering the GCM rasters
  cat('>>> ', gcm, '\n')
  bsl <- bsln.r
  ftr <- ftre.r[[grep(gcm, names(ftre.r))]]
  dfr <- dfrn[[grep(gcm, names(dfrn))]]
  
  ## To change the names
  names(bsl) <- glue('bsl_{names(bsl)}')
  names(ftr) <- glue('ftr_{names(ftr)}')
  names(ftr) <- gsub(glue('_{gcm}'), '', names(ftr))
  names(dfr) <- gsub(glue('_{gcm}'), '', names(dfr))
  
  ## Raster to table
  stk <- c(bsl, ftr)
  tbl <- terra::as.data.frame(stk, xy = T) %>% 
    as_tibble() %>% 
    mutate(gid = 1:nrow(.)) %>% 
    gather(var, value, -c(gid, x, y))
  
  ### To separate de columns
  tbl <- tbl %>% separate(data = ., col = 'var', into = c('Period', 'Bio', 'Variable'), sep = '_')
  tbl <- mutate(tbl, Period = ifelse(Period == 'bsl', 'Baseline', 'Future'), Period = factor(Period, levels = c('Baseline', 'Future')))
  tbl <- mutate(tbl, Variable = paste0('bioc_', Variable))
  vrs <- names(stk) %>% str_sub(start = 5, end = nchar(.)) %>% unique()
  
  ### To draw the map 
  ggs <- map(.x = vrs, .f = function(v){
    
    ## Grepping
    cat(v, '\n')
    trr <- stk[[grep(paste0(v, '$'), names(stk), value = F)]]
    dfn <- dfr[[grep(paste0(v, '$'), names(stk), value = F)]]
    
    ## To change the names for the raster
    names(trr) <- c('bsl', 'ftr')
    names(dfn) <- 'dfn'
    
    ## Common theme
    common_theme <- theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.key.width = unit(3, "line")
    )
    
    ## To draw the map - Raw Values
    gg1 <- ggplot() + 
      geom_spatraster(data = trr, aes(fill = bsl)) +
      scale_fill_viridis_c(na.value = 'white') +
      geom_sf(data = wrld, fill = NA, col = 'grey30') +
      coord_sf(ylim = c(-25, 25), xlim = c(-100, 170)) +
      ggtitle(label = paste0(v, ' ', gcm, ' Baseline')) +
      labs(fill = v) +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 14),
        legend.key.width = unit(3, 'line')
      )
    
    gg2 <- ggplot() + 
      geom_spatraster(data = trr, aes(fill = ftr)) +
      scale_fill_viridis_c(na.value = 'white') +
      geom_sf(data = wrld, fill = NA, col = 'grey30') +
      coord_sf(ylim = c(-25, 25), xlim = c(-100, 170)) +
      labs(fill = v) +
      ggtitle(label = paste0(v, ' ', gcm, ' Future')) +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 14),
        legend.key.width = unit(3, 'line')
      )
    
    ## To draw the map - Difference Values
    gg3 <- ggplot() + 
      geom_spatraster(data = dfn, aes(fill = dfn)) +
      scale_fill_viridis_c(na.value = 'white') +
      geom_sf(data = wrld, fill = NA, col = 'grey30') +
      coord_sf(ylim = c(-25, 25), xlim = c(-100, 170)) +
      labs(fill = v) +
      ggtitle(label = paste0(v, ' ', gcm, ' Difference')) +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 14),
        legend.key.width = unit(3, 'line')
      )
      
    ## To join both maps into only one
    gg12 <- ggpubr::ggarrange(gg1, gg2, ncol = 2, nrow = 1, common.legend = TRUE, legend = 'bottom')
    gg13 <- ggpubr::ggarrange(gg12, gg3, ncol = 1, nrow = 2, common.legend = FALSE, legend = 'bottom')
    
    ## To save the map as a png
    ggsave(plot = gg13, filename = glue('./png/maps/bio_values/{v}_{gcm}.jpg'), units = 'in', width = 11, height = 5, dpi = 300, create.dir = T)
    return(gg13)
    
    
  })
    
  
  
}

## To table 
# trr.tbl <- trr %>% terra::as.data.frame(xy = T) %>% as_tibble()
# trr.tbl <- trr.tbl %>% gather(var, value, -c(x, y)) %>% mutate(variable = v)
# trr.tbl <- trr.tbl %>% mutate(var = gsub(paste0('_', v), '', var))
# trr.tbl <- trr.tbl %>% mutate(period = ifelse(var == 'bsl', 'Baseline', 'Future'), period = factor(period, levels = c('Baseline', 'Future')))


# dfn.tbl <- dfn %>% terra::as.data.frame(xy = T) %>% as_tibble()
# dfn.tbl <- dfn.tbl %>% setNames(c('x', 'y', 'value'))
# dfn.tbl <- dfn.tbl %>% mutate(period = 'Difference')
