
## Fabio Castro - Llanos 
## Alliance Bioversity - CIAT 
## Calculate Bioclimatic variables 

## Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, sf, fs, tidyverse, glue, dismo, geodata, gtools)
pacman::p_load(terra, fs, usdm, raster, gtools, spocc, randomForest, sf, outliers, glue, CoordinateCleaner, tidyverse, rgbif, readxl, xlsx, openxlsx, rnaturalearthdata, rnaturalearth, geodata)

g <- gc(reset = T); rm(list = ls())
options(scipen = 999, warn = -1)

# Solar radiation  --------------------------------------------------------
srad <- as.character(dir_ls('//catalogue/workspace-cluster9/DATA/ET_SolRad'))
srad <- srad[-grep('info', srad)]
srad <- mixedsort(srad)
srad <- terra::rast(srad)

# Raster mask -------------------------------------------------------------
msk <- terra::rast('./common_data/atlas_hazards/cmip6/indices/ssp370_ACCESS-ESM1-5_2025_2055/PTOT/PTOT-2025-01.tif')
xtd <- terra::ext(msk)
wrl <- geodata::world(resolution = 1, path = './tmpr')

## 
srad <- terra::crop(srad, xtd)
srad <- terra::resample(srad, msk, method = 'bilinear')

# Baseline --------------------------------------------------------------
prec <- dir_ls('//catalogue/WFP_ClimateRiskPr1/1.Data/Chirps')
tmin <- dir_ls('//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmin')
tmax <- dir_ls('//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmax')

## Daily to monthly
daily.monthly <- function(fles, varb){
  
  # fles <- tmin
  # varb <- 'tmin'
  
  cat('To process!\n')
  year <- 1995:2014
  outd <- glue('./common_data/bioclimatic_variables/historical')
  
  trra <- map(.x = year, .f = function(y){
    
    cat('Year: ', y, '\n')
    
    if(!file.exists(glue('{outd}/{varb}_{y}_1-12.tif'))){
      
      if(varb == 'prec'){
        fls <- as.character(grep(paste0('0.', y, '.'), fles, value = T))
      } else {
        fls <- as.character(grep(paste0('.', y, '.'), fles, value = T))
      }
      
      ## By each month
      rstr <- map(.x = c(paste0('0', 1:9), 10:12), .f = function(m){
        
        cat('Month: ', m, '\n')
        fle <- grep(paste0('.', y, '.', m), fls, value = T)
        rst <- rast(fle)
        rst <- terra::crop(rst, xtd)
        
        if(varb == 'prec'){
          rsl <- sum(rst)
        } else {
          rsl <- mean(rst)
          rsl[rsl == -9999] <- NA
        }
        
        names(rsl) <- glue('{varb}_{m}')
        return(rsl)
        
      })
      
      ## Reduce as a stack 
      rstr <- reduce(rstr, c)
      names(rstr) <- glue('{varb}_{y}-{1:12}')
      
      ## To write the raster 
      cat('Done!\n')
      terra::writeRaster(rstr, filename = glue('{outd}/{varb}_{y}_1-12.tif'), overwrite = TRUE)
      rm(rstr); gc(reset = TRUE)
      
    }
    
  })
  
}

##
map2(list(prec, tmin, tmax), list('prec', 'tmin', 'tmax'), daily.monthly)

## Multi-monthly average 
prec <- as.character(dir_ls('./common_data/bioclimatic_variables/historical', regexp = 'prec'))
tmin <- as.character(dir_ls('./common_data/bioclimatic_variables/historical', regexp = 'tmin'))
tmax <- as.character(dir_ls('./common_data/bioclimatic_variables/historical', regexp = 'tmax'))

###
avrg.monthly <- function(fles){
  
  ## To start!
  cat('To process!\n')
  dout <- './common_data/bioclimatic_variables'
  dnme <- str_split(basename(fles), '_') %>% map_chr(1) %>% unique()
  outf <- glue('{dout}/{dnme}_hist.tif')
  
  if(!file.exists(outf)){
    rstr <- map(.x = 1:12, .f = function(i){
      
      cat('Month: ', month.abb[i], '\n')
      rst <- reduce(map(fles, rast, lyr = i), c)
      var <- str_split(basename(fles), '_') %>% map_chr(1) %>% unique()
      
      ## Condicional for prec < 0
      if(var == 'prec'){
        rst[rst < 0] <- NA
      }
      
      ## To calculate the average
      avg <- mean(rst)
      names(avg) <- glue('{var}_{i}')
      return(avg)  
      
    })
    rstr <- reduce(rstr, c)
    rstr  
    terra::writeRaster(x = rstr, filename = outf, overwrite = TRUE)
    rm(rstr); gc(reset = TRUE)
    
  }
  
  ## Finish 
  rm(dout, dnme, outf); gc(reset = TRUE)
  cat('Finish!\n')
  
}

###
map(list(prec, tmin, tmax), avrg.monthly)

## To estimate the bioclimatic variables
prec <- terra::rast('./common_data/bioclimatic_variables/prec_hist.tif')
tmin <- terra::rast('./common_data/bioclimatic_variables/tmin_hist.tif')
tmax <- terra::rast('./common_data/bioclimatic_variables/tmax_hist.tif')
tavg <- (tmax + tmin) / 2
names(tavg) <- glue('tmean_{1:12}')
etps <- terra::rast('./common_data/bioclimatic_variables/etps_hist.tif')
etps <- 0.0013 * 0.408 * srad * (tavg + 17) * (tmax - tmin - 0.0123 * prec) ^ 0.76
etps <- etps * c(31,29,31,30,31,30,31,31,30,31,30,31)
for(i in 1:12){etps[[i]][which.lyr(is.na(etps[[i]]))] <- 0}
etps <- terra::mask(etps, wrl)
terra::writeRaster(x = etps, filename = './common_data/bioclimatic_variables/etps_hist.tif')
etps <- terra::rast('./common_data/bioclimatic_variables/etps_hist.tif')

## To remove the NAs
stck <- c(prec, tmin, tmax, etps, tavg)
tble <- terra::as.data.frame(stck, xy = T, na.rm = T)
nrow(tble) == nrow(drop_na(tble))
stck <- terra::rast(tble, type = 'xyz', crs = 'EPSG:4326')
prec <- stck[[grep('prec', names(stck))]]
tmin <- stck[[grep('tmin', names(stck))]]
tmax <- stck[[grep('tmax', names(stck))]]
tavg <- stck[[grep('tmea', names(stck))]]
etps <- stck[[grep('et_s', names(stck))]]

## Bioclimatic variables 
bios <- dismo::biovars(prec = as.matrix(as.data.frame(prec)), tmin = as.matrix(as.data.frame(tmin)), tmax = as.matrix(as.data.frame(tmax)))
crds <- terra::crds(prec, na.rm = T)
bios <- cbind(crds, bios)
bios <- rast(bios, type = 'xyz', crs = 'EPSG:4326')
names(bios) <- glue('bioc_{1:19}')
terra::writeRaster(x = bios, filename = './common_data/bioclimatic_variables/bios_hist.tif', overwrite = TRUE)
bios <- terra::rast('./common_data/bioclimatic_variables/bios_hist.tif')

## ETP bioclimatic variables
source('./bioclimatic_functions.R')
names(etps) <- paste0('etp_', 1:12)
names(prec) <- paste0('prec_', 1:12)
names(tavg) <- paste0('tmean_', 1:12)

stck <- c(etps, prec, tavg, tmax, tmin)
tble <- terra::as.data.frame(stck, xy = T)
tble <- drop_na(tble)
stck <- terra::rast(tble, type = 'xyz', crs = 'EPSG:4326')
prec <- terra::as.data.frame(stck[[grep('prec', names(stck))]])
etps <- terra::as.data.frame(stck[[grep('etp_', names(stck))]])
tavg <- terra::as.data.frame(stck[[grep('tmea', names(stck))]])
tmax <- terra::as.data.frame(stck[[grep('tmax', names(stck))]])
tmin <- terra::as.data.frame(stck[[grep('tmin', names(stck))]])

nrow(tble)
nrow(drop_na(tble))

ETPAndPrec <- cbind(as.matrix(etps),as.matrix(prec),as.matrix(tavg))
etpbios    <- t(apply(ETPAndPrec, 1, etpvars))
etpbios <- cbind(crds(stck), etpbios)
etpbios <- rast(etpbios, type = 'xyz', crs = 'EPSG:4326')
names(etpbios) <- glue("bioc_{21:29}")

# Balance variables -------------------------------------------------------
prec <- rast(cbind(crds(stck), prec), type = 'xyz', crs = 'EPSG:4326')
etps <- rast(cbind(crds(stck), etps), type = 'xyz', crs = 'EPSG:4326')
tavg <- rast(cbind(crds(stck), tavg), type = 'xyz', crs = 'EPSG:4326')
tmax <- rast(cbind(crds(stck), tmax), type = 'xyz', crs = 'EPSG:4326')

defc <- prec - etps
DefAndTemp <- cbind(as.matrix(defc), as.matrix(tavg), as.matrix(tmax))
biovalues  <- t(apply(DefAndTemp, 1, cumTemp))
biovalues  <- cbind(as.data.frame(crds(defc, na.rm = F)), as.data.frame(biovalues))
biovalues  <- rast(biovalues, type = 'xyz', crs = 'EPSG:4326')
names(biovalues) <- glue('bioc_{30:33}')

terra::writeRaster(x = biovalues, filename = './common_data/bioclimatic_variables/bios-balance_hist.tif', overwrite = TRUE)
terra::writeRaster(x = etpbios, filename = './common_data/bioclimatic_variables/bios-etps_hist.tif', overwrite = TRUE)

# Bio 20 ------------------------------------------------------------------
prec <- terra::rast('./common_data/bioclimatic_variables/prec_hist.tif')
prec.rclf <- terra::rast(raster::reclassify(stack(prec), c(-Inf, 60, 1, 60, Inf, NA)))
prec.two  <- c(prec.rclf, prec.rclf)
prec.two  <- stack(prec.two)
allperiods <- stack()
for(i in 1:12){
  oneyear <- prec.two[[i:(i+11)]]
  drymonths <- cumsum(oneyear)
  maxnumber <- max(drymonths, na.rm = T)
  allperiods <- addLayer(allperiods, maxnumber) 
}
bio20 <- max(allperiods, na.rm = T)
bio20[is.na(bio20)] <- 0
bio20 <- rast(bio20) %>% crop(wrl) |> mask(wrl)
# terra::writeRaster(x = bio20, filename = './common_data/bioclimatic_variables/bioc20_100mm.tif', overwrite = TRUE)
# terra::writeRaster(x = bio20, filename = './common_data/bioclimatic_variables/bioc20_70mm.tif', overwrite = TRUE)
# terra::writeRaster(x = bio20, filename = './common_data/bioclimatic_variables/bioc20_60mm.tif', overwrite = TRUE)

# Future ------------------------------------------------------------------

dirs <- as.character(dir_ls('./common_data/atlas_hazards/cmip6/indices', type = 'directory'))

## Precipitation 
dirs.prec <- dirs %>% map(dir_ls) %>% unlist() %>% grep('PTOT', ., value = T) %>% as.character() 
prec <- map(.x = dirs.prec, .f = function(i){
  fles <- i %>% dir_ls() %>% as.character()
  gcme <- dirname(fles) %>% dirname() %>% basename() %>% str_split('_') %>% map_chr(2) %>% unique()
  rstr <- map(.x = c(paste0('0', 1:9), 10:12), .f = function(j){
    rstr <- grep(paste0('-', j, '.tif'), fles, value = T) %>% rast()
    rstr <- mean(rstr)
    return(rstr)
  })
  rstr <- reduce(rstr, c)
  names(rstr) <- glue('prec_{1:12}_{gcme}')
  return(rstr)
})
prec <- reduce(prec, c)

dirs.tmin <- dirs %>% map(dir_ls) %>% unlist() %>% grep('TMIN', ., value = T) %>% as.character() 
tmin <- map(.x = dirs.tmin, .f = function(i){
  fles <- i %>% dir_ls() %>% as.character()
  gcme <- dirname(fles) %>% dirname() %>% basename() %>% str_split('_') %>% map_chr(2) %>% unique()
  rstr <- map(.x = c(paste0('0', 1:9), 10:12), .f = function(j){
    rstr <- grep(paste0('-', j, '.tif'), fles, value = T) %>% rast()
    rstr <- mean(rstr)
    return(rstr)
  })
  rstr <- reduce(rstr, c)
  names(rstr) <- glue('tmin_{1:12}_{gcme}')
  return(rstr)
})
tmin <- reduce(tmin, c)

dirs.tmax <- dirs %>% map(dir_ls) %>% unlist() %>% grep('TMAX', ., value = T) %>% as.character() 
tmax <- map(.x = dirs.tmax, .f = function(i){
  fles <- i %>% dir_ls() %>% as.character()
  gcme <- dirname(fles) %>% dirname() %>% basename() %>% str_split('_') %>% map_chr(2) %>% unique()
  rstr <- map(.x = c(paste0('0', 1:9), 10:12), .f = function(j){
    rstr <- grep(paste0('-', j, '.tif'), fles, value = T) %>% rast()
    rstr <- mean(rstr)
    return(rstr)
  })
  rstr <- reduce(rstr, c)
  names(rstr) <- glue('tmax_{1:12}_{gcme}')
  return(rstr)
})
tmax <- reduce(tmax, c)

## To write the raster 
terra::writeRaster(x = prec, filename = './common_data/chirps_cmip6_world/Prec_GCMs-Avrg_ssp370_2025_2055.tif', overwrite = TRUE)
terra::writeRaster(x = tmax, filename = './common_data/chirts_cmip6_world/Tmax_GCMs-Avrg_ssp370_2025_2055.tif', overwrite = TRUE)
terra::writeRaster(x = tmin, filename = './common_data/chirts_cmip6_world/Tmin_GCMs-Avrg_ssp370_2025_2055.tif', overwrite = TRUE)

# To calcate bioclimatic variables ----------------------------------------
gcms <- names(prec) %>% str_split('_') %>% map_chr(3) %>% unique()
calc.bioc <- function(gcme){
  
  ## To list the files
  gcme <- gcms[5]
  cat('To process: ', gcme, '\n')
  ppt <- prec[[grep(gcme, names(prec), value = F)]]
  tmn <- tmin[[grep(gcme, names(tmin), value = F)]]  
  tmx <- tmax[[grep(gcme, names(tmax), value = F)]]
  
  ## Make a stack
  stk <- c(ppt, tmn, tmx)
  tbl <- terra::as.data.frame(stk, xy = T, na.rm = F)
  tbl <- drop_na(tbl)
  wrl <- geodata::world(resolution = 1, level = 0, path = './tmpr')
  
  ## Calc bioclim 
  bio.mtx <- dismo::biovars(
    prec = as.matrix(dplyr::select(tbl, starts_with('prec'))), 
    tmin = as.matrix(dplyr::select(tbl, starts_with('tmin'))), 
    tmax = as.matrix(dplyr::select(tbl, starts_with('tmax')))
  )
  
  bio.mtx <- cbind(tbl[,1:2], bio.mtx)
  bio <- rast(bio.mtx, type = 'xyz', crs = 'EPSG:4326')
  names(bio) <- glue('bioc_{1:19}_{gcme}')
  
  ## Bio 20
  make.bioc <- function(thr){
    # thr <- 60
    cat('To process: ', thr, '\n')
    prec.rclf <- terra::rast(raster::reclassify(stack(ppt), c(-Inf, thr, 1, thr, Inf, NA)))
    prec.two  <- c(prec.rclf, prec.rclf)
    prec.two  <- stack(prec.two)
    allperiods <- stack()
    for(i in 1:12){
      oneyear <- prec.two[[i:(i+11)]]
      drymonths <- cumsum(oneyear)
      maxnumber <- max(drymonths, na.rm = T)
      allperiods <- addLayer(allperiods, maxnumber) 
    }
    bio20 <- max(allperiods, na.rm = T)
    bio20[is.na(bio20)] <- 0
    bio20 <- rast(bio20) %>% crop(wrl) |> mask(wrl)
    return(bio20)
  }
  bio20.60 <- make.bioc(thr = 60)
  bio20.100 <- make.bioc(thr = 100) 
  
  ## ETP 
  tavg <- (tmx + tmn) / 2
  names(tavg) <- glue('tmean_{1:12}')
  etps <- 0.0013 * 0.408 * srad * (tavg + 17) * (tmx - tmn - 0.0123 * ppt) ^ 0.76
  etps <- etps * c(31,29,31,30,31,30,31,31,30,31,30,31)
  for(i in 1:12){etps[[i]][which.lyr(is.na(etps[[i]]))] <- 0}
  etps <- terra::mask(etps, wrl)
  
  ## Clean raster
  stk <- c(ppt, tmn, tmx, tavg, etps)
  tbl <- terra::as.data.frame(stk, xy = T, na.rm = T)
  
  ## BIO ETPs
  ETPAndPrec <- cbind(
    as.matrix(dplyr::select(tbl, starts_with('et_'))),
    as.matrix(dplyr::select(tbl, starts_with('prec'))),
    as.matrix(dplyr::select(tbl, starts_with('tmean')))
  )
  etpbios    <- t(apply(ETPAndPrec, 1, etpvars))
  etpbios <- cbind(tbl[,1:2], etpbios)
  etpbios <- rast(etpbios, type = 'xyz', crs = 'EPSG:4326')
  
  ## Balance
  defc <- ppt - etps
  DefAndTemp <- cbind(as.matrix(defc), as.matrix(tavg), as.matrix(tmx))
  biovalues  <- t(apply(DefAndTemp, 1, cumTemp))
  rstbiovalues  <- cbind(crds(defc, na.rm = F), biovalues)
  rstbiovalues  <- rast(rstbiovalues, type = 'xyz', crs = 'EPSG:4326')
  names(rstbiovalues) <- paste0('bioc_', 30:33)
  plot(rstbiovalues[[3]])
  
  ## Make just one stack 
  bio20.60 <- terra::resample(bio20.60, bio[[1]], method = 'near')
  bio20.100 <- terra::resample(bio20.100, bio[[1]], method = 'near')
  rstbiovalues <- terra::resample(rstbiovalues, bio[[1]], method = 'bilinear') 
  fnal <- c(bio, bio20.60, bio20.100, etpbios, rstbiovalues)
  names(fnal) <- c(glue('bioc_{1:19}_{gcme}'), glue('bioc_2060_{gcme}'), glue('bioc_20100_{gcme}'), glue('bioc_{21:33}_{gcme}'))
  terra::writeRaster(x = fnal, filename = glue('./common_data/bioclimatic_variables/bioc_all_{gcme}.tif'), overwrite = TRUE)
  rm(fnal); gc(reset = T)
}
map(gmcs, calc.bioc)


