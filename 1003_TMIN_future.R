# Compute maximum monthly temperature (TMIN)
# By: H. Achicanoy / J. Ramirez-Villegas / F. Castro-Llanos
# Alliance Bioversity-International & CIAT, 2025

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

# Arguments
# args <- commandArgs(trailingOnly = T)

# Root folder
root <- './common_data'

# Extent CHIRPS 
msk <- terra::rast('./common_data/atlas_hazards/cmip6/indices/ssp370_ACCESS-ESM1-5_2025_2055/PTOT/PTOT-2025-01.tif')
xtd <- terra::ext(msk)

# TMIN function
calc_tmin <- function(yr, mn){
  
  outfile <- paste0(out_dir,'/TMIN-',yr,'-',mn,'.tif')
  cat(outfile,'\n')
  
  if(!file.exists(outfile)){
    
    dir.create(dirname(outfile), F, T)
    
    # Sequence of dates
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    dts <- gsub('-', '\\.', dts)
    
    # Files
    # txfls <- paste0(tx_pth,'/tmin.',dts,'.tif') # Future
    txfls <- paste0(tx_pth,'/Tmin.',dts,'.tif') # Historical
    txfls <- txfls[file.exists(txfls)]
    
    ## Read daily maximum temperature data
    tmx <- terra::rast(txfls)
    tmx <- terra::crop(tmx, xtd)
    tmx[tmx == -9999] <- NA
    
    # Calculate maximum temperature
    tmax <- mean(tmx)
    terra::writeRaster(x = tmax, filename = outfile, overwrite = T)
    
    # Clean-up
    rm(tmx, tmax); gc(F, T, T)
    
  }
  
}

## Historical setup
yrs <- 1995:2014
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()
tx_pth <- '//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmin'   # Daily precipitation
out_dir <- paste0(root,'/atlas_hazards/historical/TMIN')

1:nrow(stp) %>%
  purrr::map(.f = function(i){calc_tmin(yr = stp$yrs[i], mn = stp$mns[i])})

# Future setup
# ssps  <- args[1]
# gcms  <- args[2]

# Runs
scenario <- 'historical' # historical, future
if (scenario == 'future') {
  ssps <- 'ssp370'
  yrs <- 2025:2055
} else {
  if (scenario == 'historical') {
    ssps <- 'historical'
    yrs <- 1995:2014
  }
}
gcms <- c( 'ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')

for (gcm in gcms) {
  
  for (ssp in ssps) {
    
    ## Parameters
    cmb <- paste0(ssp, '_', gcm, '_2025_2055')
    mnt <- sprintf('%02.0f',1:12)
    stp <- base::expand.grid(yrs, mnt, stringsAsFactors = F) |> setNames(c('yrs', 'mnt')) |> dplyr::arrange(yrs, mnt) |> base::as.data.frame(); rm(mnt)
    
    ## Setup in/out files
    tx_pth <- paste0(root, '/chirts_cmip6_world/Tmin_', gcm, '_', ssp, '_', prd) # Daily maximum temperatures
    out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/TMIN')
    
    
    1:nrow(stp) |> purrr::map(.f = function(i){calc_tmin(yr = stp$yrs[i], mn = stp$mnt[i])}); gc(F, T, T)
    tmpfls <- list.files(tempdir(), full.names = T)
    1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ', tmpfls[k]))})
    cat('----Finish----\n')
    
  }
  
}
