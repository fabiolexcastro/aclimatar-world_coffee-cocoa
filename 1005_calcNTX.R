# Heat stress generic crop and threshold (i.e., NTx40)
# By: H. Achicanoy, F. Castro-Llanos
# Alliance Bioversity International & CIAT, 2025

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- './common_data'

# Get CHIRPS extent
msk <- terra::rast('./common_data/atlas_hazards/cmip6/indices/ssp370_ACCESS-ESM1-5_2025_2055/PTOT/PTOT-2025-01.tif')
xtd <- terra::ext(msk)

# NTx function
calc_ntx <- function(yr, mn, thr = 40) {
  
  outfile <- paste0(out_dir,'/NTx',thr,'/NTx',thr,'-',yr,'-',mn,'.tif') 
  thr <- thr[!file.exists(outfile)]
  outfile <- outfile[!file.exists(outfile)]
  
  if (length(outfile) > 0) {
    
    cat('...processing n=', length(outfile), 'files for yr=', yr, '/ mn=', mn, '\n')
    
    # Create directories
    1:length(outfile) |> purrr::map(.f = function(j){dir.create(dirname(outfile[j]),F,T)})
    
    # Sequence of dates
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    dts <- gsub('-', '\\.', dts)
    
    # Files
    txfls <- paste0(tx_pth,'/Tmax.',dts,'.tif')
    txfls <- txfls[file.exists(txfls)]
    
    # Read daily maximum temperature data
    tmx <- terra::rast(txfls)
    tmx <- terra::crop(tmx, xtd)
    
    # Calculate heat stress generic crop
    for (j in 1:length(thr)) {
      cat('...processing threshold thr=',thr[j],'\n')
      ntx <- sum(tmx > thr[j])
      terra::writeRaster(x = ntx, filename = outfile[j], overwrite = T)
    }
    
    # Clean-up
    rm(tmx, ntx); gc(F, T, T)
    
  }
  
}

# Runs
scenario <- 'historical' # historical, future
if (scenario == 'future') {
  ssps <- c('ssp370')
  yrs <- 2025:2055
} else {
  if (scenario == 'historical') {
    ssps <- 'historical'
    yrs <- 1995:2014
  }
}


# Historical
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mnt')
stp <- stp %>% dplyr::arrange(yrs, mns) %>% base::as.data.frame()
tx_pth <- '//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmax' 
out_dir <- paste0(root,'/atlas_hazards/historical')
thr <- c(30, 35)

1:nrow(stp) |> purrr::map(.f = function(i){calc_ntx(yr = stp$yrs[i], mn = stp$mnt[i], thr = c(30, 35))}); gc(F, T, T)

# Future
gcms <- c('ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')
prd  <- '2025_2055'

for (gcm in gcms){
  
  for (ssp in ssps) {
    
    ## Parameters
    cmb <- paste0(ssp, '_', gcm, '_', prd)
    mnt <- sprintf('%02.0f',1:12)
    stp <- base::expand.grid(yrs, mnt, stringsAsFactors = F) |> base::as.data.frame() |> setNames(c('yrs', 'mnt')) |> dplyr::arrange(yrs, mnt) |> base::as.data.frame(); rm(mnt)
    
    ## Setup in/out files
    tx_pth <- paste0(root, '/chirts_cmip6_world/Tmax_', gcm, '_', ssp, '_', prd) # Daily maximum temperatures
    thr <- c(30, 35) # 35
    out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb)
    
    1:nrow(stp) |> purrr::map(.f = function(i){calc_ntx(yr = stp$yrs[i], mn = stp$mnt[i], thr = thr)}); gc(F, T, T)
    tmpfls <- list.files(tempdir(), full.names = T)
    1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ', tmpfls[k]))})
    cat('----Finish----\n')
    
  }
  
}
