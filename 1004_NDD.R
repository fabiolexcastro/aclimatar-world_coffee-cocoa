# Compute number of dry days (NDD)
# By: H. Achicanoy, F. Castro
# Alliance Bioversity International & CIAT, 2025

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- './common_data'

# Get CHIRPS extent
msk <- terra::rast('./common_data/atlas_hazards/cmip6/indices/ssp370_ACCESS-ESM1-5_2025_2055/PTOT/PTOT-2025-01.tif')
xtd <- terra::ext(msk); rm(msk)

# NDD function
calc_ndd <- function(yr, mn){
  
  outfile <- paste0(out_dir,'/NDD-',yr,'-',mn,'.tif')
  cat(outfile,'\n')
  
  if(!file.exists(outfile)){
    
    dir.create(dirname(outfile),F,T)
    
    # Sequence of dates
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    
    # Files
    fls <- paste0(pr_pth,'/chirps-v2.0.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    fls <- fls[file.exists(fls)]
    
    # Read daily precipitation data
    prc <- terra::rast(fls)
    prc <- terra::crop(prc, xtd)
    prc[prc == -9999] <- NA
    
    # Calculate number of dry days
    ndd <- sum(prc < 1)
    terra::writeRaster(x = ndd, filename = outfile, overwrite = T)
    
    # Clean-up
    rm(prc, ndd); gc(F, T, T)
    
  }
}

# Historical
yrs <- 1995:2014
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>% dplyr::arrange(yrs, mns) %>% base::as.data.frame()
pr_pth <- '//catalogue/WFP_ClimateRiskPr1/1.Data/Chirps'  # Daily precipitation
out_dir <- paste0(root,'/atlas_hazards/historical/NDD')

1:nrow(stp) %>%
  purrr::map(.f = function(i){calc_ndd(yr = stp$yrs[i], mn = stp$mns[i])})


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
gcms <- c( 'ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')

for (gcm in gcms) {
  
  for (ssp in ssps) {
    
    ## Parameters
    ssp <- 'ssp370'
    prd <- '2025_2055'
    
    cmb <- paste0(ssp,'_',gcm, '_', prd)
    mnt <- sprintf('%02.0f',1:12)
    stp <- base::expand.grid(yrs, mnt, stringsAsFactors = F) |> base::as.data.frame() |> setNames(c('yrs', 'mnt')) |> dplyr::arrange(yrs, mnt) |> base::as.data.frame(); rm(mnt)
    
    ## Setup in/out files
    pr_pth <- paste0(root,'/chirps_cmip6_world/Prec_',gcm,'_',ssp,'_',prd) 
    out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/NDD')
    
    1:nrow(stp) |> purrr::map(.f = function(i){calc_ndd(yr = stp$yrs[i], mn = stp$mnt[i])}); gc(F, T, T)
    tmpfls <- list.files(tempdir(), full.names = T)
    1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ', tmpfls[k]))})
    cat('----Finish----\n')
    
  }
  
}
