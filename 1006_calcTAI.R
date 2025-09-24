# Compute Thornthwaite's Aridity Index (TAI)
# By: H. Achicanoy, F. Castro
# Alliance Bioversity International & CIAT, 2025

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate,envirem))

root <- './common_data'

# Get CHIRPS extent
msk <- terra::rast('./common_data/chirps_cmip6_world/Prec_ACCESS-ESM1-5_ssp370_2025_2055/chirps-v2.0.2025.01.01.tif')
xtd <- terra::ext(msk); rm(msk)

# ET SRAD, load and process only once
srf <- list.dirs('//CATALOGUE/workspace-cluster9/DATA/ET_SolRad', full.names = T, recursive = F)
srf <- srf[-length(srf)]
srf <- srf |> gtools::mixedsort()
srd <- srf |> terra::rast()
names(srd) <- c(paste0('SRAD_0',1:9),paste0('SRAD_', 10:12))
gc(F, T, T)

# TAI function
calc_tai <- function(yr){
  
  outfile <- paste0(out_dir,'/TAI-',yr,'.tif')
  cat(outfile,'\n')
  
  if (!file.exists(outfile)) {
    
    dir.create(dirname(outfile),F,T)
    
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-01-01')), to = as.Date(paste0(yr,'-12-31')), by = 'day')
    
    # Files
    pr_fls <- paste0(pr_pth,'/chirps-v2.0.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    pr_fls <- pr_fls[file.exists(pr_fls)]
    tx_fls <- paste0(tx_pth,'/Tmax.',gsub('-', '\\.', dts),'.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tm_fls <- paste0(tm_pth,'/Tmin.',gsub('-', '\\.', dts),'.tif')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    
    # Read variables, and calculate monthly summaries, do by month to reduce memory consumption
    prc_ls <- tav_ls <- rng_ls <- c()
    for (j in 1:12) {
      
      cat('month=',j,'\n')
      # Load precipitation
      mth_fls <- paste0(pr_pth, '/chirps-v2.0.', yr, '.', sprintf('%02.0f', j), '.', sprintf('%02.0f', 1:31), '.tif')
      prc <- terra::rast(pr_fls[pr_fls %in% mth_fls])
      prc_month <- sum(prc)
      prc_month <- prc_month |> terra::crop(xtd); rm(prc); gc(F, T, T)
      
      # Load maximum temperature
      mth_fls <- paste0(tx_pth, '/Tmax.', yr, '.', sprintf('%02.0f', j), '.', sprintf('%02.0f', 1:31), '.tif')
      tmx <- terra::rast(tx_fls[tx_fls %in% mth_fls])
      ## Load minimum temperature
      mth_fls <- paste0(tm_pth, '/Tmin.', yr, '.', sprintf('%02.0f', j), '.', sprintf('%02.0f', 1:31), '.tif')
      tmn <- terra::rast(tm_fls[tm_fls %in% mth_fls])
      
      ## Calculate average temperature 
      tav <- (tmx + tmn) / 2
      tav_month <- mean(tav)
      tav_month <- tav_month |> terra::crop(xtd); rm(tav); gc(F, T, T)
      
      ## Calculate temperature range
      rnge <- abs(tmx - tmn)
      rng_month <- mean(rnge)
      rng_month <- rng_month |> terra::crop(xtd); rm(rnge, tmx, tmn); gc(F, T, T)
      
      # Append
      prc_ls <- c(prc_ls, prc_month)
      tav_ls <- c(tav_ls, tav_month)
      rng_ls <- c(rng_ls, rng_month)
    }; gc(F, T, T)
    
    # Assign precipitation names in envirem environment
    envirem::assignNames(solrad='SRAD_##', tmean = 'TMEAN_##', precip = 'PREC_##')
    
    TMEAN <- tav_ls |> terra::rast()
    PREC  <- prc_ls |> terra::rast()
    TRNG  <- rng_ls |> terra::rast()
    
    names(TMEAN) <- c(paste0('TMEAN_0',1:9),paste0('TMEAN_', 10:12))
    names(PREC)  <- c(paste0('PREC_0',1:9),paste0('PREC_', 10:12))
    names(TRNG)  <- c(paste0('TRNG_0',1:9),paste0('TRNG_', 10:12))
    
    ## Remove weird values
    PREC[PREC < 0] <- NA
    TMEAN[TMEAN == -9999] <- NA
    TRNG <- terra::crop(TRNG, wrld) %>% terra::mask(wrld)
    
    # Clean-up
    rm(tav_ls, prc_ls, rng_ls, tav_month, prc_month, rng_month); gc(F, T, T)
    
    # Resample extraterrestrial solar radiation to 0.25
    srd <- srd |> terra::resample(TMEAN[[1]])
    srd <- srd |> terra::crop(xtd)
    srd <- srd %>% terra::mask(wrld)
    gc(F, T, T)
    
    # Thornthwaite's Aridity Index
    PET <- envirem::monthlyPET(Tmean = TMEAN, RA = srd, TD = TRNG)
    names(PET)  <- c(paste0('PET_0',1:9), paste0('PET_', 10:12))
    TAI <- envirem::aridityIndexThornthwaite(precipStack = PREC, PETstack = PET)
    names(TAI) <- yr
    
    # Write output
    terra::writeRaster(x = TAI, filename = outfile, overwrite = T)
    
    # Clean-up
    rm(PET, TMEAN, PREC, TRNG, TAI); gc(F, T, T)
    
  }
}

# Runs
scenario <- 'historical' # historical, future
if(scenario == 'future') {
  ssps <- 'ssp370'
  yrs <- 2025:2055
} else {
  if (scenario == 'historical') {
    ssps <- 'historical'
    yrs <- 1995:2014
  }
}

# Historical
yrs <- 1995:2014

pr_pth <- '//catalogue/WFP_ClimateRiskPr1/1.Data/Chirps' # Precipitation
tm_pth <- '//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmin'# Minimum temperature
tx_pth <- '//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmax' # Maximum temperature
out_dir <- '//catalogue/aclimatar-world/aclimatar-w_v1/common_data/atlas_hazards/historical/TAI'

yrs |> purrr::map(.f = function(yr){ calc_tai(yr = yr) }); gc(F, T, T)

# Future

gcms <- c('ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')

for (gcm in gcms) {
  
  for (ssp in ssps) {
    
    cmb <- paste0(ssp, '_', gcm, '_', '2025_2055')
    cat(cmb, '\n')
    
    pr_pth <- paste0(root,'/chirps_cmip6_world/Prec_',gcm,'_',ssp,'_',prd) # Precipitation
    tm_pth <- paste0(root, '/chirts_cmip6_world/Tmin_', gcm, '_', ssp, '_', prd)# Minimum temperature
    tx_pth <- paste0(root, '/chirts_cmip6_world/Tmax_', gcm, '_', ssp, '_', prd) # Maximum temperature
    out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/TAI')
    
    yrs |> purrr::map(.f = function(yr){ calc_tai(yr = yr) }); gc(F, T, T)
    tmpfls <- list.files(tempdir(), full.names = T)
    1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ', tmpfls[k]))})
    
  }
  
}



# Check TAI bsaeline  -----------------------------------------------------
library(fs)
fls <- dir_ls('./common_data/atlas_hazards/historical/TAI')
rst <- map(fls, rast)

# Check TAI Future --------------------------------------------------------
fls <- dir_ls('./common_data/atlas_hazards/cmip6/indices/ssp370_ACCESS-ESM1-5_2025_2055/TAI')
rst <- map(fls, rast)
