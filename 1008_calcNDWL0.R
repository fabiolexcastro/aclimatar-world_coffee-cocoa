# Compute number of waterlogging days at start of saturation (NDWL0)
# By: H. Achicanoy & A. Mendez
# Alliance Bioversity International & CIAT, 2025

# R options
# args <- commandArgs(trailingOnly = T)
options(warn = -1, scipen = 999) # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate,compiler,raster,ncdf4))

peest2 <- function(srad, tmin, tmean, tmax){
  
  # Convert rasters to matrices for faster processing
  srad_m  <- terra::values(srad)
  tmin_m  <- terra::values(tmin)
  tmean_m <- terra::values(tmean)
  tmax_m  <- terra::values(tmax)
  
  # Constants
  albedo  <- 0.2
  vpd_cte <- 0.7
  a_eslope <- 611.2
  b_eslope <- 17.67
  c_eslope <- 243.5
  
  # Pre-allocate matrices
  n_cells <- ncol(srad_m)
  n_layers <- nrow(srad_m)
  rn <- matrix(0, nrow = n_layers, ncol = n_cells)
  eslope <- matrix(0, nrow = n_layers, ncol = n_cells)
  et_max <- matrix(0, nrow = n_layers, ncol = n_cells)
  
  # Optimized calculations
  rn <- (1-albedo) * srad_m
  rm(srad_m)
  
  temp_denom <- tmean_m + c_eslope
  eslope <- (a_eslope * b_eslope * c_eslope * exp(b_eslope * tmean_m/temp_denom)) / (temp_denom^2)
  rm(tmean_m, temp_denom)
  
  vpd <- vpd_cte * (0.61120 * (exp((17.67*tmax_m)/(tmax_m+243.5)) - 
                                 exp((17.67*tmin_m)/(tmin_m+243.5))))
  rm(tmin_m, tmax_m)
  
  # Constants
  pt_const <- 1.26
  pt_fact  <- 1
  vpd_ref  <- 1
  psycho   <- 62
  rho_w    <- 997
  rlat_ht  <- 2.26E6
  
  # Final calculation
  pt_coef <- 1 + (pt_fact*pt_const-1) * vpd / vpd_ref
  conversion_factor <- 1E6 * 100 / (rlat_ht * rho_w) * 10
  et_max <- pt_coef * rn * eslope/(eslope+psycho) * conversion_factor
  
  # Convert back to raster
  et_max_rast <- tmin
  terra::values(et_max_rast) <- et_max
  
  return(et_max_rast)
  
}

root <- './common_data'

# Get CHIRPS extent
msk <- terra::rast('./common_data/atlas_hazards/cmip6/indices/ssp370_ACCESS-ESM1-5_2025_2055/PTOT/PTOT-2025-01.tif')
xtd <- terra::ext(msk)

# Soil variables
scp <- terra::rast(paste0(root,'/soils/sscp_world.tif'))
sst <- terra::rast(paste0(root,'/soils/ssat_world.tif'))

# NDWL0 function
calc_ndwl0 <- function(yr, mn){
  
  outfile <- paste0(out_dir,'/NDWL0-',yr,'-',mn,'.tif')
  cat(outfile,'\n')
  
  if (!file.exists(outfile)) {
    
    dir.create(dirname(outfile),F,T)
    
    # Sequence of dates
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    
    cat('>>> Iniciando proceso:',yr,'-',mn,'\n')
    
    # Files
    pr_fls <- paste0(pr_pth,'/chirps-v2.0.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    pr_fls <- pr_fls[file.exists(pr_fls)]
    tx_fls <- paste0(tx_pth,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tm_fls <-  paste0(tm_pth,'/Tmin.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    sr_fls <- paste0(sr_pth, '/Solar-Radiation-Flux_C3S-glob-agric_AgERA5_', gsub('-', '', dts), '_final-v1.1.nc')
    sr_fls <- sr_fls[file.exists(sr_fls)]
    
    # Read variables
    prc <- terra::rast(pr_fls)
    prc <- prc |> terra::crop(xtd)
    tmx <- terra::rast(tx_fls)
    tmx <- tmx |> terra::crop(xtd)
    tmn <- terra::rast(tm_fls)
    tmn <- tmn |> terra::crop(xtd)
    tav <- (tmx + tmn)/2
    srd <- terra::rast(sr_fls)
    srd <- srd |> terra::crop(xtd)
    srd <- terra::resample(srd, prc[[1]], method = 'bilinear')
    srd <- srd / 1000000 # Just for historical
    
    ## -9999 [Historical]
    prc[prc < 0] <- NA
    tmx[tmx == -9999] <- NA
    tmn[tmn == -9999] <- NA
    tav[tav == -9999] <- NA
    
    # Maximum evapotranspiration
    ETMAX <<- peest2(srad = srd, tmin = tmn, tmean = tav, tmax = tmx)
    rm(list = c('tmn','tmx','tav','srd'))
    
    scp <- scp |> terra::resample(ETMAX[[1]])
    scp <- scp |> terra::crop(xtd)
    sst <- sst |> terra::resample(ETMAX[[1]])
    sst <- sst |> terra::crop(xtd)
    
    # Compute water balance model
    date <- paste0(yr,'-',mn)
    if(date %in% c('1995-01','2025-01')){
      AVAIL <<- ETMAX[[1]] * 0
      AVAIL[!is.na(AVAIL)] <- 0
    } else {
      avail_fl <- list.files(path = dirname(outfile), pattern = 'AVAIL-')
      avail_fl <- avail_fl[grep(pattern = '\\.tif', avail_fl)]
      avail_fl <- avail_fl[length(avail_fl)]
      AVAIL <<- terra::rast(paste0(dirname(outfile),'/',avail_fl))
      AVAIL <<- AVAIL[[terra::nlyr(AVAIL)]]
    }
    eabyep_calc <- function(soilcp = scp, soilsat = ssat, avail = AVAIL,rain = prc[[1]], evap = ETMAX[[1]]){
      
      avail   <- min(avail, soilcp)
      
      # ERATIO
      percwt <- min(avail/soilcp*100, 100)
      percwt <- max(percwt, 1)
      eratio <- min(percwt/(97-3.868*sqrt(soilcp)), 1)
      
      demand  <- eratio * evap
      result  <- avail + rain - demand
      logging <- result - soilcp
      logging <- max(logging, 0)
      logging <- min(logging, soilsat)
      # runoff  <- result - logging + soilcp
      avail   <- min(soilcp, result)
      avail   <- max(avail, 0)
      # runoff  <- max(runoff, 0)
      
      out     <- list(Availability = c(AVAIL, avail),
                      # Demand       = demand,
                      # Eratio       = eratio
                      Logging      = logging
      )
      return(out)
    }
    ceabyep_calc <- compiler::cmpfun(eabyep_calc)
    
    watbal <- vector('list', terra::nlyr(ETMAX))
    for(j in 1:terra::nlyr(ETMAX)){
      water_balance <- ceabyep_calc(soilcp  = scp,
                                    soilsat = sst,
                                    avail   = AVAIL[[terra::nlyr(AVAIL)]],
                                    rain    = prc[[j]],
                                    evap    = ETMAX[[j]])
      # Update AVAIL with deep copy to avoid memory leaks
      AVAIL <- terra::deepcopy(water_balance$Availability)
      # Store result and clean temporary objects
      watbal[[j]] <- water_balance
      rm(water_balance)
    }
    
    LOGGING <- watbal |> purrr::map('Logging') |> terra::rast()
    
    # Calculate number of soil waterlogging days (if logging is above 0)
    # Note NDWL50 uses ssat * 0.5 (so this means soil is at 50% toward saturation)
    # here it suffices if the soil is above field capacity
    NDWL0  <- sum(LOGGING > 0)
    terra::writeRaster(NDWL0, outfile)
    terra::writeRaster(AVAIL, paste0(dirname(outfile),'/AVAIL-',yr,'-',mn,'.tif'))
    
    ## Clean up
    rm(list = c('prc','ETMAX','AVAIL','watbal','ERATIO','LOGGING','NDWL0'))
    gc(F, T, T)
  }
}

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


## Historical
yrs <- 1995:2014
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>% dplyr::arrange(yrs, mns) %>% base::as.data.frame()

pr_pth <- '//catalogue/WFP_ClimateRiskPr1/1.Data/Chirps'
tx_pth <- '//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmax' 
tm_pth <- '//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmin' 
sr_pth <- '//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5/solar_radiation_flux'
out_dir <- paste0(root,'/atlas_hazards/historical/NDWL0')

1:nrow(stp) %>%
  purrr::map(.f = function(i){calc_ndwl0(yr = stp$yrs[i], mn = stp$mns[i])})

## Future
gcms <- c('ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')
gcm <- gcms[1]
ssp <- ssps[1]

for (gcm in gcms) {
  
  for (ssp in ssps) {
    
    prd <- '2025_2055'
    cmb <- paste0(ssp,'_',gcm, '_', prd)
    cat('To process -----> ', cmb, '\n')
    mns <- sprintf('%02.0f',1:12)
    stp <- base::expand.grid(yrs, mns, stringsAsFactors = F) |> base::as.data.frame(); rm(mns)
    names(stp) <- c('yrs','mns')
    stp <- stp |> dplyr::arrange(yrs, mns) |> base::as.data.frame()
    
    
    pr_pth <- paste0(root,'/chirps_cmip6_world/Prec_',gcm,'_',ssp,'_',prd) # Precipitation
    tm_pth <- paste0(root, '/chirts_cmip6_world/Tmin_', gcm, '_', ssp, '_', prd)# Minimum temperature
    tx_pth <- paste0(root, '/chirts_cmip6_world/Tmax_', gcm, '_', ssp, '_', prd) # Maximum temperature
    sr_pth <- paste0(root, '/agera5_cmip6_africa/solar_radiation_flux_', gcm, '_', ssp, '_', prd) # Solar radiation
    out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/NDWL0') # Output directory
    
    1:nrow(stp) |>
      purrr::map(.f = function(i){
        calc_ndwl0(yr = stp$yrs[i], mn = stp$mns[i]); gc(F, T, T)
        if (i%%5 == 0) {
          tmpfls <- list.files(tempdir(), full.names = T)
          1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ',tmpfls[k]))})
        }
      })
    
  }
  
}

r <- rast('./common_data/atlas_hazards/cmip6/indices/ssp370_MRI-ESM2-0_2025_2055/NDWL0/NDWL0-2055-10.tif')
