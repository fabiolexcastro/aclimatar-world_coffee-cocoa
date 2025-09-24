## Get daily future data
## By: H. Achicanoy
## December, 2022

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,fs,lubridate,furrr))

root <- './common_data' # root <- '/home/jovyan/common_data'

# Africa reference raster
ref <- terra::rast(paste0('//catalogue/WFP_ClimateRiskPr1/1.Data/chirps/chirps-v2.0.1981.01.01.tif'))
ext <- ext(rast('./common_data/esfg_metagrid/intermediate/interpolated_mthly_anomaly/ACCESS-ESM1-5_ssp370_r1i1p1f1_hurs_Africa_monthly_intp_anomaly_2025_2055_fixed.tif'))
ref <- terra::crop(ref, ext)

# Interpolated monthly anomalies directory
anm_pth <- paste0(root,'/esfg_metagrid/intermediate/interpolated_mthly_anomaly')

# Parameters 
gcms <- c('ACCESS-ESM1-5', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')
ssp  <- 'ssp370'
vars <- c('pr', 'hurs', 'rsds', 'tasmin', 'tasmax')
mpgd <- tibble(Baseline = 1984:2014, Future = 2025:2055)

## Parameters
gcm <- gcms[1]; ssp <- ssp; var <- 'pr'; prd <- '2025_2055'

# Read monthly deltas
get_daily_future_data <- function(gcm, ssp, var, prd){
  cat(paste0('processing ',var,'_',gcm,'_',ssp,'_',prd,'\n'))
  prd <- as.character(prd)
  file <- paste0(gcm,'_',ssp,'_r1i1p1f1_',var,'_Africa_monthly_intp_anomaly_',prd,'_fixed.tif')
  # Load deltas
  dlts <- terra::rast(paste0(anm_pth,'/',file))
  dlts <- dlts |> terra::resample(ref)
  dlts <- dlts |> terra::crop(terra::ext(ref)) |> terra::mask(ref)
  prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
  his_yrs <- 1984:2014
  fut_yrs <- prd_num[1]:prd_num[2]
  
  # Temporal mapping
  Baseline = seq(from = as.Date('1984-01-01'),
                 to   = as.Date('2014-12-31'),
                 by   = 'day')
  Future   = seq(from = as.Date(paste0(prd_num[1],'-01-01')),
                 to   = as.Date(paste0(prd_num[2],'-12-31')),
                 by   = 'day')
  # Remove feb-29 from all leap years (not coincidence between years)
  Baseline <- Baseline[!(format(Baseline,'%m') == '02' & format(Baseline,'%d') == '29'), drop = F]
  Future   <- Future[!(format(Future,'%m') == '02' & format(Future,'%d') == '29'), drop = F]
  mpg <- data.frame(Baseline, Future)
  mpg$year  <- lubridate::year(mpg$Baseline)
  mpg$year_fut <- lubridate::year(mpg$Future)
  mpg$month <- lubridate::month(mpg$Baseline)
  # mpg <- mpg |> dplyr::filter(month == as.numeric(mn)) # OFF
  
  if(var %in% c('pr','rsds')){
    # Paths
    if (var == 'pr') {
      his_pth <- paste0('//catalogue/WFP_ClimateRiskPr1/1.Data/Chirps')
      fut_pth <- paste0(root,'/chirps_cmip6_world/Prec_',gcm,'_',ssp,'_',prd); dir.create(fut_pth, F, T)
    } else if (var == 'rsds') {
      his_pth <- paste0('//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5/solar_radiation_flux')
      fut_pth <- paste0(root,'/agera5_cmip6_africa/solar_radiation_flux_',gcm,'_',ssp,'_',prd)
      dir.create(fut_pth, F, T)
    }
    #if(length(list.files(fut_pth)) < 7300){
    # File structure
    if (var == 'pr') {
      his_str <- paste0('chirps-v2.0.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
      fut_str <- paste0('chirps-v2.0.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
    } else if (var == 'rsds') {
      his_str <- paste0('Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=mpg$Baseline, fixed=T),'_final-v1.0.nc')
      fut_str <- paste0('Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=mpg$Future, fixed=T),'_final-v1.0.tif')
    }
    # Split by months
    his_lst <- split(his_str, mpg$month)
    fut_lst <- split(fut_str, mpg$month)
    1:length(his_lst) |>
      purrr::map(.f = function(j){
        delta <- dlts[[j]]
        his_daily <- his_lst[[j]]
        fut_daily <- fut_lst[[j]]
        future::plan(future::multisession, workers = 5)
        1:length(his_daily) |>
          purrr::map(.f = function(k){ # furrr::future_map
            outfile <- paste0(fut_pth,'/',fut_daily[k])
            r <- terra::rast(paste0(his_pth,'/',his_daily[k]))
            if (var == 'pr') {
              r <- r |> terra::crop(terra::ext(ref)) #|> terra::mask(ref)
              r <- terra::classify(r, rcl = cbind(-Inf, -9990, NA))
            } else if (var == 'rsds') {
              r <- r |> terra::crop(terra::ext(ref)) |> terra::resample(ref) #|> terra::mask(ref)
              r <- r * 1e-6
            }
            r <- r * (1 + delta)
            r <- terra::classify(r, cbind(-Inf, 0, 0))
            terra::writeRaster(r, outfile, overwrite = T)
          })
        # future::plan(future::sequential); gc(reset = T)
      })
    #}
  }
  if(var %in% c('tasmax','tasmin','hurs')){
    # Paths
    his_pth <- ifelse(var == 'tasmax',
                      paste0('//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmax'),
                      ifelse(var == 'tasmin',
                             paste0('//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/Tmin'), 
                             paste0('//catalogue/WFP_ClimateRiskPr1/1.Data/chirts_global/RHum')))
    fut_pth <- ifelse(var == 'tasmax',
                      paste0(root,'/chirts_cmip6_world/Tmax_',gcm,'_',ssp,'_',prd),
                      ifelse(var == 'tasmin',
                             paste0(root,'/chirts_cmip6_world/Tmin_',gcm,'_',ssp,'_',prd),
                             paste0(root,'/chirts_cmip6_world/RHum_',gcm,'_',ssp,'_',prd)))
    if(length(list.files(fut_pth)) < 11315){
      # File structure
      if(var == 'tasmax'){
        his_str <- paste0('Tmax.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('Tmax.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      } else if (var == 'tasmin'){
        his_str <- paste0('Tmin.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('Tmin.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      } else {
        his_str <- paste0('RH.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('RH.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      }
      yrs_str   <- mpg$year
      yrs_f_str <- mpg$year_fut
      # Split by months
      his_lst   <- split(his_str, mpg$month)
      fut_lst   <- split(fut_str, mpg$month)
      yrs_lst   <- split(yrs_str, mpg$month)
      yrs_f_lst <- split(yrs_f_str, mpg$month)
      1:length(his_lst) |>
        purrr::map(.f = function(j){
          delta <- dlts[[j]]
          his_daily   <- his_lst[[j]]
          fut_daily   <- fut_lst[[j]]
          yrs_daily   <- yrs_lst[[j]]
          yrs_f_daily <- yrs_f_lst[[j]]
          future::plan(future::multisession, workers = 5)
          1:length(his_daily) |>
            furrr::future_map(.f = function(k){
              outfile <- paste0(fut_pth,'/',yrs_f_daily[k],'/',fut_daily[k]); dir.create(dirname(outfile),F,T)
              if(!file.exists(outfile)){
                r <- terra::rast(paste0(his_pth,'/',yrs_daily[k],'/',his_daily[k]))
                r <- r |> terra::crop(terra::ext(ref)) #|> terra::mask(ref)
                r <- terra::classify(r, rcl = cbind(-Inf, -9990, NA))
                r <- r + delta
                #for hurs limit min to 0, max to 100
                if (var == 'hurs') {
                  r <- min(r, 100)
                  r <- max(r, 0)
                }
                terra::writeRaster(r, outfile)
              }
            })
          future::plan(future::sequential); gc(reset = T)
        })
    }
  }
  
  return(cat(paste0(var,'_',gcm,'_',ssp,'_',prd,' ready.\n')))
  
}

## Manually
get_daily_future_data(gcm = gcm, ssp = ssp, var = 'pr', prd = prd)
gc(verbose = F, full = T, reset = T)

## Loop
map(.x = 1:length(gcms), .f = function(g){
  get_daily_future_data(gcm = gcms[g], ssp = ssp, var = 'pr', prd = prd)
})
