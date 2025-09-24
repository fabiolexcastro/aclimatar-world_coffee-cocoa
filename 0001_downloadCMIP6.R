

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1, timeout = 36000)

# Functions to use --------------------------------------------------------

my.download <- function(gcme, urls){
  
  gcme <- 'MPI-ESM1-2-HR'
  urls <- urls
  
  # varb <- 'tasmin'
  # urls <- glue('http://esgf3.dkrz.de/thredds/dodsC/cmip6/ScenarioMIP/EC-Earth-Consortium/EC-Earth3/ssp370/r1i1p1f1/day/{varb}/gr/v20200310/{varb}_day_EC-Earth3_ssp370_r1i1p1f1_gr_{year}0101-{year}1231.nc') %>% as.character()
  
  
  dir.create(oute, '/', gcme)
  
  for(i in 1:length(urls)){
    
    try(
      expr = {
        
        ## Tidy output
        urle <- urls[[i]]
        # urle <- 'http://esgf3.dkrz.de/thredds/dodsC/cmip6/ScenarioMIP/EC-Earth-Consortium/EC-Earth3/ssp370/r1i1p1f1/day/tasmin/gr/v20200310/tasmin_day_EC-Earth3_ssp370_r1i1p1f1_gr_20250101-20251231.nc'
        cat(basename(urle), '\n')
        outp <- paste0(oute, '/', gcme, '/', basename(urle))
        urle_descarga <- gsub("dodsC", "fileServer", urle)
        
        ## To download
        if(!file.exists(outp)){
          cat('To download\n')
          download.file(url = urle_descarga, destfile = outp, mode = 'wb')
        }
        
        ## To read the results
        rstr <- terra::rast(outp)
        rm(rstr)
        Sys.sleep(1)
        
      }
    )
    
  }
  
  dir_ls(paste0(oute, '/', 'EC-Earth3'))
  
  
}
my.check <- function(gcme, nmes){
  
  # gcme <- 'ACCESS-ESM1-5'
  # nmes <- map_chr(urls.accs, basename)
  
  cat('To list the files\n')
  fles <- as.character(dir_ls(glue('./common_data/esfg_metagrid/raw_world/{gcme}')))
  
  cat('To make the plots\n')
  plts <- map(.x = 1:length(fles), .f = function(x){
    print(basename(fles[x]))
    rstr <- rast(fles[x])
    plot(rstr[[1]], main = basename(fles[x]))
    rm(rstr)
    cat('Done!\n')
  })
  cat('Done!\n')
  
}

# Parameters --------------------------------------------------------------

## GCMs / Periods
gcms <- c('ACCESS-CM-2', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0')
prds <- c('1991-2014', '2021-2055') 
ssps <- c('historical', 'ssp370')
vars <- c('pr', 'tasmin', 'tasmax', 'hurs', 'rsds')

## Output directories 
dir_create('./common_data/esfg_metagrid/raw_world')
oute <- './common_data/esfg_metagrid/raw_world'

# ACCESS-ESM1-5 -----------------------------------------------------------

## Historical / SSP 370
root.accs <- 'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/CSIRO/ACCESS-ESM1-5' # 'ssp370/r1i1p1f1/day'
urls.accs <- c(
  paste0(root.accs, '/ssp370/r1i1p1f1/day/', vars, '/', 'gn/v20191115/', vars, '_day_ACCESS-ESM1-5_ssp370_r1i1p1f1_gn_20150101-20641231.nc'), # SSP 370
  paste0(root.accs, '/historical/r1i1p1f1/day/', vars, '/', 'gn/v20191115/', vars, '_day_ACCESS-ESM1-5_historical_r1i1p1f1_gn_20150101-20641231.nc') # Historical
)

# Download
my.download(gcme = 'ACCESS-ESM1-5', urls = urls.accs)

## Check 
my.check(gcme = 'ACCESS-ESM1-5', nmes = basename(urls.accs))

# EC-Earth3 ---------------------------------------------------------------

### Historical / SSP 370
root.ecea <- 'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3'
vars <- c('pr', 'tasmin', 'tasmax', 'rsds', 'hurs')

#### Historical
urls.ecea <- list()
for(i in 1:length(vars)){
  
  urls.ecea[[i]] <- paste0(
    root.ecea, '/historical/r1i1p1f1/day/', vars[i], '/', 'gr/v20200310/', vars[i], '_day_EC-Earth3_historical_r1i1p1f1_gr_', 1995:2014, '0101-', 1995:2014, '1231.nc'
  ) 
  
}
urls.ecea <- unlist(urls.ecea)

#### SSP 370
urls.ecea <- list()
for(i in 1:length(vars)){
  urls.ecea[[i]] <- glue(
    'http://esgf3.dkrz.de/thredds/dodsC/cmip6/ScenarioMIP/EC-Earth-Consortium/EC-Earth3/ssp370/r1i1p1f1/day/{vars[[i]]}/gr/v20200310/{vars[[i]]}_day_EC-Earth3_ssp370_r1i1p1f1_gr_{2025:2055}0101-{2025:2055}1231.nc'
  ) 
}
urls.ecea <- unlist(urls.ecea)

# Download
my.download(gcme = 'EC-Earth3', urls = urls.ecea)

# Check 
my.check(gcme = 'EC-Earth3', nmes = basename(urls.ecea))

# INM-CM5-0 ---------------------------------------------------------------

### Historical /SSP 370
root.inmc <- 'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/INM/INM-CM5-0'
urls.inmc <- list()
for(i in 1:length(vars)){
  urls.inmc[[i]] <- c(
    paste0(root.inmc, '/historical/r1i1p1f1/day/', vars[i], '/', 'gr1/v20200310/', vars[i], '_day_INM-CM5-0_historical_r1i1p1f1_gr1_', c('1950', '2000'), '0101-', c('1999', '2014'), '1231.nc'),
    paste0('https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/INM/INM-CM5-0', '/ssp370/r1i1p1f1/day/', vars[i], '/', 'gr1/v20190618/', vars[i], '_day_INM-CM5-0_ssp370_r1i1p1f1_gr1_', c('2015'), '0101-', '2064', '1231.nc')
  )
}
urls.inmc <- unlist(urls.inmc)

# Download
urls.inmc <- grep('rsds', urls.inmc, value = T)
urls.inmc <- grep('historical', urls.inmc, value = T)
my.download(gcme = 'INM-CM5-0', urls = urls.inmc)

# Check
my.check(gcme = 'INM-CM5-0', nmes = basename(urls.inmc))


# MPI-ESM1-2-HR -----------------------------------------------------------
urls.prec <- c(
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/pr/gn/v20190710/pr_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_19950101-19991231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/pr/gn/v20190710/pr_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20000101-20041231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/pr/gn/v20190710/pr_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20050101-20091231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/pr/gn/v20190710/pr_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20100101-20141231.nc'
)

urls.tmin <- c(
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/tasmin/gn/v20190710/tasmin_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_19950101-19991231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/tasmin/gn/v20190710/tasmin_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20000101-20041231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/tasmin/gn/v20190710/tasmin_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20050101-20091231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/tasmin/gn/v20190710/tasmin_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20100101-20141231.nc'
)

urls.tmax <- c(
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/tasmax/gn/v20190710/tasmax_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_19950101-19991231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/tasmax/gn/v20190710/tasmax_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20000101-20041231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/tasmax/gn/v20190710/tasmax_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20050101-20091231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/tasmax/gn/v20190710/tasmax_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20100101-20141231.nc'
)

urls.rsds <- c(
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/rsds/gn/v20190710/rsds_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_19950101-19991231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/rsds/gn/v20190710/rsds_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20000101-20041231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/rsds/gn/v20190710/rsds_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20050101-20091231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/rsds/gn/v20190710/rsds_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20100101-20141231.nc'
)

urls.hurs <- c(
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/hurs/gn/v20190710/hurs_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_19950101-19991231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/hurs/gn/v20190710/hurs_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20000101-20041231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/hurs/gn/v20190710/hurs_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20050101-20091231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/historical/r1i1p1f1/day/hurs/gn/v20190710/hurs_day_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_20100101-20141231.nc'
)



for(i in 1:length(urls.prec)){
  
  try(
    expr = {
      
      ## Tidy output
      urle <- urls.prec[[i]]
      outp <- paste0(oute, '/MPI-ESM1-2-HR/', basename(urle))
      urle_descarga <- gsub("dodsC", "fileServer", urle)
      
      ## To download
      download.file(url = urle_descarga, destfile = outp, mode = 'wb')
      
      ## To read the results
      rstr <- terra::rast(outp)
      rm(rstr)
      
    }
  )
  
}

## SSP 370
urls.prec <- glue('https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp370/r1i1p1f1/day/pr/gn/v20190710/pr_day_MPI-ESM1-2-HR_ssp370_r1i1p1f1_gn_{c(2025,2030,2035,2040,2045,2050,2055)}0101-{c(2029,2034,2039,2044,2049,2054,2059)}1231.nc')
urls.tmin <- glue('https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp370/r1i1p1f1/day/tasmin/gn/v20190710/tasmin_day_MPI-ESM1-2-HR_ssp370_r1i1p1f1_gn_{c(2025,2030,2035,2040,2045,2050,2055)}0101-{c(2029,2034,2039,2044,2049,2054,2059)}1231.nc')
urls.tmax <- glue('https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp370/r1i1p1f1/day/tasmax/gn/v20190710/tasmax_day_MPI-ESM1-2-HR_ssp370_r1i1p1f1_gn_{c(2025,2030,2035,2040,2045,2050,2055)}0101-{c(2029,2034,2039,2044,2049,2054,2059)}1231.nc')
urls.rsds <- glue('https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp370/r1i1p1f1/day/hurs/gn/v20190710/hurs_day_MPI-ESM1-2-HR_ssp370_r1i1p1f1_gn_{c(2025,2030,2035,2040,2045,2050,2055)}0101-{c(2029,2034,2039,2044,2049,2054,2059)}1231.nc')
urls.hurs <- glue('https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/DKRZ/MPI-ESM1-2-HR/ssp370/r1i1p1f1/day/rsds/gn/v20190710/rsds_day_MPI-ESM1-2-HR_ssp370_r1i1p1f1_gn_{c(2025,2030,2035,2040,2045,2050,2055)}0101-{c(2029,2034,2039,2044,2049,2054,2059)}1231.nc')

urls <- list(urls.prec, urls.tmin, urls.tmax, urls.rsds, urls.hurs)
urls <- unlist(urls)

# Download
my.download(gcme = 'INM-CM5-0', urls = as.character(urls.hurs))


# MRI-ESM2-0 --------------------------------------------------------------

## Historical 
urls.mrie <- list(
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/day/hurs/gn/v20190603/hurs_day_MRI-ESM2-0_historical_r1i1p1f1_gn_19500101-19991231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/day/hurs/v20190603/hurs_day_MRI-ESM2-0_historical_r1i1p1f1_gn_20000101-20141231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/day/pr/gn/v20190603/pr_day_MRI-ESM2-0_historical_r1i1p1f1_gn_19500101-19991231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/day/pr/gn/v20190603/pr_day_MRI-ESM2-0_historical_r1i1p1f1_gn_20000101-20141231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/day/tasmin/gn/v20190603/tasmin_day_MRI-ESM2-0_historical_r1i1p1f1_gn_19500101-19991231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/day/tasmin/gn/v20190603/tasmin_day_MRI-ESM2-0_historical_r1i1p1f1_gn_20000101-20141231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/day/tasmax/gn/v20190603/tasmax_day_MRI-ESM2-0_historical_r1i1p1f1_gn_19500101-19991231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/day/tasmax/gn/v20190603/tasmax_day_MRI-ESM2-0_historical_r1i1p1f1_gn_20000101-20141231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/day/rsds/gn/v20190603/rsds_day_MRI-ESM2-0_historical_r1i1p1f1_gn_19500101-19991231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/day/rsds/gn/v20190603/rsds_day_MRI-ESM2-0_historical_r1i1p1f1_gn_20000101-20141231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/MRI/MRI-ESM2-0/ssp370/r1i1p1f1/day/hurs/gn/v20190603/hurs_day_MRI-ESM2-0_ssp370_r1i1p1f1_gn_20150101-20641231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/MRI/MRI-ESM2-0/ssp370/r1i1p1f1/day/tasmin/gn/v20190603/tasmin_day_MRI-ESM2-0_ssp370_r1i1p1f1_gn_20150101-20641231.nc', 
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/MRI/MRI-ESM2-0/ssp370/r1i1p1f1/day/tasmax/gn/v20190603/tasmax_day_MRI-ESM2-0_ssp370_r1i1p1f1_gn_20150101-20641231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/MRI/MRI-ESM2-0/ssp370/r1i1p1f1/day/pr/gn/v20190603/pr_day_MRI-ESM2-0_ssp370_r1i1p1f1_gn_20150101-20641231.nc',
  'https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/ScenarioMIP/MRI/MRI-ESM2-0/ssp370/r1i1p1f1/day/rsds/gn/v20190603/rsds_day_MRI-ESM2-0_ssp370_r1i1p1f1_gn_20150101-20641231.nc'
  
)


