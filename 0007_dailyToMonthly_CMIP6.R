


# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, geodata, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
root <- './common_data'
dirs <- dir_ls(root) %>% as.character()
gcms <- c('ACCESS-CM2', 'EC-Earth3', 'INM-CM5-0', 'MPI-ESM1-2-HR')

# Function ----------------------------------------------------------------
day2mnt <- function(dire){
  
  dire <- glue('{dirs[2]}')
  
  ## To list the files
  cat('To starth the process!\n')
  fles <- dir_ls(dire) %>% as.character()
  seqn <- seq(
    from = as.Date('2025-01-01'), 
    to = as.Date('2055-12-31'), 
    by = 'day'
  )
  varb <- basename(dire) %>% str_split('_') %>% map_chr(1)
  
  ## Year / Month
  mnts <- expand.grid(yr = 2025:2055, mn = c(paste0('0', 1:9), 10:12)) %>% 
    mutate(yearmonth = paste0(yr, '-', mn)) %>% 
    pull(3)
  
  if(varb == 'Prec'){
    mnts <- gsub('-', '\\.', mnts)
  }
  
  ## To aggregate
  rstr <- map(.x = mnts, .f = function(i){
    
    ## Grepping the files
    cat(i, '\t')
    fls <- grep(i, fles, value = T)
    
    ## To aggregate
    if(varb == 'Prec'){
      rst <- rast(fls)
      rst <- sum(rst)
      names(rst) <- glue('{varb}_{i}')
    }
    
    ## Finish 
    cat('Done!\n')
    return(rst)
    
  })
  
  ## Reduce yearly 
  year <- unique(str_sub(mnts, 1, 4))
  
  
  
  
  
  
}



