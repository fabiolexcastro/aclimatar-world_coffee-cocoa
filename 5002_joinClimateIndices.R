
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, glue, geodata, ggspatial, hrbrthemes, RColorBrewer)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Raster
fles <- dir_ls('./common_data/stack_final', regexp = '.tif$')
clma <- rast(grep('stack_climate.tif', fles, value = T))
# indx <- rast(grep('indices', fles, value = T)[2:3])

## ISOS
isos <- c("BRA","VNM","IDN","COL","ETH","HND","UGA","PER","IND","CAF","GTM","GIN","MEX","LAO","NIC","CHN","CHN","CIV","CRI","TZA","COD","VEN","MDG","KEN","PNG","SLV","YEM","PHL","RWA","CMR")

isos <- c("BRA","VNM","IDN","COL","ETH","HND","UGA","PER","IND","GTM","MEX",
          "NIC","CHN","CHN","CIV","CRI","TZA",,"VEN","MDG","KEN","SLV","YEM"
          ,"PHL","RWA","CMR")

## World
wrld <- world(resolution = 1, level = 0, path = './tmpr')
zone <- wrld[wrld$GID_0 %in% isos,]

# Mosaicking --------------------------------------------------------------
indx <- as.character(grep('indices', fles, value = T))
indx <- indx[-grep('ara-rob', (indx), value = F)]

indx.1 <- indx[-grep('oth', indx, value = F)]
indx.2 <- indx[grep('oth', indx, value = F)]
indx.1 <- rast(indx.1)
indx.2 <- rast(indx.2)

names(indx.1)
names(indx.2) <- c('ndd_ftre', 'n30_ftre', 'n35_ftre', 'tai_ftre', 'ndw_ftre', 'hsh_ftre', 'ndd_hist', 'n30_hist', 'n35_hist', 'hsh_hist', 'tai_hist', 'ndw_hist')

indx <- mosaic(indx.1, indx.2)
terra::writeRaster(x = indx, filename = './common_data/stack_final/stack_indices_hist-ftre.tif', overwrite = TRUE)

# Check duplicated --------------------------------------------------------
dupv <- duplicated(names(clma))
table(dupv)

# Mask for the index ------------------------------------------------------
mskr <- indx[[1]] * 0 + 1
pols <- as.polygons(mskr)

# To extract by mask ------------------------------------------------------
clma <- terra::crop(clma, pols)
clma <- terra::mask(clma, pols)

indx <- terra::crop(indx, pols)
indx <- terra::mask(indx, pols)

# To make the stack  ------------------------------------------------------
stck <- c(clma, indx)
terra::writeRaster(x = stck, filename = './common_data/stack_final/stack_climate-indices_coffee_ara-rob_countries.tif', overwrite = TRUE)

stck <- terra::rast('./common_data/stack_final/stack_climate-indices_coffee_ara-rob_countries.tif')
rraw <- terra::rast('./common_data/stack_final/stack_indices-hist.tif')

# Raster to table ---------------------------------------------------------
map(.x = 1:length(isos), .f = function(i){
  
  try(
    expr = {
      
      ## Filtering the adminsitrative zone
      cat('To process: ', isos[i], '\n')
      dir <- './common_data/table_climate-indices/byCountry'
      out <- glue('{dir}/{isos[i]}_climate-indices_coffee.csv')
      
      if(!file.exists(out)){
        
        iso <- isos[i]
        zne <- wrld[wrld$GID_0 == iso,]
        
        ## To extract by mask 
        stk <- terra::crop(stck, zne)
        stk <- terra::mask(stk,  zne)  
        
        ## Raster to tabhle 
        tbl <- terra::as.data.frame(stk, xy = T)
        
        ## To remove the NAs
        tbl.2 <- drop_na(tbl)
        
        ## Table to raster
        rst <- terra::rast(tbl.2, type = 'xyz', crs = 'EPSG:4326')
        
        ## Raster to table again
        tbl <- terra::as.data.frame(tbl, xy = T)
        tbl <- mutate(tbl, country = iso, .before = 'x')
        head(tbl)
        
        ## To write the table 
        write.csv(tbl, glue('{dir}/{iso}_climate-indices_coffee.csv'), row.names = FALSE)
        cat('Done!\n')
      }
      
      
    }
    
  )
     
})


# To make some maps -------------------------------------------------------
fles <- as.character(dir_ls('./common_data/table_climate-indices/byCountry'))
isos <- str_split(basename(fles), '_') %>% map_chr(1) %>% unique()
iso  <- isos[1]

make.maps <- function(iso){
  
  ## To read the table
  cat('To process: ', iso, '\n')
  tble <- read_csv(grep(iso, fles, value = T), show_col_types = FALSE)
  
  ## Tidy the table
  tble <- tble %>% gather(var, value, -c(country, x, y))
  vars <- unique(tble$var)
  indx <- c('ndd', 'n30', 'n35', 'tai', 'ndw', 'hsh')
  tble.indx <- tble[grep(paste0(indx, collapse = '|'), tble$var, value = F),]
  tble.clma <- tble[-grep(paste0(indx, collapse = '|'), tble$var, value = F),]
  
  ## Country 
  shpe <- geodata::gadm(country = iso, level = 1, path = './tmpr')
  shpe <- st_as_sf(shpe)
  
  ## Colors
  clrs <- rev(RColorBrewer::brewer.pal(n = 3, name = 'RdYlGn'))
  names(clrs) <- c('Low', 'Medium', 'High')
  
  ## Index 
  tble.indx <- tble.indx %>% separate(data = ., col = 'var', into = c('index', 'period'), sep = '_')
  tble.indx <- tble.indx %>% mutate(value = ifelse(value == 1, 'Low', ifelse(value == 2, 'Medium', ifelse(value == 3, 'High', 'Very High'))))
  tble.indx <- tble.indx %>% mutate(value = factor(value, levels = c('Low', 'Medium', 'High', 'Very High')))
  tble.indx <- tble.indx %>% mutate(period = factor(period, levels = c('hist', 'ftre')))
  
  g.indx <- ggplot() +
    geom_tile(data = tble.indx, aes(x = x, y = y, fill = value)) +
    facet_wrap(.~index + period) +
    scale_fill_manual(values = clrs, name = '') +
    ggtitle(label = iso) +
    geom_sf(data= shpe, fill = NA, col = 'grey30') +
    coord_sf() +
    labs(x = 'Lon', y = 'Lat', fill = 'Class') +
    theme_bw() +
    theme(
      plot.title = element_text(face = 'bold', hjust = 0.5),
      legend.position = 'bottom'
    ) 
  
  
  ## Finish 
  cat('Done!\n')
  return(g.indx)
  
}

ggin <- map(isos, make.maps)

# Save the map as a PDF file ----------------------------------------------
pdf("maps_indices_todos2.pdf", width = 10, height = 10)
walk(ggin, print)
dev.off()
