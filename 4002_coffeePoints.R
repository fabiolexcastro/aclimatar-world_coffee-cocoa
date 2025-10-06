

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, outliers, readxl, hablar, glue, readxl, geodata, RColorBrewer)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Vector -----------------------------
wrld <- geodata::world(resolution = 1, path = './tmpr')
plot(wrld)

## Tabular data --------------------
cffe <- read_csv('./tbl/points/coffee_v3_2022.csv')
spcs <- sort(unique(cffe$Species))
oth1 <- read_excel('./tbl/points/Coordenadas Anserma y Ecuador.xlsx', 1) %>% dplyr::select(Lon, Lat) %>% mutate(Specie = 'Arabica')
oth2 <- read_excel('./tbl/points/Coordenadas Anserma y Ecuador.xlsx', 2) %>% dplyr::select(Lon, Lat) %>% mutate(Specie = 'Arabica')

rbst <- filter(cffe, Species %in% c('robusta', 'Robusta'))
arab <- filter(cffe, !Species %in% c('robusta', 'Robusta'))

### Tidy the tables
rbst <- dplyr::select(rbst, Lon = Longitude, Lat = Latitude, Specie = Species)
arab <- dplyr::select(arab, Lon = Longitude, Lat = Latitude, Specie = Species)
arab <- rbind(arab, oth1, oth2)

## Raster data -----------------------
bioc <- terra::rast('./common_data/input_bios/bioc_hist.tif')
mskr <- bioc[[1]] * 0 + 1

# Remove duplicated by cell -----------------------------------------------
rmve.dupv <- function(tble){
  
  cat('To start the process!\n')
  clls <- terra::extract(mskr, tble[,c('Lon', 'Lat')], cell = T)$cell
  tble <- mutate(tble, cell = clls)
  clls <- unique(clls)
  fnal <- map(.x = clls, .f = function(i){tble %>% filter(cell == i) %>% slice(1)})
  fnal <- bind_rows(fnal)
  return(fnal)
  
}
rbst <- rmve.dupv(tble = rbst)
arab <- rmve.dupv(tble = arab)

## 
rbst <- mutate(rbst, Specie = 'robusta')
arab <- mutate(arab, Specie = 'arabica')

## To write the table 
cffe <- rbind(arab, rbst)
cffe <- dplyr::select(cffe, -cell)
write.csv(cffe, './tbl/points/clean/arabica_robusta.csv', row.names = FALSE)

# Sample with data --------------------------------------------------------
cffe <- read_csv('./tbl/points/clean/arabica_robusta.csv')
cffe <- as_tibble(cbind(cffe, terra::extract(bioc, cffe[,c('Lon', 'Lat')])))
cffe <- dplyr::select(cffe, -ID)
cffe <- dplyr::select(cffe, -bioc_20_100)
cffe <- drop_na(cffe)

# Remove outliers ---------------------------------------------------------

## Function
rmvOutliers <- function(pnts){
  norm <- scores(pnts[,4:ncol(pnts)], 'z')
  norm_na <- norm
  norm_na[abs(norm_na) > 3.5] <- NA
  normpoints <- cbind(pnts[,c('Lon', 'Lat')], norm_na) %>% 
    na.omit() %>% 
    as_data_frame()
  print('Done...!')
  normpoints <- normpoints[,c('Lon', 'Lat')]
  return(normpoints)
}
arab <- rmvOutliers(pnts = filter(cffe, Specie == 'arabica')) # 6828 --> 6413
robu <- rmvOutliers(pnts = filter(cffe, Specie == 'robusta')) # 645 --> 604
arab <- mutate(arab, Specie = 'arabica')
robu <- mutate(robu, Specie = 'robusta')

## Join both tables into only one
cffe <- rbind(arab, robu)
write.csv(cffe, './tbl/points/clean/arabica_robusta_otl.csv', row.names = FALSE)

# Sample with data again --------------------------------------------------
cffe <- read_csv('./tbl/points/clean/arabica_robusta_otl.csv')
cffe <- as_tibble(cbind(cffe, terra::extract(bioc, cffe[,c('Lon', 'Lat')])))
write.csv(cffe, './tbl/points/clean/arabica_robusta_otl_swd.csv', row.names = FALSE)

# Clustering using random forest ------------------------------------------



