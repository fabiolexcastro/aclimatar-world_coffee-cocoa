

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, outliers, readxl, hablar, glue, readxl, geodata, RColorBrewer)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Vector -----------------------------
wrld <- geodata::world(resolution = 1, path = './tmpr')
vtnm <- wrld[wrld$GID_0 == 'VNM',]
ethi <- wrld[wrld$GID_0 == 'ETH',]
plot(wrld)

## Tabular data --------------------
cffe <- read_csv('./tbl/points/coffee_v3_2022.csv')
spcs <- sort(unique(cffe$Species))
oth1 <- read_excel('./tbl/points/Coordenadas Anserma y Ecuador.xlsx', 1) %>% dplyr::select(Lon, Lat) %>% mutate(Specie = 'Arabica')
oth2 <- read_excel('./tbl/points/Coordenadas Anserma y Ecuador.xlsx', 2) %>% dplyr::select(Lon, Lat) %>% mutate(Specie = 'Arabica')

rbst <- filter(cffe, Species %in% c('robusta', 'Robusta', 'Coffee Robusta'))
arab <- filter(cffe, !Species %in% c('robusta', 'Robusta', 'Coffee Robusta'))

### Tidy the tables
rbst <- dplyr::select(rbst, Lon = Longitude, Lat = Latitude, Specie = Species)
arab <- dplyr::select(arab, Lon = Longitude, Lat = Latitude, Specie = Species)
arab <- rbind(arab, oth1, oth2)

# Check Vietnam -----------------------------------------------------------
cffe.vnm <- cffe %>% filter(Country == 'Vietnam') %>% group_by(Country, Species, Source, Year) %>% reframe(count = n()) 
cffe.vnm <- cffe %>% filter(Country == 'Vietnam')
cffe.vnm <- cffe.vnm %>% mutate(Species = ifelse(Species %in% c('Coffee Robusta', 'Robusta'), 'Robusta', Species))
cffe.vnm <- as_tibble(cbind(cffe.vnm, terra::extract(vtnm, cffe.vnm[,c('Longitude', 'Latitude')])))
cffe.vnm <- drop_na(cffe.vnm, GID_0)
g.vnm <- ggplot() + geom_sf(data = st_as_sf(vtnm), fill = NA, col = 'grey30') + geom_point(data = cffe.vnm, aes(x = Longitude, y = Latitude, col = Species)) + coord_sf() + theme_bw() + theme(axis.text = element_text(size = 6)) +  guides(color = guide_legend(override.aes = list(size = 5)))
g.vnm
ggsave(plot = g.vnm, filename = './png/maps/points/coffee-points_arabica-robusta_vnm.jpg', units = 'in')

# Add Perú ----------------------------------------------------------------
pnts.peru <- read_csv('./tbl/points/Coffee_data_peru.csv', show_col_types = FALSE) %>% mutate(Specie = 'Arabica') %>% dplyr::select(Specie, Lon, Lat)
arab <- rbind(arab, pnts.peru)

# Add Ethiopia ------------------------------------------------------------
pnts.ethi <- read_csv('//catalogue/workspace-cluster9/2024/UgaEthTha/tble/points_eth_v2.csv', show_col_types = FALSE)
g.eth <- ggplot() + geom_sf(data = st_as_sf(ethi), fill = NA, col = 'grey30') + geom_point(data = pnts.ethi, aes(x = Longitude, y = Latitude, col = Species)) + coord_sf() + theme_bw() + theme(axis.text = element_text(size = 6)) +  guides(color = guide_legend(override.aes = list(size = 5)))
pnts.ethi <- dplyr::select(pnts.ethi, Lon = Longitude, Lat = Latitude, Specie = Species)
arab <- rbind(arab, pnts.ethi)
unique(arab$Specie)

## Raster data -----------------------
bioc <- terra::rast('./common_data/input_bios/bioc_hist.tif')
mskr <- bioc[[1]] * 0 + 1

indx <- terra::rast('./common_data/atlas_hazards/hist_indices2.tif')
indx <- indx[[grep('hist', names(indx))]]
bioc <- c(indx, bioc)

# Mask --------------------------------------------------------------------
mskr <- rast(ext(bioc[[1]]), crs = 'EPSG:4326', res = 0.45)
mskr[] <- 1; mskr <- terra::mask(mskr, wrld)

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
write.csv(cffe, './tbl/points/clean/v2/arabica_robusta.csv', row.names = FALSE)

# Sample with data --------------------------------------------------------
cffe <- read_csv('./tbl/points/clean/v2/arabica_robusta.csv')
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
write.csv(cffe, './tbl/points/clean/v2/arabica_robusta_otl.csv', row.names = FALSE)

# Sample with data again --------------------------------------------------
cffe <- read_csv('./tbl/points/clean/v2/arabica_robusta_otl.csv')
cffe <- as_tibble(cbind(cffe, terra::extract(bioc, cffe[,c('Lon', 'Lat')])))
write.csv(cffe, './tbl/points/clean/v2/arabica_robusta_otl_swd.csv', row.names = FALSE)


# Mapping -----------------------------------------------------------------

cffe
g.pnt <- ggplot() + 
  geom_point(data = cffe, aes(x = Lon, y = Lat, col = Specie), size = 0.5) +
  geom_sf(data = wrld, fill = NA, col = 'grey40') +
  coord_sf(ylim = c(-30, 30), xlim = c(-110, 150)) +
  theme_bw() +
  labs(caption = 'Arabica: 1034\nRobusta: 259') +
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6, angle = 90),
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(plot = g.pnt, filename = './png/maps/points/v2/point_ara-rob.jpg', units = 'in', width = 9, height = 4, dpi = 300, create.dir = T)

# Add new points and check Vietnam -------------------------------------------
# cffe <- read_csv('./tbl/points/clean/arabica_robusta_otl_swd.csv')
# read_csv('./tbl/points/Coffee_data_peru.csv', show_col_types = FALSE)




