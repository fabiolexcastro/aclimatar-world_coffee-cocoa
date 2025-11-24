

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
unique(arab$Specie) %>% mixedsort()

write.csv(arab, './tbl/points/all_points_coffee-arabica.csv', row.names = FALSE)
arab <- read_csv('./tbl/points/all_points_coffee-arabica.csv', show_col_types = FALSE)

# Check Ecuador -----------------------------------------------------------
cffe.ecu <- cffe %>% filter(Country == 'Ecuador') #%>% pull(Source) %>% table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% kableExtra::kable()
write.csv(cffe.ecu, './tbl/points/coffee_ecuador.csv', row.names = FALSE)
ecu1 <- geodata::gadm(country = 'ECU', level = 1, path = './tmpr') %>% st_as_sf()
g.ecu <- ggplot() + geom_sf(data = ecu1, fill = NA, col = 'grey30') + geom_point(data = cffe.ecu, aes(x = Longitude, y = Latitude, col = Source), size = 0.5) + coord_sf(xlim = c(-81, -75)) + theme_minimal() + theme()
g.ecu
ggsave(plot = g.ecu, filename = './workspace_coffee/png/mapa_ecu.jpg', units = 'in', width = 8, height = 7, dpi = 300, create.dir = T)
cffe.ecu <- st_as_sf(cffe.ecu, coords = c('Longitude', 'Latitude'), crs = st_crs(4326))
st_write(cffe.ecu, './workspace_coffee/shp/coffee_ecuador.kml')

cffe.ecu <- as_tibble(cbind(cffe.ecu, terra::extract(bioc, cffe.ecu[,c('Longitude', 'Latitude')])))
cffe.ecu <- dplyr::select(cffe.ecu, Species, Longitude, Latitude, Country, Source, all_of(vars))
cffe.ecu <- cffe.ecu %>% gather(var, value, -c(Species, Longitude, Latitude, Country, Source))
cffe.ecu <- cffe.ecu %>% mutate(var = factor(var, levels = vars))

g.box <- ggplot(data = cffe.ecu, aes(y = value)) + 
  geom_boxplot() +
  facet_wrap(.~var, scales = 'free_y') +
  labs(x = '', y = 'Value', title = 'Coffee arabica - Ecuador [bioclimatic variables]') +
  theme_minimal() +
  theme(strip.text = element_text(face = 'bold', hjust = 0.5), axis.text.x = element_blank())
ggsave(plot = g.box, filename = './png/graphs/boxplot_ecuador_arabica-raw.jpg', units = 'in', width = 8, height = 7, dpi = 300, create.dir = T)


# Add Robusta -------------------------------------------------------------
rbst.gbif <- read_csv('./tbl/points/coffee_robusta-canephora_rgbif.csv', show_col_types = FALSE)
rbst.gbif <- rbst.gbif %>% dplyr::select(Lon = decimalLongitude, Lat = decimalLatitude) %>% mutate(Specie = 'robusta')

## Join with the main table 
rbst <- rbind(rbst, rbst.gbif)

## Raster data -----------------------
bioc <- terra::rast('./common_data/input_bios/bioc_hist.tif')
mskr <- bioc[[1]] * 0 + 1

indx <- terra::rast('./common_data/atlas_hazards/hist_indices2.tif')
indx <- indx[[grep('hist', names(indx))]]
bioc <- c(indx, bioc)

# Mask --------------------------------------------------------------------
mskr <- rast(ext(bioc[[1]]), crs = 'EPSG:4326', res = 0.45)
# mskr <- rast(ext(bioc[[1]]), crs = 'EPSG:4326', res = 0.90)
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
nrow(rbst)
nrow(arab)

## 
rbst <- mutate(rbst, Specie = 'robusta')
arab <- mutate(arab, Specie = 'arabica')

## To write the table 
cffe <- rbind(arab, rbst)
cffe <- dplyr::select(cffe, -cell)
dir.create('./tbl/points/coffee')
write.csv(cffe, './tbl/points/coffee/arabica_robusta.csv', row.names = FALSE)
isos <- terra::extract(wrld, cffe[,c('Lon', 'Lat')])$GID_0

# Sample with data --------------------------------------------------------
cffe <- read_csv('./tbl/points/coffee/arabica_robusta.csv', show_col_types = FALSE)
cffe <- as_tibble(cbind(cffe, terra::extract(bioc, cffe[,c('Lon', 'Lat')])))
cffe <- dplyr::select(cffe, -ID)
cffe <- dplyr::select(cffe, -bioc_20_100)
cffe <- drop_na(cffe)


# Remove outliers ---------------------------------------------------------

## Function
rmvOutliers <- function(pnts){
  norm <- scores(pnts[,6:ncol(pnts)], 'z')
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
write.csv(cffe, './tbl/points/coffee/arabica_robusta_otl.csv', row.names = FALSE)

# Sample with data again --------------------------------------------------
cffe <- read_csv('./tbl/points/coffee/arabica_robusta_otl.csv', show_col_types = FALSE)
cffe <- as_tibble(cbind(cffe, terra::extract(bioc, cffe[,c('Lon', 'Lat')])))
write.csv(cffe, './tbl/points/coffee/arabica_robusta_otl_swd.csv', row.names = FALSE)

# Mapping -----------------------------------------------------------------
arab <- filter(cffe, Specie == 'arabica')
robu <- filter(cffe, Specie == 'robusta')

cffe
g.pnt <- ggplot() + 
  geom_point(data = cffe, aes(x = Lon, y = Lat, col = Specie), size = 0.5) +
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey40') +
  coord_sf(ylim = c(-30, 30), xlim = c(-110, 150)) +
  theme_bw() +
  labs(caption = glue::glue('Arabica: {nrow(arab)} Robusta: {nrow(robu)}'), col = '') +
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6, angle = 90),
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 5)))

g.pnt
ggsave(plot = g.pnt, filename = './png/maps/coffee/point_ara-rob.jpg', units = 'in', width = 9, height = 4, dpi = 300, create.dir = T)

# Add new points and check Vietnam -------------------------------------------
# cffe <- read_csv('./tbl/points/clean/arabica_robusta_otl_swd.csv')
# read_csv('./tbl/points/Coffee_data_peru.csv', show_col_types = FALSE)

# Add robusta points  -----------------------------------------------------


