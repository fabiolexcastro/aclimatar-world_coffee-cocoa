
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, outliers, readxl, hablar, glue, readxl, geodata, RColorBrewer)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Tabular data
pnts.1 <- read_csv('./tbl/points/points_eric-removeDupCell_v2.csv', show_col_types = FALSE)
pnts.2 <- read_csv('./tbl/points/occ_cacao.csv')
pnts.3 <- read_excel('./tbl/points/Sitios_Cacao_ClimaLoCa.xlsx')

## Vector data
wrld <- geodata::world(resolution = 1, path = './tmpr')

## Mask
bioc <- terra::rast('./common_data/stack_final/stack_climate-indices_coffee_ara-rob_countries.tif')
mskr <- rast(ext(bioc[[1]]), crs = 'EPSG:4326', res = 0.45)
mskr[] <- 1; mskr <- terra::mask(mskr, wrld)

# Join points into only one -----------------------------------------------

## Points 1 / 2 / 3
pnts.1 <- pnts.1 %>% dplyr::select(Lon, Lat, country)
pnts.2 <- pnts.2 %>% dplyr::select(Lon = Longitude, Lat = Latitude)
pnts.3 <- pnts.3 %>% dplyr::select(Lon = `Longitud (WGS84)`, Lat = `Latitud (WGS84)`) 

## Extract the country
pnts.1 <- mutate(pnts.1, country = terra::extract(wrld, pnts.1[,c('Lon', 'Lat')])[,2])
pnts.2 <- mutate(pnts.2, country = terra::extract(wrld, pnts.2[,c('Lon', 'Lat')])[,2])
pnts.3 <- mutate(pnts.3, country = 'COL')

## Join the points 
pnts <- rbind(pnts.1, pnts.2, pnts.3)
pnts
write.csv(pnts, './tbl/points/cocoa/points_cocoa_all.csv', row.names = FALSE)

# Remove duplicated by cell -----------------------------------------------
clls <- terra::extract(mskr, pnts[,c('Lon', 'Lat')], cell = T)$cell
unqc <- unique(clls)
pnts <- mutate(pnts, cell = clls)

pnts.cln <- map(.x = 1:length(unqc), .f = function(i){
  pnt <- filter(pnts, cell == unqc[i])
  pnt <- slice(pnt, 1)
  return(pnt)
})
pnts.cln <- bind_rows(pnts.cln)
write.csv(pnts, './tbl/points/cocoa/points_cocoa_all_rmv50.csv', row.names = FALSE)

# Draw the map ------------------------------------------------------------
gg <- ggplot() +
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') + 
  geom_point(data = pnts.cln, aes(x = Lon, y = Lat), size = 0.1, fill = NA, col = 'darkred') +
  coord_sf(xlim = c(-120, 140), ylim = c(-30, 30)) +
  labs(x = 'Lon', y = 'Lat') +
  ggtitle(label = 'Cocoa location - Remove dup 50 km') +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6), 
    axis.text.y = element_text(size = 6)
  )

ggsave(plot = gg, filename = './png/maps/cocoa/cocoa_location.jpg', units = 'in', width = 12, height = 5, dpi = 300, create.dir = T)


# Sample with data --------------------------------------------------------

## Raster
bioc <- terra::rast('./common_data/input_bios/bioc_hist.tif')
mskr <- bioc[[1]] * 0 + 1
indx <- terra::rast('./common_data/atlas_hazards/hist_indices2.tif')
indx <- indx[[grep('hist', names(indx))]]
bioc <- c(indx, bioc)

## To extract the values
pnts <- as_tibble(cbind(pnts.cln, terra::extract(bioc, pnts.cln[,c('Lon', 'Lat')])))
pnts <- dplyr::select(pnts, -cell)
pnts <- drop_na(pnts)

# Remove outliers ---------------------------------------------------------

## Function
rmvOutliers <- function(pnts){
  norm <- scores(pnts[,5:ncol(pnts)], 'z')
  norm_na <- norm
  norm_na[abs(norm_na) > 3.5] <- NA
  normpoints <- cbind(pnts[,c('Lon', 'Lat')], norm_na) %>% 
    na.omit() %>% 
    as_data_frame()
  print('Done...!')
  normpoints <- normpoints[,c('Lon', 'Lat')]
  return(normpoints)
}

## To apply the function 
pnts.otl <- rmvOutliers(pnts = pnts) # 1247 --> 1144
pnts.otl
write.csv(pnts.otl, './tbl/points/cocoa/points_cocoa_all_rmv50_NoOTL.csv', row.names = FALSE)
