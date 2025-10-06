


# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, hablar, glue, readxl, geodata, RColorBrewer)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)


# Load data ---------------------------------------------------------------

## Vector -----------------------------
wrld <- geodata::world(resolution = 1, path = './tmpr')

## Raster -----------------------------
mskr <- geodata::worldclim_global(var = 'prec', res = 2.5, path = './tmpr')
mskr <- mskr[[1]]
mskr <- mskr * 0

## Points -----------------------------
fles <- as.character(dir_ls('./tbl/points', regexp = '.xlsx$'))
fles.txt <- dir_ls('./tbl/points', regexp = '.txt$') 

# To process the points ---------------------------------------------------

### Table 1
pnt1 <- read_excel("./tbl/points/CDI_baseline_GPS_cleaned.xlsx")
pnt1 <- dplyr::select(pnt1, -ID, -starts_with('Alt'), -starts_with('Acc'))
pnt1 <- mutate(pnt1, across(everything(), as.numeric))
pnt1 <- map(.x = 1:8, .f = function(i){
  pn <- pnt1 %>% dplyr::select(matches(as.character(i), '$'))
  colnames(pn) <- gsub(i, '', colnames(pn))
  pn <- drop_na(pn)
  return(pn)
})
pnt1 <- bind_rows(pnt1)
pnt1 <- mutate(pnt1, country = terra::extract(wrld, pnt1[,c('Lon', 'Lat')])$GID_0)
pnt1 <- as_tibble(pnt1)
pnt1 <- dplyr::select(pnt1, country, Lon, Lat)

### Table 2
pnt2 <- read_excel('./tbl/points/CDI_baseline_GPS_Point_cleaned.xlsx')
pnt2 <- dplyr::select(pnt2, Lon = Lon1, Lat = Lat1)
pnt2 <- mutate(pnt2, across(everything(), as.numeric))
pnt2 <- drop_na(pnt2)
pnt2 <- mutate(pnt2, country = terra::extract(wrld, pnt2[,c('Lon', 'Lat')])$GID_0)
pnt2 <- dplyr::select(pnt2, country, Lon, Lat)

### Table 3
pnt3 <- read_excel("./tbl/points/CocoaMapping - all versions - labels - 2020-03-18-17-12-58.xlsx") 
pnt3 <- pnt3 %>% dplyr::select(matches('atitude'), matches('ongitude')) %>% dplyr::select(matches('center')) %>% mutate(across(everything(), as.numeric))
pnt3 <- pnt3 %>% setNames(c('Lat', 'Lon')) %>% mutate(country = terra::extract(wrld, .[,c('Lon', 'Lat')])$GID_0)
pnt3 <- dplyr::select(pnt3, Lon, Lat, country)

### Table 4
pnt4 <- read_excel("./tbl/points/CocoaSoils_Baseline_Ghana_cleaned.xlsx")
pnt4 <- dplyr::select(pnt4, -ID, -starts_with('Alt'), -starts_with('Acc'))
pnt4 <- mutate(pnt4, across(everything(), as.numeric))
pnt4 <- map(.x = 1:8, .f = function(i){
  pn <- pnt4 %>% dplyr::select(matches(as.character(i), '$'))
  colnames(pn) <- gsub(i, '', colnames(pn))
  pn <- drop_na(pn)
  return(pn)
})
pnt4 <- bind_rows(pnt4)
pnt4 <- drop_na(pnt4)
pnt4 <- mutate(pnt4, country = terra::extract(wrld, pnt4[,c('Long', 'Lat')])$GID_0)
pnt4 <- drop_na(pnt4)
pnt4 <- as_tibble(pnt4)
pnt4 <- dplyr::select(pnt4, Lon = Long, Lat, country)

### Table 5
pnt5 <- read_excel("./tbl/points/NIGERIA_GEOTRACE_PLOT_LIST_cleaned.xlsx")
pnt5 <- pnt5 %>% dplyr::select(starts_with('Lon'), starts_with('Lat'))
pnt5 <- mutate(pnt5, across(everything(), as.numeric))
pnt5 <- map(.x = 1:8, .f = function(i){
  pn <- pnt5 %>% dplyr::select(matches(as.character(i), '$'))
  colnames(pn) <- gsub(i, '', colnames(pn))
  pn <- drop_na(pn)
  return(pn)
})
pnt5 <- bind_rows(pnt5)
pnt5 <- mutate(pnt5, country = terra::extract(wrld, pnt5[,c('Long', 'Lat')])$GID_0)
pnt5 <- dplyr::select(pnt5, country, Lon = Long, Lat)

### TXT file 
tble <- read.table(fles.txt, sep = ';')
tble <- as_tibble(tble)
tble <- dplyr::select(tble, Lon = coords.x1, Lat = coords.x2)
tble <- mutate(tble, country = terra::extract(wrld, tble[,c(1,2)])$GID_0)
tble <- dplyr::select(tble, country, Lon, Lat)

### WA Cocoa GPS
shps <- dir_ls('./tbl/points/WA_Cocoa_GPS') %>% map(dir_ls) %>% unlist() %>% as.character() %>% grep('.shp$', ., value = T)
shps <- map(shps, st_read)
crds <- map(1:length(shps), function(i){
  shps[[i]] %>% 
    st_coordinates() %>% 
    as_tibble() %>% 
    setNames(c('Lon', 'Lat'))
})
crds <- bind_rows(crds)
crds <- mutate(crds, country = terra::extract(wrld, crds[,1:2])$GID_0)
crds <- dplyr::select(crds, country, Lon, Lat)

### Join all the points into only one
pnts <- bind_rows(list(pnt1, pnt2, pnt3, pnt4, pnt5, crds, tble))
write.csv(pnts, './tbl/points/points_eric.csv', row.names = FALSE)

# To remove the duplicated by cell ----------------------------------------
pnts <- as_tibble(cbind(pnts, cell = terra::extract(mskr, pnts[,c('Lon', 'Lat')], cell = T)$cell))
pnts <- drop_na(pnts)
clls <- unique(pnts$cell)

pnts.clean <- map(.x = 1:length(clls), .f = function(i){
  pnt <- filter(pnts, cell == clls[i])
  pnt <- slice(pnt, 1)
  return(pnt)
})
pnts.clean <- bind_rows(pnts.clean)
pnts.clean <- filter(pnts.clean, country != 'PNG')
pnts.clean <- mutate(pnts.clean, croop = 'Cocoa')
write.csv(pnts.clean, './tbl/points/points_eric-removeDupCell.csv', row.names = FALSE)

# Draw the map ------------------------------------------------------------

cnts <- unique(pnts.clean$country)
wafr <- wrld[wrld$GID_0 %in% c('GHA', 'CIV', 'CMR'),]
idns <- wrld[wrld$GID_0 %in% c('IDN'),]

g.waf <- ggplot() +
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey50') + 
  geom_point(data = pnts.clean, aes(x = Lon, y = Lat)) +
  geom_sf_text(data = st_as_sf(wrld), aes(label = NAME_0), col = 'grey30') +
  coord_sf(xlim = ext(wafr)[1:2], ylim = ext(wafr)[3:4]) +
  theme_bw() +
  ggtitle(label = 'Cocoa points at: CIV, GHA, NGA, CMR') +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5)
  )

g.waf

g.idn <- ggplot() +
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey50') + 
  geom_point(data = pnts.clean, aes(x = Lon, y = Lat)) +
  geom_sf_text(data = st_as_sf(wrld), aes(label = NAME_0), col = 'grey30') +
  coord_sf(xlim = ext(idns)[1:2], ylim = ext(idns)[3:4]) +
  theme_bw() +
  ggtitle(label = 'Cocoa points at: IDN') +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5)
  )

g.idn

## To save the map
ggsave(plot = g.waf, filename = './png/maps/points/westafrica_points_eric.jpg', units = 'in', width = 9, height = 5, dpi = 300, create.dir = T)
ggsave(plot = g.idn, filename = './png/maps/points/idn-png_points_eric.jpg', units = 'in', width = 9, height = 4, dpi = 300, create.dir = T)


## CIV Farm points 
shp1 <- dir_ls('./tbl/points/CIV_CRAFT_farm_points', regexp = '.shp$')
shp1 <- st_read(shp1)
shp1 <- st_transform(shp1, st_crs(4326))

### Get the coordinates as a table
shp1 <- shp1 %>% st_coordinates() %>% as_tibble() %>% mutate(country = 'CIV')
shp1 <- dplyr::select(shp1, country, Lon = X, Lat = Y)

## GHA Farm points
shp2 <- dir_ls('./tbl/points/GHA_Farm_Points', regexp = '.shp$')
shp2 <- st_read(shp2)
shp2 <- st_transform(shp2, st_crs(4326))

### Get the coordinates as a table 
shp2 <- shp2 %>% st_coordinates() %>% as_tibble() %>% mutate(country = 'GHA') 
shp2 <- dplyr::select(shp2, country, Lon = X, Lat = Y)

## Join both tables into only one
shps <- rbind(shp1, shp2)
gha.civ <- wafr[wafr$GID_0 %in% c('CIV', 'GHA'),]

g.gha.civ <- ggplot() +
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey50') + 
  geom_point(data = shps, aes(x = Lon, y = Lat)) +
  geom_sf_text(data = st_as_sf(wrld), aes(label = NAME_0), col = 'grey30') +
  coord_sf(xlim = ext(gha.civ)[1:2], ylim = ext(gha.civ)[3:4]) +
  theme_bw() +
  ggtitle(label = 'Cocoa points at: CIV - GHA') +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5)
  )

g.gha.civ
ggsave(plot = g.gha.civ, filename = './png/maps/points/gha-civ_points_eric.jpg', units = 'in', width = 9, height = 7, dpi = 300, create.dir = T)


## Add to the main table
shps <- mutate(shps, croop = 'Cocoa', cell = NA)
allp <- bind_rows(pnts.clean, shps)
allp <- mutate(allp, cell = terra::extract(mskr, allp[,c('Lon', 'Lat')], cell = T)$cell)
clls <- unique(allp$cell)
allp <- map_dfr(.x = clls, .f = function(i){
  allp %>% filter(cell == i) %>% slice(1)
})
write.csv(allp, './tbl/points/points_eric-removeDupCell_v2.csv', row.names = FALSE)


