


# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, glue, rgbif)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
spce <- 'Coffea canephora'
occr <- occ_search(scientificName = spce, limit = 200000, hasCoordinate = TRUE, hasGeospatialIssue = FALSE)
wrld <- geodata::world(resolution = 1, path = './tmpr')

## 
occr <- occr$data

# A simple map ------------------------------------------------------------
g1 <- ggplot() + 
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') + 
  geom_point(data = occr, aes(x = decimalLongitude, y = decimalLatitude)) + coord_sf() + theme_bw() + theme()

# Remove duplicated by each cell ------------------------------------------
mskr <- rast(ext(wrld), crs = 'EPSG:4326', res = 0.45)
mskr[] <- 1
mskr <- terra::mask(mskr, wrld)
clls <- terra::extract(mskr, occr[,c('decimalLongitude', 'decimalLatitude')], cell = T)$cell
occr <- mutate(occr, cell = clls)
clls <- unique(clls)

occr.cln <- map(.x = 1:length(clls), .f = function(i){
  occ <- filter(occr, cell == clls[i])
  occ <- slice(occ, 1)
  return(occ)
})
occr.cln <- bind_rows(occr.cln)
write.csv(occr.cln, './tbl/points/coffee_robusta-canephora_rgbif.csv', row.names = FALSE)
occr.cln <- read_csv('./tbl/points/coffee_robusta-canephora_rgbif.csv', show_col_types = FALSE)

# A simple map ------------------------------------------------------------
g2 <- ggplot() + 
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') + 
  geom_point(data = occr.cln, aes(x = decimalLongitude, y = decimalLatitude), col = 'red', size = 0.5) + 
  coord_sf(ylim = c(-30, 30), xlim = c(-120, 140)) +
  labs(x = '', y = '') + theme_bw() + theme(axis.text.y = element_text(angle = 90, hjust = 0.5))

ggsave(plot = g2, filename = './png/maps/points/c canephora robusta rgbif.jpg', units = 'in', width = 10, height = 3, dpi = 300)
nrow(occr.cln)

# Read model results ------------------------------------------------------
crnt <- terra::rast('./rf/output/run_4/results/process/rf-mixed_bsl.tif')
vles <- terra::extract(crnt, occr.cln[,c('decimalLongitude', 'decimalLatitude')])
vles <- pull(vles, 2)
frqn <- table(vles) %>% as.data.frame() %>% rownames_to_column()
frqn <- frqn %>% arrange(desc(Freq))
frqn <- mutate(frqn, vles = factor(vles, levels = vles))

clrs <- c('#FFFFFF',  '#3E5C23', '#C6E01B',  '#D9D9D9', '#FFFFB4')
names(clrs) <- c('Unsuitable', 'Robusta', 'Arabica', 'Limitations', 'Mixed')

gfrqn <- ggplot(data = frqn, aes(x = vles, y = Freq, fill = vles)) + 
  geom_col()  + 
  scale_fill_manual(values = clrs) + 
  labs(fill = '', x = '', y = 'Frequence') +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5)
  )

ggsave(plot = gfrqn, filename = './png/graphs/v2/validationRobusta_RGBIF.jpg', units = 'in', width = 9, height = 7, dpi = 300)

kableExtra::kable(frqn)
frqn %>% mutate(porc = Freq / sum(Freq) * 100, porc = round(porc, digits = 0))
