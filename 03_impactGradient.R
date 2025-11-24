
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, tidyterra, usdm, sp, raster, RColorBrewer, randomForest, outliers, readxl, hablar, glue, readxl, geodata, RColorBrewer)
pacman::p_load(tidyverse, raster, ggpubr, cptcity, corrplot, cclust, dismo, gtools, sp, FactoMineR, pROC, randomForest, Hmisc)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Raster
crnt <- raster('./workspace_coffee/rf/output/run_2/results/process/rf_mixed-bsl.tif')
ftre <- raster('./workspace_coffee/rf/output/run_2/results/process/rf-mix_ftr-gcms.tif')

## Tabular 
all_options <- read_csv('./tbl/classesImpGraLimMix_2clusters.csv', show_col_types = FALSE)
labelss <- tibble(value = 0:5, category = c('Unsuit', 'cope', 'adjust', 'transform', 'opportunity', 'resilience'))
all_options %>% View()

## Vector
wrld <- geodata::world(resolution = 1, path = './tmpr')

# Impact gradient ---------------------------------------------------------

## Function
make.impact <- function(crn, ftr){
  
  # crn <- crnt
  # ftr <- ftre
  
  cat('To start the analysis!\n')
  msk <- crn * 0
  crd_df <- coordinates(crn)
  
  x <- raster::extract(crn, crd_df, cellnumbers = TRUE) %>% as_data_frame()
  ncell <- dplyr::select(x, cells)
  x <- select_(x, names(crn))
  colnames(x) <- 'current'
  
  y <- raster::extract(ftr, crd_df[,c('x', 'y')], cellnumbers = TRUE) %>% as_data_frame()
  y <- select_(y, names(ftr))
  colnames(y) <- 'future'
  
  z <- data.frame(x, y, ncell) %>% as_tibble()
  
  print('To Results')
  rslts <- left_join(z, all_options, by = c('current', 'future'))
  labls <- as_tibble(labelss) %>% mutate(category = as.character(category))
  
  final <- full_join(rslts, labls, by = 'category') %>% dplyr::select(value) %>% pull(1)
  final <- left_join(rslts, labls, by = 'category') %>% dplyr::select(value) %>% pull(1)
  
  length(final)
  length(msk)
  hist(final)
  
  rst <- raster::setValues(msk, final)
  rst <- rast(rst)
  return(rst)
  
}

## To apply the function
impg <- make.impact(crn = crnt, ftr = ftre)
plot(impg)
terra::writeRaster(x = impg, filename = './workspace_coffee/rf/output/run_2/results/process/rf-impactGradient_mdl.tif', overwrite = TRUE)

## As a factor table
levels(impg)[[1]] <- data.frame(ID = 0:5, clase = c('Unsuitable', 'Incremental adaptation', 'Systemic adaptation', 'Transform', 'Opportunities', 'Systemic resilience'))

## Colors for the map
clrs <- c('white', '#80bc8d', '#e5d51c', '#c95855', '#5c8a4d', '#e2971a')
names(clrs) <- c('Unsuitable', 'Incremental adaptation', 'Systemic adaptation', 'Transform', 'Opportunities', 'Systemic resilience')

# To draw the map ---------------------------------------------------------

## Mapping
gmap <- ggplot() +
  geom_spatraster(data = impg, aes(fill = clase)) + 
  scale_fill_manual(values = clrs, na.value = 'transparent', na.translate = FALSE) +
  geom_sf(data = wrld, fill = NA, col = 'grey40') +
  coord_sf(xlim = c(-110, 150), ylim = c(-30, 30)) +
  labs(x = '', y = '', fill = '') +
  ggtitle(label = 'Impact gradient - Coffee arabica & robusta') +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, size = 14),
    legend.position = 'bottom',
    axis.text.x = element_text(size = 6), 
    axis.text.y = element_text(size = 6, hjust = 0.5, angle = 90)
  )

gmap
ggsave(plot = gmap, filename = './workspace_coffee/png/maps/run_2/rf_impactGradient.jpg', units = 'in', width = 7, height = 3.5, dpi = 300, create.dir = TRUE)


levels(impg)[[1]]
