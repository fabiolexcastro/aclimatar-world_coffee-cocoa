
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(climateStability,cclust,corrplot,cptcity,dismo,FactoMineR,fs,ggpubr,glue,geodata,gtools,hablar,Hmisc,openxlsx,outliers,parallelDist,pROC,randomForest,raster,RColorBrewer,rnaturalearth,rnaturalearthdata,Rfast,scales,sf,sp,spatialEco,spatstat,terra,tidyverse,usdm)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ----------------------------------------------------------------

## Vector data
wrld <- geodata::world(resolution = 1, path = './tmpr')

## Models
mdls <- dir_ls('./workspace_coffee/rf/output/run_1/models', regexp = '.rdata$')
mdls <- mixedsort(mdls); for(i in 1:length(mdls)){load(mdls[i])}
mdls <- as.character(mdls)
mdl  <- list()

for(i in 1:length(mdls)){
  load(mdls[i])
  mdl[[i]] <- rfmodel 
}

rff <- do.call(randomForest::combine, mdl)

## Future data
clma <- rast('./common_data/input_bios/bioc_ftre-gcms.tif')
vars <- c('ndd_hist', 'n30_hist', 'n35_hist', 'ndw_hist', 'bioc_3', 'bioc_4', 'bioc_6',
          'bioc_13', 'bioc_14', 'bioc_18', 'bioc_19', 'bioc_20_60', 'bioc_22', 'bioc_23',
          'bioc_24', 'bioc_27', 'bioc_29', 'bioc_31', 'bioc_32')
gcms <- names(clma) %>% str_split('_') %>% map_chr(3) %>% unique(); gcme <- gcms[1]
indx <- terra::rast('./common_data/atlas_hazards/ftre-gcms_indices2.tif')
names(clma) <- gsub('bioc_2060', 'bioc_20_60', names(clma))

vars

## Grepping Climate and Index
clma <- clma[[grep(paste0(paste0(vars, '_'), collapse = '|'), names(clma))]]
vrs <- map_chr(str_split(vars, '_'), 1)
vrs.idx <- vrs[-grep('bioc', vrs)]

## Presences and thresholds
load('./workspace_coffee/rData/run_2/clustereddata.rData')
thrp <- readRDS(file = './workspace_coffee/rds/run_2/threshold_prob_coffee.rds')
thru <- readRDS(file = './workspace_coffee/rds/run_2/threshold_uncr_coffee.rds')

## AEZ
lbls <- tibble(value = 1:6, class = c('Unsuitable',  'Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed'))

# To predict --------------------------------------------------------------
make.predict <- function(gcme){
  
  cat('To process: ', gcme, '\n')
  
  ## Grepping Climate and Index
  clm <- clma[[grep(gcme, names(clma))]]
  idx <- indx[[grep(gcme, names(indx))]]
  idx <- idx[[grep(paste0(vrs.idx, collapse = '|'), names(idx))]]
  
  ## To make the stack (Climate + Index)
  stk <- c(clm, idx)
  names(stk) <- gsub(paste0('_', gcme), '', names(stk))
  names(stk) <- gsub('bioc_20100', 'bioc_20_100', names(stk))
  
  ## Change the names for the stacking
  names(stk) <- gsub('ndd', 'ndd_hist', names(stk))
  names(stk) <- gsub('n30', 'n30_hist', names(stk))
  names(stk) <- gsub('n35', 'n35_hist', names(stk))
  names(stk) <- gsub('ndw', 'ndw_hist', names(stk))
  
  names(stk) <- gsub(paste0('_', gcme), '', names(stk))
  
  vls <- values(stk)
  lyr <- stk
  
  ## To predict
  rasterProbs <- predict(rff, vls, type = 'prob')
  
  rasterRF <- rowSums(rasterProbs[,c(3:(2+2))])
  uncertainty <- apply(rasterProbs, 1, max)  
  
  rasterRFprob <- lyr[[1]]
  values(rasterRFprob) <- rasterRF 
  
  rasterRFuncertainty <- lyr[[1]]
  values(rasterRFuncertainty) <- uncertainty 
  
  rasterRF <- max.col(rasterProbs, 'first')
  rasterRFclass <- lyr[[1]]
  values(rasterRFclass) <- rasterRF
  
  ## Plots 
  plot(rasterRFclass)
  plot(rasterRFprob)
  plot(rasterRFuncertainty)
  
  ## Return the raster 
  names(rasterRFclass) <- glue('rf_{gcme}')
  names(rasterRFprob) <- glue('prob_{gcme}')
  names(rasterRFuncertainty) <- glue('uncr_{gcme}')
  
  ## Finish 
  cat("Done!\n")
  return(c(rasterRFclass, rasterRFprob, rasterRFuncertainty))
  
  
}
prdc <- map(gcms, make.predict)
prdc <- reduce(prdc, c)

dir_create('./workspace_coffee/rf/output/run_2/results/process')
terra::writeRaster(x = prdc, filename = './workspace_coffee/rf/output/run_2/results/process/rf-raw_ftr-gcms.tif', overwrite = TRUE)

# To add Limitations ------------------------------------------------------
make.uncr <- function(gcme){
  
  ## To filter 
  cat('To process: ', gcme, '\n')
  prd <- prdc[[grep(gcme, names(prdc), value = F)]]
  
  ## Prob / Cluster / Uncertainty
  prb <- prd[[grep('prob', names(prd))]]
  rfr <- prd[[grep('rf', names(prd))]]
  unc <- prd[[grep('uncr', names(prd))]]
  
  ## To classify prob and uncertainty
  prb.bin <- terra::ifel(prb < thrp, 0, 2)
  cls.bin <- terra::ifel(rfr < 2.1,  0, 1)
  
  ### Make the difference / Limitations
  dfr <- prb.bin - cls.bin
  rsl <- rfr
  rsl[dfr == -1] <- 5
  rsl[dfr ==  2] <- 5
  
  ### Make mixed zones
  fnl <- rsl
  fnl[unc < thru & prb > thrp] <- 6
  fnl <- as.factor(fnl)
  
  ## Change the names 
  names(fnl) <- glue('mix_{gcme}')
  return(fnl)
  
}
mixd <- map(gcms, make.uncr)
mixd <- reduce(mixd, c)
names(mixd) <- glue('mixed_{gcms}')
terra::writeRaster(x = mixd, filename = './workspace_coffee/rf/output/run_2/results/process/rf-mix_ftr-gcms.tif', overwrite = TRUE)

# To draw the maps --------------------------------------------------------

## Tidy the table
tble <- as_tibble(terra::as.data.frame(mixd, xy = T)) %>% gather(var, value, -c(x, y ))
tble <- mutate(tble, gcme = gsub('mix_', '', var), value = as.numeric(value))
tble <- inner_join(tble, lbls, by = 'value')
tble <- mutate(tble, class = factor(class, levels = unique(lbls$class)))

## Tidy the raster
mixd <- as.factor(mixd)
mixd <- map(.x = 1:nlyr(mixd), .f = function(i){
  # levels(mixd)[[i]] <- data.frame(ID = 1:9, clase = c('Unsuitable', 'Unsuitable', 'Very hot - Wet', 'Very cold - Very dry', 'Cold - Very dry', 'Hot - Very wet', 'Hot - Dry', 'Limitations', 'Mixed'))
  levels(mixd)[[i]] <- data.frame(ID = 1:6, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed'))
  return(mixd[[i]])
})
mixd <- reduce(mixd, c)

## Colors
# clrs <- c('#FFFFFF', '#F9A03F', '#6BAED6', '#2171B5', '#238B45', '#D94801', '#D9D9D9', '#FFFFB4')
clrs <- c('#FFFFFF',  '#3E5C23', '#C6E01B',  '#D9D9D9', '#FFFFB4')
# names(clrs) <- c('Unsuitable', 'Very hot - Wet', 'Very cold - Very dry', 'Cold - Very dry', 'Hot - Very wet', 'Hot - Dry', 'Limitations', 'Mixed')
names(clrs) <- c('Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed')

## To make the map 
ggps <- map(.x = 1:nlyr(mixd), .f = function(i){
  
  mix <- mixd[[i]]
  gcm <- gcms[i]
  
  cat('To make the map: ', gcm, '\n')
  gft <- ggplot() +
    geom_spatraster(data = mix) + 
    scale_fill_manual(values = clrs, na.value = 'white', na.translate = FALSE) +
    geom_sf(data = st_as_sf((wrld)), fill = NA, col = 'grey30') +
    ggtitle(label = gcm) +
    coord_sf(xlim = c(-120, 140), ylim = c(-30, 30)) +
    labs(fill = '', x = '', y = '') +
    theme_bw() +
    theme(legend.position = 'bottom', axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold'), axis.text.x = element_blank(), axis.text.y = element_blank()) # axis.text.y = element_text(angle = 90, hjust = 0.5), 
  
  cat('Done!\n')
  return(gft)
  
})

## Join all the maps into only one
ggps.all <- ggpubr::ggarrange(ggps[[1]], ggps[[2]], ggps[[3]], ggps[[4]], ggps[[5]], ncol = 3, nrow = 2, common.legend = TRUE, legend = 'bottom')
ggps.all
ggsave(plot = ggps.all, filename = './workspace_coffee/png/maps/run_2/rf_cluster5_ara-rob_gcms.jpg', units = 'in', width = 12, height = 4, dpi = 300)

## Modal
mdal <- terra::modal(mixd)
# levels(mdal)[[1]] <- data.frame(ID = 1:9, clase = c('Unsuitable', 'Unsuitable', 'Very hot - Wet', 'Very cold - Very dry', 'Cold - Very dry', 'Hot - Very wet', 'Hot - Dry', 'Limitations', 'Mixed'))
levels(mdal)[[1]] <- data.frame(ID = 1:6, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed'))

### Draw the map for the modal
gmdl <- ggplot() +
  geom_spatraster(data = mdal) + 
  scale_fill_manual(values = clrs, na.value = 'white', na.translate = FALSE) +
  geom_sf(data = st_as_sf((wrld)), fill = NA, col = 'grey30') +
  ggtitle(label = 'Future [modal]') +
  coord_sf(xlim = c(-120, 140), ylim = c(-30, 30)) +
  labs(fill = '', x = '', y = '') +
  theme_bw() +
  theme(legend.position = 'bottom', axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold'), axis.text.x = element_blank(), axis.text.y = element_blank()) # axis.text.y = element_text(angle = 90, hjust = 0.5), 

ggps.mdl <- ggpubr::ggarrange(ggps[[1]], ggps[[2]], ggps[[3]], ggps[[4]], ggps[[5]], gmdl, ncol = 3, nrow = 2, common.legend = TRUE, legend = 'bottom')
ggsave(plot = ggps.mdl, filename = './png/maps/coffee/rf_cluster5_ara-rob_gcms-modal.jpg', units = 'in', width = 12, height = 4, dpi = 300)
ggsave(plot = gmdl, filename = './workspace_coffee/png/maps/run_2/rf_cluster5_ara-rob_modal.jpg', units = 'in', width = 7, height = 3, dpi = 300)

## Join both 
ggall <- ggpubr::ggarrange(ggps.mdl, gmdl, ncol = 1, nrow = 2)
ggsave(plot = ggall, filename = './png/maps/coffee/rf_cluster5_ara-rob_modal_gcms.jpg', units = 'in', width = 15, height = 9.5, dpi = 300)

# To write the rasters ----------------------------------------------------
terra::writeRaster(x = mdal, filename = './workspace_coffee/rf/output/run_2/results/process/rf-mix_ftr-modal.tif', overwrite = TRUE)
terra::writeRaster(x = mixd, filename = './workspace_coffee/rf/output/run_2/results/process/rf-mix_ftr-gcms.tif', overwrite = TRUE)

# Join current and future map  --------------------------------------------
rasterRFclass <- rast('./workspace_coffee/rf/output/run_2/results/process/rf_mixed-bsl.tif')
rasterRFclass <- as.factor(rasterRFclass)
# levels(rasterRFclass)[[1]] <- data.frame(ID = 1:7, clase = c('Unsuitable', 'Unsuitable', paste0('Type ', 1:5)))
levels(rasterRFclass)[[1]] <- data.frame(ID = 1:6, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed'))

clrs <- c('white', '#3E5C23', '#C6E01B', '#D9D9D9', '#FFFFB4')
names(clrs) <- c('Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed')

gg <- ggplot() + 
  geom_spatraster(data = rasterRFclass) + 
  # scale_fill_manual(values = c('white', '#F9A03F', '#6BAED6', '#2171B5', '#238B45', '#D94801'), na.value = 'transparent', na.translate = F) +
  scale_fill_manual(values = clrs, na.value = 'transparent', na.translate = F) +
  geom_sf(data = st_as_sf((wrld)), fill = NA, col = 'grey30') +
  coord_sf(xlim = c(-115, 170), ylim = c(-30, 30)) +
  ggtitle(label = 'Current') +
  labs(fill = '') +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    plot.title = element_text(hjust = 0.5, face = 'bold'),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 6),
    axis.text.x = element_text(hjust = 0.5, size = 6)
  )

gg.two <- ggpubr::ggarrange(gg, gmdl, ncol = 1, nrow = 2, common.legend = TRUE, legend = 'bottom')
ggsave(plot = gg.two, filename = './workspace_coffee/png/maps/run_2/rf_cluster5_ara-rob_crnt-modal.jpg', units = 'in', width = 13, height = 7.5, dpi = 300)
