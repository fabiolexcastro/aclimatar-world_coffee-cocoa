

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, usdm, sp, raster, RColorBrewer, randomForest, outliers, readxl, hablar, glue, readxl, geodata, RColorBrewer)
pacman::p_load(tidyverse, raster, ggpubr, cptcity, corrplot, cclust, dismo, gtools, sp, FactoMineR, pROC, randomForest, Hmisc)
p_load(terra, sf, fs, tidyverse, spatstat, RColorBrewer, parallelDist, Rfast, glue, outliers, spatialEco, climateStability, scales, glue, rnaturalearthdata, rnaturalearth, openxlsx)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

source('https://raw.githubusercontent.com/fabiolexcastro/rfSalvador/refs/heads/master/FunctionsRFclustering.R')

# Function ----------------------------------------------------------------
rf.clust <- function(occ, nforest, ntrees, nVars, nclasses){
  
  datRF_presences <- occ[,4:ncol(occ)]
  print(nrow(datRF))
  
  attach(datRF_presences)
  no.forests <- nforest
  no.trees <- ntrees
  distRF_presences <- RFdist(datRF_presences, mtry1 = nVars, no.trees, no.forests, addcl1 = T, addcl2 = F, imp = T, oob.prox1 = T)
  no.presencesclasses <- nclasses
  labelRF <- pamNew(distRF_presences$cl1, no.presencesclasses)
  print(table(labelRF))
  clusterdata <- hclust(as.dist(distRF_presences$cl1), method = 'ward.D2')
  
  return(list(labelRF, clusterdata))
  
}

# Load data ---------------------------------------------------------------

## Vector -----------------------------
wrld <- geodata::world(resolution = 1, path = './tmpr')
cntn <- aggregate(wrld)
vtnm <- wrld[wrld$GID_0 == 'VNM',]
plot(wrld)

## Tabular data ----------------------
# cffe <- read_csv('./tbl/points/clean/v2/arabica_robusta_otl_swd.csv', show_col_types = FALSE)
cffe <- read_csv('./tbl/points/coffee/arabica_robusta_otl_swd.csv', show_col_types = FALSE)

cffe.ara <- filter(cffe, Specie == 'arabica')
cffe.rob <- filter(cffe, Specie == 'robusta')
 
## Sampling arabica ## Add arabica and robusta into only one table
# cffe.ara <- cffe.ara %>% sample_n(size = nrow(cffe.rob), replace = FALSE, tbl = .)
# cffe <- rbind(cffe.ara, cffe.rob)
# write.csv(cffe, './tbl/points/clean/v2/arabica_robusta_otl_swd_sampling.csv', row.names = FALSE)

# Map ---------------------------------------------------------------------
g.smp <- ggplot() + 
  geom_point(data = cffe, aes(x = Lon, y = Lat, col = Specie), size = 0.5) +
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey40') +
  coord_sf(ylim = c(-30, 30), xlim = c(-110, 150)) +
  theme_bw() +
  labs(caption = 'Arabica: 259\nRobusta: 259') +
  theme(axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6, angle = 90),
        legend.position = 'bottom') +
  guides(color = guide_legend(override.aes = list(size = 5)))

g.smp
ggsave(plot = g.smp, filename = './png/maps/points/v2/point_ara-rob_sampling.jpg', units = 'in', width = 9, height = 4, dpi = 300, create.dir = T)
 
## Presences first model --------------------------------------------------------
# load('./rData/run_3_arabica-robusta/presences.rData')
# clusteredpresdata
# occ <- clusteredpresdata %>% dplyr::select(-ID)
# cffe <- occ


# New presences -----------------------------------------------------------
cffe

## Raster data -----------------------
bioc <- terra::rast('./common_data/input_bios/bioc_hist.tif')
indx <- terra::rast('./common_data/atlas_hazards/hist_indices2.tif')
indx <- indx[[grep('hist', names(indx))]]
bioc <- c(bioc, indx)
# bioc <- bioc[[-grep('100', names(bioc))]]

## Boxplot indices 
# gg.box <- cffe %>% dplyr::select(Lon, Lat, Specie, ends_with('hist')) %>% gather(var, value, -c(Lon, Lat, Specie)) %>% separate(col = 'var', into = c('index', 'period'), sep = '_') %>% ggplot(data = ., aes(x = Specie, y = value, fill = Specie)) + geom_boxplot() + facet_wrap(.~index, scales = 'free_y') + labs(title = 'Boxplot Indices - Baseline', y = 'Value', x = 'Specie') + theme_bw() + theme(plot.title = element_text(face = 'bold', hjus = 0.5), strip.text = element_text(face = 'bold', hjust = 0.5))
# ggsave(plot = gg.box, filename = './png/graphs/v2/boxplot_indices_arabica-robusta.jpg', units = 'in', width = 8, height = 6, dpi = 300, create.dir = TRUE)


# Function for the VIF ----------------------------------------------------
occ <- cffe %>% drop_na()

vars <- map(.x = c(2.5, 5, 10, 15), .f = function(i){
  vif.res <- vifstep(x = as.data.frame(occ)[5:ncol(occ)], th = i, keep = 'bioc_6')
  vrs <- vif.res@results$Variables %>% as.character()
  rsl <- tibble(threshold = i, variable = vrs)
  return(rsl)
})
vars <- bind_rows(vars)
write.csv(vars, './tbl/vif_coffee-arabica-robusta_keepBioc6.csv', row.names = FALSE)

# VIF analysis ------------------------------------------------------------
occ <- cffe
occ <- drop_na(occ)
vif.res <- vif(x = as.data.frame(occ)[,5:ncol(occ)])
vif.step <- vifstep(x = as.data.frame(occ)[,5:ncol(occ)], th = 10)
vrs <- vif.step@results$Variables %>% as.character()
vrs <- vrs[-grep('_100', vrs)]
dir.create('./rds/v4_coffee')
saveRDS(vrs, './rds/v4_coffee/vars_coffee_arabica-robusta.rds')
vrs <- readRDS('./rds/v4_coffee/vars_coffee_arabica-robusta.rds')
# vrs
# [1] "ndd_hist"   "n30_hist"   "n35_hist"   "ndw_hist"   "bioc_3"     "bioc_4"     "bioc_8"    
# [8] "bioc_13"    "bioc_14"    "bioc_18"    "bioc_19"    "bioc_20_60" "bioc_22"    "bioc_25"   
# [15] "bioc_26"    "bioc_27"    "bioc_29"    "bioc_31"    "bioc_32"

## Corrplot
mtrx <- occ[,5:ncol(occ)]
corr <- as.matrix(mtrx)
dir.create('./png/graphs/coffee')
png(filename = './png/graphs/coffee/corplot_arabica-robusta.jpg', units = 'in', width = 9, height = 8, res = 300)
corrplot::corrplot(cor(corr))
dev.off()

# Mapping the variables ---------------------------------------------------
bioc.vars <- bioc[[grep(paste0(paste0(vrs, '$'), collapse = '|'), names(bioc), value = F)]]
bioc.vars
plot(bioc.vars)

bioc.vars <- bioc.vars %>% terra::as.data.frame(xy = T) %>% as_tibble() %>% gather(var, value, -c(x, y))
unique(bioc.vars$var)
bioc.vars.tasm <- bioc.vars %>% filter(var %in% c('bioc_2', 'bioc_3', 'bioc_4', 'bioc_8')) ## Temperature 
bioc.vars.prec <- bioc.vars %>% filter(var %in% c('bioc_13', 'bioc_14', 'bioc_18', 'bioc_19')) ## Precipitation
bioc.vars.etps <- bioc.vars %>% filter(var %in% c('bioc_22', 'bioc_27', 'bioc_29'))
bioc.vars.baln <- bioc.vars %>% filter(var %in% c('bioc_30', 'bioc_31'))
bioc.vars.idts <- bioc.vars %>% filter(var %in% c('n30_hist', 'n35_hist', 'ndd_hist', 'ndw_hist'))

## Mapping variables
make.mapp <- function(vari, clrs){
  
  # vari <- 'bioc_19'; clrs <- (brewer.pal(n = 9, name = 'BrBG'))
  
  ## To filter the table
  cat('To process: ', vari, '\n')
  tbl <- filter(bioc.vars, var == vari)
  
  ## Treshold 98
  hist(tbl$value)
  if(vari %in% c('bioc_13', 'bioc_14', 'bioc_18', 'bioc_19', 'bioc_22', 'bioc_27', 'bioc_29')){
    tbl <- mutate(tbl, value = ifelse(value >= quantile(tbl$value, 0.99, na.rm = T), quantile(tbl$value, 0.99, na.rm = T), value))  
  }

  ## To draw the map 
  ggm <- ggplot() + 
    geom_tile(data = tbl, aes(x = x, y = y, fill = value)) + 
    scale_fill_gradientn(colors = clrs, name = gsub('_', ' ', vari)) + 
    geom_sf(data = st_as_sf(cntn), fill = NA, col = 'grey30') +
    coord_sf(xlim = c(-130, 150), ylim = c(-25, 25)) +
    labs(title = vari, x = '', y = '') + 
    theme_bw() +
    theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5), legend.key.width = unit(3, 'line'), axis.text = element_text(size = 6), legend.title.position = 'top', legend.title = element_text(hjust = 0.5))

  ## To return the map
  cat('Done!\n')
  return(ggm)
  
  
}
maps.tasm <- map2(.x = c('bioc_2', 'bioc_3', 'bioc_8', 'bioc_8'), .y = replicate(4, brewer.pal(n = 9, name = "YlOrRd"), simplify = FALSE), .f = make.mapp)
maps.prec <- map2(.x = c('bioc_13', 'bioc_14', 'bioc_18', 'bioc_19'), .y = replicate(4, brewer.pal(n = 9, name = "BrBG"), simplify = FALSE), .f = make.mapp)
maps.etps <- map2(.x = c('bioc_22', 'bioc_27', 'bioc_29'), .y = replicate(3, rev(brewer.pal(n = 9, name = "BrBG")), simplify = FALSE), .f = make.mapp)
maps.baln <- map2(.x = c('bioc_30', 'bioc_31'), .y = replicate(2, rev(brewer.pal(n = 9, name = "BrBG")), simplify = FALSE), .f = make.mapp)
maps.ints <- map2(.x = c('n30_hist', 'n35_hist'), .y = replicate(2, (brewer.pal(n = 9, name = "YlOrRd")), simplify = FALSE), .f = make.mapp)
maps.ndds <- map2(.x = c('ndd_hist'), .y = replicate(1, rev(brewer.pal(n = 9, name = "BrBG")), simplify = FALSE), .f = make.mapp)
maps.ndws <- map2(.x = c('ndw_hist'), .y = replicate(1, (brewer.pal(n = 9, name = "BrBG")), simplify = FALSE), .f = make.mapp)
maps
gmap <- cowplot::plot_grid(plotlist = c(maps.tasm, maps.prec, maps.etps, maps.baln, maps.ints, maps.ndds, maps.ndws), ncol = 5, nrow = 4, common.legend = FALSE)
ggplot2::ggsave(plot = gmap, filename = './png/maps/v2/bios_bsl_vars_ara-rob.jpg', width = 17, height = 11,  dpi = 300, bg = "white")

# To select the variables -------------------------------------------------
occ <- occ %>% dplyr::select(Lon, Lat, Specie, vrs)
nrw.ara <- occ %>% filter(Specie == 'arabica') %>% nrow()
nrw.rob <- occ %>% filter(Specie == 'robusta') %>% nrow()

# Draw mapping ------------------------------------------------------------
g.cff.ara <- ggplot() +
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') +
  geom_point(data = occ, aes(x = Lon, y = Lat, col = Specie), size = 0.3) +
  coord_sf(ylim = c(-30, 30), xlim = c(-110, 135)) +
  ggtitle(label = 'C. arabica - C. robusta') +
  labs(caption = paste0('C arabica: ', nrw.ara, '\n', 'C Robusta: ', nrw.rob)) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        legend.position = 'bottom',
        text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
g.cff.ara
ggsave(plot = g.cff.ara, filename = './png/maps/points/arabica-robusta_world.jpg', units = 'in', width = 8, height = 3, dpi = 300, create.dir = T)

# Sampling for arabica ----------------------------------------------------
# write.csv(occ, './tbl/points/clean/arabica-robusta_otl_swd_sampling.csv', row.names = FALSE)

# Clustering --------------------------------------------------------------
env_values <- as.matrix(occ[,4:ncol(occ)]); nrow(env_values)
datRF <- as.data.frame(occ[,4:ncol(occ)]); nrow(datRF)
d <- dist(datRF, method = "euclidean")  
rfClust <- rf.clust(occ = occ, nforest = 25, ntrees = 100, nVars = 8, nclasses = 5)
labelRF <- rfClust[[1]]
clusterdata <- rfClust[[2]]
classdata <- cbind(pb = as.factor(labelRF), occ[,3:ncol(occ)])
clusteredpresdata <- cbind(occ, cluster = labelRF) %>% na.omit() %>% tbl_df()
no.clusters <- 5
clusteredpresdata <- clusteredpresdata %>% mutate(cluster = ifelse(Specie == 'arabica', 1, 2))

dir_ls('./rData')
dir.create('./rData/run_7_arabica-robusta')
run <- 'run_7_arabica-robusta'
save(datRF, file = paste0('./rData/', run, '/datRF.rData'))
save(clusterdata, file = paste0('./rData/', run, '/clusterdata.rData'))
save(occ, clusteredpresdata, no.clusters, labelRF, file = paste0('./rData/', run, '/clustereddata.rData'))
save(clusteredpresdata, file = './rData/run_7_arabica-robusta/presences.rData')

# Mapping clustering ------------------------------------------------------
clusteredpresdata %>% distinct(Specie, cluster)

g.clst <- ggplot() + 
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') + 
  geom_point(data = clusteredpresdata, aes(x = Lon, y = Lat, col = factor(cluster))) + 
  scale_color_manual(labels = c('arabica', 'robusta'), values = c('tomato3', 'slateblue2')) +
  # facet_wrap(.~cluster) +
  coord_sf(ylim = c(-30, 30), xlim = c(-130, 140)) +
  theme_bw() +
  theme(
    legend.position = 'bottom'
  )

g.clst
ggsave(plot = g.clst, filename = './png/maps/coffee/rf_cluster5_ara-rob.jpg', units = 'in', width = 9, height = 5, dpi = 300)

## Check arabica vs robusta -----------------------------------------------
clst <- clusteredpresdata %>% dplyr::select(Lon, Lat, Specie, cluster)
clst <- clst %>% group_by(Specie, cluster) %>% reframe(count = n()) 
clst <- clst %>% mutate(cluster = factor(cluster, levels = 1:5))
clst %>% spread(cluster, count) %>% kableExtra::kable()

gcls <- ggplot(data = clst, aes(x = cluster, y = count, fill = Specie)) + 
  geom_col(position = 'dodge') +
  labs(x = 'Cluster', y = 'Frequence') +
  ggtitle(label = 'Count coffee presences (ara / rob) by each cluster [RF]') +
  theme_bw() +
  theme(
    legend.position = 'bottom'
  )

gcls
ggsave(plot = gcls, filename = './png/maps/coffee/count_rf-cluster_arabica-robusta.jpg', units = 'in', width = 6, height = 5, dpi = 300, create.dir = T)

# Bias background sampling ------------------------------------------------

## Harvested area map
faos <- read_csv('./tbl/faostat/FAOSTAT_data_coffee.csv', show_col_types = FALSE)
faos <- faos %>% group_by(Área) %>% reframe(Value = mean(Valor, na.rm = T)) %>% ungroup() %>% drop_na() %>% mutate(Value = round(Value, digits = 1)) %>% arrange(desc(Value))
# write.csv(faos, './tbl/faostat/faostat_average_2000-2023_coffee.csv', row.names = FALSE)
faos <- read_csv('./tbl/faostat/faostat_coffee_iso3.csv', show_col_types = FALSE)[,c('ISO3', 'Value')] 
faos <- full_join(st_as_sf(wrld), faos, by = c('GID_0' = 'ISO3'))


## Mapping
g.faos <- ggplot() +
  geom_sf(data = faos, aes(fill = Value)) + 
  scale_fill_gradientn(labels = scales::comma, na.value = 'transparent', colors = RColorBrewer::brewer.pal(n = 9, name = 'YlOrBr'), breaks = seq(0, 3000000, by = 500000)) +
  geom_sf(data = st_as_sf(cntn), fill = NA, col = 'grey40') +
  coord_sf(xlim = c(-125, 150), ylim = c(-30, 30)) +
  ggtitle(label = 'Harvested area - Coffee [Faostat]') +
  labs(fill = 'Hectareas', caption = 'Average 2000-2023') +
  theme_minimal() +
  theme(legend.title = element_text(hjust = 0.5), legend.title.position = 'top',
        legend.position = 'bottom', legend.key.width = unit(5, 'line'), legend.key.height = unit(0.5, 'line'))
g.faos
ggsave(plot = g.faos, filename = './png/maps/coffee/harvested_area_faostat.jpg', units = 'in', width = 8, height = 4, dpi = 300, create.dir = T)

## Presences
clusteredpresdata
frqn <- terra::extract(wrld, clusteredpresdata[,c('Lon', 'Lat')])$GID_0
frqn <- as.data.frame(table(frqn)) %>% arrange(desc(Freq)) %>% setNames(c('iso', 'frqn'))
faos <- full_join(faos, frqn, by = c('GID_0' = 'iso'))
faos <- faos %>% mutate(prob1 = frqn / Value, prob2 = Value / frqn)
faos <- faos %>% drop_na(prob1)
# faos <- faos %>% mutate(prob1 = ifelse(is.na(prob1), 0.00001, prob1))

## Calc CDF
dns <- density(faos$prob1)
cdf <- CDF(dns)
cnv <- cdf(faos$prob1)
faos <- mutate(faos, prob1 = cnv)

# faos <- mutate(faos, prob1_2 = prob1 / sum(prob1, na.rm = T)) # Version 2
# faos %>% drop_na(prob1)
# faos <- mutate(faos, prob1_2 = ifelse(is.na(prob1_2), 0.00001, prob1_2)) # min(faos$prob1_2, na.rm = T)

## Raster bias background sampling
mskr <- bioc[[1]] * 0 + 1
mskr <- terra::rasterize(x = faos, y = mskr, field = 'prob1') #prob1_2
bckg <- as_tibble(dismo::randomPoints(mask = raster::raster(mskr), n = nrow(clusteredpresdata), prob = TRUE))
plot(mskr)

## Count the presences by each country 
bckg <- mutate(bckg, iso = terra::extract(wrld, bckg[,c(1, 2)])$GID_0)
frqn.bckg <- bckg %>% group_by(iso) %>% reframe(count = n()) %>% arrange(desc(count))
faos <- full_join(faos, frqn.bckg, by = c('GID_0' = 'iso')) %>% rename(count_bckg = count, count_occr = frqn)

## Mapping the pseudo-absences
gg.faos <- g.faos
gg.mskr <- ggplot() + ggtitle(label = 'Probability for the background sampling') + geom_spatraster(data = mskr, aes(fill = prob1)) + scale_fill_viridis_c(na.value = 'transparent') + geom_sf(data = faos, fill = NA, col = 'grey30') + labs(fill = 'Prob bckg') + coord_sf(xlim = c(-120, 140), ylim = c(-30, 30)) + theme_minimal() + theme(legend.title = element_text(hjust = 0.5), legend.title.position = 'top', legend.position = 'bottom', legend.key.width = unit(5, 'line'), legend.key.height = unit(0.5, 'line'))
gg.bckg <- ggplot() + geom_point(data = bckg, aes(x = x, y = y), col = 'grey30', size = 0.1) + geom_sf(data = faos, fill = NA, col = 'grey30') + labs(x = '', y = '') + ggtitle(label = 'Background sampling') + coord_sf(xlim = c(-120, 140), ylim = c(-30, 30)) + theme_minimal() + theme(legend.title = element_text(hjust = 0.5), legend.title.position = 'top', legend.position = 'bottom', legend.key.width = unit(5, 'line'), legend.key.height = unit(0.5, 'line'))
gg.frbc <- ggplot() + ggtitle(label = 'Count background by country') + geom_sf(data = faos, aes(fill = count_bckg), col = 'grey30') + scale_fill_viridis_c() +coord_sf(xlim = c(-120, 140), ylim = c(-30, 30)) + labs(fill = 'Count background') + theme_minimal() + theme(legend.position = 'bottom', legend.key.width = unit(3, 'line'), legend.title.position = 'top', legend.title = element_text(hjust = 0.5))
gg.allp <- ggpubr::ggarrange(gg.faos, gg.mskr, gg.bckg, gg.frbc, ncol = 2, nrow = 2)
ggsave(plot = gg.allp, filename = './png/maps/coffee/background_sampling_v2.jpg', units = 'in', width = 14, height = 7, dpi = 300, create.dir = T)

# Background --------------------------------------------------------------
SPspecies <- SpatialPoints(occ[,1:2]) 
crs(SPspecies) <- crs(bioc)
back_raster <- bioc[[1]]
speciescell <- raster::extract(raster(bioc[[1]]), SPspecies, cellnumber = TRUE)
back_raster[speciescell[,1]]  <- NA #remove the cell with presences
samplesize <- round(min(summary(as.factor(clusteredpresdata$cluster))) / 2, 0) 
NumberOfClusters <- max(clusteredpresdata$cluster) 
ratio <- NumberOfClusters/1
numberofpresences <- nrow(clusteredpresdata) 
crs(back_raster) <- crs(bioc)
back <- randomPoints(raster(back_raster), 1*numberofpresences) %>% as_data_frame() # change 1 by 2
back_sp <- back
coordinates(back) <- ~ x + y
back_swd <- terra::extract(bioc, back_sp) %>% cbind(back_sp, .)
back_swd <- as_tibble(back_swd)
back_swd <- dplyr::select(back_swd, x, y, vrs)
nrow(back_swd)
nrow(drop_na(back_swd))
back_swd <- drop_na(back_swd)
dir.create('./tbl/points/clean/v3')
back_swd %>% write.csv('./tbl/points/clean/v3/arabica-robusta_background_swd-vars.csv', row.names = FALSE)
back_swd

## Back swd
bckr
back <- bckg[,1:2]
back_swd <- terra::extract(bioc, back[,1:2]) %>% cbind(back, .)
back_swd <- dplyr::select(back_swd, x, y, vrs)
back_swd <- drop_na(back_swd)
back_swd %>% write.csv('./tbl/points/coffee/arabica-robusta_background_swd-vars.csv', row.names = FALSE)
back_swd <- mutate(back_swd, Specie = 'back', .after = 'y')

# Cluster analysis to pseudoabsences
bckclust <- rf.clust(occ = back_swd, nforest = 50, ntrees = 500, nVars = 8, nclasses = 2)
datRF <- as.data.frame(back_swd[,3:ncol(back_swd)])
attach(datRF)
no.forests <- 50#raw = 25
no.trees <- 500
distRF <- RFdist(datRF, mtry1 = 8, no.trees, no.forests, addcl1 = T, addcl2 = F, imp =T, oob.prox1 = T)# mtry1 = 4 raw  # es la cantidad de variables a utilizar en cada no
no.absenceclasses <- 2
labelRF <- pamNew(distRF$cl1, no.absenceclasses)
labelRF <- c(rep(1, 500), rep(2, 958))
detach(datRF)
classdata <- cbind(pb = as.factor(labelRF), back_swd[,3:ncol(back_swd)])

presvalue_swd  <- clusteredpresdata[,4:ncol(clusteredpresdata)] %>%
  cbind(pb = (clusteredpresdata$cluster + no.absenceclasses), .) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(cluster = cluster + no.absenceclasses)
presvalue_swd <- dplyr::select(presvalue_swd, pb, vrs)
presvalue_swd <- mutate(presvalue_swd, pb = as.factor(pb))
classdata_2 <- cbind(pb = as.data.frame(classdata)$pb, classdata[,2:ncol(classdata)]) # Background
classdata_2 <- dplyr::select(classdata_2, -Specie)

dim(classdata_2); dim(presvalue_swd)

allclasses_swd <- rbind(classdata_2, presvalue_swd[,1:ncol(classdata_2)])
unique(allclasses_swd$pb)
dir_create('./tbl/rf/points')
write.csv(allclasses_swd, './tbl/rf/points/run-7_all_classes_swd.csv', row.names = FALSE)

# To make the random forest analysis --------------------------------------
model1 <- as.formula(paste('factor(pb) ~', paste(paste(vrs), collapse = '+', sep =' ')))
rflist <- vector('list', 50) 
auc <- vector('list', 50)

allclasses_swd %>% pull(pb) %>% table()

for(repe in 1:50){ # 50 bosques
  
  print(repe)
  pressample <- list()
  
  for (i in 1:(NumberOfClusters+no.absenceclasses)){
    
    print(i)  
    if(any(i==c(1:no.absenceclasses))) { 
      
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), 
                     size = samplesize*NumberOfClusters/2/no.absenceclasses,replace=F)
    } else {
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), size=samplesize, replace=F)
    }
    pressample[[i]] <- allclasses_swd[rows,] 
  }
  
  species <- na.omit(do.call(rbind, pressample)) 
  head(species)
  Samplesplit <- sample(rownames(species)) 
  
  envtrain <- species[Samplesplit[1:(0.8*nrow(species))],] 
  envtest <- species[Samplesplit[(0.8*nrow(species)):nrow(species)],] 
  
  rfmodel <- randomForest(model1, data = envtrain, ntree = 500, na.action = na.omit, nodesize = 2) 
  
  save(rfmodel, file = paste('./rf/output/run_8_coffee/models/', NumberOfClusters, 'Prob_' , 'rep_' ,repe, '.rdata' ,sep=''))
  rflist[[repe]] <- rfmodel
  
  # AUC 
  predicted <- as.numeric(predict(rfmodel, envtest))
  observed <- as.vector(envtest[,'pb'])
  auc[[repe]] <- auc(as.numeric(observed), as.numeric(predicted))
  rm(rfmodel)
  
  cat(auc[[repe]] ,'\n')
  
}

auc <- unlist(auc)
boxplot(auc)
rff <- do.call(randomForest::combine, rflist)
importance <- as.data.frame(rff$importance)

saveRDS(object = rff, file = './rds/v6_coffee/rff.rds')

# Predict modell
lyr <- bioc[[grep(paste0(paste0(vrs, '$'), collapse = '|'), names(bioc), value = F)]]
climatevalues  <- data.frame(values(lyr))
NumberOfClusters <- 2

rasterProbs <- predict(rff, climatevalues, type = 'prob') # proximity = T
rasterProbs_na <- na.omit(rasterProbs)
sum_rasterProbs_na <- apply(rasterProbs_na, 1, sum)

rasterRF <- rowSums(rasterProbs[,c(3:(NumberOfClusters+2))])
uncertainty <- apply(rasterProbs, 1, max)  

rasterRFprob <- lyr[[1]]
values(rasterRFprob) <- rasterRF 

rasterRFuncertainty <- lyr[[1]]
values(rasterRFuncertainty) <- uncertainty 

rasterRF <- max.col(rasterProbs, 'first')
rasterRFclass <- lyr[[1]]
values(rasterRFclass) <- rasterRF

plot(rasterRFclass)
plot(rasterRFprob)
plot(rasterRFuncertainty)

## To write the rasters
writeRaster(rasterRFclass, paste0('./rf/output/run_8_coffee/results/raw/RF_2Clust_current.tif'), overwrite = T)
writeRaster(rasterRFprob, paste0('./rf/output/run_8_coffee/results/raw/RF_2Prob_current.tif'), overwrite = T)
writeRaster(rasterRFuncertainty, paste0('./rf/output/run_8_coffee/results/raw/RF_2Unc_current.tif'), overwrite = T)


# Mapping -----------------------------------------------------------------


## Clustering
rasterRFclass <- as.factor(rasterRFclass)
# levels(rasterRFclass)[[1]] <- data.frame(ID = 1:7, clase = c('Unsuitable', 'Unsuitable', paste0('Type ', 1:5)))
levels(rasterRFclass)[[1]] <- data.frame(ID = 1:4, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta'))

gg <- ggplot() + 
  geom_spatraster(data = rasterRFclass) + 
  # scale_fill_manual(values = c('white', '#F9A03F', '#6BAED6', '#2171B5', '#238B45', '#D94801'), na.value = 'transparent', na.translate = F) +
  scale_fill_manual(values = c('white', 'white', '#3E5C23', '#C6E01B'), na.value = 'transparent', na.translate = F) +
  geom_sf(data = st_as_sf(terra::aggregate(wrld)), fill = NA, col = 'grey30') +
  coord_sf(xlim = c(-115, 170), ylim = c(-30, 30)) +
  labs(fill = '') +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 6),
    axis.text.x = element_text(hjust = 0.5, size = 6)
  )

gg
ggsave(plot = gg, filename = './png/maps/clustering/v3/v3_rf_map_ara-rob.jpg', units = 'in', width = 9, height =4, dpi = 300, create.dir = T)

## Check Vietnam
r <- rasterRFclass %>% terra::crop(vtnm) %>% terra::mask(vtnm)
plot(r, col = c('white', 'blue', 'yellow', 'green')); plot(vtnm, add = T, border = 'red')

## Probability
gp <- ggplot() + 
  geom_spatraster(data = rasterRFprob, aes(fill = ndd_hist )) + 
  scale_fill_gradientn(colors = cpt(pal = 'imagej_gyr_centre', n = 10, rev = TRUE), na.value = 'transparent') +
  geom_sf(data = st_as_sf(terra::aggregate(wrld)), fill = NA, col = 'grey30') +
  coord_sf(xlim = c(-115, 170), ylim = c(-30, 30)) +
  ggtitle(label = 'Suitability score [prob]') +
  labs(fill = '') +
  theme_bw() +
  theme(
    legend.key.width = unit(3, 'line'), legend.position = 'bottom', axis.text.y = element_text(angle = 90, hjust = 0.5, size = 6), axis.text.x = element_text(hjust = 0.5, size = 6)
  )
ggsave(plot = gp, filename = './png/maps/coffee/prob_map_ara-rob.jpg', units = 'in', width = 7, height = 4, dpi = 300, create.dir = T)

# Add limitations ---------------------------------------------------------
names(rasterRFprob) <- 'prob'
occ
vles <- as_tibble(terra::extract(rasterRFprob, occ[,c('Lon', 'Lat')]))
qnt5 <- as.data.frame(quantile(vles$prob, seq(0, 1, 0.01)))
thrp <- subset(qnt5, rownames(qnt5) == '5%') %>% as.numeric()

prob.rcl <- terra::ifel(rasterRFprob > thrp, 2, 0)
clst.rcl <- terra::ifel(as.numeric(rasterRFclass) %in% c(1, 2), 0, 1) %>% terra::crop(wrld) %>% terra::mask(wrld)
diff <- prob.rcl - clst.rcl
rslt <- rasterRFclass
rslt[diff == -1] <- 5
rslt[diff ==  2] <- 5
# levels(rslt)[[1]] <- data.frame(ID = 1:8, clase = c('Unsuitable', 'Unsuitable', paste0('Type ', 1:5), 'Limitations'))

plot(rslt)
levels(rslt)[[1]] <- data.frame(ID = 1:5, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta', 'Limitations'))

# Add mixed ---------------------------------------------------------------
vles <- as_tibble(terra::extract(rasterRFuncertainty, occ[,c('Lon', 'Lat')]))
colnames(vles)[2] <- 'prob'
qnt1 <- as.data.frame(quantile(vles$prob, seq(0, 1, 0.01)))
thru <- subset(qnt5, rownames(qnt1) == '1%') %>% as.numeric()

fnal <- rslt
fnal[rasterRFuncertainty < thru & rasterRFprob > thrp] <- 6
# levels(fnal)[[1]] <- data.frame(ID = 1:9, clase = c('Unsuitable', 'Unsuitable', paste0('Type ', 1:5), 'Limitations', 'Mixed'))
levels(fnal)[[1]] <- data.frame(ID = 1:6, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed'))
plot(fnal)

# Save thresholds ---------------------------------------------------------
dir.create('./rds/v6_coffee')
saveRDS(object = thrp, file = './rds/v6_coffee/threshol_prob_coffee.rds') #0.609926, # 0.661486
saveRDS(object = thru, file = './rds/v6_coffee/threshol_uncr_coffee.rds') #0.400889, # 0.5488172

## 0.524912
## 0.358360

# Save rasters ------------------------------------------------------------
dir_create('./rf/output/run_8_coffee/results/process')
terra::writeRaster(x = fnal, filename = './rf/output/run_8_coffee/results/process/rf-mixed_bsl.tif')

# Boxplot for the agroclimatic zones --------------------------------------
occ.cls <- clusteredpresdata %>% mutate(cluster = factor(cluster, levels = 1:5)) %>% gather(var, value, -c(Lon, Lat, Specie, cluster))
occ.cls <- filter(occ.cls, var != 'bioc_20_100')
occ.cls <- mutate(occ.cls, var = factor(var, levels = vrs))
g.box <- ggplot(data = occ.cls, aes(x = cluster, y = value)) + geom_boxplot() + facet_wrap(.~var, scales = 'free_y') + labs(x = 'Cluster', y = 'Value') + theme_bw() + theme(strip.text = element_text(face = 'bold'), axis.text.y = element_text(angle = 90, hjust = 0.5, size = 6), axis.text.x = element_text(size = 6))
ggsave(plot = g.box, filename = './png/graphs/coffee/boxplot_variables-cluster_arabica-robusta.jpg', units = 'in', width = 9, height = 7, dpi = 300, create.dir = T)

# Mapping AEZ -------------------------------------------------------------
levels(fnal)
# clrs <- c('#FFFFFF', '#F9A03F', '#6BAED6', '#2171B5', '#238B45', '#D94801', '#D9D9D9', '#FFFFB4')
clrs <- c('#FFFFFF',  '#3E5C23', '#C6E01B',  '#D9D9D9', '#FFFFB4')
# names(clrs) <- c('Unsuitable', 'Very hot - Wet', 'Very cold - Very dry', 'Cold - Very dry', 'Hot - Very wet', 'Hot - Dry', 'Limitations', 'Mixed')
names(clrs) <- c('Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed')
# levels(fnal)[[1]] <- data.frame(ID = 1:9, clase = c('Unsuitable', 'Unsuitable', 'Very hot - Wet', 'Very cold - Very dry', 'Cold - Very dry', 'Hot - Very wet', 'Hot - Dry', 'Limitations', 'Mixed'))
levels(fnal)[[1]] <- data.frame(ID = 1:6, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed'))

gg.aez <- ggplot() +
  geom_spatraster(data = fnal, aes(fill = clase)) + 
  scale_fill_manual(values = clrs, na.value = 'white', na.translate = FALSE) +
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') +
  coord_sf(ylim = c(-30, 30), xlim = c(-120, 130)) +
  ggtitle(label = 'Agroclimatic zones - Arabica / Robusta') +
  labs(fill = '') +
  theme_bw() +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, face = 'bold'), axis.text.y = element_text(angle = 90, hjust = 0.5))
gg.aez
ggsave(plot = gg.aez, filename = './png/maps/coffee/map_rf-mixed_bsl.jpg', units = 'in', width = 8, height = 4, dpi = 300, create.dir = TRUE)


# Presences and pseudo-absences -------------------------------------------

bckg <- read_csv('./tbl/points/coffee/arabica-robusta_background_swd-vars.csv')
bckg <- bckg %>% dplyr::select(x, y) %>% mutate(type = 'Background')
cffe <- cffe %>% dplyr::select(x = Lon, y = Lat, Specie) %>% rename(type = Specie)
allp <- rbind(cffe, bckg)
allp <- mutate(allp, type = factor(type, levels = c('arabica', 'robusta', 'Background')))

g.pnt <- ggplot() + 
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') + 
  geom_point(data = allp, aes(x = x, y = y, col = type), size = 0.2) + 
  coord_sf(xlim = c(-140, 150), y = c(-30, 30)) +
  labs(x = 'Lon', y = 'Lat') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6), 
    axis.text.y = element_text(size = 6)
  )

ggsave(plot = g.pnt, filename = './png/maps/coffee/arabica-robusta-background.jpg', units = 'in', width = 10, height = 4, dpi = 300)


vles <- terra::extract(bioc[[6]], allp[,c('x', 'y')])[,2]
allp <- mutate(allp, bioc_6 = vles)

g.b06 <- ggplot(data = allp, aes(x = type, y = bioc_6)) + 
  geom_boxplot() + 
  labs(x = '', y = 'Bioc 6') + 
  ggtitle(label = 'Bioc 6') + 
  theme_bw() + 
  theme()
g.b06
ggsave(plot = g.b06, filename = './png/maps/coffee/boxplot_bioc6-ara-rob-bck.jpg', units = 'in', width = 8, height = 7, dpi = 300, create.dir = T)
