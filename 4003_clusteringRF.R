

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, usdm, sp, raster, randomForest, outliers, readxl, hablar, glue, readxl, geodata, RColorBrewer)
pacman::p_load(tidyverse, raster, ggpubr, rgdal, corrplot, cclust, dismo, gtools, sp, FactoMineR, pROC, randomForest, Hmisc)

mg <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

source('https://raw.githubusercontent.com/fabiolexcastro/rfSalvador/refs/heads/master/FunctionsRFclustering.R')

# Function ----------------------------------------------------------------
rf.clust <- function(occ, nforest, ntrees, nVars, nclasses){
  
  datRF_presences <- occ[,3:ncol(occ)]
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
plot(wrld)

## Tabular data ----------------------
cffe <- read_csv('./tbl/points/clean/arabica_robusta_otl_swd.csv', show_col_types = FALSE)

cffe.ara <- filter(cffe, Specie == 'arabica')
cffe.rob <- filter(cffe, Specie == 'robusta')
 
## Sampling arabica ## Add arabica and robusta into only one table
cffe.ara <- cffe.ara %>% sample_n(size = 600, replace = FALSE, tbl = .)
cffe <- rbind(cffe.ara, cffe.rob)

## Raster data -----------------------
bioc <- terra::rast('./common_data/input_bios/bioc_hist.tif')
indx <- terra::rast('./common_data/atlas_hazards/hist_indices2.tif')
indx <- indx[[grep('hist', names(indx))]]
bioc <- c(bioc, indx)
bioc <- bioc[[-grep('100', names(bioc))]]

# Add the indices variables -----------------------------------------------
cffe <- cbind(
  cffe, 
  terra::extract(indx, cffe[,c('Lon', 'Lat')])[,-1]
) %>% 
  as_tibble()

## Boxplot indices 
gg.box <- cffe %>% dplyr::select(Lon, Lat, Specie, ends_with('hist')) %>% gather(var, value, -c(Lon, Lat, Specie)) %>% separate(col = 'var', into = c('index', 'period'), sep = '_') %>% 
  ggplot(data = ., aes(x = Specie, y = value, fill = Specie)) + 
  geom_boxplot() + 
  facet_wrap(.~index, scales = 'free_y') + 
  labs(title = 'Boxplot Indices - Baseline', y = 'Value', x = 'Specie') +
  theme_bw() +
  theme(
    plot.title = element_text(face = 'bold', hjus = 0.5),
    strip.text = element_text(face = 'bold', hjust = 0.5)
  )

ggsave(plot = gg.box, filename = './png/graphs/boxplot_indices_arabica-robusta.jpg', units = 'in', width = 8, height = 6, dpi = 300, create.dir = TRUE)

# VIF analysis ------------------------------------------------------------
occ <- cffe
occ <- drop_na(occ)
vif.res <- vif(x = as.data.frame(occ)[,5:ncol(occ)])
vif.step <- vifstep(x = as.data.frame(occ)[,5:ncol(occ)], th = 10)
vrs <- vif.step@results$Variables %>% as.character()
saveRDS(vrs, './rds/vars_coffee_arabica-robusta.rds')

##
mtrx <- occ[,5:ncol(occ)]
corr <- as.matrix(mtrx)
png(filename = './png/graphs/corplot_arabica-robusta.jpg', units = 'in', width = 9, height = 8, res = 300)
corrplot::corrplot(cor(corr))
dev.off()

# Mapping the variables ---------------------------------------------------
bioc.vars <- bioc[[grep(paste0(paste0(vrs, '$'), collapse = '|'), names(bioc), value = F)]]
bioc.vars
plot(bioc.vars)

bioc.vars <- bioc.vars %>% terra::as.data.frame(xy = T) %>% as_tibble() %>% gather(var, value, -c(x, y))
unique(bioc.vars$var)

bioc.vars.tasm <- bioc.vars %>% filter(var %in% c('bioc_2', 'bioc_3', 'bioc_8')) ## Temperature 
bioc.vars.prec <- bioc.vars %>% filter(var %in% c('bioc_13', 'bioc_14', 'bioc_18', 'bioc_19'))
bioc.vars.etps <- bioc.vars %>% filter(var %in% c('bioc_22', 'bioc_27', 'bioc_29'))
bioc.vars.baln <- bioc.vars %>% filter(var %in% c('bioc_30', 'bioc_31'))
bioc.vars.idts <- bioc.vars %>% filter(var %in% c('n30_hist', 'n35_hist'))

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
maps.tasm <- map2(.x = c('bioc_2', 'bioc_3', 'bioc_8'), .y = replicate(3, brewer.pal(n = 9, name = "YlOrRd"), simplify = FALSE), .f = make.mapp)
maps.prec <- map2(.x = c('bioc_13', 'bioc_14', 'bioc_18', 'bioc_19'), .y = replicate(4, brewer.pal(n = 9, name = "BrBG"), simplify = FALSE), .f = make.mapp)
maps.etps <- map2(.x = c('bioc_22', 'bioc_27', 'bioc_29'), .y = replicate(3, rev(brewer.pal(n = 9, name = "BrBG")), simplify = FALSE), .f = make.mapp)
maps.baln <- map2(.x = c('bioc_30', 'bioc_31'), .y = replicate(2, rev(brewer.pal(n = 9, name = "BrBG")), simplify = FALSE), .f = make.mapp)
maps.ints <- map2(.x = c('n30_hist', 'n35_hist'), .y = replicate(2, (brewer.pal(n = 9, name = "YlOrRd")), simplify = FALSE), .f = make.mapp)
maps.ndds <- map2(.x = c('ndd_hist'), .y = replicate(1, rev(brewer.pal(n = 9, name = "BrBG")), simplify = FALSE), .f = make.mapp)
maps
gmap <- cowplot::plot_grid(plotlist = c(maps.tasm, maps.prec, maps.etps, maps.baln, maps.ints, maps.ndds), ncol = 5, nrow = 4, common.legend = FALSE)
ggplot2::ggsave(plot = gmap, filename = './png/maps/bios_bsl_vars_ara-rob.jpg', width = 17, height = 11,  dpi = 300, bg = "white")

# To select the variables -------------------------------------------------
occ <- occ %>% dplyr::select(Lon, Lat, Specie, vrs)

# Draw mapping ------------------------------------------------------------
g.cff.ara <- ggplot() +
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') +
  geom_point(data = occ, aes(x = Lon, y = Lat, col = Specie), size = 0.3) +
  coord_sf(ylim = c(-30, 30), xlim = c(-110, 135)) +
  ggtitle(label = 'C. arabica') +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        legend.position = 'bottom',
        text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
g.cff.ara
ggsave(plot = g.cff.ara, filename = './png/maps/points/arabica-robusta_world.jpg', units = 'in', width = 8, height = 3, dpi = 300, create.dir = T)



# Sampling for arabica ----------------------------------------------------
write.csv(occ, './tbl/points/clean/arabica-robusta_otl_swd_sampling.csv', row.names = FALSE)

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

dir.create('./rData/run_1_arabica-robusta')
run <- 'run_1_arabica-robusta'
save(datRF, file = paste0('./rData/', run, '/datRF.rData'))
save(clusterdata, file = paste0('./rData/', run, '/clusterdata.rData'))
save(occ, clusteredpresdata, no.clusters, labelRF, file = paste0('./rData/', run, '/clustereddata.rData'))
save(clusteredpresdata, file = './rData/run_1_arabica/presences.rData')


# Mapping clustering ------------------------------------------------------
clusteredpresdata
g.clst <- ggplot() + 
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') + 
  geom_point(data = clusteredpresdata, aes(x = Lon, y = Lat, col = factor(cluster))) + 
  facet_wrap(.~cluster) +
  coord_sf(ylim = c(-30, 30), xlim = c(-130, 140)) +
  theme_bw() +
  theme(
    legend.position = 'bottom'
  )

ggsave(plot = g.clst, filename = './png/maps/clustering/rf_cluster5_ara-rob.jpg', units = 'in', width = 9, height = 5, dpi = 300)




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
ggsave(plot = gcls, filename = './png/graphs/count_rf-cluster_arabica-robusta.jpg', units = 'in', width = 6, height = 5, dpi = 300, create.dir = T)


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
back <- randomPoints(raster(back_raster), 1*numberofpresences) %>% as_data_frame()
back_sp <- back
coordinates(back) <- ~ x + y
back_swd <- terra::extract(bioc, back_sp) %>% cbind(back_sp, .)
back_swd <- as_tibble(back_swd)
back_swd <- dplyr::select(back_swd, x, y, vrs)
back_swd %>% write.csv('./tbl/points/clean/arabica-robusta_background_swd-vars.csv', row.names = FALSE)

# Cluster analysis to pseudoabsences
bckclust <- rf.clust(occ = back_swd, nforest = 50, ntrees = 500, nVars = 8, nclasses = 2)
datRF <- as.data.frame(back_swd[,3:ncol(back_swd)])
attach(datRF)
no.forests <- 50#raw = 25
no.trees <- 500
distRF <- RFdist(datRF, mtry1 = 8, no.trees, no.forests, addcl1 = T, addcl2 = F, imp =T, oob.prox1 = T)# mtry1 = 4 raw  # es la cantidad de variables a utilizar en cada no
no.absenceclasses <- 2
labelRF <- pamNew(distRF$cl1, no.absenceclasses)
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

dim(classdata_2); dim(presvalue_swd)
# presvalue_swd <- presvalue_swd %>% dplyr::select(-cluster)
