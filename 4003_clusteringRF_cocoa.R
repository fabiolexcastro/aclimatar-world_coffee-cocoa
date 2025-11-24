
# Load libraries ----------------------------------------------------------
require(pacman)
p_load(terra, fs, sf, tidyverse, usdm, sp, raster, RColorBrewer, randomForest, outliers, readxl, hablar, glue, readxl, geodata, RColorBrewer)
pacman::p_load(tidyverse, raster, ggpubr, cptcity, rgdal, corrplot, cclust, dismo, gtools, sp, FactoMineR, pROC, randomForest, Hmisc)

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
make.mapp <- function(vari, clrs){
  
  # vari <- 'bioc_2'; clrs <- (brewer.pal(n = 9, name = 'YlOrRd'))
  
  ## To filter the table
  cat('To process: ', vari, '\n')
  bio <- bioc[[grep(paste0(vari, '$'), names(bioc))]]
  tbl <- terra::as.data.frame(bio, xy = T, na.rm = T)
  colnames(tbl)[3] <- 'value'
  
  ## Treshold 98
  hist(tbl$value)
  if(vari %in% c('bioc_13', 'bioc_14', 'bioc_18', 'bioc_19', 'bioc_22', 'bioc_27', 'bioc_29')){
    tbl <- mutate(tbl, value = ifelse(value >= quantile(tbl$value, 0.99, na.rm = T), quantile(tbl$value, 0.99, na.rm = T), value))  
  }
  
  ## To draw the map 
  ggm <- ggplot() + 
    geom_tile(data = tbl, aes(x = x, y = y, fill = value)) + 
    scale_fill_gradientn(colors = clrs, name = '') + 
    geom_sf(data = st_as_sf(cntn), fill = NA, col = 'grey30') +
    coord_sf(xlim = c(-130, 150), ylim = c(-25, 25)) +
    labs(title = vari, x = '', y = '', fill = '') + 
    theme_bw() +
    theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5), legend.key.width = unit(3, 'line'), axis.text = element_text(size = 6), legend.title.position = 'top', legend.title = element_text(hjust = 0.5))
  
  ## To return the map
  cat('Done!\n')
  return(ggm)
  
  
}

# Load data ---------------------------------------------------------------

## Vector -----------------------------
wrld <- geodata::world(resolution = 1, path = './tmpr')
cntn <- aggregate(wrld)

## Tabular data ----------------------
pnts <- read_csv('./tbl/points/cocoa/points_cocoa_all_rmv50_NoOTL.csv', show_col_types = FALSE)

## Raster data -----------------------
bioc <- terra::rast('./common_data/input_bios/bioc_hist.tif')
indx <- terra::rast('./common_data/atlas_hazards/hist_indices2.tif')
indx <- indx[[grep('hist', names(indx))]]
bioc <- c(bioc, indx)
bioc <- bioc[[-grep('bioc_20_60', names(bioc))]]

# Extract the values for the presences ------------------------------------
pnts <- as_tibble(cbind(pnts, terra::extract(bioc, pnts[,c('Lon', 'Lat')]))) 

# VIF analysis ------------------------------------------------------------
occ <- pnts
vif.res <- vif(x = as.data.frame(occ)[,4:ncol(occ)])
vif.step <- vifstep(x = as.data.frame(occ)[,4:ncol(occ)], th = 10)
vrs <- vif.step@results$Variables %>% as.character()
saveRDS(vrs, './rds/v3_cocoa/vars_cocoa.rds') #c("bioc_2", "bioc_3", "bioc_4", "bioc_13", "bioc_14", "bioc_18", "bioc_19" ,"bioc_20_100", "bioc_20_60", "bioc_22", "bioc_27", "bioc_31", "bioc_32", "ndd_hist", "n30_hist", "n35_hist", "ndw_hist")

## Corrplot
mtrx <- occ[,4:ncol(occ)]
corr <- as.matrix(mtrx)
png(filename = './png/graphs/v3_cocoa/corplot_cocoa.jpg', units = 'in', width = 9, height = 8, res = 300)
corrplot::corrplot(cor(corr))
dev.off()

# Mapping the variables ---------------------------------------------------
bioc.vars <- bioc[[grep(paste0(paste0(vrs, '$'), collapse = '|'), names(bioc), value = F)]]
names(bioc.vars)

maps.tasm <- map2(.x = c('bioc_2', 'bioc_3', 'bioc_4', 'bioc_32'), .y = replicate(4, brewer.pal(n = 9, name = "YlOrRd"), simplify = FALSE), .f = make.mapp)
maps.prec <- map2(.x = c('bioc_13', 'bioc_14', 'bioc_18', 'bioc_19'), .y = replicate(4, brewer.pal(n = 9, name = "BrBG"), simplify = FALSE), .f = make.mapp)
maps.etps <- map2(.x = c('bioc_20_100', 'bioc_22', 'bioc_27', 'bioc_30', 'bioc_31'), .y = replicate(5, rev(brewer.pal(n = 9, name = "BrBG")), simplify = FALSE), .f = make.mapp)
maps.ints <- map2(.x = c('n30_hist', 'n35_hist'), .y = replicate(2, (brewer.pal(n = 9, name = "YlOrRd")), simplify = FALSE), .f = make.mapp)
maps.ndws <- map2(.x = c('ndw_hist'), .y = replicate(1, (brewer.pal(n = 9, name = "BrBG")), simplify = FALSE), .f = make.mapp)
gmap <- cowplot::plot_grid(plotlist = c(maps.tasm, maps.prec, maps.etps, maps.ints, maps.ndws), ncol = 5, nrow = 4, common.legend = FALSE)
ggplot2::ggsave(plot = gmap, filename = './png/maps/cocoa/bios_bsl_vars_cocoa.jpg', width = 17, height = 11,  dpi = 300, bg = "white")

# To select the variables -------------------------------------------------
occ <- occ %>% dplyr::select(Lon, Lat, ID, vrs)

# Clustering --------------------------------------------------------------
env_values <- as.matrix(occ[,4:ncol(occ)]); nrow(env_values)
datRF <- as.data.frame(occ[,4:ncol(occ)]); nrow(datRF)
d <- dist(datRF, method = "euclidean")  
rfClust <- rf.clust(occ = occ, nforest = 25, ntrees = 100, nVars = 8, nclasses = 5)
labelRF <- rfClust[[1]]
clusterdata <- rfClust[[2]]
classdata <- cbind(pb = as.factor(labelRF), occ[,4:ncol(occ)])
clusteredpresdata <- cbind(occ, cluster = labelRF) %>% na.omit() %>% tbl_df()
no.clusters <- 5

dir.create('./rData/run_6_cocoa')
run <- 'run_6_cocoa'
save(datRF, file = paste0('./rData/', run, '/datRF.rData'))
save(clusterdata, file = paste0('./rData/', run, '/clusterdata.rData'))
save(occ, clusteredpresdata, no.clusters, labelRF, file = paste0('./rData/', run, '/clustereddata.rData'))
save(clusteredpresdata, file = './rData/run_1_arabica/presences.rData')

# Mapping clustering ------------------------------------------------------
clusteredpresdata
g.clst <- ggplot() + 
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') + 
  geom_point(data = clusteredpresdata, aes(x = Lon, y = Lat, col = factor(cluster))) + 
  # facet_wrap(.~cluster) +
  coord_sf(ylim = c(-30, 30), xlim = c(-130, 140)) +
  theme_bw() +
  theme(
    legend.position = 'bottom'
  )
g.clst
ggsave(plot = g.clst, filename = './png/maps/cocoa/rf_cluster5_cocoa.jpg', units = 'in', width = 9, height = 5, dpi = 300)

# Bias background sampling ------------------------------------------------

## Harvested area
faos <- read_csv('./tbl/faostat/FAOSTAT_data_en_11-11-2025.csv', show_col_types = FALSE)
faos <- faos %>% group_by(Area) %>% reframe(Value = mean(Value, na.rm = T)) %>% ungroup() %>% drop_na() %>% mutate(Value = round(Value, digits = 1)) %>% arrange(desc(Value))
write.csv(faos, './tbl/faostat/faostat_average_2000-2023.csv', row.names = FALSE)
faos <- read_csv('./tbl/faostat/faostat_average_2000-2023_iso.csv', show_col_types = FALSE)
faos <- full_join(st_as_sf(wrld), faos, by = c('GID_0' = 'ISO3'))

## Mapping Harvested Area
g1 <- ggplot() +
  geom_sf(data = faos, aes(fill = Value)) + 
  scale_fill_gradientn(labels = scales::comma, colors = RColorBrewer::brewer.pal(n = 9, name = 'YlOrBr'), breaks = seq(0, 3000000, by = 500000)) +
  geom_sf(data = st_as_sf(cntn), fill = NA, col = 'grey40') +
  coord_sf(xlim = c(-125, 150), ylim = c(-30, 30)) +
  ggtitle(label = 'Harvested area - Cocoa [Faostat]') +
  labs(fill = 'Hectareas', caption = 'Average 2000-2023') +
  theme_minimal() +
  theme(legend.title = element_text(hjust = 0.5), legend.title.position = 'top',
        legend.position = 'bottom', legend.key.width = unit(5, 'line'), legend.key.height = unit(0.5, 'line'))

ggsave(plot = g1, filename = './png/maps/cocoa/harvested_area_faostat.jpg', units = 'in', width = 8, height = 4, dpi = 300, create.dir = T)

## Presences
clusteredpresdata
frqn <- terra::extract(wrld, clusteredpresdata[,c('Lon', 'Lat')])$GID_0
frqn <- as.data.frame(table(frqn)) %>% arrange(desc(Freq)) %>% setNames(c('iso', 'frqn'))
faos <- full_join(faos, frqn, by = c('GID_0' = 'iso'))
faos <- faos %>% mutate(prob1 = frqn / Value, prob2 = Value / frqn)

# faos <- mutate(faos, prob1_2 = (prob1 / sum(prob1, na.rm = T))) # Version 1
faos <- mutate(faos, prob1_2 = prob1) # Version 2
faos %>% drop_na(prob1)
faos <- mutate(faos, prob1_2 = ifelse(is.na(prob1_2), 0.00001, prob1_2)) # min(faos$prob1_2, na.rm = T)

## Raster bias background sampling
mskr <- bioc[[1]] * 0 + 1
mskr <- terra::rasterize(x = faos, y = mskr, field = 'prob1_2')
bckg <- as_tibble(dismo::randomPoints(mask = raster::raster(mskr), n = nrow(clusteredpresdata), prob = TRUE))
plot(mskr)

## Count the presences by each country 
bckg <- mutate(bckg, iso = terra::extract(wrld, bckg[,c(1, 2)])$GID_0)
frqn.bckg <- bckg %>% group_by(iso) %>% reframe(count = n()) %>% arrange(desc(count))
faos <- full_join(faos, frqn.bckg, by = c('GID_0' = 'iso')) %>% rename(count_bckg = count, count_occr = frqn)

## Mapping the pseudo-absences
gg.faos <- g1
gg.mskr <- ggplot() + ggtitle(label = 'Probability for the background sampling') + geom_spatraster(data = mskr, aes(fill = prob1_2)) + scale_fill_viridis_c(na.value = 'transparent') + geom_sf(data = faos, fill = NA, col = 'grey30') + labs(fill = 'Prob bckg') + coord_sf(xlim = c(-120, 140), ylim = c(-30, 30)) + theme_minimal() + theme(legend.title = element_text(hjust = 0.5), legend.title.position = 'top', legend.position = 'bottom', legend.key.width = unit(5, 'line'), legend.key.height = unit(0.5, 'line'))
gg.bckg <- ggplot() + geom_point(data = bckg, aes(x = x, y = y), col = 'grey30', size = 0.1) + geom_sf(data = faos, fill = NA, col = 'grey30') + labs(x = '', y = '') + ggtitle(label = 'Background sampling') + coord_sf(xlim = c(-120, 140), ylim = c(-30, 30)) + theme_minimal() + theme(legend.title = element_text(hjust = 0.5), legend.title.position = 'top', legend.position = 'bottom', legend.key.width = unit(5, 'line'), legend.key.height = unit(0.5, 'line'))
gg.frbc <- ggplot() + ggtitle(label = 'Count background by country') + geom_sf(data = faos, aes(fill = count_bckg), col = 'grey30') + scale_fill_viridis_c() +coord_sf(xlim = c(-120, 140), ylim = c(-30, 30)) + labs(fill = 'Count background') + theme_minimal() + theme(legend.position = 'bottom', legend.key.width = unit(3, 'line'), legend.title.position = 'top', legend.title = element_text(hjust = 0.5))
gg.allp <- ggpubr::ggarrange(gg.faos, gg.mskr, gg.bckg, gg.frbc, ncol = 2, nrow = 2)
ggsave(plot = gg.allp, filename = './png/maps/cocoa/background_sampling_v2.jpg', units = 'in', width = 14, height = 7, dpi = 300, create.dir = T)

faos %>% 
  st_drop_geometry() %>% 
  as_tibble() %>% 
  dplyr::select(GID_0, NAME_0, Area, Value, count_occr, prob1_2, count_bckg) %>% 
  write.csv(., './tbl/points/cocoa/count_background-presences.csv', row.names = FALSE)

write.csv(bckg, './tbl/points/cocoa/background_v2.csv', row.names = FALSE)

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

##
back <- bckg
back <- back[,1:2]
back_sp <- back
coordinates(back) <- ~ x + y

##
back_swd <- terra::extract(bioc, back_sp[,1:2]) %>% cbind(back_sp, .)
back_swd <- as_tibble(back_swd)
back_swd <- dplyr::select(back_swd, x, y, vrs)
back_swd <- drop_na(back_swd)
back_swd %>% write.csv('./tbl/points/cocoa/cocoa_background_swd-vars_v2.csv', row.names = FALSE)

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

allclasses_swd <- rbind(classdata_2, presvalue_swd[,1:ncol(classdata_2)])
unique(allclasses_swd$pb)
dir_create('./tbl/rf/points_cocoa')
write.csv(allclasses_swd, './tbl/rf/points_cocoa/run-5_all_classes_swd_cocoa.csv', row.names = FALSE)

# To make the random forest analysis --------------------------------------
model1 <- as.formula(paste('factor(pb) ~', paste(paste(vrs), collapse = '+', sep =' ')))
rflist <- vector('list', 50) 
auc <- vector('list', 50)

allclasses_swd %>% pull(pb) %>% table()
NumberOfClusters <- 5

for(repe in 1:50){ # 50 bosques
  
  print(repe)
  pressample <- list()
  
  for (i in 1:(NumberOfClusters+no.absenceclasses)){
    
    if(any(i==c(1:no.absenceclasses))) { 
      
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), 
                     size = samplesize*NumberOfClusters/2/no.absenceclasses)
    } else {
      rows <- sample(rownames(allclasses_swd[allclasses_swd$pb==i,]), size=samplesize)
    }
    pressample[[i]] <- allclasses_swd[rows,] 
  }
  
  species <- na.omit(do.call(rbind, pressample)) 
  head(species)
  Samplesplit <- sample(rownames(species)) 
  
  envtrain <- species[Samplesplit[1:(0.8*nrow(species))],] 
  envtest <- species[Samplesplit[(0.8*nrow(species)):nrow(species)],] 
  
  rfmodel <- randomForest(model1, data = envtrain, ntree = 500, na.action = na.omit, nodesize = 2) 
  
  save(rfmodel, file = paste('./rf/output/run_6_cocoa/models/', NumberOfClusters, 'Prob_' , 'rep_' ,repe, '.rdata' ,sep=''))
  rflist[[repe]] <- rfmodel
  
  # AUC 
  predicted <- as.numeric(predict(rfmodel, envtest))
  observed <- as.vector(envtest[,'pb'])
  auc[[repe]] <- auc(observed, predicted) 
  rm(rfmodel)
  
  cat(auc[[repe]] ,'\n')
  
}

auc <- unlist(auc)
boxplot(auc)
rff <- do.call(randomForest::combine, rflist)
importance <- as.data.frame(rff$importance)

# Predict modell
lyr <- bioc[[grep(paste0(paste0(vrs, '$'), collapse = '|'), names(bioc), value = F)]]
climatevalues  <- data.frame(values(lyr))
NumberOfClusters <- 5

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
dir_create('./rf/output/run_7_cocoa/results/raw')
writeRaster(rasterRFclass, paste0('./rf/output/run_7_cocoa/results/raw/RF_5Clust_current.tif'), overwrite = T)
writeRaster(rasterRFprob, paste0('./rf/output/run_7_cocoa/results/raw/RF_5Prob_current.tif'), overwrite = T)
writeRaster(rasterRFuncertainty, paste0('./rf/output/run_7_cocoa/results/raw/RF_5Unc_current.tif'), overwrite = T)

# Mapping -----------------------------------------------------------------

## Clustering
rasterRFclass <- as.factor(rasterRFclass)
levels(rasterRFclass)[[1]] <- data.frame(ID = 1:7, clase = c('Unsuitable', 'Unsuitable', paste0('Type ', 1:5)))
# levels(rasterRFclass)[[1]] <- data.frame(ID = 1:4, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta'))

gg <- ggplot() + 
  geom_spatraster(data = rasterRFclass) + 
  scale_fill_manual(values = c('white', '#F9A03F', '#6BAED6', '#2171B5', '#238B45', '#D94801'), na.value = 'transparent', na.translate = F) +
  # scale_fill_manual(values = c('white', 'white', '#3E5C23', '#C6E01B'), na.value = 'transparent', na.translate = F) +
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
ggsave(plot = gg, filename = './png/maps/cocoa/v2/rf_map_cocoa_raw.jpg', units = 'in', width = 9, height =4, dpi = 300, create.dir = T)

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
rslt[diff == -1] <- 8
rslt[diff ==  2] <- 8
levels(rslt)[[1]] <- data.frame(ID = 1:8, clase = c('Unsuitable', 'Unsuitable', paste0('Type ', 1:5), 'Limitations'))
# levels(rslt)[[1]] <- data.frame(ID = 1:5, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta', 'Limitations'))

# Add mixed ---------------------------------------------------------------
vles <- as_tibble(terra::extract(rasterRFuncertainty, occ[,c('Lon', 'Lat')]))
colnames(vles)[2] <- 'prob'
qnt1 <- as.data.frame(quantile(vles$prob, seq(0, 1, 0.01)))
thru <- subset(qnt5, rownames(qnt1) == '1%') %>% as.numeric()

fnal <- rslt
fnal[rasterRFuncertainty < thru & rasterRFprob > thrp] <- 9
levels(fnal)[[1]] <- data.frame(ID = 1:9, clase = c('Unsuitable', 'Unsuitable', paste0('Type ', 1:5), 'Limitations', 'Mixed'))
# levels(fnal)[[1]] <- data.frame(ID = 1:6, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed'))
plot(fnal)

# Save thresholds ---------------------------------------------------------
dir.create('./rds/v5_cocoa')
saveRDS(object = thrp, file = './rds/v5_cocoa/threshol_prob_cocoa.rds') # 0.539626
saveRDS(object = thru, file = './rds/v5_cocoa/threshol_uncr_cocoa.rds') # 0.4141716

# Cocoa

# Save rasters ------------------------------------------------------------
dir_create('./rf/output/run_7_cocoa/results/process')
terra::writeRaster(x = fnal, filename = './rf/output/run_7_cocoa/results/process/rf-mixed_bsl.tif', overwrite = TRUE)

# Boxplot for the agroclimatic zones --------------------------------------
occ.cls <- clusteredpresdata %>% mutate(cluster = factor(cluster, levels = 1:5)) %>% 
  gather(var, value, -c(Lon, Lat, ID, cluster))
occ.cls <- mutate(occ.cls, var = factor(var, levels = vrs))
g.box <- ggplot(data = occ.cls, aes(x = cluster, y = value)) + geom_boxplot() + facet_wrap(.~var, scales = 'free_y') + labs(x = 'Cluster', y = 'Value') + theme_bw() + theme(strip.text = element_text(face = 'bold'), axis.text.y = element_text(angle = 90, hjust = 0.5, size = 6), axis.text.x = element_text(size = 6))
ggsave(plot = g.box, filename = './png/graphs/cocoa/boxplot_variables-cluster_cocoa.jpg', units = 'in', width = 9, height = 7, dpi = 300, create.dir = T)

# Mapping AEZ -------------------------------------------------------------
levels(fnal)
# clrs <- c('#FFFFFF', '#F9A03F', '#6BAED6', '#2171B5', '#238B45', '#D94801', '#D9D9D9', '#FFFFB4')
# clrs <- c('#FFFFFF',  '#3E5C23', '#C6E01B',  '#D9D9D9', '#FFFFB4')
clrs <- c('#FFFFFF', "#F4BFBF", "#A3C4F3", "#FFD6A5", "#BDE0FE", "#C8E6C9", '#D9D9D9', '#FFFFB4')
names(clrs) <- c('Unsuitable', 'Constant – Very dry', 'Cold - Wet', 'Very hot - Wet', 'Hot - Very wet', 'Very cold - Very dry', 'Limitations', 'Mixed')


clrs <- c('#FFFFFF', '#F7A399', '#F9C784', '#A8E6CF', '#A0C4FF', '#CDB4DB',  '#D9D9D9', '#FFFFB4')
names(clrs) <- c('Unsuitable', 'Hot - Very Dry', ' Very hot - Dry', 'Hot - Wet', 'Cold - Very Wet', 'Very cold - Dry', 'Limitations', 'Mixed')
# names(clrs) <- c('Unsuitable', 'Very hot - Wet', 'Very cold - Very dry', 'Cold - Very dry', 'Hot - Very wet', 'Hot - Dry', 'Limitations', 'Mixed')
# names(clrs) <- c('Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed')
# levels(fnal)[[1]] <- data.frame(ID = 1:9, clase = c('Unsuitable', 'Unsuitable', 'Very hot - Wet', 'Very cold - Very dry', 'Cold - Very dry', 'Hot - Very wet', 'Hot - Dry', 'Limitations', 'Mixed'))
levels(fnal)[[1]] <- data.frame(ID = 1:9, clase = c('Unsuitable', 'Unsuitable', 'Hot - Very Dry', ' Very hot - Dry', 'Hot - Wet', 'Cold - Very Wet', 'Very cold - Dry', 'Limitations', 'Mixed'))
levels(fnal)[[1]] <- data.frame(ID = 1:9, clase = c('Unsuitable', 'Unsuitable', 'Constant – Very dry', 'Cold - Wet', 'Very hot - Wet', 'Hot - Very wet', 'Very cold - Very dry', 'Limitations', 'Mixed'))

gg.aez <- ggplot() +
  geom_spatraster(data = fnal, aes(fill = clase)) + 
  scale_fill_manual(values = clrs, na.value = 'white', na.translate = FALSE) +
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') +
  coord_sf(ylim = c(-30, 30), xlim = c(-120, 150)) +
  ggtitle(label = 'Agroclimatic zones - Cocoa') +
  labs(fill = '') +
  theme_bw() +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, face = 'bold'), axis.text.y = element_text(angle = 90, hjust = 0.5))
gg.aez
ggsave(plot = gg.aez, filename = './png/maps/cocoa/v2/map_rf-mixed_bsl.jpg', units = 'in', width = 8, height = 4, dpi = 300, create.dir = TRUE)


