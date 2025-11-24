

# Load libraries ----------------------------------------------------------
require(pacman)
p_load(climateStability,cclust,corrplot,cptcity,dismo,FactoMineR,fs,ggpubr,glue,geodata,gtools,hablar,Hmisc,openxlsx,outliers,parallelDist,pROC,randomForest,raster,RColorBrewer,rnaturalearth,rnaturalearthdata,Rfast,scales,sf,sp,spatialEco,spatstat,terra,tidyverse,usdm)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)


# Function ----------------------------------------------------------------
make.map.var <- function(vari, clrs){
  
  # vari <- 'bioc_3'; clrs <- brewer.pal(n = 9, name = 'YlOrRd')
  
  ## Filtering 
  cat('To process: ', vari, '\n')
  rst <- bioc[[grep(paste0(vari, '$'), names(bioc))]]
  
  ## To mapping 
  ggm <- ggplot() + 
    geom_spatraster(data = rst, aes_string(fill = vari)) + 
    scale_fill_gradientn(colors = clrs, na.value = 'transparent') +
    geom_sf(data = wrld, fill = NA, col = 'grey30') +
    coord_sf(xlim = ext(rst)[1:2], ylim = ext(rst)[3:4]) +
    labs(fill = '', title = vari) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 6), plot.title = element_text(hjust = 0.5, face = 'bold'), legend.title = element_text(hjust = 0.5), legend.title.position = 'top', legend.position = 'bottom', legend.key.width = unit(3, 'line'), legend.key.height = unit(0.5, 'line')
    )
  
  ## Finish 
  return(ggm)
  
}
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
source('https://raw.githubusercontent.com/fabiolexcastro/rfSalvador/refs/heads/master/FunctionsRFclustering.R')

# Load data ---------------------------------------------------------------

## Vector -----------------------------
wrld <- geodata::world(resolution = 1, path = './tmpr')
cntn <- aggregate(wrld)
peru <- geodata::gadm(country = 'PER', level = 1, path = './tmpr')
mexi <- geodata::gadm(country = 'MEX', level = 1, path = './tmpr')

## Tabular ----------------------------
cffe <- read_csv('./tbl/points/coffee/arabica_robusta_otl_swd.csv', show_col_types = FALSE)
cffe <- dplyr::select(cffe, -bioc_20_60)

## Raster data ------------------------
bioc <- terra::rast('./common_data/input_bios/bioc_hist.tif')
indx <- terra::rast('./common_data/atlas_hazards/hist_indices2.tif')
indx <- indx[[grep('hist', names(indx))]]
bioc <- c(bioc, indx)
bioc <- bioc[[-grep('bioc_20_100', names(bioc))]]

# Extract the values ------------------------------------------------------
cffe <- dplyr::select(cffe, Lon, Lat, Specie)
cffe <- cbind(cffe, terra::extract(bioc, cffe[,c('Lon', 'Lat')]))
cffe <- as_tibble(cffe)

# Check points  -----------------------------------------------------------
cffe <- mutate(cffe, iso = terra::extract(wrld, cffe[,c('Lon', 'Lat')])$NAME_0)
cffe.per <- filter(cffe, iso == 'Peru')
cffe.mex <- filter(cffe, iso == 'Mexico')

cffe.per <- mutate(cffe.per, name_1 = terra::extract(peru, cffe.per[,c('Lon', 'Lat')])$NAME_1)
cffe.mex <- mutate(cffe.mex, name_1 = terra::extract(mexi, cffe.mex[,c('Lon', 'Lat')])$NAME_1)

## Remove from the main table
cffe.mex <- filter(cffe.mex, name_1 != 'Yucatán')
cffe.per <- filter(cffe.per, !name_1 %in% c('Loreto', 'Lima', 'Lima Province'))
cffe.oth <- filter(cffe, !iso %in% c('Peru', 'Mexico'))

cffe.mex <- dplyr::select(cffe.mex, Lon, Lat, Specie)
cffe.per <- dplyr::select(cffe.per, Lon, Lat, Specie)
cffe.oth <- dplyr::select(cffe.oth, Lon, Lat, Specie)

cffe.fnl <- rbind(cffe.oth, cffe.mex, cffe.per)

## Adding
gg <- ggplot() +
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') + 
  geom_point(data = cffe.fnl, aes(x = Lon, y = Lat, col = Specie), size = 0.05) +
  coord_sf(xlim = c(-120, 140), ylim = c(-30, 30)) +
  theme_minimal() +
  ggtitle(label = 'Location coffee points') +
  theme(
    legend.position = 'bottom'
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

gg
ggsave(plot = gg, filename = './workspace_coffee/png/maps/run_2/coffee_points.jpg', units = 'in', width = 7, height = 4, dpi = 300, create.dir = T)
write.csv(cffe.fnl, './workspace_coffee/tbl/points_v2.csv', row.names = FALSE)

cffe <- cffe.fnl

# Extract value for the points --------------------------------------------
cffe <- cbind(cffe, terra::extract(bioc, cffe[,c('Lon', 'Lat')]))
write.csv(cffe, './workspace_coffee/tbl/points_v2.csv', row.names = FALSE)

# VIF ---------------------------------------------------------------------
occ <- cffe %>% drop_na()
vars <- map_dfr(.x = c(2.5, 5, 10, 15), .f = function(i){
  vif.res <- vifstep(x = as.data.frame(occ)[5:ncol(occ)], th = i, keep = 'bioc_6')
  vrs <- vif.res@results$Variables %>% as.character()
  rsl <- tibble(threshold = i, variable = vrs)
  return(rsl)
}) %>% filter(threshold == 10) %>% pull(2)
vars <- c('ndd_hist', 'n30_hist', 'n35_hist', 'ndw_hist', 'bioc_3', 'bioc_4', 'bioc_6',
          'bioc_13', 'bioc_14', 'bioc_18', 'bioc_19', 'bioc_20_60', 'bioc_22', 'bioc_23',
          'bioc_24', 'bioc_27', 'bioc_29', 'bioc_31', 'bioc_32')
dir_create('./workspace_coffee/rds/run_1')
# saveRDS(object = vars, file = './workspace_coffee/rds/run_1/vars.rds')
# vars <- readRDS(file = './workspace_coffee/rds/run_1/vars.rds')

# Mapping variables -------------------------------------------------------
bioc <- bioc[[grep(paste0(paste0(vars, '$'), collapse = '|'), names(bioc), value = F)]]

## Temperature
ggts <- map2(.x = c('bioc_3', 'bioc_4', 'bioc_6', 'bioc_32', 'n30_hist', 'n35_hist'), .y = replicate(6, brewer.pal(n = 9, name = "YlOrRd"), simplify = F), .f = make.map.var)
ggts <- ggpubr::ggarrange(ggts[[1]], ggts[[2]], ggts[[3]], ggts[[4]], ggts[[5]], ggts[[6]], ncol = 3, nrow = 2)
ggsave(plot = ggts, filename = './workspace_coffee/png/maps/vars_tas_run1.jpg', units = 'in', width = 11, height = 4, dpi = 300, create.dir = T)

## Precipitation
ggpp <- map2(.x = c('bioc_13', 'bioc_14', 'bioc_18', 'bioc_19'), .y = replicate(4, brewer.pal(n = 9, name = "BrBG"), simplify = F),  .f = make.map.var)
ggpp <- ggpubr::ggarrange(ggpp[[1]], ggpp[[2]], ggpp[[3]], ggpp[[4]], ncol = 2, nrow = 2)
ggsave(plot = ggpp, filename = './workspace_coffee/png/maps/vars_ppt_run1.jpg', units = 'in', width = 11, height = 4, dpi = 300, create.dir = T)

## Evapotranspiration
gget <- map2(.x = c('bioc_20_60', 'bioc_22', 'bioc_23', 'bioc_24', 'bioc_27', 'bioc_29', 'bioc_31', 'ndw_hist', 'ndd_hist'), .y = replicate(9, brewer.pal(n = 9, name = "BrBG"), simplify = F), .f = make.map.var)
gget <- ggpubr::ggarrange(gget[[1]], gget[[2]], gget[[3]], gget[[4]], gget[[5]], gget[[6]], gget[[7]], gget[[8]], gget[[9]], ncol = 4, nrow = 3)
ggsave(plot = gget, filename = './workspace_coffee/png/maps/vars_etp_run1.jpg', units = 'in', width = 12, height = 7, dpi = 300, create.dir = T)

# Clustering --------------------------------------------------------------

##
occ <- occ %>% dplyr::select(Lon, Lat, Specie, ID, vars)

##
env_values <- as.matrix(occ[,5:ncol(occ)]); nrow(env_values)
datRF <- as.data.frame(occ[,5:ncol(occ)]); nrow(datRF)
d <- dist(datRF, method = "euclidean")  
rfClust <- rf.clust(occ = occ, nforest = 25, ntrees = 100, nVars = 8, nclasses = 5)
labelRF <- rfClust[[1]]
clusterdata <- rfClust[[2]]
classdata <- cbind(pb = as.factor(labelRF), occ[,3:ncol(occ)])
clusteredpresdata <- cbind(occ, cluster = labelRF) %>% na.omit() %>% tbl_df()
no.clusters <- 5
clusteredpresdata <- clusteredpresdata %>% mutate(cluster = ifelse(Specie == 'arabica', 1, 2))

dir_ls('./workspace_coffee/rData')
dir_create('./workspace_coffee/rData/run_2')
run <- 'run_2'
save(datRF, file = paste0('./workspace_coffee/rData/', run, '/datRF.rData'))
save(clusterdata, file = paste0('./workspace_coffee/rData/', run, '/clusterdata.rData'))
save(occ, clusteredpresdata, no.clusters, labelRF, file = paste0('./workspace_coffee/rData/', run, '/clustereddata.rData'))
save(clusteredpresdata, file = './workspace_coffee/rData/run_1/presences.rData')

# Mapping clustering data -------------------------------------------------
g.clst <- ggplot() + 
  geom_sf(data = st_as_sf(wrld), fill = NA, col = 'grey30') + 
  geom_point(data = clusteredpresdata, aes(x = Lon, y = Lat, col = factor(cluster)), size = 0.2) + 
  scale_color_manual(labels = c('arabica', 'robusta'), values = c('#3E5C23', '#C6E01B')) +
  labs(fill = '', col = '') +
  ggtitle(label = 'Clustering presences - Random Forest') +
  # facet_wrap(.~cluster) +
  coord_sf(ylim = c(-30, 30), xlim = c(-130, 140)) +
  theme_bw() +
  theme(
    legend.position = 'bottom'
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(plot = g.clst, filename = './workspace_coffee/png/maps/run_2/rf_cluster_ara-rob.jpg', units = 'in', width = 9, height = 3, dpi = 300)

# Bias - background sampling ----------------------------------------------
bckg <- read_csv('./tbl/points/coffee/arabica-robusta_background_swd-vars.csv', show_col_types = FALSE)[,c('x', 'y')]
bckg <- cbind(bckg, terra::extract(bioc, bckg[,1:2]))
bckg <- dplyr::select(bckg, -ID)
bckg <- as_tibble(bckg)
back_swd <- bckg

samplesize <- round(min(summary(as.factor(clusteredpresdata$cluster))) / 2, 0) 

# Cluster analysis to pseudoabsences
bckclust <- rf.clust(occ = back_swd, nforest = 50, ntrees = 500, nVars = 8, nclasses = 2)
datRF <- as.data.frame(back_swd[,3:ncol(back_swd)])
attach(datRF)
no.forests <- 50#raw = 25
no.trees <- 500
distRF <- RFdist(datRF, mtry1 = 8, no.trees, no.forests, addcl1 = T, addcl2 = F, imp =T, oob.prox1 = T)# mtry1 = 4 raw  # es la cantidad de variables a utilizar en cada no
no.absenceclasses <- 2
labelRF <- pamNew(distRF$cl1, no.absenceclasses)
# labelRF <- c(rep(1, 500), rep(2, 958))
detach(datRF)
classdata <- cbind(pb = as.factor(labelRF), back_swd[,3:ncol(back_swd)])

presvalue_swd  <- clusteredpresdata[,4:ncol(clusteredpresdata)] %>%
  cbind(pb = (clusteredpresdata$cluster + no.absenceclasses), .) %>%
  na.omit() %>%
  as.data.frame() %>%
  mutate(cluster = cluster + no.absenceclasses)
presvalue_swd <- dplyr::select(presvalue_swd, pb, vars)
presvalue_swd <- mutate(presvalue_swd, pb = as.factor(pb))
classdata_2 <- cbind(pb = as.data.frame(classdata)$pb, classdata[,2:ncol(classdata)]) # Background
# classdata_2 <- dplyr::select(classdata_2, -Specie)
dim(classdata_2); dim(presvalue_swd)

allclasses_swd <- rbind(classdata_2, presvalue_swd[,1:ncol(classdata_2)])
unique(allclasses_swd$pb)
dir_create('./workspace_coffee/tbl/rf/points')
write.csv(allclasses_swd, './workspace_coffee/tbl/rf/points/run-2_all_classes_swd.csv', row.names = FALSE)
write.csv(back_swd, './workspace_coffee/tbl/rf/points/run-2_background.csv', row.names = FALSE)
write.csv(clusteredpresdata, './workspace_coffee/tbl/rf/points/run-2_presences.csv', row.names = FALSE)

# To make the random forest analysis --------------------------------------
model1 <- as.formula(paste('factor(pb) ~', paste(paste(vars), collapse = '+', sep =' ')))
rflist <- vector('list', 50) 
auc <- vector('list', 50)

allclasses_swd %>% pull(pb) %>% table()

for(repe in 1:50){ # 50 bosques
  
  print(repe)
  pressample <- list()
  NumberOfClusters <- 2
  
  for (i in 1:(NumberOfClusters+no.absenceclasses)){
    
    cat(i, ' model', '\n')  
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
  
  # dir_create('./workspace_coffee/rf/output/run_2/models')
  NumberOfClusters <- 2
  save(rfmodel, file = paste('./workspace_coffee/rf/output/run_2/models/', NumberOfClusters, 'Prob_' , 'rep_' ,repe, '.rdata' ,sep=''))
  rflist[[repe]] <- rfmodel
  
  # AUC 
  predicted <- as.numeric(predict(rfmodel, envtest))
  observed <- as.vector(envtest[,'pb'])
  auc[[repe]] <- pROC::auc(as.numeric(observed), as.numeric(predicted))
  rm(rfmodel)
  
  cat(auc[[repe]] ,'\n')
  
}

auc <- unlist(auc)
boxplot(auc)
rff <- do.call(randomForest::combine, rflist)
importance <- as.data.frame(rff$importance)
dir_create('./workspace_coffee/rds/run_2')
saveRDS(object = rff, file = './workspace_coffee/rds/run_2/rff.rds')

# Predict modell
lyr <- bioc[[grep(paste0(paste0(vars, '$'), collapse = '|'), names(bioc), value = F)]]
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
dir_create('./workspace_coffee/rf/output/run_2/results/raw')
writeRaster(rasterRFclass, paste0('./workspace_coffee/rf/output/run_2/results/raw/RF_2Clust_current.tif'), overwrite = T)
writeRaster(rasterRFprob, paste0('./workspace_coffee/rf/output/run_2/results/raw/RF_2Prob_current.tif'), overwrite = T)
writeRaster(rasterRFuncertainty, paste0('./workspace_coffee/rf/output/run_2/results/raw/RF_2Unc_current.tif'), overwrite = T)

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
levels(rslt)[[1]] <- data.frame(ID = 1:5, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta', 'Limitations'))

plot(rslt)

# Add mixed ---------------------------------------------------------------
vles <- as_tibble(terra::extract(rasterRFuncertainty, occ[,c('Lon', 'Lat')]))
colnames(vles)[2] <- 'prob'
qnt1 <- as.data.frame(quantile(vles$prob, seq(0, 1, 0.01)))
thru <- subset(qnt5, rownames(qnt1) == '1%') %>% as.numeric()

fnal <- rslt
fnal[rasterRFuncertainty < thru & rasterRFprob > thrp] <- 6
levels(fnal)[[1]] <- data.frame(ID = 1:6, clase = c('Unsuitable', 'Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed'))
plot(fnal)

# Save thresholds / raster -------------------------------------------------
dir.create('./workspace_coffee/rds/run_2')
saveRDS(object = thrp, file = './workspace_coffee/rds/run_2/threshold_prob_coffee.rds') #0.609926, # 0.661486
saveRDS(object = thru, file = './workspace_coffee/rds/run_2/threshold_uncr_coffee.rds') #0.400889, # 0.5488172
terra::writeRaster(fnal, './workspace_coffee/rf/output/run_2/results/process/rf_mixed-bsl.tif', overwrite = TRUE)
saveRDS(clusteredpresdata, file = './workspace_coffee/rds/run_2/occ.rds')

# Boxplot for the agroclimatic zones --------------------------------------
occ.cls <- clusteredpresdata %>% mutate(cluster = factor(cluster, levels = 1:2)) %>% gather(var, value, -c(Lon, Lat, ID, Specie, cluster))
occ.cls <- filter(occ.cls, var != 'bioc_20_100')
occ.cls <- mutate(occ.cls, var = factor(var, levels = vars))
g.box <- ggplot(data = occ.cls, aes(x = cluster, y = value)) + geom_boxplot() + facet_wrap(.~var, scales = 'free_y') + labs(x = 'Cluster', y = 'Value') + theme_bw() + theme(strip.text = element_text(face = 'bold'), axis.text.y = element_text(angle = 90, hjust = 0.5, size = 6), axis.text.x = element_text(size = 6))
# ggsave(plot = g.box, filename = './workspace_coffee/png/graphs/boxplot_variables-cluster_arabica-robusta.jpg', units = 'in', width = 9, height = 7, dpi = 300, create.dir = T)

# Mapping AEZ -------------------------------------------------------------
levels(fnal)
clrs <- c('#FFFFFF',  '#3E5C23', '#C6E01B',  '#D9D9D9', '#FFFFB4')
names(clrs) <- c('Unsuitable', 'Arabica', 'Robusta', 'Limitations', 'Mixed')
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
ggsave(plot = gg.aez, filename = './workspace_coffee/png/maps/run_2/map_rf-mixed_bsl_run2.jpg', units = 'in', width = 8, height = 4, dpi = 300, create.dir = TRUE)


clusteredpresdata
sftr <- st_as_sf(clusteredpresdata, coords = c('Lon', 'Lat'), crs = st_crs(4326))
dir.create('./workspace_coffee/shp')
st_write(sftr, './workspace_coffee/shp/arabica_robusta.shp')
