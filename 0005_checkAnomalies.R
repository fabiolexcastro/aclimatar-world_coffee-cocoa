
# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, cowplot, ggpubr, tidyverse, glue)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1, timeout = 36000)

# Load data ---------------------------------------------------------------
root <- './common_data/esfg_metagrid/intermediate/interpolated_mthly_anomaly'
fles <- as.character(dir_ls(root, regexp = '.tif$'))
gcms <- c('ACCESS-ESM1-5', 'EC-Earth3', 'MPI-ESM1-2-HR', 'INM-CM5-0', 'MRI-ESM2-0')
vars <- c('pr', 'tasmin', 'tasmax', 'rsds', 'hurs')
wrld <- geodata::world(resolution = 1, level = 0, path = './tmpr')
prmt <- aggregate(wrld)


# Function ----------------------------------------------------------------
makeMap <- function(gcme, varb){
  
  # gcme <- gcms[2]
  # varb <- 'pr'
  
  ## To filter
  cat('To process: ', gcme, ' ', varb, '\n')
  fls <- grep(gcme, fles, value = T)
  fls <- grep(varb, fls, value = T)  
  rstr <- rast(fls)
  
  ## Raster to table
  tble <- rstr %>% 
    terra::as.data.frame(xy = T) %>% 
    as_tibble() %>% 
    gather(var, value, -c(x, y)) %>% 
    inner_join(tibble(var = c(paste0('0', 1:9), 10:12), month = month.abb), by = c('var')) %>% 
    mutate(month = factor(month, levels = c(month.abb)))
  
  ## To draw the map 
  ggmp <- ggplot() + 
    geom_tile(data = tble, aes(x = x, y = y, fill = value)) +
    facet_wrap(.~month) +
    scale_fill_viridis_c() +
    # geom_sf(data = st_as_sf(prmt), fill = NA, col = 'grey30') + 
    coord_sf() +
    ggtitle(label = glue(gcme, ' ', varb)) +
    labs(x = '', y = ' ', fill = 'Anomaly') +
    theme_minimal() +
    theme(
      plot.title = element_text(face = 'bold', hjust = 0.5),
      legend.title.position = 'top',
      legend.title = element_text(face = 'bold', hjust = 0.5),
      legend.position = 'bottom', 
      legend.key.width = unit(3, 'line'),
      strip.text = element_text(face = 'bold')
    )
  
  ## Output map 
  dout <- glue('./png/maps/plot_anomaly-inter')
  ggsave(plot = ggmp, filename = glue('{dout}/{gcme}_{varb}.jpg'), units = 'in', width = 8, height = 4, dpi = 300, create.dir = TRUE)
  return(ggmp)
  
}

# Apply the function ------------------------------------------------------
gg1 <- map(.x = 1:length(vars), .f = function(i){makeMap(gcme = gcms[1], varb = vars[i])})
gg2 <- map(.x = 1:length(vars), .f = function(i){makeMap(gcme = gcms[2], varb = vars[i])})
gg3 <- map(.x = 1:length(vars), .f = function(i){makeMap(gcme = gcms[3], varb = vars[i])})
gg4 <- map(.x = 1:length(vars), .f = function(i){makeMap(gcme = gcms[4], varb = vars[i])})
gg5 <- map(.x = 1:length(vars), .f = function(i){makeMap(gcme = gcms[5], varb = vars[i])})

ggs <- list(gg1, gg2, gg3, gg4, gg5)

for(i in 1:length(ggs)){
  print(i)
  gg <- ggarrange(ggs[[i]][[1]], ggs[[i]][[2]], ggs[[i]][[3]], ggs[[i]][[4]], ggs[[i]][[5]], ncol = 2, nrow = 3)
  ou <- glue('./png/maps/plot_anomaly-inter/{gcms[i]}_vars.jpg')
  ggsave(plot = gg, filename = ou, units = 'in', width = 15, height = 10, dpi = 300, create.dir = T)  
}


gg1.all <- ggarrange(gg1[[1]], gg1[[2]], gg1[[3]], gg1[[4]], gg1[[5]], ncol = 2, nrow = 3)
ggsave(plot = gg1.all, filename = './png/maps/plot_anomaly-inter/ACCESS-ESM1-5_vars.jpg', units = 'in', width = 15, height = 10, dpi = 300, create.dir = T)

gg2.all <- ggarrange(gg2[[1]], gg1[[2]], gg1[[3]], gg1[[4]], gg1[[5]], ncol = 2, nrow = 3)
ggsave(plot = gg1.all, filename = './png/maps/plot_anomaly-inter/ACCESS-ESM1-5_vars.jpg', units = 'in', width = 15, height = 10, dpi = 300, create.dir = T)

gg3.all <- ggarrange(gg1[[1]], gg1[[2]], gg1[[3]], gg1[[4]], gg1[[5]], ncol = 2, nrow = 3)
ggsave(plot = gg1.all, filename = './png/maps/plot_anomaly-inter/ACCESS-ESM1-5_vars.jpg', units = 'in', width = 15, height = 10, dpi = 300, create.dir = T)

gg4.all <- ggarrange(gg1[[1]], gg1[[2]], gg1[[3]], gg1[[4]], gg1[[5]], ncol = 2, nrow = 3)
ggsave(plot = gg1.all, filename = './png/maps/plot_anomaly-inter/ACCESS-ESM1-5_vars.jpg', units = 'in', width = 15, height = 10, dpi = 300, create.dir = T)

gg5.all <- ggarrange(gg1[[1]], gg1[[2]], gg1[[3]], gg1[[4]], gg1[[5]], ncol = 2, nrow = 3)
ggsave(plot = gg1.all, filename = './png/maps/plot_anomaly-inter/ACCESS-ESM1-5_vars.jpg', units = 'in', width = 15, height = 10, dpi = 300, create.dir = T)


