##### Librarys ####
library(terra)
library(dplyr)

# Set working directory to source file location 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Load all disturbance rasters
disturbance.year.map <- terra::rast("../data/DisturbanceYear/disturbance_year_austria.tif")
disturbance.agent.map <- terra::rast("../data/DisturbanceAgent/attributed_agents_V3.tif")
disturbance.severity.map <- terra::rast("../data/DisturbanceSeverity/disturbance_severity_austria.tif")
high.severity.map <- terra::rast("../data/DisturbanceSeverity/high_severity_map.tif")
disturbance.patches.natural <- terra::extend(terra::rast("../data/DisturbancePatches/Centroids/natural_patches.tif"), disturbance.year.map)
disturbance.patches.barkbeetle <- terra::extend(terra::rast("../data/DisturbancePatches/Centroids/barkbeetle_patches.tif"), disturbance.year.map)
disturbance.patches.wind <- terra::extend(terra::rast("../data/DisturbancePatches/Centroids/wind_patches.tif"), disturbance.year.map)

forest.cover <- terra::rast("../data/ForestCover/forestcover_austria.tif")

# Load site raster
site_raster <- terra::rast("../data/Betriebe/site_raster.tif")

# Include rasters in stack
disturbance_stack <- c(forest.cover, disturbance.year.map, disturbance.agent.map, disturbance.severity.map, high.severity.map, 
                       disturbance.patches.natural, disturbance.patches.barkbeetle, disturbance.patches.wind)
names(disturbance_stack) <- c("forest.cover", "disturbance.year", "disturbance.agent", "severity", "high.severity", 
                              "patch.size.natural", "patch.size.barkbeetle", "patch.size.wind")
disturbance_stack <- terra::mask(disturbance_stack, forest.cover)

# Create Research Stack
Research_disturbance_stack <- terra::mask(disturbance_stack, site_raster)

# Create Control Stack
Control_disturbance_stack <- terra::mask(disturbance_stack, site_raster, inverse = TRUE)

rm(disturbance_stack, forest.cover, disturbance.year.map, disturbance.agent.map, disturbance.severity.map, high.severity.map, 
   disturbance.patches.natural, disturbance.patches.barkbeetle, disturbance.patches.wind, site_raster)


# function to extract disturbance_information of sublandscape
get_disturbance_metrics <- function(x, y, Research){
  if (Research == 1){
    temp_stack <- Research_disturbance_stack
  } else{
    temp_stack <- Control_disturbance_stack
  }
  cell  <- terra::cellFromXY(temp_stack, cbind(x, y)) #Testcell: x <- 4604600; y <- 2855400
  
  # Define Neighborhood matrix
  neighborhood_matrix <- rbind(cbind(matrix(1,nrow=20,ncol=20), 0), 0)
  
  cellnumbers <- as.vector(terra::adjacent(temp_stack, cell, directions = neighborhood_matrix))
  disturbance_information <- as.data.frame(temp_stack[cellnumbers, drop = FALSE], xy=TRUE)
  
  #plot(temp_stack[cellnumbers, drop = FALSE])
  
  cells_in_sublandscape <- nrow(disturbance_information)
  
  #get numbers of disturbed cells per year
  cells.per.year <- disturbance_information %>%
    dplyr::filter(complete.cases(disturbance.year),
                  !is.na(disturbance.agent)) %>%
    dplyr::group_by(disturbance.year, disturbance.agent) %>%
    dplyr::summarise(ncells.per.year = n()) %>%
    rbind(disturbance_information %>%
            dplyr::filter(complete.cases(disturbance.year),
                          disturbance.agent != 3) %>%
            dplyr::group_by(disturbance.year) %>%
            dplyr::summarise(ncells.per.year = n()) %>%
            dplyr::mutate(disturbance.agent = 0)) %>%
    dplyr::right_join(data.frame(disturbance.year = rep(1986:2020, times=4),
                                 disturbance.agent = rep(0:3, each=35),
                                 default.cells.per.year = 0), by = c("disturbance.year", "disturbance.agent")) %>%
    dplyr::mutate(ncells.per.year = pmax(ncells.per.year, default.cells.per.year, na.rm=TRUE)) %>%
    dplyr::arrange(disturbance.year) %>%
    dplyr::mutate(disturbance.year = paste0("ncells.", dplyr::case_when(disturbance.agent == 0 ~ "natural", 
                                                                        disturbance.agent == 1 ~ "barkbeetle", 
                                                                        disturbance.agent == 2 ~ "wind", 
                                                                        disturbance.agent == 3 ~ "harvest"), disturbance.year)) %>%
    dplyr::select(-c(default.cells.per.year, disturbance.agent)) %>% 
    tidyr::spread(key=disturbance.year, value = ncells.per.year) %>%
    as.matrix()
  
  
  #get number of severe disturbed cells per year
  severity.cells.per.year <- disturbance_information %>%
    dplyr::filter(high.severity == 1,
                  disturbance.agent != 3) %>%
    dplyr::group_by(disturbance.year, disturbance.agent) %>%
    dplyr::summarise(ncells.per.year = n()) %>%
    rbind(disturbance_information %>%
            dplyr::filter(high.severity == 1,
                          disturbance.agent != 3) %>%
            dplyr::group_by(disturbance.year) %>%
            dplyr::summarise(ncells.per.year = n()) %>%
            dplyr::mutate(disturbance.agent = 0))%>%
    dplyr::right_join(data.frame(disturbance.year = rep(1986:2020, times=3),
                                 disturbance.agent = rep(0:2, each=35),
                                 default.cells.per.year = 0), by = c("disturbance.year", "disturbance.agent")) %>%
    dplyr::mutate(ncells.per.year = pmax(ncells.per.year, default.cells.per.year, na.rm=TRUE)) %>%
    dplyr::arrange(disturbance.year) %>%
    dplyr::mutate(disturbance.year = paste0("ncells.severe.", dplyr::case_when(disturbance.agent == 0 ~ "natural", 
                                                                               disturbance.agent == 1 ~ "barkbeetle", 
                                                                               disturbance.agent == 2 ~ "wind"), disturbance.year)) %>%
    dplyr::select(-c(default.cells.per.year, disturbance.agent)) %>%
    tidyr::spread(key=disturbance.year, value = ncells.per.year) %>%
    as.matrix()
  
  #get patch.size per year
  patch.size.per.year <- disturbance_information %>%
    dplyr::filter(disturbance.agent %in% c(1,2)) %>%
    dplyr::group_by(disturbance.year) %>%
    dplyr::summarise(disturbance.agent = 0,
                     patch.size = max(patch.size.natural, 0, na.rm=TRUE)) %>%
    rbind(disturbance_information %>%
            dplyr::filter(disturbance.agent == 1) %>%
            dplyr::group_by(disturbance.year, disturbance.agent) %>%
            dplyr::summarise(patch.size = max(patch.size.barkbeetle, 0, na.rm=TRUE))) %>%
    rbind(disturbance_information %>%
            dplyr::filter(disturbance.agent == 2) %>%
            dplyr::group_by(disturbance.year, disturbance.agent) %>%
            dplyr::summarise(patch.size = max(patch.size.wind, 0, na.rm=TRUE))) %>%
    dplyr::right_join(data.frame(disturbance.year = rep(1986:2020, times=3),
                                 disturbance.agent = rep(0:2, each=35),
                                 default.patch.size = 0), by = c("disturbance.year", "disturbance.agent")) %>%
    dplyr::mutate(patch.size = pmax(patch.size, default.patch.size, na.rm=TRUE)) %>%
    dplyr::arrange(disturbance.year) %>%
    dplyr::mutate(disturbance.year = paste0("patch.size.", dplyr::case_when(disturbance.agent == 0 ~ "natural", 
                                                                            disturbance.agent == 1 ~ "barkbeetle", 
                                                                            disturbance.agent == 2 ~ "wind"), disturbance.year)) %>%
    dplyr::select(-c(default.patch.size, disturbance.agent)) %>%
    tidyr::spread(key=disturbance.year, value = patch.size) %>%
    as.matrix()
  
  return(list("ncells" = cells_in_sublandscape,
              "cells.per.year" = cells.per.year,
              "severity.cells.per.year" = severity.cells.per.year,
              "patch.size.per.year" = patch.size.per.year))
}

options(dplyr.summarise.inform = FALSE)

##### --------------------- 5 km -------------------- #####
load(file = paste0("../output/matched_Sublandscapes/matched_Sublandscapes_before_wrangling_5km.RData"))
export_df <- matched_Sublandscapes_before_wrangling_5km
rm(matched_Sublandscapes_before_wrangling_5km)

# calculations
pb <- dplyr::progress_estimated(nrow(export_df))
for (i in 1:nrow(export_df)){
  sublandscape <- export_df[i,]
  disturbance_informations <- get_disturbance_metrics(x = sublandscape$x, y = sublandscape$y, Research = sublandscape$Research)
  if (disturbance_informations$ncells != sublandscape$ncells) print(paste("Fehler! Anzahl Zellen ist unterschiedlich! Reihe i:", i))
  
  export_df[i, paste0("disturbed.cells.", rep(c("barkbeetle", "harvest", "natural", "wind"), each=35), 1986:2020)] <- disturbance_informations$cells.per.year
  export_df[i, paste0("severe.disturbed.cells.", rep(c("barkbeetle", "natural", "wind"), each=35), 1986:2020)] <- disturbance_informations$severity.cells.per.year
  export_df[i, paste0("patch.size.", rep(c("barkbeetle", "natural", "wind"), each=35), 1986:2020)] <- disturbance_informations$patch.size.per.year
  pb$tick()$print()
  
}

head(export_df)
assign("Data_Sublandscapes_before_wrangling_5km", export_df)
save(list = "Data_Sublandscapes_before_wrangling_5km", file = "../output_clean/data_Sublandscapes/Data_Sublandscapes_before_wrangling_5km.RData")

temp <- as.data.frame(MatchIt::mahalanobis_dist(Research ~ elevation + slope + topex + north.westerness + coniferous.rate, 
                                                data = export_df %>% dplyr::arrange(subclass)))

export_df_including_dist <- export_df %>%
  dplyr::arrange(subclass) %>%
  dplyr::mutate(match.dist = rep(diag(as.matrix(temp)), each=2))

save(export_df_including_dist, file="../output_clean/data_Sublandscapes/Data_Sublandscapes_incl_mahabDist_5km.RData")


##### --------------------- 10 km -------------------- #####
load(file = paste0("../output/matched_Sublandscapes/matched_Sublandscapes_before_wrangling_10km.RData"))
export_df <- matched_Sublandscapes_before_wrangling_10km
rm(matched_Sublandscapes_before_wrangling_10km)

# calculations
pb <- dplyr::progress_estimated(nrow(export_df))
for (i in 1:nrow(export_df)){
  sublandscape <- export_df[i,]
  disturbance_informations <- get_disturbance_metrics(x = sublandscape$x, y = sublandscape$y, Research = sublandscape$Research)
  if (disturbance_informations$ncells != sublandscape$ncells) print(paste("Fehler! Anzahl Zellen ist unterschiedlich! Reihe i:", i))
  
  export_df[i, paste0("disturbed.cells.", rep(c("barkbeetle", "harvest", "natural", "wind"), each=35), 1986:2020)] <- disturbance_informations$cells.per.year
  export_df[i, paste0("severe.disturbed.cells.", rep(c("barkbeetle", "natural", "wind"), each=35), 1986:2020)] <- disturbance_informations$severity.cells.per.year
  export_df[i, paste0("patch.size.", rep(c("barkbeetle", "natural", "wind"), each=35), 1986:2020)] <- disturbance_informations$patch.size.per.year
  pb$tick()$print()
  
}

head(export_df)
assign("Data_Sublandscapes_before_wrangling_10km", export_df)
save(list = "Data_Sublandscapes_before_wrangling_10km", file = "../output_clean/data_Sublandscapes/Data_Sublandscapes_before_wrangling_10km.RData")

temp <- as.data.frame(MatchIt::mahalanobis_dist(Research ~ elevation + slope + topex + north.westerness + coniferous.rate,  
                                                data = export_df %>% dplyr::arrange(subclass)))

export_df_including_dist <- export_df %>%
  dplyr::arrange(subclass) %>%
  dplyr::mutate(match.dist = rep(diag(as.matrix(temp)), each=2))

save(export_df_including_dist, file="../output_clean/data_Sublandscapes/Data_Sublandscapes_incl_mahabDist_10km.RData")


##### --------------------- 15 km -------------------- #####
load(file = paste0("../output/matched_Sublandscapes/matched_Sublandscapes_before_wrangling_15km.RData"))
export_df <- matched_Sublandscapes_before_wrangling_15km
rm(matched_Sublandscapes_before_wrangling_15km)

# calculations
pb <- dplyr::progress_estimated(nrow(export_df))
for (i in 1:nrow(export_df)){
  sublandscape <- export_df[i,]
  disturbance_informations <- get_disturbance_metrics(x = sublandscape$x, y = sublandscape$y, Research = sublandscape$Research)
  if (disturbance_informations$ncells != sublandscape$ncells) print(paste("Fehler! Anzahl Zellen ist unterschiedlich! Reihe i:", i))
  
  export_df[i, paste0("disturbed.cells.", rep(c("barkbeetle", "harvest", "natural", "wind"), each=35), 1986:2020)] <- disturbance_informations$cells.per.year
  export_df[i, paste0("severe.disturbed.cells.", rep(c("barkbeetle", "natural", "wind"), each=35), 1986:2020)] <- disturbance_informations$severity.cells.per.year
  export_df[i, paste0("patch.size.", rep(c("barkbeetle", "natural", "wind"), each=35), 1986:2020)] <- disturbance_informations$patch.size.per.year
  pb$tick()$print()
  
}

head(export_df)
assign("Data_Sublandscapes_before_wrangling_15km", export_df)
save(list = "Data_Sublandscapes_before_wrangling_15km", file = "../output_clean/data_Sublandscapes/Data_Sublandscapes_before_wrangling_15km.RData")

temp <- as.data.frame(MatchIt::mahalanobis_dist(Research ~ elevation + slope + topex + north.westerness + coniferous.rate, 
                                                data = export_df %>% dplyr::arrange(subclass)))

export_df_including_dist <- export_df %>%
  dplyr::arrange(subclass) %>%
  dplyr::mutate(match.dist = rep(diag(as.matrix(temp)), each=2))

save(export_df_including_dist, file="../output_clean/data_Sublandscapes/Data_Sublandscapes_incl_mahabDist_15km.RData")
