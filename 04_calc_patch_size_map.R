##### Librarys ####
library(terra)
library(dplyr)

# Set working directory to source file location 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
site.raster <- terra::rast("../data/Betriebe/site_raster.tif")
buffer <- terra::buffer(site.raster, width = 20000) # to reduce computation time

disturbance.year.map <- terra::mask(terra::crop(terra::rast("../data/DisturbanceYear/disturbance_year_austria.tif"), buffer), buffer, maskvalues=0)
disturbance.agent.map <- terra::mask(terra::crop(terra::rast("../data/DisturbanceAgent/attributed_agents_V3.tif"), buffer), buffer, maskvalues=0)

rm(buffer)

barkbeetle.map <- disturbance.agent.map == 1
wind.map <- disturbance.agent.map == 2
natural.dist.map <- disturbance.agent.map %in% c(1,2)


years <- 1986:2020

natural.patches.df <- data.frame()
barkbeetle.patch.df <- data.frame()
wind.patch.df <- data.frame()


for (year in years){
  print(paste("----------------------", year, "---------------------"))
  temp.rast <- disturbance.year.map 
  temp.rast[!temp.rast %in% c(year, year + 1)] <- NA
  temp.rast[!is.na(temp.rast)] <- year
  
  # natural Disturbance patches
  natural.patch.map <- terra::patches(terra::mask(natural.dist.map, temp.rast), directions=8, zeroAsNA=TRUE) #queens case
  
  
  if (!is.null(unique(natural.patch.map))){
    natural.patches <- terra::as.data.frame(natural.patch.map, xy=TRUE) %>%
      dplyr::group_by(patches) %>%
      dplyr::summarise(n = n(),
                       centroid_x = floor((mean(x) - terra::ext(natural.patch.map)[1])/30)*30 + terra::ext(natural.patch.map)[1] + 15,
                       centroid_y = floor((mean(y) - terra::ext(natural.patch.map)[3])/30)*30 + terra::ext(natural.patch.map)[3] + 15) %>%
      as.data.frame(.) %>%
      dplyr::rename(x = centroid_x, y = centroid_y) %>%
      dplyr::mutate(coords = paste(x, y)) %>%
      dplyr::group_by(coords) %>%
      dplyr::summarise(natural.patch.size = max(n)) %>%
      as.data.frame(.) %>%
      dplyr::mutate(Year = year, Agent = "Natural") %>%
      dplyr::select(coords, Year, Agent, natural.patch.size)
    
    write.csv(natural.patches, paste0("../data/DisturbancePatches/CSVs/natural_patches", year, ".csv"))
    
    natural.patches.df <- rbind(natural.patches.df, natural.patches)
    
    rm(natural.patches)
  }
  rm(natural.patch.map)
  print("Natural done!")
  
  # Agent 1
  barkbeetle.patch.map <- terra::patches(terra::mask(barkbeetle.map, temp.rast), directions=8, zeroAsNA=TRUE)
  
  
  if (!is.null(unique(barkbeetle.patch.map))){
    barkbeetle.patches <- as.data.frame(barkbeetle.patch.map, xy=TRUE) %>%
      dplyr::group_by(patches) %>%
      dplyr::summarise(n = n(),
                       centroid_x = floor((mean(x) - terra::ext(barkbeetle.patch.map)[1])/30)*30 + terra::ext(barkbeetle.patch.map)[1] + 15,
                       centroid_y = floor((mean(y) - terra::ext(barkbeetle.patch.map)[3])/30)*30 + terra::ext(barkbeetle.patch.map)[3] + 15) %>%
      as.data.frame(.) %>%
      dplyr::rename(x = centroid_x, y = centroid_y) %>%
      dplyr::mutate(coords = paste(x, y)) %>%
      dplyr::group_by(coords) %>%
      dplyr::summarise(barkbeetle.patch.size = max(n)) %>%
      as.data.frame(.) %>%
      dplyr::mutate(Year = year, Agent = "Barkbeetle") %>%
      dplyr::select(coords, Year, Agent, barkbeetle.patch.size)
    
    write.csv(barkbeetle.patches, paste0("../data/DisturbancePatches/CSVs/barkbeetle_patches", year, ".csv"))
    
    barkbeetle.patch.df <- rbind(barkbeetle.patch.df, barkbeetle.patches)
    
    rm(barkbeetle.patches)
  }
  rm(barkbeetle.patch.map)
  print("Barkbeetle done!")
  
  # Agent 2
  wind.patch.map <- terra::patches(terra::mask(wind.map, temp.rast), directions=8, zeroAsNA=TRUE)
  
  
  if (!is.null(unique(wind.patch.map))) {
    wind.patches <- as.data.frame(wind.patch.map, xy=TRUE) %>%
      dplyr::group_by(patches) %>%
      dplyr::summarise(n = n(),
                       centroid_x = floor((mean(x) - terra::ext(wind.patch.map)[1])/30)*30 + terra::ext(wind.patch.map)[1] + 15,
                       centroid_y = floor((mean(y) - terra::ext(wind.patch.map)[3])/30)*30 + terra::ext(wind.patch.map)[3] + 15) %>%
      as.data.frame(.) %>%
      dplyr::rename(x = centroid_x, y = centroid_y) %>%
      dplyr::mutate(coords = paste(x, y)) %>%
      dplyr::group_by(coords) %>%
      dplyr::summarise(wind.patch.size = max(n)) %>%
      as.data.frame(.) %>%
      dplyr::mutate(Year = year, Agent = "Wind") %>%
      dplyr::select(coords, Year, Agent, wind.patch.size)
    
    write.csv(wind.patches, paste0("../data/DisturbancePatches/CSVs/wind_patches", year, ".csv"))
    
    wind.patch.df <- rbind(wind.patch.df, wind.patches)
    
    rm(wind.patches)
  }
  rm(wind.patch.map)
  print("Wind done!")
}

patches.df <- dplyr::select(natural.patches.df, -Agent) %>%
  dplyr::full_join(dplyr::select(barkbeetle.patch.df, -Agent)) %>%
  dplyr::full_join(dplyr::select(wind.patch.df, -Agent)) %>%
  dplyr::group_by(coords, Year) %>%
  dplyr::summarise(natural.patch.size = max(natural.patch.size), 
                   barkbeetle.patch.size = max(barkbeetle.patch.size), 
                   wind.patch.size = max(wind.patch.size)) %>%
  tidyr::separate(coords, sep=" ", into=c("x", "y")) %>%
  as.data.frame() %>%
  dplyr::mutate(x=as.numeric(x), y=as.numeric(y))

natural.patches.raster <- terra::rast(dplyr::select(patches.df, x, y, natural.patch.size), type="xyz", crs=terra::crs(disturbance.year.map))
origin(natural.patches.raster) <- origin(disturbance.year.map)
terra::writeRaster(natural.patches.raster, "../data/DisturbancePatches/Centroids/natural_patches_centroids.tif", overwrite = TRUE)

natural.patches.raster <- terra::rast(dplyr::select(patches.df, x, y, barkbeetle.patch.size), type="xyz", crs=terra::crs(disturbance.year.map))
origin(natural.patches.raster) <- origin(disturbance.year.map)
terra::writeRaster(natural.patches.raster, "../data/DisturbancePatches/Centroids/barkbeetle_patches_centroids.tif", overwrite = TRUE)

natural.patches.raster <- terra::rast(dplyr::select(patches.df, x, y, wind.patch.size), type="xyz", crs=terra::crs(disturbance.year.map))
origin(natural.patches.raster) <- origin(disturbance.year.map)
terra::writeRaster(natural.patches.raster, "../data/DisturbancePatches/Centroids/wind_patches_centroids.tif", overwrite = TRUE)

