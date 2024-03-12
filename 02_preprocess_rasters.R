##### Librarys ####
library(terra)
library(dplyr)

# Set working directory to source file location 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Import data
# Forest Data
forest.cover.map <- terra::rast("../data/ForestCover/forestcover_austria.tif")
coniferous.rate.map <- terra::rast("../data/DominantLeafType/coniferous_rate_raster.tif")

# Elevation Data
elevation.map <- terra::rast("../data/Elevation/elevation_map.tif")
slope.map <- terra::rast("../data/Elevation/slope_map.tif")
topex.map <- terra::rast("../data/Elevation/topex_map.tif")
north.westerness.map <- terra::rast("../data/Elevation/north.westerness_map.tif")


# load forest maps
max.extent <- terra::ext(4480000, 4820000, 260000, 2890000)
sites <- terra::crop(terra::rast("../data/Betriebe/site_raster.tif"), max.extent)
buffer <- terra::buffer(sites, width = 20000, background=NA)


matching_stack <- c(elevation.map, slope.map, topex.map, north.westerness.map, coniferous.rate.map) 
rm(elevation.map, slope.map, topex.map, north.westerness.map, coniferous.rate.map) 

research_mask <- terra::crop(forest.cover.map, sites, mask=TRUE)
control_mask <-  terra::mask(terra::crop(forest.cover.map, buffer, mask=TRUE), sites, inverse = TRUE)


research_matching_stack <- terra::crop(matching_stack, research_mask, mask=TRUE)
control_matching_stack <- terra::crop(matching_stack, control_mask, mask=TRUE)

rm(matching_stack)


#### Funktionen für focal ####
most <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    vector <- as.double(vector)
    tbl <- table(vector, useNA = "no")
    return(as.double(names(tbl[which.max(tbl)])))
  }
}

binned.most <- function(vector, digits=1){
  if (all(is.na(vector))) {
    return(NA)
  } else { 
    vector <- round(as.double(vector), digits=digits)
    tbl <- table(vector, useNA = "no")
    return(as.double(names(tbl[which.max(tbl)]))) 
  }
}

mean_values <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    return(mean(as.double(vector), na.rm=TRUE))
  }
}

ncells <- function(vector) {
  if (all(is.na(vector))) {
    return(NA)
  } else {
    return(sum(!is.na(as.double(vector))))
  }
}



#### Preprocess Data for matching ####
neighborhood_matrix <- rbind(cbind(matrix(1,nrow=20,ncol=20), NA), NA)

# Preprocess Research matching stack
preprocessed_research_matching_stack <- NULL

for (i in 1:terra::nlyr(research_matching_stack)){
  layer <- research_matching_stack[[i]]
  if (names(layer) == "north.westerness"){
    temp <- terra::focal(x=layer, w=neighborhood_matrix, fun=binned.most)
  } else {
    temp <- terra::focal(x=layer, w=neighborhood_matrix, fun=mean_values)
  }
  names(temp) <- names(layer)
  preprocessed_research_matching_stack <- c(preprocessed_research_matching_stack, temp)
}
ncell_layer <- terra::focal(research_matching_stack$elevation, neighborhood_matrix, fun=ncells)
names(ncell_layer) <- "ncells"
preprocessed_research_matching_stack_all <- c(preprocessed_research_matching_stack, ncell_layer)

preprocessed_research_matching_stack_all <- terra::mask(terra::rast(preprocessed_research_matching_stack_all), research_mask)
writeRaster(preprocessed_research_matching_stack_all, "../data/Preprocessed/preprocessed_research_matching_stack.tif", filetype = "GTiff", overwrite = TRUE)


# Preprocess Control matching stack
preprocessed_control_matching_stack <- NULL
for (i in 1:terra::nlyr(control_matching_stack)){
  layer <- control_matching_stack[[i]]
  if (names(layer) == "north.westerness"){
    temp <- terra::focal(x=layer, w=neighborhood_matrix, fun=binned.most)
  } else {
    temp <- terra::focal(x=layer, w=neighborhood_matrix, fun=mean_values)
  }
  names(temp) <- names(layer)
  preprocessed_control_matching_stack <- c(preprocessed_control_matching_stack, temp)
}
ncell_layer <- terra::focal(control_matching_stack$elevation, neighborhood_matrix, fun=ncells)
names(ncell_layer) <- "ncells"
preprocessed_control_matching_stack_all <- c(preprocessed_control_matching_stack, ncell_layer)

preprocessed_control_matching_stack_all <- terra::mask(terra::rast(preprocessed_control_matching_stack_all), control_mask)
writeRaster(preprocessed_control_matching_stack_all, "../data/Preprocessed/preprocessed_control_matching_stack.tif", filetype = "GTiff", overwrite = TRUE)

