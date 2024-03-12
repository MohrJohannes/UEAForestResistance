##### Librarys ####
library(terra)
library(dplyr)

# Set working directory to source file location 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Load Data ####
proj_string_data <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +datum=WGS84"


# Disturbance data
disturbance.agent.map <- terra::rast("../data/DisturbanceAgent/attributed_agents_V3.tif")
disturbance.severity.map <- terra::rast("../data/DisturbanceSeverity/disturbance_severity_austria.tif")
disturbance.year.map <- terra::rast("../data/DisturbanceYear/disturbance_year_austria.tif")


# Calculate high severity map
high.severity <- disturbance.severity.map >= 0.8
names(high.severity) <- "high.severity"
terra::writeRaster(high.severity, "high_severity_map.tif", filetype = "GTiff", overwrite = TRUE)

# Forest Data
forest.cover.map <- terra::rast("../data/ForestCover/forestcover_austria.tif")

# Create coniferous_rate_raster
# Bind all dominant.leaf.type.data
dominant.leaf.data.list <- grep(list.files("../data/DominantLeafType/DATA", pattern="tif", full.names = TRUE), pattern="dbf", invert=TRUE, value=TRUE)
coniferous.rate.map <- terra::classify(terra::init(terra::rast(dominant.leaf.data.list[1]), terra::values(terra::rast(dominant.leaf.data.list[1]))), matrix(c(0:2, NA, 0, 1), ncol=2))

for (dlt.file in dominant.leaf.data.list[-1]) {
  coniferous.rate.map <- terra::merge(coniferous.rate.map,
                                      terra::classify(terra::init(terra::rast(dlt.file), terra::values(terra::rast(dlt.file))), matrix(c(0:2, NA, 0, 1), ncol=2))
  )
}

coniferous.rate.map <- terra::mask(terra::resample(coniferous.rate.map, forest.cover.map, method = "average"), forest.cover.map)
names(temp) <- "coniferous.rate"
terra::writeRaster(temp, "../data/DominantLeafType/coniferous_rate_raster.tif", filetype = "GTiff")


#### Load and Create Elevation Data ####
elevation.map <- terra::rast("../data/Elevation/dhm_at_lamb_10m_2018.tif")
elevation.map <- terra::project(elevation.map, forest.cover.map)
names(elevation.map) <- "elevation"
terra::writeRaster(elevation.map, "../data/Elevation/elevation_map.tif", filetype = "GTiff", overwrite = TRUE)

# Create slope map
slope.map <- terra::terrain(elevation.map, v = "slope", unit = "degrees")
terra::writeRaster(slope.map, "../data/Elevation/slope_map.tif", filetype = "GTiff")

# Create Northerness map
north.westerness.map <- cos(terra::terrain(elevation.map, v = "aspect", unit="radians") + pi/4)
names(north.westerness.map) <- "north.westerness"
terra::writeRaster(north.westerness.map, "../data/Elevation/north.westerness_map.tif", filetype = "GTiff", overwrite = TRUE)

# Create TOPEX map

slope_radians_map <- terra::terrain(terra::rast(elevation.map), v = "slope", unit="radians")
aspect_radians_map <- terra::terrain(terra::rast(elevation.map), v = "aspect", unit="radians")


hillshade_N <- terra::shade(slope_radians_map, aspect_radians_map, angle=5, direction=0)
terra::values(hillshade_N) <- ceiling(values(hillshade_N))

hillshade_NE <- terra::shade(slope_radians_map, aspect_radians_map, angle=5, direction=45)
terra::values(hillshade_NE) <- ceiling(values(hillshade_NE))

hillshade_E <- terra::shade(slope_radians_map, aspect_radians_map, angle=5, direction=90)
terra::values(hillshade_E) <- ceiling(values(hillshade_E))

hillshade_SE <- terra::shade(slope_radians_map, aspect_radians_map, angle=5, direction=135)
terra::values(hillshade_SE) <- ceiling(values(hillshade_SE))

hillshade_S <- terra::shade(slope_radians_map, aspect_radians_map, angle=5, direction=180)
terra::values(hillshade_S) <- ceiling(values(hillshade_S))

hillshade_SW <- terra::shade(slope_radians_map, aspect_radians_map, angle=5, direction=225)
terra::values(hillshade_SW) <- ceiling(values(hillshade_SW))

hillshade_W <- terra::shade(slope_radians_map, aspect_radians_map, angle=5, direction=270)
terra::values(hillshade_W) <- ceiling(values(hillshade_W))

hillshade_NW <- terra::shade(slope_radians_map, aspect_radians_map, angle=5, direction=315)
terra::values(hillshade_NW) <- ceiling(values(hillshade_NW))

hillshades <- c(hillshade_N,
                hillshade_NE,
                hillshade_E,
                hillshade_SE,
                hillshade_S,
                hillshade_SW,
                hillshade_W,
                hillshade_NW)

topex <- terra::app(hillshades, sum)
names(topex) <- "topex"
terra::writeRaster(topex, "../data/Elevation/topex_map.tif", filetype = "GTiff")