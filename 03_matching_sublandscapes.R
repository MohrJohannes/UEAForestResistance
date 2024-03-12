##### Librarys ####
library(terra)
library(sp)
library(dplyr)
library(MatchIt)

# Set working directory to source file location 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#### Import #### 
# Import preprocessed stacks
preprocessed_research_matching_stack <- terra::rast("../data/Preprocessed/preprocessed_research_matching_stack.tif")
names(preprocessed_research_matching_stack) <- c("elevation", "slope", "topex", "north.westerness", "coniferous.rate", "ncells")


preprocessed_control_matching_stack <- terra::rast("../data/Preprocessed/preprocessed_control_matching_stack.tif")
names(preprocessed_control_matching_stack) <- c("elevation", "slope", "topex", "north.westerness", "coniferous.rate", "ncells")



# Import Disturbance information rasters
proj_data <- terra::crs(preprocessed_control_matching_stack) 

# load disturbance map 
disturbance_year <- terra::rast("../data/DisturbanceYear/disturbance_year_austria.tif")
terra::crs(disturbance_year) <- proj_data


# import reseach shapes
# load forest maps
site_raster <- terra::rast("../data/Betriebe/site_raster.tif") 



# add site_raster to preprocessed_stacks
preprocessed_research_matching_stack <- c(preprocessed_research_matching_stack, site_raster)
preprocessed_control_matching_stack <- c(preprocessed_control_matching_stack, site_raster)


#### Create searching Dataframe ####
set.seed(42) 

pre_research_df <- terra::as.data.frame(preprocessed_research_matching_stack, xy=TRUE) 
research_df <- pre_research_df %>% 
  .[complete.cases(.),] %>%
  dplyr::filter(ncells >= 334) %>% #only use sublandscapes with at least 30ha of forest inside
  dplyr::group_by(site) %>%
  dplyr::sample_n(10000/4, replace = FALSE) %>% #sample 10`000 areas in total with 2`500 in each site
  dplyr::mutate(Research = 1) %>%
  as.data.frame(.)
head(research_df)

pre_control_df <- as.data.frame(preprocessed_control_matching_stack, xy=TRUE) 
control_df <- pre_control_df %>% 
  .[complete.cases(.[c("elevation", "slope", "topex", "north.westerness", "coniferous.rate", "ncells")]),] %>%
  dplyr::filter(ncells >= 334) %>% #only use sublandscapes with at least 30ha of forest inside
  dplyr::mutate(Research = 0) %>%
  as.data.frame(.)
head(control_df)

searching_df <- rbind(control_df, research_df) 

head(searching_df)

assign("searching_df", searching_df)
save(list = "searching_df", file = "../output/searching_df/searching_df.RData")

load(file = "../output/searching_df/searching_df.RData")

searching_distance <- 5000

matchout <- MatchIt::matchit(Research ~ elevation + slope + topex + north.westerness + coniferous.rate, 
                             data = searching_df, 
                             method = "nearest", 
                             distance = "mahalanobis",
                             caliper = c("x" = searching_distance,
                                         "y" = searching_distance),
                             std.caliper = FALSE,
                             replace = FALSE,
                             verbose = TRUE)


match_data <- MatchIt::match.data(matchout)
head(match_data)

assign("matched_Sublandscapes_before_wrangling_5km", match_data)
save(list = "matched_Sublandscapes_before_wrangling_5km", file = "../output/matched_Sublandscapes/matched_Sublandscapes_before_wrangling_5km.RData")


searching_distance <- 10000

matchout <- MatchIt::matchit(Research ~ elevation + slope + topex + north.westerness + coniferous.rate, 
                             data = searching_df, 
                             method = "nearest", 
                             distance = "mahalanobis",
                             caliper = c("x" = searching_distance,
                                         "y" = searching_distance),
                             std.caliper = FALSE,
                             replace = FALSE,
                             verbose = TRUE)


match_data <- MatchIt::match.data(matchout)
head(match_data)

assign("matched_Sublandscapes_before_wrangling_10km", match_data)
save(list = "matched_Sublandscapes_before_wrangling_10km", file = "../output/matched_Sublandscapes/matched_Sublandscapes_before_wrangling_10km.RData")


rm(match_data, matched_Sublandscapes_before_wrangling_10km, matchout)
searching_distance <- 15000

matchout <- MatchIt::matchit(Research ~ elevation + slope + topex + north.westerness + coniferous.rate, 
                             data = searching_df, 
                             method = "nearest", 
                             distance = "mahalanobis",
                             caliper = c("x" = searching_distance,
                                         "y" = searching_distance),
                             std.caliper = FALSE,
                             replace = FALSE,
                             verbose = TRUE)


match_data <- MatchIt::match.data(matchout)
head(match_data)

assign("matched_Sublandscapes_before_wrangling_15km", match_data)
save(list = "matched_Sublandscapes_before_wrangling_15km", file = "../output/matched_Sublandscapes/matched_Sublandscapes_before_wrangling_15km.RData")

