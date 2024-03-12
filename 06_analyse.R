# Librarys
library(dplyr)
library(terra)

library(rnaturalearth)
library(ggplot2)
library(ggpubr)

# Set working directory to source file location 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load Main Data
load(file = "../../output/data_Sublandscapes_noTRI/Data_Sublandscapes_incl_mahabDist_10km.RData")
data <- export_df_including_dist; rm(export_df_including_dist)


max.dist <- max(data$match.dist)
data <- data %>% 
  dplyr::full_join(dplyr::select(dplyr::filter(data, Research==1), site, subclass), by="subclass") %>%
  dplyr::mutate(rel.dist = (1 - match.dist/max.dist) * 100,
                Site = site.y,
                Research = dplyr::case_when(Research == 1 ~ "Uneven-aged", TRUE ~ "Even-aged")) %>% 
  dplyr::relocate(Site, .before = elevation) %>%
  dplyr::select(-c(site.x, site.y))

data <- data %>% dplyr::filter(match.dist <= 1); rm(max.dist)
data$Research <- factor(data$Research, levels = c("Uneven-aged", "Even-aged"))


# Resample data, so that each site has same number of samples
set.seed(25) #5
subclass.data <- data %>%
  dplyr::filter(Research == "Uneven-aged") %>%
  dplyr::group_by(Site) %>%
  dplyr::sample_n(1250, replace=FALSE) %>%
  dplyr::select(subclass)

data <- dplyr::filter(data, subclass %in% subclass.data$subclass)
rm(subclass.data)

managament.colors <- c("#f0b030", "#cccccc")
metric.colors <- c("#E69F00", "#0072B2", "#CC79A7", "#D55E00")


#### -------------------------------- Calc area studied ------------------------------------ ####
all.rast <- data %>%
  dplyr::select(x, y, Research) %>%
  dplyr::mutate(Research = Research == "Uneven-aged") %>% 
  as.data.frame() %>%
  terra::rast(type="xyz")

terra::plot(all.rast)


# Define a custom neighborhood matrix with a size of 21 by 21 cells
neighborhood_matrix <- rbind(cbind(matrix(1,nrow=20,ncol=20), NA), NA)

# Use focal function to apply the neighborhood operation
result_rast <- terra::focal(all.rast, w = neighborhood_matrix, fun = function(x) {
  # Check if there are any non-NA values in the neighborhood
  if (any(!is.na(x))) {
    # Assign the value of the central cell to all non-NA cells in the neighborhood
    return(max(x, na.rm=TRUE))
  } else {
    # If all values in the neighborhood are NA, return NA
    return(NA)
  }
})

# Plot the result if needed
terra::plot(result_rast)
table(terra::values(result_rast))*30*30/10000

rm(all.rast, neighborhood_matrix, result_rast)

#### -------------------------------- Figure 1 - site map  -------------------------------- ####
site.raster <- terra::rast("../data/Betriebe/site_raster.tif")
buffer <- terra::buffer(site.raster, width = 10000)

AOI <- rnaturalearth::ne_countries(scale = "medium", continent = "Europe")
AOI <- terra::crop(terra::project(terra::vect(AOI), buffer), terra::ext(4220000, 4900000, 2550000, 2950000))

elevation.map <- terra::rast("../data/Elevation/elevation_map.tif")
elevation.map <- terra::extend(elevation.map, AOI)
elevation.colors <- colorRampPalette(c("#70C820", "#8FB92A", "#008200", "#006400", "#B47E34", "#BF8171", "#EAD6D1", "#FFFFFF"), bias=1.5)(200)

buffer <- terra::extend(buffer, AOI)

buffer <- terra::mask(buffer, terra::crop(elevation.map, buffer))
buffer[buffer==0] <- NA

country.names <- terra::init(terra::aggregate(elevation.map, 100), NA)
Country.df <- data.frame(ID = 1:9,
                         Name=c("DE", "CZ", "SK", "AT", "HU", "CH", "IT", "SI", "HR"),
                         x=c(4370000, 4720000, 4880000, 4650000, 4880000, 4250000, 4450000, 4700000, 4820000), 
                         y=c(2850000, 2910000, 2830000, 2720000, 2680000, 2640000, 2600000, 2570000, 2570000)) 

for (cntr in 1:nrow(Country.df)){
  cell <- terra::cellFromXY(country.names, Country.df[cntr, c("x", "y")])
  country.names[cell] <- Country.df[cntr, "ID"]
}
levels(country.names) <- Country.df[,c("ID", "Name")]

png("../graphs/site_map.png", res=1000, units="cm",  width=13, height=8)
terra::plot(elevation.map, col=elevation.colors, type="continuous", legend=TRUE, box=FALSE, axes=FALSE, plg=list(title="Elevation [m]", title.cex=0.9))
terra::plot(buffer, col="black", alpha=0.6, add=TRUE, legend=FALSE, box=FALSE, axes=FALSE)
terra::plot(AOI, add=TRUE)
terra::text(country.names, digits=2, cex=1)
terra::north(xy=c(4300000, 2900000), type=2, cex=0.7) 
terra::sbar(xy=c(4320000, 2900000), type="bar", divs=4, scaleby=1000, below="km", cex=0.5) 
dev.off()

rm(site.raster, buffer, AOI, elevation.map, elevation.colors, country.names, Country.df, cell, cntr); gc()


#### -------------------------------- Table 1 -------------------------------- ####
site.raster <- terra::rast("../data/Betriebe/site_raster.tif")
elevation.raster <- terra::rast("../data/Elevation/elevation_map.tif")
slope.raster <- terra::rast("../data/Elevation/slope_map.tif")
conifer.raster <- terra::rast("../data/DominantLeafType/coniferous_rate_raster.tif")

# Area
Area.data <- terra::values(site.raster, dataframe=TRUE, na.rm=TRUE) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(Area_ha = n()*30*30/10000) 

# Site and Stand conditions
tmp.data <- terra::values(c(site.raster, elevation.raster, slope.raster, conifer.raster), dataframe=TRUE, na.rm=TRUE) %>%
  dplyr::mutate(Site = ifelse(site=="aaa", "Site A",                            # site names hidden
                              ifelse(site=="bbb", "Site B",
                                     ifelse(site=="ccc", "Site C", 
                                            ifelse(site=="ddd", "Site D", "Else"))))) %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(min.elevation = min(elevation),
                   mean.elevaton = mean(elevation),
                   max.elevation = max(elevation),
                   min.slope = min(slope),
                   mean.slope = mean(slope),
                   max.slope = max(slope),
                   min.conifer = min(coniferous.rate), 
                   mean.conifer = mean(coniferous.rate),
                   max.conifer = max(coniferous.rate))

rm(site.raster, elevation.raster, slope.raster, conifer.raster, Area.data, tmp.data)


#### -------------------------------- Table 2 -------------------------------- ####
# Data for table 1
UEA_data <- dplyr::filter(data, Research=="Uneven-aged")
round(c(quantile(UEA_data$elevation, probs=0.5), quantile(UEA_data$elevation, probs=0.05), quantile(UEA_data$elevation, probs=0.95)),1)
round(c(quantile(UEA_data$slope, probs=0.5), quantile(UEA_data$slope, probs=0.05), quantile(UEA_data$slope, probs=0.95)),1)
round(c(quantile(UEA_data$topex, probs=0.5), quantile(UEA_data$topex, probs=0.05), quantile(UEA_data$topex, probs=0.95)),1)
round(c(quantile(UEA_data$north.westerness, probs=0.5), quantile(UEA_data$north.westerness, probs=0.05), quantile(UEA_data$north.westerness, probs=0.95)),1)
round(100*c(quantile(UEA_data$coniferous.rate, probs=0.5), quantile(UEA_data$coniferous.rate, probs=0.05), quantile(UEA_data$coniferous.rate, probs=0.95)),1)


EA_data <- dplyr::filter(data, Research=="Even-aged")
round(c(quantile(EA_data$elevation, probs=0.5), quantile(EA_data$elevation, probs=0.05), quantile(EA_data$elevation, probs=0.95)))
round(c(quantile(EA_data$slope, probs=0.5), quantile(EA_data$slope, probs=0.05), quantile(EA_data$slope, probs=0.95)),1)
round(c(quantile(EA_data$topex, probs=0.5), quantile(EA_data$topex, probs=0.05), quantile(EA_data$topex, probs=0.95)),1)
round(c(quantile(EA_data$north.westerness, probs=0.5), quantile(EA_data$north.westerness, probs=0.05), quantile(EA_data$north.westerness, probs=0.95)),1)
round(100*c(quantile(EA_data$coniferous.rate, probs=0.5), quantile(EA_data$coniferous.rate, probs=0.05), quantile(EA_data$coniferous.rate, probs=0.95)),1)

rm(UEA_data, EA_data)

##### -------------------------------- Main Plots -------------------------------- #####
#### -------------------------------- Figure 3 - boxplot comparison -------------------------------- ####
temp.data <- data %>% 
  dplyr::select(Site, Research, subclass, ncells, disturbed.cells.barkbeetle1986:patch.size.wind2020) %>%
  tidyr::gather(-c(1,2,3, 4), key="key", value="value") %>%
  dplyr::mutate(DisturbanceParameter = dplyr::case_when(as.logical(grepl("severe", key)) ~ "ncells.severe.disturbed",
                                                        as.logical(grepl("size", key)) ~ "Patch.Size",
                                                        TRUE ~ "ncells.disturbed"),
                agent = dplyr::case_when(grepl("wind", key) ~ "wind", 
                                         grepl("barkbeetle", key) ~ "barkbeetle", 
                                         grepl("natural", key) ~ "natural"),
                year = as.numeric(gsub("barkbeetle", "", 
                                       gsub("wind", "",  
                                            gsub("natural", "", 
                                                 gsub("disturbed.cells.", "", 
                                                      gsub("severe.disturbed.cells.", "",
                                                           gsub("patch.size.", "", key)))))))) %>%
  dplyr::select(-key) %>%
  tidyr::spread(DisturbanceParameter, value) %>%
  dplyr::mutate(DistEvent = ncells.disturbed > 0)

temp.data <- temp.data %>% 
  dplyr::right_join(temp.data %>% 
                      dplyr::select(Site:ncells.disturbed) %>%
                      tidyr::spread(key=agent, value=ncells.disturbed) %>%
                      dplyr::mutate(Wind.ratio = wind/(wind+barkbeetle)) %>%
                      dplyr::select(-c(barkbeetle, natural, wind)), by=c("Site", "Research", "subclass", "year", "ncells")) %>% 
  dplyr::group_by(Site, Research, agent, subclass, ncells) %>%
  dplyr::summarise(ncells.disturbed = sum(ncells.disturbed),
                   nEvents = sum(DistEvent),
                   ncells.severe.disturbed =sum(ncells.severe.disturbed),
                   Patch.Size = max(Patch.Size, na.rm=TRUE),
                   Wind.ratio = mean(Wind.ratio, na.rm=TRUE)) %>% 
  dplyr::mutate(annual.disturbance.rate = ncells.disturbed/ncells/35,
                Severity.rate = ncells.severe.disturbed/ncells.disturbed) %>%
  dplyr::select(-c(ncells, ncells.disturbed, ncells.severe.disturbed)) %>%
  tidyr::gather(5:9, key="DisturbanceMetric", value="Value") 
temp.data$Research <- factor(temp.data$Research, levels = c("Uneven-aged", "Even-aged"))
temp.data$DisturbanceMetricName <- dplyr::case_when(temp.data$DisturbanceMetric == "annual.disturbance.rate" ~ "Annual Disturbance Rate",
                                                    temp.data$DisturbanceMetric == "nEvents" ~ "Annual Number of Events",
                                                    temp.data$DisturbanceMetric == "Severity.rate" ~ "Severity Rate",
                                                    temp.data$DisturbanceMetric == "Patch.Size" ~ "Patch Size")

rate.plot <- ggplot(dplyr::filter(temp.data, DisturbanceMetric == "annual.disturbance.rate", agent=="natural"), aes(x=Research, y=Value*100, fill=Research))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(data = dplyr::filter(temp.data, DisturbanceMetric == "annual.disturbance.rate", agent=="natural") %>%
               dplyr::group_by(Research) %>%
               dplyr::summarise(mean=mean(Value)), aes(x=Research, y=mean*100))+
  coord_cartesian(ylim = c(0, 0.6))+
  scale_x_discrete(labels = NULL)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.15))+ 
  scale_fill_manual(values=managament.colors)+
  labs(y=expression(paste('Rate [% year' ^'-1', ']')), x="", fill="Management")+
  theme_bw()

n.plot <- ggplot(dplyr::filter(temp.data, DisturbanceMetric == "nEvents", agent=="natural"), aes(x=Research, y=Value/35, fill=Research))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(data = dplyr::filter(temp.data, DisturbanceMetric == "nEvents", agent=="natural") %>%
               dplyr::group_by(Research) %>%
               dplyr::summarise(mean=mean(Value)), aes(x=Research, y=mean/35))+
  coord_cartesian(ylim = c(0, 0.2))+
  scale_x_discrete(labels = NULL)+
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05))+
  scale_fill_manual(values=managament.colors)+
  labs(x="", y=expression(paste('Frequency [year' ^'-1', ']')), fill="Management")+
  theme_bw()

severity.plot <- ggplot(dplyr::filter(temp.data, DisturbanceMetric == "Severity.rate", agent=="natural"), aes(x=Research, y=Value*100, fill=Research))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(data = dplyr::filter(temp.data, DisturbanceMetric == "Severity.rate", agent=="natural") %>%
               dplyr::group_by(Research) %>%
               dplyr::summarise(mean=mean(Value, na.rm=TRUE)), aes(x=Research, y=mean*100))+
  scale_x_discrete(labels = NULL)+
  scale_y_continuous(breaks = seq(0, 100, by = 25))+ 
  scale_fill_manual(values=managament.colors)+
  labs(x="", y="Severity [%]", fill="Management")+
  theme_bw()

size.plot <- ggplot(dplyr::filter(temp.data, DisturbanceMetric == "Patch.Size", agent=="natural"), aes(x=Research, y=Value*30*30/10000, fill=Research))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(data = dplyr::filter(temp.data, DisturbanceMetric == "Patch.Size", agent=="natural") %>%
               dplyr::group_by(Research) %>%
               dplyr::summarise(mean=mean(Value)), aes(x=Research, y=mean*30*30/10000))+
  coord_cartesian(ylim = c(0, 3))+
  scale_x_discrete(labels = NULL)+
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), labels = paste0("   ", seq(0, 3, by = 0.5)))+
  scale_fill_manual(values=managament.colors)+
  labs(x="", y="Size [ha]", fill="Management")+
  theme_bw()


png("../graphs/Figure_3.png", res = 1000, width=16, height=5, units="cm")
ggarrange(rate.plot, n.plot, size.plot, severity.plot, label.x = c(0.42, 0.42, 0.405, 0.365), label.y = 0.98,
          common.legend = TRUE, legend = "top", ncol=4, nrow=1, widths = c(1, 1, 0.98, 0.94), labels = "auto") 
dev.off()

png("../graphs/boxplots_Comparison_GraphicalAbstract.png", res = 1000, width=15, height=7, units="cm", bg="#EDD9BA")
ggarrange(rate.plot, n.plot, size.plot, severity.plot, label.x = c(0.42, 0.42, 0.405, 0.365), label.y = 0.98,
          common.legend = TRUE, legend = "top", ncol=4, nrow=1, labels="auto", widths = c(1, 1, 0.98, 0.94))
dev.off()
rm(temp.data, rate.plot, n.plot, severity.plot, size.plot)



#### -------------------------------- Table 3 - t.test -------------------------------- ####
all.data <- data %>%
  dplyr::select(elevation:coniferous.rate, Research, subclass, ncells, disturbed.cells.natural1986:disturbed.cells.natural2020) %>%
  tidyr::gather(c(9:43), key="key", value="ncells.disturbed") %>%
  dplyr::mutate(Year = as.numeric(gsub("disturbed.cells.natural", "", key))) %>%
  dplyr::select(subclass, Research, elevation:coniferous.rate, ncells, Year, ncells.disturbed) %>%
  dplyr::right_join(data %>%
                      dplyr::select(elevation:coniferous.rate, Research, subclass, ncells, severe.disturbed.cells.natural1986:severe.disturbed.cells.natural2020) %>%
                      tidyr::gather(c(9:43), key="key", value="ncells.severe.disturbed") %>%
                      dplyr::mutate(Year = as.numeric(gsub("severe.disturbed.cells.natural", "", key))) %>%
                      dplyr::select(subclass, Research, elevation:coniferous.rate, ncells, Year, ncells.severe.disturbed)) %>%
  dplyr::right_join(data %>%
                      dplyr::select(elevation:coniferous.rate, Research, subclass, ncells, patch.size.natural1986:patch.size.natural2020) %>%
                      tidyr::gather(c(9:43), key="key", value="patch.size") %>%
                      dplyr::mutate(Year = as.numeric(gsub("patch.size.natural", "", key))) %>%
                      dplyr::select(subclass, Research, elevation:coniferous.rate, ncells, Year, patch.size)) %>%
  dplyr::group_by(subclass, Research, elevation, slope, topex, north.westerness, coniferous.rate, ncells) %>%
  dplyr::summarise(Frequency = sum(ncells.disturbed>0)/35,
                   Rate = sum(ncells.disturbed)/35,
                   Size = max(patch.size),
                   Severity = sum(ncells.severe.disturbed)/sum(ncells.disturbed)) %>%
  dplyr::mutate(Rate = Rate/ncells) 

diff.data <- all.data %>%
  dplyr::group_by(subclass) %>%
  dplyr::summarise(Frequency.diff = sum((-1)**ifelse(Research == "Uneven-aged", 0, 1) * Frequency),
                   Rate.diff = sum((-1)**ifelse(Research == "Uneven-aged", 0, 1) * Rate),
                   Size.diff = sum((-1)**ifelse(Research == "Uneven-aged", 0, 1) * Size),
                   Severity.diff = sum((-1)**ifelse(Research == "Uneven-aged", 0, 1) * Severity)) %>%
  tidyr::gather(-1, key="DisturbanceMetric", value="Value")


png("../graphs/Supplement/Figure_S6.png", res = 1000, width=15, height=12, units="cm")
ggplot(diff.data, aes(x=Value, fill=gsub(".diff", "", DisturbanceMetric)))+
  geom_histogram(bins=50, alpha=0.7)+
  geom_vline(xintercept = 0, linetype="dashed", col="grey20")+
  geom_vline(data=diff.data %>% dplyr::group_by(DisturbanceMetric) %>% dplyr::summarise(mean.value=mean(Value, na.rm=TRUE)),
             aes(xintercept = mean.value))+
  facet_wrap(~gsub(".diff", "", DisturbanceMetric), scales="free")+
  scale_fill_manual(values = metric.colors)+
  labs(x="Difference within sublandscape pairs", y="Count", fill="Disturbance metric")+
  theme_bw()
dev.off()

t.test(dplyr::filter(diff.data, DisturbanceMetric == "Frequency.diff")$Value, alternative = "two.sided")
t.test(dplyr::filter(diff.data, DisturbanceMetric == "Rate.diff")$Value, alternative = "two.sided")
t.test(dplyr::filter(diff.data, DisturbanceMetric == "Severity.diff")$Value, alternative = "two.sided")
t.test(dplyr::filter(diff.data, DisturbanceMetric == "Size.diff")$Value*30*30/10000, alternative = "two.sided")


#### -------------------------------- Figure 4 - disturbance rate over time -------------------------------- ####
png("../graphs/Figure_4.png", res = 1000, width=15, height=9, units="cm")
F1.data <- data %>% 
  dplyr::select(ncells, Site, Research, disturbed.cells.natural1986:disturbed.cells.natural2020) %>%
  tidyr::gather(-c(1,2,3), key="key", value="ncells.disturbed") %>%
  dplyr::mutate(year = as.numeric(gsub("disturbed.cells.natural", "", key)),
                DisturbanceRate = ncells.disturbed/ncells) %>% 
  dplyr::group_by(Site, Research, year) %>%
  dplyr::summarise(mean = mean(DisturbanceRate),
                   low = quantile(DisturbanceRate, probs = 0.10),
                   high = quantile(DisturbanceRate, probs = 0.90)) %>% 
  dplyr::mutate(Site = dplyr::case_when(Site == "aaa" ~ "Site C",
                                        Site == "bbb" ~ "Site D",
                                        Site == "ccc" ~ "Site A",
                                        Site == "ddd" ~ "Site B")) 
F1.data$Research <- factor(F1.data$Research, levels=c("Uneven-aged", "Even-aged"))

ggplot(F1.data, aes(x=year, col=Research))+
  geom_ribbon(data=F1.data, aes(ymin=low*100, ymax=high*100, fill=Research), alpha=0.3, linewidth=0.2)+
  geom_line(aes(y=mean*100, linetype=Research), linewidth=0.7)+ 
  scale_x_continuous(breaks=seq(1980, 2020, by=10))+
  scale_y_continuous(trans = "log1p", breaks=c(seq(0,1, 0.2), seq(2, 10, 2), 20), labels=c(0, rep("", 4), 1, rep("", 4), 10, 20))+ 
  scale_colour_manual(values=managament.colors)+
  scale_fill_manual(values=managament.colors)+
  theme_bw()+
  theme(legend.position = "top", strip.background = element_rect(fill="white"))+
  labs(x="Year", y=expression(paste('disturbance rate  [% year' ^'-1', ']')), 
       col="Management", fill="Management", linetype="Management") +
  facet_wrap(~Site)

dev.off()
rm(F1.data)




#### -------------------------------- Figure 5 - gams -------------------------------- ####
all.data <- data %>%
  dplyr::select(Site, elevation:coniferous.rate, Research, subclass, ncells, disturbed.cells.natural1986:disturbed.cells.natural2020) %>%
  tidyr::gather(c(10:44), key="key", value="ncells.disturbed") %>%
  dplyr::mutate(Year = as.numeric(gsub("disturbed.cells.natural", "", key))) %>%
  dplyr::select(Site, subclass, Research, elevation:coniferous.rate, ncells, Year, ncells.disturbed) %>%
  dplyr::right_join(data %>%
                      dplyr::select(elevation:coniferous.rate, Research, subclass, ncells, severe.disturbed.cells.natural1986:severe.disturbed.cells.natural2020) %>%
                      tidyr::gather(c(9:43), key="key", value="ncells.severe.disturbed") %>%
                      dplyr::mutate(Year = as.numeric(gsub("severe.disturbed.cells.natural", "", key))) %>%
                      dplyr::select(subclass, Research, elevation:coniferous.rate, ncells, Year, ncells.severe.disturbed)) %>%
  dplyr::right_join(data %>%
                      dplyr::select(elevation:coniferous.rate, Research, subclass, ncells, patch.size.natural1986:patch.size.natural2020) %>%
                      tidyr::gather(c(9:43), key="key", value="patch.size") %>%
                      dplyr::mutate(Year = as.numeric(gsub("patch.size.natural", "", key))) %>%
                      dplyr::select(subclass, Research, elevation:coniferous.rate, ncells, Year, patch.size)) %>%
  dplyr::group_by(Site, subclass, Research, elevation, slope, topex, north.westerness, coniferous.rate, ncells) %>%
  dplyr::summarise(Frequency = sum(ncells.disturbed>0)/35,
                   Rate = sum(ncells.disturbed)/35*100,
                   Size = max(patch.size)*30*30/10000,
                   Severity = sum(ncells.severe.disturbed)/sum(ncells.disturbed)*100) %>%
  dplyr::mutate(Rate = Rate/ncells) 
all.data$Research <- factor(all.data$Research, levels=c("Uneven-aged", "Even-aged"))

all.data <- as.data.frame(all.data)


denisty.data <- all.data %>%
  dplyr::select(Research:coniferous.rate) %>%
  tidyr::gather(2:6, key="Predictor", value="Values")


# Data collect function
collect.data <- function(Target, gam.model, mult.factor, inverse.link, smooth.data.frame, predictors, all.data) {
  for (predictor in predictors){
    df <- data.frame(elevation = rep(median(all.data$elevation), times=400),
                     slope = rep(median(all.data$slope), times=400),
                     north.westerness = rep(median(all.data$north.westerness), times=400),
                     coniferous.rate = rep(median(all.data$coniferous.rate), times=400),
                     Research = rep(c("Uneven-aged", "Even-aged"), each=200))
    
    df[,predictor] <- rep(seq(min(all.data[, predictor]), max(all.data[, predictor]), length=200), times=2)
    
    pred <- predict(gam.model, df, type = "link", se.fit = TRUE)
    pred <- cbind(pred, df)
    pred <- transform(pred,
                      fitted = inverse.link(fit)*mult.factor, 
                      lwr_ci = inverse.link(fit-se.fit)*mult.factor,
                      upr_ci = inverse.link(fit+se.fit)*mult.factor)
    
    
    tmp <- data.frame(Target=Target, 
                      Predictor=predictor,
                      x=df[,predictor],
                      Research=df$Research,
                      fit=pred$fitted,
                      lower=pred$lwr_ci,
                      upper=pred$upr_ci)
    smooth.data.frame <- rbind(smooth.data.frame, tmp)
  }
  return(smooth.data.frame)
}

predictors <- c("elevation", "slope", "north.westerness", "coniferous.rate")


smooth.data.frame <- data.frame()

# Rate 
gam.model <- mgcv::gam(I(Rate/100) ~ s(elevation, k=3, by=Research) + s(slope, k=3, by=Research) + s(north.westerness, k=3, by=Research) + s(coniferous.rate, k=3, by=Research), 
                       data=all.data, method = "REML", family = stats::quasibinomial())

inverse.link <- family(gam.model)$linkinv
smooth.data.frame <- collect.data("Rate", gam.model, 100, inverse.link, smooth.data.frame, predictors, all.data)



# Frequency
gam.model <- mgcv::gam(Frequency ~ s(elevation, k=3, by=Research) + s(slope, k=3, by=Research) + s(north.westerness, k=3, by=Research) + s(coniferous.rate, k=3, by=Research), 
                       data=all.data, method = "REML", family=quasibinomial())

inverse.link <- family(gam.model)$linkinv
smooth.data.frame <- collect.data("Frequency", gam.model, 1, inverse.link, smooth.data.frame, predictors, all.data)



# Size 
gam.model <- mgcv::gam(I(Size/50) ~ s(elevation, k=3, by=Research) + s(slope, k=3, by=Research) + s(north.westerness, k=3, by=Research) + s(coniferous.rate, k=3, by=Research), 
                       data=all.data, method = "REML", family=quasibinomial())

inverse.link <- family(gam.model)$linkinv
smooth.data.frame <- collect.data("Size", gam.model, 50, inverse.link, smooth.data.frame, predictors, all.data)



# Severity
gam.model <- mgcv::gam(Severity/100 ~ s(elevation, k=3, by=Research) + s(slope, k=3, by=Research) + s(north.westerness, k=3, by=Research) + s(coniferous.rate, k=3, by=Research), 
                       data=all.data, method = "REML", family=quasibinomial())

inverse.link <- family(gam.model)$linkinv
smooth.data.frame <- collect.data("Severity", gam.model, 100, inverse.link, smooth.data.frame, predictors, all.data)


smooth.data.frame$Target <- factor(smooth.data.frame$Target, levels = c("Rate", "Frequency", "Size", "Severity"))
smooth.data.frame$Predictor <- ifelse(smooth.data.frame$Predictor=="elevation", "Elevation", 
                                      ifelse(smooth.data.frame$Predictor=="slope", "Slope", 
                                             ifelse(smooth.data.frame$Predictor=="north.westerness", "North-Westerness", "Coniferous Rate")))

smooth.data.frame$Predictor <- factor(smooth.data.frame$Predictor, levels = c("Elevation", "Slope", "North-Westerness", "Coniferous Rate"))


denisty.data$Predictor <- ifelse(denisty.data$Predictor=="elevation", "Elevation", 
                                 ifelse(denisty.data$Predictor=="slope", "Slope", 
                                        ifelse(denisty.data$Predictor=="north.westerness", "North-Westerness", "Coniferous Rate")))


LocalParameters <- c("Elevation", "Slope", "North-Westerness", "Coniferous Rate"); LocalParameter <- LocalParameters[1]
DisturbanceMetrics <-c("Rate", "Frequency", "Size", "Severity"); DisturbanceMetric <- DisturbanceMetrics[1]

limits <- vector("list"); breaks <- vector("list"); names <- vector("list")
limits["Elevation"] <- list(c(200, 1800)); breaks["Elevation"] <- list(0:4*500); names["Elevation"] <- list(paste("", unlist(breaks["Elevation"]), "m"))
limits["Slope"] <- list(c(0, 35)); breaks["Slope"] <- list(0:4*10); names["Slope"] <- list(paste("    ", unlist(breaks["Slope"]), "°"))
limits["North-Westerness"] <- list(c(-1, 1)); breaks["North-Westerness"] <- list(-2:2/2); names["North-Westerness"] <- list(paste("   ", unlist(breaks["North-Westerness"]), ""))
limits["Coniferous Rate"] <- list(c(0.35, 1)); breaks["Coniferous Rate"] <- list(2:5*0.20); names["Coniferous Rate"] <- list(paste(" ", unlist(breaks["Coniferous Rate"])*100, "%"))

limits["Rate"] <- list(c(0, 0.6)); breaks["Rate"] <- list(0:3/5); names["Rate"] <- list(paste(unlist(breaks["Rate"]), "% p.a."))
limits["Frequency"] <- list(c(0, 0.3)); breaks["Frequency"] <- list(0:3/10); names["Frequency"] <- list(paste("   ", unlist(breaks["Frequency"]), "p.a."))
limits["Size"] <- list(c(0, 4)); breaks["Size"] <- list(0:4); names["Size"] <- list(paste("        ", unlist(breaks["Size"]), "ha"))
limits["Severity"] <- list(c(0, 80)); breaks["Severity"] <- list(0:3*25); names["Severity"] <- list(paste("       ", unlist(breaks["Severity"]), "%"))

plot.list <- vector("list")
i <- 1
for (DisturbanceMetric in DisturbanceMetrics){
  for (LocalParameter in LocalParameters){
    if (DisturbanceMetric == "Severity"){
      xlab <- LocalParameter
      xlabels <- unlist(names[LocalParameter])
    } else {
      xlab <- ""
      xlabels <- rep("", length(unlist(names[LocalParameter])))
    }
    
    if (LocalParameter == "Elevation"){
      ylab <- DisturbanceMetric
      ylabels <- unlist(names[DisturbanceMetric])
    } else {
      ylab <- ""
      ylabels <- rep("", length(unlist(names[DisturbanceMetric])))
    }
    
    xlimits <- c(unlist(limits[LocalParameter])[c(1, 2)])
    ylimits <- c(unlist(limits[DisturbanceMetric])[c(1, 2)])
    
    xbreaks <- unlist(breaks[LocalParameter])
    ybreaks <- unlist(breaks[DisturbanceMetric])
    
    
    plot.df <- smooth.data.frame %>%
      as.data.frame() %>%
      dplyr::filter(Target == DisturbanceMetric, Predictor == LocalParameter)
    
    plot <- ggplot(plot.df, aes(x=x, y=fit, col=Research, fill=Research))+
      geom_hline(yintercept = 0, linetype="dashed")+
      geom_ribbon(aes(ymin = lower, ymax = upper, y = NULL),
                  alpha = 0.3)+
      geom_line()+
      scale_x_continuous(breaks = xbreaks, labels=xlabels)+ 
      scale_y_continuous(breaks = ybreaks, labels=ylabels)+
      scale_colour_manual(values=rev(managament.colors))+
      scale_fill_manual(values=rev(managament.colors))+
      coord_cartesian(xlim=xlimits, ylim=ylimits)+
      labs(x=xlab, y=ylab, col="Management", fill="Management")+
      theme_bw()+
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
    
    if (DisturbanceMetric == "Rate") {
      tmp.data <- denisty.data %>%
        dplyr::filter(Predictor == LocalParameter)
      
      if (LocalParameter == "Elevation") {
        density.plot <- ggplot(tmp.data, aes(x=Values, col=Research, fill=Research))+
          geom_density(alpha=0.4, show.legend = FALSE)+
          scale_colour_manual(values=managament.colors)+
          scale_fill_manual(values=managament.colors)+
          labs(x="", y="")+
          theme_bw()+
          theme(strip.background = element_blank(),
                strip.text = element_blank(),
                axis.text = element_blank(),
                panel.border = element_blank(),
                axis.ticks = element_blank(),
                plot.margin = unit(c(5.5, 5.5, 5.5, 50.5), "pt"))+
          coord_cartesian(xlim=xlimits)
      } else {
        density.plot <- ggplot(tmp.data, aes(x=Values, col=Research, fill=Research))+
          geom_density(alpha=0.4, show.legend = FALSE)+
          scale_colour_manual(values=managament.colors)+
          scale_fill_manual(values=managament.colors)+
          labs(x="", y="")+
          theme_bw()+
          theme(strip.background = element_blank(),
                strip.text = element_blank(),
                axis.text = element_blank(),
                panel.border = element_blank(),
                axis.ticks = element_blank())+
          coord_cartesian(xlim=xlimits)
      }
      
      plot <- ggpubr::ggarrange(density.plot, plot,
                                ncol=1, nrow=2,
                                legend="none",
                                heights = c(1,2)) 
    }
    plot.list[[i]] <- plot
    i <- i+1
  }
}

png("../graphs/Figure_5.png", res = 1000, width=14, height=18, units="cm")
ggpubr::ggarrange(plotlist=plot.list[-c(4, 8, 12, 16)],
                  ncol=3, nrow=4,
                  common.legend = TRUE, legend="bottom",
                  widths = c(1.3,1,1),
                  heights = c(1.5,1,1,1.15)) 
dev.off()

##### Supplement #####

#Load Data
load(file = "../output/data_Sublandscapes_noTRI/Data_Sublandscapes_incl_mahabDist_10km.RData")
data <- export_df_including_dist; rm(export_df_including_dist)


max.dist <- max(data$match.dist)
data <- data %>% 
  dplyr::full_join(dplyr::select(dplyr::filter(data, Research==1), site, subclass), by="subclass") %>%
  dplyr::mutate(rel.dist = (1 - match.dist/max.dist) * 100,
                Site = site.y,
                Research = dplyr::case_when(Research == 1 ~ "Uneven-aged", TRUE ~ "Even-aged")) %>% 
  dplyr::relocate(Site, .before = elevation) %>%
  dplyr::select(-c(site.x, site.y))

# make super data by filter for matched sublandscapes with match.dist <= 1
data <- data %>% dplyr::filter(match.dist <= 1); rm(max.dist)


# Resample data, so that each site has same number of samples
set.seed(25) 
subclass.data <- data %>%
  dplyr::filter(Research == "Uneven-aged") %>%
  dplyr::group_by(Site) %>%
  dplyr::sample_n(1250, replace=FALSE) %>%
  dplyr::select(subclass)

data <- dplyr::filter(data, subclass %in% subclass.data$subclass)

# set colors
managament.colors <- c("#f0b030", "#cccccc")
metric.colors <- c("#E69F00", "#0072B2", "#D55E00", "#CC79A7")
even.aged.colors <- c("#f0b030", "#aaaaaa", "#cccccc", "#eeeeee") 



#### -------------------------------- Figure S1 -------------------------------- ####
site.raster <- terra::rast("../data/Betriebe/site_raster.tif")
elevation.raster <- terra::rast("../data/Elevation/elevation_map.tif")
slope.raster <- terra::rast("../data/Elevation/slope_map.tif")
conifer.raster <- terra::rast("../data/DominantLeafType/coniferous_rate_raster.tif")

# bind data
tmp.data <- terra::values(c(site.raster, elevation.raster, slope.raster, conifer.raster), dataframe=TRUE, na.rm=TRUE)%>%
  dplyr::mutate(Site = ifelse(site=="aaa", "Site A", 
                              ifelse(site=="bbb", "Site B",
                                     ifelse(site=="ccc", "Site C", 
                                            ifelse(site=="ddd", "Site D", "Else")))))


png("../graphs/Supplement/Figure_S1.png", res = 1000, width=15, height=15, units="cm")
ggplot(tmp.data, aes(x=slope, y=elevation, col=Site))+
  geom_vline(xintercept = 20)+
  geom_hline(yintercept = 1000)+
  geom_density_2d(show.legend = FALSE)+
  facet_wrap(~Site)+
  coord_cartesian(xlim = c(0,40), ylim=c(0,2000))+
  scale_x_continuous(breaks = seq(0,40,10), labels = paste0(seq(0,40,10), "°"))+
  scale_y_continuous(breaks = seq(0,2000,500), labels = paste0(seq(0,2000,500), "m"))+
  labs(x="Slope", y="Elevation")+
  theme_bw()
dev.off()

rm(site.raster, elevation.raster, slope.raster, conifer.raster, tmp.data)

#### -------------------------------- Figure S2 -------------------------------- ####
density.plot.data <- data %>%
  dplyr::select(elevation:coniferous.rate, Research, subclass, ncells) %>%
  tidyr::gather(c(1:5), key="LocalParameter", value="LocalParameterValue") %>%
  dplyr::mutate(LocalParameter = dplyr::case_when(LocalParameter == "coniferous.rate" ~ "Conifer Rate", 
                                                  LocalParameter == "north.westerness" ~ "North-Westerness",  
                                                  LocalParameter == "topex" ~ "Topographic Exposure",  
                                                  LocalParameter == "slope" ~ "Slope",  
                                                  LocalParameter == "elevation" ~ "Elevation"))
density.plot.data$LocalParameter <- factor(density.plot.data$LocalParameter, levels = c("Elevation", "Slope", "North-Westerness", "Topographic Exposure", "Conifer Rate"))
density.plot.data$Research <- factor(density.plot.data$Research, levels=c("Uneven-aged", "Even-aged"))

png("../graphs/Supplement/Figure_S2.png", res = 1000, width=15, height=10, units="cm")
elevation.plot <- ggplot(dplyr::filter(density.plot.data, LocalParameter == "Elevation"), aes(x=LocalParameterValue, fill=Research))+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=managament.colors)+
  scale_x_continuous(breaks = seq(0, 2000, by=500), labels = paste0(seq(0, 2000, by=500), "m"))+
  facet_wrap(~LocalParameter, scales="free")+
  labs(x="", fill="Management")+  
  theme_bw()+
  theme(legend.position = "top")

slope.plot <- ggplot(dplyr::filter(density.plot.data, LocalParameter == "Slope"), aes(x=LocalParameterValue, fill=Research))+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=managament.colors)+
  scale_x_continuous(breaks = seq(0, 40, by=10), labels = paste0(seq(0, 40, by=10), "°"))+
  facet_wrap(~LocalParameter, scales="free")+
  labs(x="", y= "", fill="Management")+  
  theme_bw()+
  theme(legend.position = "top")

northerness.plot <- ggplot(dplyr::filter(density.plot.data, LocalParameter == "North-Westerness"), aes(x=LocalParameterValue, fill=Research))+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=managament.colors)+
  scale_x_continuous(breaks = seq(-0.5, 0.5, by=0.5))+
  facet_wrap(~LocalParameter, scales="free")+
  labs(x="", y= "", fill="Management")+  
  theme_bw()+
  theme(legend.position = "top")

topex.plot <- ggplot(dplyr::filter(density.plot.data, LocalParameter == "Topographic Exposure"), aes(x=LocalParameterValue, fill=Research))+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=managament.colors)+
  scale_x_continuous(breaks = 4:8)+
  facet_wrap(~LocalParameter, scales="free")+
  labs(x="", fill="Management")+  
  theme_bw()+
  theme(legend.position = "top")

conifer.plot <- ggplot(dplyr::filter(density.plot.data, LocalParameter == "Conifer Rate"), aes(x=LocalParameterValue, fill=Research))+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=managament.colors)+
  scale_x_continuous(breaks = seq(0, 1, by=0.25), labels = paste0(seq(0, 100, by=25), "%"))+
  facet_wrap(~LocalParameter, scales="free")+
  labs(x="", y= "", fill="Management")+  
  theme_bw()+
  theme(legend.position = "top")

ggpubr::ggarrange(elevation.plot, slope.plot, northerness.plot, topex.plot, conifer.plot,
                  nrow=2, ncol=3, common.legend = TRUE, legend = "top")
dev.off()

rm(density.plot.data, elevation.plot, slope.plot, northerness.plot, topex.plot, conifer.plot)



#### -------------------------------- Figure S3 -------------------------------- ####
temp.data <- data %>% 
  dplyr::select(Site, Research, subclass, ncells, disturbed.cells.barkbeetle1986:patch.size.wind2020) %>%
  tidyr::gather(-c(1,2,3, 4), key="key", value="value") %>%
  dplyr::mutate(DisturbanceParameter = dplyr::case_when(as.logical(grepl("severe", key)) ~ "ncells.severe.disturbed",
                                                        as.logical(grepl("size", key)) ~ "Patch.Size",
                                                        TRUE ~ "ncells.disturbed"),
                agent = dplyr::case_when(grepl("wind", key) ~ "wind", 
                                         grepl("barkbeetle", key) ~ "barkbeetle", 
                                         grepl("natural", key) ~ "natural"),
                year = as.numeric(gsub("barkbeetle", "", 
                                       gsub("wind", "",  
                                            gsub("natural", "", 
                                                 gsub("disturbed.cells.", "", 
                                                      gsub("severe.disturbed.cells.", "",
                                                           gsub("patch.size.", "", key)))))))) %>%
  dplyr::select(-key) %>%
  tidyr::spread(DisturbanceParameter, value) %>%
  dplyr::mutate(DistEvent = ncells.disturbed > 0)

temp.data <- temp.data %>% 
  dplyr::right_join(temp.data %>% 
                      dplyr::select(Site:ncells.disturbed) %>%
                      tidyr::spread(key=agent, value=ncells.disturbed) %>%
                      dplyr::mutate(Wind.ratio = wind/(wind+barkbeetle)) %>%
                      dplyr::select(-c(barkbeetle, natural, wind)), by=c("Site", "Research", "subclass", "year", "ncells")) %>% 
  dplyr::group_by(Site, Research, agent, subclass, ncells) %>%
  dplyr::summarise(ncells.disturbed = sum(ncells.disturbed),
                   nEvents = sum(DistEvent),
                   ncells.severe.disturbed =sum(ncells.severe.disturbed),
                   Patch.Size = max(Patch.Size, na.rm=TRUE),
                   Wind.ratio = mean(Wind.ratio, na.rm=TRUE)) %>% 
  dplyr::mutate(annual.disturbance.rate = ncells.disturbed/ncells/35,
                Severity.rate = ncells.severe.disturbed/ncells.disturbed) %>%
  dplyr::select(-c(ncells, ncells.disturbed, ncells.severe.disturbed)) %>%
  tidyr::gather(5:9, key="DisturbanceMetric", value="Value") 
temp.data$Research <- factor(temp.data$Research, levels = c("Uneven-aged", "Even-aged"))
temp.data$DisturbanceMetricName <- dplyr::case_when(temp.data$DisturbanceMetric == "annual.disturbance.rate" ~ "Annual Disturbance Rate",
                                                    temp.data$DisturbanceMetric == "nEvents" ~ "Annual Number of Events",
                                                    temp.data$DisturbanceMetric == "Severity.rate" ~ "Severity Rate",
                                                    temp.data$DisturbanceMetric == "Patch.Size" ~ "Patch Size")



png("../graphs/Figure_S3.png", res = 1000, width=15, height=9, units="cm")
ggplot(dplyr::filter(temp.data, DisturbanceMetric == "Wind.ratio", agent=="natural"), aes(x=Research, y=Value*100, fill=Research))+
  geom_violin(draw_quantiles = 0.5, show.legend=FALSE)+
  geom_point(data=dplyr::filter(temp.data, DisturbanceMetric == "Wind.ratio", agent=="natural") %>%
               dplyr::group_by(Research) %>%
               dplyr::summarise(mean=mean(Value, na.rm=TRUE)), aes(x=Research, y=mean*100), show.legend=FALSE)+
  scale_y_continuous(breaks = seq(0, 100, by = 25), labels = paste0("   ", seq(0, 100, by = 25), "%"))+
  scale_fill_manual(values=managament.colors)+
  labs(y="Proportion of wind disturbance [%]", fill="Management")+
  theme_bw()
dev.off()

quantile(dplyr::filter(temp.data, DisturbanceMetric == "Wind.ratio", agent=="natural", Research == "Uneven-aged")$Value*100, probs=c(0.025, 0.1, 0.5, 0.9, 0.975), na.rm=TRUE)
quantile(dplyr::filter(temp.data, DisturbanceMetric == "Wind.ratio", agent=="natural", Research == "Even-aged")$Value*100, probs=c(0.025, 0.1, 0.5, 0.9, 0.975), na.rm=TRUE)

mean(dplyr::filter(temp.data, DisturbanceMetric == "Wind.ratio", agent=="natural", Research == "Uneven-aged")$Value*100, na.rm=TRUE)
mean(dplyr::filter(temp.data, DisturbanceMetric == "Wind.ratio", agent=="natural", Research == "Even-aged")$Value*100, na.rm=TRUE)



#### -------------------------------- Figure S4 - compare searching distances -------------------------------- ####
# Load Data 
data <- NULL
load(file = "../output/data_Sublandscapes_noTRI/Data_Sublandscapes_incl_mahabDist_5km.RData")
data <- rbind(data, cbind(searching.distance = 10, export_df_including_dist)); rm(export_df_including_dist)
load(file = "../output/data_Sublandscapes_noTRI/Data_Sublandscapes_incl_mahabDist_10km.RData")
data <- rbind(data, cbind(searching.distance = 20, export_df_including_dist)); rm(export_df_including_dist)
load(file = "../output/data_Sublandscapes_noTRI/Data_Sublandscapes_incl_mahabDist_15km.RData")
data <- rbind(data, cbind(searching.distance = 30, export_df_including_dist)); rm(export_df_including_dist)

data$searching.distance <- factor(data$searching.distance, levels = c(30, 20, 10))

png("../graphs/Supplement/Figure_S4_A.png", res = 1000, width=14, height=6, units="cm")
ggplot(data, aes(x=match.dist, y=as.factor(searching.distance), fill=as.factor(searching.distance)))+
  geom_boxplot(outlier.shape = NA, alpha=0.9)+
  geom_point(data = data %>%
               dplyr::group_by(searching.distance) %>%
               dplyr::summarise(mean=mean(match.dist, na.rm=TRUE)), aes(x=mean, y=as.factor(searching.distance)), position = position_dodge2(0.75))+
  coord_cartesian(xlim = c(NA, 2.5))+
  scale_x_continuous(breaks = seq(0,5,by=0.5))+
  labs(x="Similarity",
       y="Search window [km]",
       fill="Search window")+
  theme_bw()+
  geom_vline(xintercept = 1)+
  scale_fill_manual(breaks = c(10, 20, 30), values=even.aged.colors[2:4]) 
dev.off()



# Compare difference in disturbance metrics for different searching distances
png("../graphs/Supplement/Figure_S4_B.png", res = 1000, width=15, height=10, units="cm")
set.seed(25) #5
subclass.data <- data %>%
  dplyr::filter(match.dist <= 1) %>%
  dplyr::filter(Research == 1) %>%
  dplyr::group_by(searching.distance, site) %>% 
  dplyr::sample_n(1250, replace=FALSE) %>%
  dplyr::select(subclass) %>%
  dplyr::mutate(Use = 1)

comparison_data <- dplyr::right_join(data, subclass.data, by=c("subclass", "searching.distance")) %>%
  dplyr::filter(!is.na(Use))%>%
  dplyr::mutate(Site = site.y) %>%
  dplyr::select(-site.x, -site.y, -Use) %>%
  dplyr::filter(Research == 0 | searching.distance == 20) %>%
  dplyr::mutate(Research = ifelse(Research == 1, "uneven-aged", paste0("surrounding_", searching.distance, "km"))) %>%
  dplyr::select(Research, ncells, subclass:patch.size.wind2020)
comparison_data$Research <- factor(comparison_data$Research, levels = c("uneven-aged", "surrounding_10km", "surrounding_20km", "surrounding_30km"))

comparison_data <- comparison_data %>% 
  dplyr::select(Research, subclass, ncells) %>%
  dplyr::right_join(comparison_data %>%
                      dplyr::select(Research:subclass, disturbed.cells.natural1986:disturbed.cells.natural2020) %>%
                      tidyr::gather(-c(1:3), key="key", value="value") %>% 
                      dplyr::mutate(year=as.numeric(gsub("disturbed.cells.natural", "", key))) %>% 
                      dplyr::group_by(Research, subclass) %>%
                      dplyr::summarise(ncells.disturbed = sum(value, na.rm=TRUE))) %>%
  dplyr::right_join(comparison_data %>%
                      dplyr::select(Research:subclass, disturbed.cells.natural1986:disturbed.cells.natural2020) %>%
                      tidyr::gather(-c(1:3), key="key", value="value") %>% 
                      dplyr::mutate(year=as.numeric(gsub("disturbed.cells.natural", "", key))) %>% 
                      dplyr::group_by(Research, subclass) %>%
                      dplyr::summarise(n.disturbance.events = sum(value>0, na.rm=TRUE))) %>%
  dplyr::right_join(comparison_data %>%
                      dplyr::select(Research:subclass, severe.disturbed.cells.natural1986:severe.disturbed.cells.natural2020) %>%
                      tidyr::gather(-c(1:3), key="key", value="value") %>% 
                      dplyr::mutate(year=as.numeric(gsub("severe.disturbed.cells.natural", "", key))) %>%
                      dplyr::group_by(Research, subclass) %>%
                      dplyr::summarise(ncells.severe.disturbed = sum(value, na.rm=TRUE))) %>%
  dplyr::right_join(comparison_data %>%
                      dplyr::select(Research:subclass, patch.size.natural1986:patch.size.natural2020) %>%
                      tidyr::gather(-c(1:3), key="key", value="value") %>% 
                      dplyr::mutate(year=as.numeric(gsub("patch.size.natural", "", key))) %>%
                      dplyr::group_by(Research, subclass) %>%
                      dplyr::summarise(patch.size = max(value, na.rm=TRUE))) %>%
  dplyr::mutate(disturbance.rate = ncells.disturbed/(35*ncells)*100,
                high.severity.share = ncells.severe.disturbed/ncells.disturbed*100,
                patch.size = patch.size * 30*30/10000) %>%
  dplyr::select(Research, subclass, disturbance.rate, n.disturbance.events, high.severity.share, patch.size)

comparison_data$main <- ifelse(comparison_data$Research %in% c("uneven-aged", "surrounding_20km"), "Yes", "No")
comparison_data$main <- factor(comparison_data$main, levels = c("Yes", "No"))

comparison1 <- ggplot(comparison_data, aes(x=Research, y=disturbance.rate, fill=Research))+
  geom_boxplot(outlier.shape = NA, alpha=0.7, aes(col=main))+
  geom_point(data = comparison_data %>%
               dplyr::group_by(Research) %>%
               dplyr::summarise(mean=mean(disturbance.rate, na.rm=TRUE)), aes(x=Research, y=mean), position = position_dodge2(0.75))+
  scale_x_discrete(labels=NULL)+
  coord_cartesian(ylim = c(0, 0.6))+
  scale_y_continuous(breaks = seq(0,1,by=0.1), labels = paste(seq(0,1,by=0.1), "%"))+
  scale_fill_manual(values=even.aged.colors)+ 
  scale_color_manual(values = c("#000000", "#888888"))+
  labs(x="",
       y="Rate [% p.a.]",
       col="Used in study",
       fill="Search window")+
  theme_bw()+
  theme(legend.position = c(0.8, 0.8))

comparison2 <- ggplot(comparison_data, aes(x=Research, y=patch.size, fill=Research))+
  geom_boxplot(outlier.shape = NA, alpha=0.7, aes(col=main))+
  geom_point(data = comparison_data %>%
               dplyr::group_by(Research) %>%
               dplyr::summarise(mean=mean(patch.size, na.rm=TRUE)), aes(x=Research, y=mean), position = position_dodge2(0.75))+
  scale_x_discrete(labels=NULL)+
  coord_cartesian(ylim = c(0, 3))+
  scale_y_continuous(breaks = seq(0,10,by=0.5), labels = paste(seq(0,10,by=0.5), "ha"))+
  scale_fill_manual(values=even.aged.colors)+
  scale_color_manual(values = c("#000000", "#888888"))+
  labs(x="",
       y="Size [ha]",
       col="Used in study",
       fill="Search window")+
  theme_bw()+
  theme(legend.position = c(0.8, 0.8))

comparison3 <- ggplot(comparison_data, aes(x=Research, y=n.disturbance.events/35, fill=Research))+
  geom_boxplot(outlier.shape = NA, alpha=0.7, aes(col=main))+
  geom_point(data = comparison_data %>%
               dplyr::group_by(Research) %>%
               dplyr::summarise(mean=mean(n.disturbance.events, na.rm=TRUE)), aes(x=Research, y=mean/35), position = position_dodge2(0.75))+
  scale_x_discrete(labels=NULL)+
  coord_cartesian(ylim = c(0, 0.2))+
  scale_y_continuous(breaks = seq(0, 0.3, by=0.1))+
  scale_fill_manual(values=even.aged.colors)+ 
  scale_color_manual(values = c("#000000", "#888888"))+
  labs(x="",
       y="Frequency [1/year]",
       col="Used in study",
       fill="Search window")+
  theme_bw()+
  theme(legend.position = c(0.8, 0.8))

comparison4 <- ggplot(comparison_data, aes(x=Research, y=high.severity.share, fill=Research))+
  geom_boxplot(outlier.shape = NA, alpha=0.7, aes(col=main))+
  geom_point(data = comparison_data %>%
               dplyr::group_by(Research) %>%
               dplyr::summarise(mean=mean(high.severity.share, na.rm=TRUE)), aes(x=Research, y=mean), position = position_dodge2(0.75))+
  scale_x_discrete(labels=NULL)+
  coord_cartesian(ylim = c(0, 100))+
  scale_y_continuous(breaks = seq(0,100,by=25), labels = paste(seq(0,100,by=25), "%"))+
  scale_fill_manual(values=even.aged.colors)+ 
  scale_color_manual(values = c("#000000", "#888888"))+
  labs(x="",
       y="Severity [%]",
       col="Used in study",
       fill="Search window")+
  theme_bw()+
  theme(legend.position = c(0.8, 0.8))

ggpubr::ggarrange(comparison1, comparison2, comparison3, comparison4, nrow=2, ncol=2, common.legend = TRUE, legend="right")
dev.off()
rm(comparison_data, comparison1, comparison2, comparison3, comparison4)



#### -------------------------------- Figure S5 - compare mahalanobis threshold -------------------------------- ####
load(file = "../output/data_Sublandscapes_noTRI/Data_Sublandscapes_incl_mahabDist_10km.RData")
data <- export_df_including_dist; rm(export_df_including_dist)

data <- data %>% 
  dplyr::full_join(dplyr::select(dplyr::filter(data, Research==1), site, subclass), by="subclass") %>%
  dplyr::mutate(Site = site.y,
                Research = dplyr::case_when(Research == 1 ~ "Uneven-aged", TRUE ~ "Surrounding")) %>% 
  dplyr::relocate(Site, .before = elevation) %>%
  dplyr::select(-c(site.x, site.y))

data <- rbind(data.frame(threshold = 1, 
                         dplyr::filter(data, match.dist <= 1)),
              data.frame(threshold = 0.2, 
                         dplyr::filter(data, match.dist <= 0.2)),
              data.frame(threshold = 0.4, 
                         dplyr::filter(data, match.dist <= 0.4)),
              data.frame(threshold = 0.6, 
                         dplyr::filter(data, match.dist <= 0.6)),
              data.frame(threshold = 0.8, 
                         dplyr::filter(data, match.dist <= 0.8)),
              data.frame(threshold = 1.2, 
                         dplyr::filter(data, match.dist <= 1.2)),
              data.frame(threshold = 1.4, 
                         dplyr::filter(data, match.dist <= 1.4)),
              data.frame(threshold = 1.6, 
                         dplyr::filter(data, match.dist <= 1.6)),
              data.frame(threshold = 1.8, 
                         dplyr::filter(data, match.dist <= 1.8)),
              data.frame(threshold = 2, 
                         dplyr::filter(data, match.dist <= 2)))


comparison_data <- data %>% 
  dplyr::select(Research, subclass, ncells, threshold) %>%
  dplyr::right_join(data %>%
                      dplyr::select(Research:subclass, threshold, disturbed.cells.natural1986:disturbed.cells.natural2020) %>%
                      tidyr::gather(-c(1:4), key="key", value="value") %>% 
                      dplyr::mutate(year=as.numeric(gsub("disturbed.cells.natural", "", key))) %>% 
                      dplyr::group_by(Research, subclass, threshold) %>%
                      dplyr::summarise(ncells.disturbed = sum(value, na.rm=TRUE))) %>%
  dplyr::right_join(data %>%
                      dplyr::select(Research:subclass, threshold, disturbed.cells.natural1986:disturbed.cells.natural2020) %>%
                      tidyr::gather(-c(1:4), key="key", value="value") %>% 
                      dplyr::mutate(year=as.numeric(gsub("disturbed.cells.natural", "", key))) %>% 
                      dplyr::group_by(Research, subclass, threshold) %>%
                      dplyr::summarise(n.disturbance.events = sum(value>0, na.rm=TRUE))) %>%
  dplyr::right_join(data %>%
                      dplyr::select(Research:subclass, threshold, severe.disturbed.cells.natural1986:severe.disturbed.cells.natural2020) %>%
                      tidyr::gather(-c(1:4), key="key", value="value") %>% 
                      dplyr::mutate(year=as.numeric(gsub("severe.disturbed.cells.natural", "", key))) %>%
                      dplyr::group_by(Research, subclass, threshold) %>%
                      dplyr::summarise(ncells.severe.disturbed = sum(value, na.rm=TRUE))) %>%
  dplyr::right_join(data %>%
                      dplyr::select(Research:subclass, threshold, patch.size.natural1986:patch.size.natural2020) %>%
                      tidyr::gather(-c(1:4), key="key", value="value") %>% 
                      dplyr::mutate(year=as.numeric(gsub("patch.size.natural", "", key))) %>%
                      dplyr::group_by(Research, subclass, threshold) %>%
                      dplyr::summarise(patch.size = max(value, na.rm=TRUE))) %>%
  dplyr::mutate(disturbance.rate = ncells.disturbed/(35*ncells)*100,
                n.Events.per.year = n.disturbance.events/35,
                high.severity.share = ncells.severe.disturbed/ncells.disturbed*100,
                patch.size = patch.size * 30*30/10000) %>%
  dplyr::select(Research, subclass, threshold, disturbance.rate, n.Events.per.year, high.severity.share, patch.size)

comparison_data$main <- ifelse(comparison_data$threshold == 1, "Yes", "No")
comparison_data$main <- factor(comparison_data$main, levels = c("Yes", "No"))
comparison_data$threshold <- factor(comparison_data$threshold, levels = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2))
comparison_data$Research <- factor(comparison_data$Research, levels = c("Uneven-aged", "Surrounding"))

comparison_data <- comparison_data %>%
  tidyr::gather(4:7, key="DisturbanceMetric", value="Value")


rate.plot <- ggplot(dplyr::filter(comparison_data, DisturbanceMetric=="disturbance.rate", threshold %in% c(0.2, 0.6, 1, 1.4, 1.8)), aes(x=threshold, y=Value, fill=Research))+
  geom_boxplot(outlier.shape = NA, aes(col=main))+
  geom_point(data = dplyr::filter(comparison_data, DisturbanceMetric=="disturbance.rate", threshold %in% c(0.2, 0.6, 1, 1.4, 1.8)) %>%
               dplyr::group_by(Research, threshold) %>%
               dplyr::summarise(mean=mean(Value)), aes(x=threshold, y=mean), position = position_dodge2(0.75))+
  coord_cartesian(ylim = c(0,0.6))+
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = paste0(seq(0, 1, 0.2), "%"))+
  scale_color_manual(values = c("#000000", "#888888"))+
  scale_fill_manual(values=managament.colors)+ 
  labs(x="", y="Rate [% p.a.]", col="Used in study", fill="Management")+
  theme_bw()

n.plot <- ggplot(dplyr::filter(comparison_data, DisturbanceMetric=="n.Events.per.year", threshold %in% c(0.2, 0.6, 1, 1.4, 1.8)), aes(x=threshold, y=Value, fill=Research))+
  geom_boxplot(outlier.shape = NA, aes(col=main))+
  geom_point(data = dplyr::filter(comparison_data, DisturbanceMetric=="n.Events.per.year", threshold %in% c(0.2, 0.6, 1, 1.4, 1.8)) %>%
               dplyr::group_by(Research, threshold) %>%
               dplyr::summarise(mean=mean(Value)), aes(x=threshold, y=mean), position = position_dodge2(0.75))+
  coord_cartesian(ylim = c(0,0.2))+
  scale_color_manual(values = c("#000000", "#888888"))+
  scale_fill_manual(values=managament.colors)+ 
  labs(x="", y="Frequency [1/year]", col="Used in study", fill="Management")+
  theme_bw()

size.plot <- ggplot(dplyr::filter(comparison_data, DisturbanceMetric=="patch.size", threshold %in% c(0.2, 0.6, 1, 1.4, 1.8)), aes(x=threshold, y=Value, fill=Research))+
  geom_boxplot(outlier.shape = NA, aes(col=main))+
  geom_point(data = dplyr::filter(comparison_data, DisturbanceMetric=="patch.size", threshold %in% c(0.2, 0.6, 1, 1.4, 1.8)) %>%
               dplyr::group_by(Research, threshold) %>%
               dplyr::summarise(mean=mean(Value)), aes(x=threshold, y=mean), position = position_dodge2(0.75))+
  coord_cartesian(ylim = c(0,3))+
  scale_y_continuous(breaks = 0:4, labels = paste0(0:4, " ha"))+
  scale_color_manual(values = c("#000000", "#888888"))+
  scale_fill_manual(values=managament.colors)+ 
  labs(x = "", y="Size [ha]", col="Used in study", fill="Management")+
  theme_bw()

severity.plot <- ggplot(dplyr::filter(comparison_data, DisturbanceMetric=="high.severity.share", threshold %in% c(0.2, 0.6, 1, 1.4, 1.8)), aes(x=threshold, y=Value, fill=Research))+
  geom_boxplot(outlier.shape = NA, aes(col=main))+
  geom_point(data = dplyr::filter(comparison_data, DisturbanceMetric=="high.severity.share", threshold %in% c(0.2, 0.6, 1, 1.4, 1.8)) %>%
               dplyr::group_by(Research, threshold) %>%
               dplyr::summarise(mean=mean(Value, na.rm=TRUE)), aes(x=threshold, y=mean), position = position_dodge2(0.75))+
  coord_cartesian(ylim = c(0,100))+
  scale_y_continuous(breaks = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"))+
  scale_color_manual(values = c("#000000", "#888888"))+
  scale_fill_manual(values=managament.colors)+ 
  labs(x="Mahalanobis threshold", y="Severity [%]", col="Used in study", fill="Management")+
  theme_bw()


png("../graphs/Supplement/Figure_S5.png", res = 1000, width=15, height=15, units="cm")
ggpubr::ggarrange(rate.plot, n.plot, size.plot, severity.plot, 
                  nrow=4, ncol=1, common.legend = TRUE)
dev.off()



#### -------------------------------- Figure S7 -------------------------------- ####
load(file = "../output/data_Sublandscapes_noTRI/Data_Sublandscapes_incl_mahabDist_10km.RData")
data <- export_df_including_dist; rm(export_df_including_dist)


max.dist <- max(data$match.dist)
data <- data %>% 
  dplyr::full_join(dplyr::select(dplyr::filter(data, Research==1), site, subclass), by="subclass") %>%
  dplyr::mutate(rel.dist = (1 - match.dist/max.dist) * 100,
                Site = site.y,
                Research = dplyr::case_when(Research == 1 ~ "Uneven-aged", TRUE ~ "Even-aged")) %>% 
  dplyr::relocate(Site, .before = elevation) %>%
  dplyr::select(-c(site.x, site.y))

# make super data by filter for matched sublandscapes with match.dist <= 1
data <- data %>% dplyr::filter(match.dist <= 1); rm(max.dist)


# Resample data, so that each site has same number of samples
set.seed(25) 
subclass.data <- data %>%
  dplyr::filter(Research == "Uneven-aged") %>%
  dplyr::group_by(Site) %>%
  dplyr::sample_n(1250, replace=FALSE) %>%
  dplyr::select(subclass)

data <- dplyr::filter(data, subclass %in% subclass.data$subclass)

all.data <- data %>%
  dplyr::select(Site, elevation:coniferous.rate, Research, subclass, ncells, disturbed.cells.natural1986:disturbed.cells.natural2020) %>%
  tidyr::gather(c(10:44), key="key", value="ncells.disturbed") %>%
  dplyr::mutate(Year = as.numeric(gsub("disturbed.cells.natural", "", key))) %>%
  dplyr::select(Site, subclass, Research, elevation:coniferous.rate, ncells, Year, ncells.disturbed) %>%
  dplyr::right_join(data %>%
                      dplyr::select(elevation:coniferous.rate, Research, subclass, ncells, severe.disturbed.cells.natural1986:severe.disturbed.cells.natural2020) %>%
                      tidyr::gather(c(9:43), key="key", value="ncells.severe.disturbed") %>%
                      dplyr::mutate(Year = as.numeric(gsub("severe.disturbed.cells.natural", "", key))) %>%
                      dplyr::select(subclass, Research, elevation:coniferous.rate, ncells, Year, ncells.severe.disturbed)) %>%
  dplyr::right_join(data %>%
                      dplyr::select(elevation:coniferous.rate, Research, subclass, ncells, patch.size.natural1986:patch.size.natural2020) %>%
                      tidyr::gather(c(9:43), key="key", value="patch.size") %>%
                      dplyr::mutate(Year = as.numeric(gsub("patch.size.natural", "", key))) %>%
                      dplyr::select(subclass, Research, elevation:coniferous.rate, ncells, Year, patch.size)) %>%
  dplyr::group_by(Site, subclass, Research, elevation, slope, topex, north.westerness, coniferous.rate, ncells) %>%
  dplyr::summarise(Frequency = sum(ncells.disturbed>0)/35,
                   Rate = sum(ncells.disturbed)/35*100,
                   Size = max(patch.size)*30*30/10000,
                   Severity = sum(ncells.severe.disturbed)/sum(ncells.disturbed)*100) %>%
  dplyr::mutate(Rate = Rate/ncells) 
all.data$Research <- factor(all.data$Research, levels=c("Uneven-aged", "Even-aged"))

all.data <- as.data.frame(all.data)

png("../graphs/Supplement/Figure_S7.png", res = 1000, width=15, height=10, units="cm")
par(mfrow=c(2,4))

gam.model <- mgcv::gam(I(Rate/100) ~ s(elevation, k=3, by=Research) + s(slope, k=3, by=Research) + s(north.westerness, k=3, by=Research) + s(coniferous.rate, k=3, by=Research), 
                       data=all.data, method = "REML", family = stats::quasibinomial()); summary(gam.model)
hist(mgcv::residuals.gam(gam.model), main="Histogram of residuals \n Rate", xlab="Residuals", nclass=10)
mgcv::qq.gam(gam.model, main="Normal Q-Q Plot \n Rate")

gam.model <- mgcv::gam(Frequency ~ s(elevation, k=3, by=Research) + s(slope, k=3, by=Research) + s(north.westerness, k=3, by=Research) + s(coniferous.rate, k=3, by=Research), 
                       data=all.data, method = "REML", family=quasibinomial()); summary(gam.model)
hist(mgcv::residuals.gam(gam.model), main="Histogram of residuals \n Frequency", xlab="Residuals", nclass=10)
mgcv::qq.gam(gam.model, main="Normal Q-Q Plot \n Frequency")

gam.model <- mgcv::gam(I(Size/50) ~ s(elevation, k=3, by=Research) + s(slope, k=3, by=Research) + s(north.westerness, k=3, by=Research) + s(coniferous.rate, k=3, by=Research), 
                       data=all.data, method = "REML", family=quasibinomial()); summary(gam.model)
hist(mgcv::residuals.gam(gam.model), main="Histogram of residuals \n Size", xlab="Residuals", nclass=10)
mgcv::qq.gam(gam.model, main="Normal Q-Q Plot \n Size")

gam.model <- mgcv::gam(Severity/100 ~ s(elevation, k=3, by=Research) + s(slope, k=3, by=Research) + s(north.westerness, k=3, by=Research) + s(coniferous.rate, k=3, by=Research), 
                       data=all.data, method = "REML", family=quasibinomial()); summary(gam.model)
hist(mgcv::residuals.gam(gam.model), main="Histogram of residuals \n Severity", xlab="Residuals", nclass=10)
mgcv::qq.gam(gam.model, main="Normal Q-Q Plot \n Severity")

dev.off()


#### -------------------------------- Figure S8 - MonteCarlo Simulation -------------------------------- ####
load(file = "../output_clean/data_Sublandscapes/Data_Sublandscapes_incl_mahabDist_10km.RData")
data <- export_df_including_dist; rm(export_df_including_dist)


max.dist <- max(data$match.dist)
data <- data %>% 
  dplyr::full_join(dplyr::select(dplyr::filter(data, Research==1), site, subclass), by="subclass") %>%
  dplyr::mutate(rel.dist = (1 - match.dist/max.dist) * 100,
                Site = site.y,
                Research = dplyr::case_when(Research == 1 ~ "Uneven-aged", TRUE ~ "Even-aged")) %>% 
  dplyr::relocate(Site, .before = elevation) %>%
  dplyr::select(-c(site.x, site.y))

data <- data %>% dplyr::filter(match.dist <= 1); rm(max.dist)
data$Research <- factor(data$Research, levels = c("Uneven-aged", "Even-aged"))


# Resample data, so that each site has same number of samples
set.seed(25) #5
subclass.data <- data %>%
  dplyr::filter(Research == "Uneven-aged") %>%
  dplyr::group_by(Site) %>%
  dplyr::sample_n(1250, replace=FALSE) %>%
  dplyr::select(subclass)

data <- dplyr::filter(data, subclass %in% subclass.data$subclass)

names(data)

wrong.positive <- 1-0.740 # predicted natural is actually harvest (data from Sebald et al. 2021)
wrong.negative <- 1-0.676 # predicted harvest is actually natural (data from Sebald et al. 2021)

options(dplyr.summarise.inform = FALSE)

runs <- 1000
pb <- dplyr::progress_estimated(runs)
MonteCarlo.df <- data.frame()
for (run in 1:runs) {
  set.seed(run)
  tmp.data <- data
  for (sublandscape in 1:nrow(data)){
    # Check wrong positive
    for (column in 83:117){ 
      n.positive <- tmp.data[sublandscape,column]
      if (n.positive != 0) {
        n.wrong.positive <- sum(runif(n.positive, 0, 1) <= wrong.positive)
        tmp.data[sublandscape,column] <- n.positive - n.wrong.positive 
      } 
    }
    
    # Check wrong negative
    for (column in 48:82){ 
      n.negative <- tmp.data[sublandscape,column]
      if (n.negative != 0) {
        n.wrong.negative <- sum(runif(n.negative, 0, 1) <= wrong.negative)
        tmp.data[sublandscape,column+35] <- tmp.data[sublandscape,column+35] + n.wrong.negative 
      } 
    }
  }
  
  tmp.data <- tmp.data %>% 
    dplyr::select(Site, Research, subclass, ncells, disturbed.cells.natural1986:disturbed.cells.natural2020) %>%
    tidyr::gather(-c(1,2,3, 4), key="key", value="value") %>%
    dplyr::mutate(DisturbanceParameter = "ncells.disturbed",
                  agent =  "natural",
                  year = as.numeric(gsub("disturbed.cells.natural", "", key))) %>%
    dplyr::select(-key) %>%
    tidyr::spread(DisturbanceParameter, value) %>%
    dplyr::mutate(DistEvent = ncells.disturbed > 0) %>% 
    dplyr::group_by(Site, Research, agent, subclass, ncells) %>%
    dplyr::summarise(ncells.disturbed = sum(ncells.disturbed),
                     nEvents = sum(DistEvent)) %>% 
    dplyr::mutate(annual.disturbance.rate = ncells.disturbed/ncells/35) %>%
    dplyr::select(-c(ncells, ncells.disturbed)) 
  
  MonteCarlo.df <- rbind(MonteCarlo.df,
                         tmp.data %>%
                           dplyr::group_by(Research) %>%
                           dplyr::summarise(DisturbanceRate = mean(annual.disturbance.rate),
                                            Frequency = mean(nEvents)/35) %>%
                           dplyr::mutate(MCRun = run) %>%
                           dplyr::relocate(MCRun))
  pb$tick()$print()
}


managament.colors <- c("#f0b030", "#cccccc")

png("../graphs/Supplement/Figure_S8.png", res = 1000, width=15, height=9, units="cm")
MonteCarlo.df %>%
  dplyr::group_by(Research) %>%
  dplyr::summarise(ymin = min(DisturbanceRate),
                   ymean = mean(DisturbanceRate),
                   ymax = max(DisturbanceRate)) %>%
  dplyr::mutate(label=paste0("variation (max-min): \n", round((ymax-ymin)*100, 4), "% points")) %>%
  ggplot(aes(x=Research, y=ymean*100, ymin=ymin*100, ymax=ymax*100, col=Research, label=label))+
  geom_errorbar(width=0.2)+
  geom_point()+
  geom_label(nudge_y = -0.03, show.legend = FALSE)+
  coord_cartesian(ylim=c(0, 0.25))+ 
  scale_colour_manual(values=managament.colors)+
  labs(x="", y="Disturbance Rate [% p.a.]")+
  theme_bw()
dev.off()
