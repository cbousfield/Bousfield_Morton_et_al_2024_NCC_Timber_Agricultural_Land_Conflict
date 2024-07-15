#### Code for Bousfield, Morton et al. 2024, Nature Climate Change ####

### Key steps in the analysis are as follows ###

### 1. Data prep - load in agricultural suitability layers from Zabel et al. and forestry data from Curtis et al.
###              - edit Curtis forestry data to only include areas with >10% tree cover in 2000, and not mapped as current cropland, urban settlements, bare rock or water in the ESA land cover data, then reproject to match Zabel et al. projection/resolution
### 
### 2. Spatial overlays - overlay forestry land rasters with agricultural change rasters from Zabel et al. to get predictions of future change in agricultural suitability in current forestry land
###
### 3. Global suitability area - extract area for each suitability category globally through time and RCPs
###
### 4. Global gain/loss/persist area - extract the area that has become or lost productivity globally through time and RCPs
###
### 5. Travel time to population centers summary - extract the area weighted travel time for forestry and non forestry land.
###
### 6. Distance to agriculture summary - extract the area weighted travel distance to current agriculture for forestry and non forestry land.
### 
### 7. Country classified area - extract area for each land cover.
###
### 8. All countries suitability - extract area for each suitability category per country through time and RCPs
###
### 9. All countries productive gain/loss/persist - extract the area that has became or lost productivity globally through time and RCPs
###
### 10. All countries change - extract the area improving or decreasing per country per land cover per rcp
###
### 11. Country frontier area - calculate the area of newly productive land per country
###
###Data sources:  Zabel et al. - https://zenodo.org/records/5982577), 
###                Curtis et al. - https://data.globalforestwatch.org/documents/gfw::tree-cover-loss-by-dominant-driver-2022/about
###                Forest cover - https://storage.googleapis.com/earthenginepartners-hansen/GFC-2022-v1.10/download.html
###                ESA land cover - https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-land-cover?tab=form
###                Travel time to urban area - https://figshare.com/articles/dataset/Travel_time_to_cities_and_ports_in_the_year_2015/7638134/3
###                Distance to current cropland - same as ESA land cover


library(terra)
library(dplyr)
library(sf)
sf_use_s2(F)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(viridis)
library(tidyterra)
library(rnaturalearth)
library(terra)



#### Data Prep ####
## load in agricultural suitability layers from Zabel et al. ##

ag.suitability.current<- rast("Data/overall_suitability_subset_1to17.present.bil")
#export as a tif file
writeRaster(ag.suitability.current, "Data/overall_suitability_subset_1to17.present.tif")
ag.suitability.2010.2039.rcp2p6<- rast("Data/overall_suitability_subset_1to17.2010.2039.rcp2p6.bil")
ag.suitability.2010.2039.rcp8p5<- rast("Data/overall_suitability_subset_1to17.2010.2039.rcp8p5.bil")
ag.suitability.2040.2069.rcp2p6<- rast("Data/overall_suitability_subset_1to17.2040.2069.rcp2p6.bil")
ag.suitability.2040.2069.rcp8p5<- rast("Data/overall_suitability_subset_1to17.2040.2069.rcp8p5.bil")
ag.suitability.2070.2099.rcp2p6<- rast("Data/overall_suitability_subset_1to17.2070.2099.rcp2p6.bil")
ag.suitability.2070.2099.rcp8p5<- rast("Data/overall_suitability_subset_1to17.2070.2099.rcp8p5.bil")
plot(ag.suitability.2070.2099.rcp8p5)
expanse(ag.suitability.2070.2099.rcp8p5, unit = "ha")

## load in curtis forestry layer ##

curtis.forest.loss<- rast("Data/TCL_DD_2022_20230407.tif")
expanse(curtis.forest.loss, unit = "ha", byValue = T)
#retain only forestry forest loss, code 3 #
curtis.forestry<- subst(curtis.forest.loss, 3,1,others = NA)
curtis.forestry.new<- project(curtis.forestry, "epsg:4326")
expanse(curtis.forestry.new, unit = "ha")
plot(curtis.forestry)
# reproject to match crs and resolution of agriculture layer #
curtis.forestry.reprojected<- project(curtis.forestry, ag.suitability.current, filename = "Data/curtis.forestry.reprojected.tif", overwrite = T)
curtis.forestry.reprojected<- rast("Data/curtis.forestry.reprojected.tif")
plot(curtis.forestry.reprojected)
expanse(curtis.forestry.reprojected, unit = "ha")

## load in Hansen tree cover data - 2000 % tree cover aggregated to 1km cells ##

hansen.tree.cover<- rast("Data/hansen_tree_cover_1km.tif")
plot(hansen.tree.cover)

## filter Curtis layer for 1km cells with tree cover > 10% ##

hansen.reprojected.tree.cover<- project(hansen.tree.cover, curtis.forestry.reprojected, filename = "Data/hansen.reprojected.tree.cover.1km.tif")
plot(hansen.reprojected.tree.cover)

curtis.forestry.ten.percent<- mask(hansen.reprojected.tree.cover, curtis.forestry.reprojected)
plot(curtis.forestry.ten.percent)

# keep only cells with min tree cover 10%
curtis.forestry.ten.percent<- clamp(curtis.forestry.ten.percent, lower = 10, values = F, filename = "curtis.forestry.ten.percent.tif")
plot(curtis.forestry.ten.percent)


## Load in ESA land cover - retain only cells not mapped as current cropland, urban settlements, bare rock or water ##
esa.land.cover<- rast("Data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc")
esa.land.cover<- subset(esa.land.cover, 1, filename = "Data/esa.land.cover.original.resolution.tif")
esa.land.cover<- rast( "Data/esa.land.cover.original.resolution.tif")
plot(esa.land.cover)
# aggregate data to closest possible value to forestry layer (900m ) before reprojecting to this resolution (1km)
esa.land.cover.aggregated<- aggregate(esa.land.cover, fact = 3, fun = "modal", filename = "Data/esa.land.cover.900m.resolution.tif")
esa.land.cover.reprojected<- resample(esa.land.cover.aggregated, no.forestry.ag.suitability.current, method = "near",filename = "Data/esa.land.cover.original.aggregated.and.resampled.1km.tif")

# get agriculture only layer - used for later analysis to calculate distance to current cropland #
existing.agriculture<- subst(esa.land.cover.reprojected, 10:40, 1, others = NA, filename = "Data/esa.existing.agriculture.1km.tif")
existing.agriculture<- rast("Data/esa.existing.agriculture.1km.tif")
plot(existing.agriculture, col = "darkgreen")
#layer with agriculture, settlements, bare land and water removed
no.go.areas.removed<- subst(esa.land.cover.reprojected, c(10:40, 190:220), NA, filename = "Data/esa.no.go.areas.1km.tif", overwrite = T)
plot(no.go.areas.removed)
freq(esa.land.cover.reprojected)
freq(no.go.areas.removed)

### Final Forestry Layer ###

# remove agriculture, urban, bare and water from curtis forestry layer  #
curtis.forestry.ten.percent<- rast("curtis.forestry.ten.percent.tif")
curtis.forestry.final.area<- mask(curtis.forestry.ten.percent, no.go.areas.removed, filename = "Data/curtis.forestry.final.area.tif")

### Final Non-Forestry Layer ###
no.forestry.area<- mask(ag.suitability.current,curtis.forestry.final.area, inverse = T,  filename = "Data/no.curtis.area.tif", overwrite = T )
no.forestry.available.area<- mask(no.forestry.area, no.go.areas.removed, filename = "Data/no.forestry.available.area.tif")



#### Spatial Overlays ####

## overlay forestry land rasters with agricultural change rasters from Zabel et al. to get predicitons of future change in agricultural suitability in current forestry land ##

curtis.forestry.final.area<- rast("Data/curtis.forestry.final.area.tif")

## Current conditions ##

forestry.ag.suitability.current<- mask(ag.suitability.current,curtis.forestry.final.area, filename = "Data/suitability.rasters/curtis.forestry.current.ag.suitability.tif", overwrite = T )
plot(forestry.ag.suitability.current)

forestry.ag.suitability.current<- rast("Data/suitability.rasters/curtis.forestry.current.ag.suitability.tif")


## 2010-2039 ##

## RCP 2.6 ##
forestry.ag.suitability.2010.2039.rcp2p6<- mask(ag.suitability.2010.2039.rcp2p6,curtis.forestry.final.area, filename = "Data/suitability.rasters/curtis.forestry.2010.2039.rcp2p6.ag.suitability.tif", overwrite = T )
plot(forestry.ag.suitability.2010.2039.rcp2p6)

## RCP 8.5 ##
forestry.ag.suitability.2010.2039.rcp8p5<- mask(ag.suitability.2010.2039.rcp8p5,curtis.forestry.final.area, filename = "Data/suitability.rasters/curtis.forestry.2010.2039.rcp8p5.ag.suitability.tif", overwrite = T )
plot(forestry.ag.suitability.2010.2039.rcp8p5)

## 2040-2069 ##

## RCP 2.6 ##
forestry.ag.suitability.2040.2069.rcp2p6<- mask(ag.suitability.2040.2069.rcp2p6,curtis.forestry.final.area, filename = "Data/suitability.rasters/curtis.forestry.2040.2069.rcp2p6.ag.suitability.tif", overwrite = T )
plot(forestry.ag.suitability.2040.2069.rcp2p6)

## RCP 8.5 ##
forestry.ag.suitability.2040.2069.rcp8p5<- mask(ag.suitability.2040.2069.rcp8p5,curtis.forestry.final.area, filename = "Data/suitability.rasters/curtis.forestry.2040.2069.rcp8p5.ag.suitability.tif", overwrite = T )
plot(forestry.ag.suitability.2040.2069.rcp8p5)

## 2070-2099 ##

## RCP 2.6 ##
forestry.ag.suitability.2070.2099.rcp2p6<- mask(ag.suitability.2070.2099.rcp2p6,curtis.forestry.final.area, filename = "Data/suitability.rasters/curtis.forestry.2070.2099.rcp2p6.ag.suitability.tif", overwrite = T )
plot(forestry.ag.suitability.2070.2099.rcp2p6)

## RCP 8.5 ##
forestry.ag.suitability.2070.2099.rcp8p5<- mask(ag.suitability.2070.2099.rcp8p5,curtis.forestry.final.area, filename = "Data/suitability.rasters/curtis.forestry.2070.2099.rcp8p5.ag.suitability.tif", overwrite = T )
plot(forestry.ag.suitability.2070.2099.rcp8p5)


### Categorise into Marginal, Medium and High suitability - as per Zabel et al ###

m <- c(0, 0, 1,
       1, 32, 2,
       33, 74, 3,
       75, 100, 4)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
### Forestry Areas ###

# current agricultural suitability #
forestry.ag.suitability.current.classified <- classify(forestry.ag.suitability.current, rclmat, right = NA, filename = "Data/classified.rasters/curtis.forestry.current.ag.suitability.classified.tif", overwrite = T)
forestry.ag.suitability.current.classified<- rast("Data/classified.rasters/curtis.forestry.current.ag.suitability.classified.tif")

plot(forestry.ag.suitability.current.classified)

## RCP 2.6 ##

## 2010-2039 ##
forestry.ag.suitability.2010.2039.rcp2p6.classified <- classify(forestry.ag.suitability.2010.2039.rcp2p6, rclmat, right = NA, filename = "Data/classified.rasters/curtis.forestry.2010.2039.rcp2p6.ag.suitability.classified.tif", overwrite = T)
plot(forestry.ag.suitability.2010.2039.rcp2p6.classified)
## 2040-2069 ##
forestry.ag.suitability.2040.2069.rcp2p6.classified <- classify(forestry.ag.suitability.2040.2069.rcp2p6, rclmat, right = NA, filename = "Data/classified.rasters/curtis.forestry.2040.2069.rcp2p6.ag.suitability.classified.tif", overwrite = T)
plot(forestry.ag.suitability.2040.2069.rcp2p6.classified)
## 2070-2099 ##
forestry.ag.suitability.2070.2099.rcp2p6.classified <- classify(forestry.ag.suitability.2070.2099.rcp2p6, rclmat, right = NA, filename = "Data/classified.rasters/curtis.forestry.2070.2099.rcp2p6.ag.suitability.classified.tif", overwrite = T)
plot(forestry.ag.suitability.2070.2099.rcp2p6.classified)

## RCP 8.5 ##

## 2010-2039 ##
forestry.ag.suitability.2010.2039.rcp8p5.classified <- classify(forestry.ag.suitability.2010.2039.rcp8p5, rclmat, right = NA, filename = "Data/classified.rasters/curtis.forestry.2010.2039.rcp8p5.ag.suitability.classified.tif", overwrite = T)
plot(forestry.ag.suitability.2010.2039.rcp8p5.classified)
## 2040-2069 ##
forestry.ag.suitability.2040.2069.rcp8p5.classified <- classify(forestry.ag.suitability.2040.2069.rcp8p5, rclmat, right = NA, filename = "Data/classified.rasters/curtis.forestry.2040.2069.rcp8p5.ag.suitability.classified.tif", overwrite = T)
plot(forestry.ag.suitability.2040.2069.rcp8p5.classified)
## 2070-2099 ##
forestry.ag.suitability.2070.2099.rcp8p5.classified <- classify(forestry.ag.suitability.2070.2099.rcp8p5, rclmat, right = NA, filename = "Data/classified.rasters/curtis.forestry.2070.2099.rcp8p5.ag.suitability.classified.tif", overwrite = T)
forestry.ag.suitability.2070.2099.rcp8p5.classified<- rast("Data/classified.rasters/curtis.forestry.2070.2099.rcp8p5.ag.suitability.classified.tif")
plot(forestry.ag.suitability.2070.2099.rcp8p5.classified)
freq(forestry.ag.suitability.2070.2099.rcp8p5.classified)

### Get Change in suitability from present to future ###

### Forestry Areas ###

#current
forestry.ag.suitability.current<- rast("Data/suitability.rasters/curtis.forestry.current.ag.suitability.tif")
plot(forestry.ag.suitability.current)
## RCP 2.6 ##

## 2010 - 2039 ##
forestry.ag.suitability.2010.2039.rcp2p6<- rast("Data/suitability.rasters/curtis.forestry.2010.2039.rcp2p6.ag.suitability.tif")
plot(forestry.ag.suitability.2010.2039.rcp2p6)
# get difference #
forestry.ag.suitability.2010.2039.rcp2p6.change<- forestry.ag.suitability.2010.2039.rcp2p6 - forestry.ag.suitability.current
plot(forestry.ag.suitability.2010.2039.rcp2p6.change)
writeRaster(forestry.ag.suitability.2010.2039.rcp2p6.change, "Data/change.rasters/curtis.forestry.2010.2039.rcp2p6.ag.suitability.change.tif", overwrite = T)

## 2040 - 2069 ##
forestry.ag.suitability.2040.2069.rcp2p6<- rast("Data/suitability.rasters/curtis.forestry.2040.2069.rcp2p6.ag.suitability.tif")
plot(forestry.ag.suitability.2040.2069.rcp2p6)
# get difference #
forestry.ag.suitability.2040.2069.rcp2p6.change<- forestry.ag.suitability.2040.2069.rcp2p6 - forestry.ag.suitability.current
plot(forestry.ag.suitability.2040.2069.rcp2p6.change)
writeRaster(forestry.ag.suitability.2040.2069.rcp2p6.change, "Data/change.rasters/curtis.forestry.2040.2069.rcp2p6.ag.suitability.change.tif", overwrite = T)

## 2070 - 2099 ##
forestry.ag.suitability.2070.2099.rcp2p6<- rast("Data/suitability.rasters/curtis.forestry.2070.2099.rcp2p6.ag.suitability.tif")
plot(forestry.ag.suitability.2070.2099.rcp2p6)
# get difference #
forestry.ag.suitability.2070.2099.rcp2p6.change<- forestry.ag.suitability.2070.2099.rcp2p6 - forestry.ag.suitability.current
plot(forestry.ag.suitability.2070.2099.rcp2p6.change)
writeRaster(forestry.ag.suitability.2070.2099.rcp2p6.change, "Data/change.rasters/curtis.forestry.2070.2099.rcp2p6.ag.suitability.change.tif", overwrite = T)

## RCP 8.5 ##

## 2010 - 2039 ##
forestry.ag.suitability.2010.2039.rcp8p5<- rast("Data/suitability.rasters/curtis.forestry.2010.2039.rcp8p5.ag.suitability.tif")
plot(forestry.ag.suitability.2010.2039.rcp8p5)
# get difference #
forestry.ag.suitability.2010.2039.rcp8p5.change<- forestry.ag.suitability.2010.2039.rcp8p5 - forestry.ag.suitability.current
plot(forestry.ag.suitability.2010.2039.rcp8p5.change)
writeRaster(forestry.ag.suitability.2010.2039.rcp8p5.change, "Data/change.rasters/curtis.forestry.2010.2039.rcp8p5.ag.suitability.change.tif", overwrite = T)

## 2040 - 2069 ##
forestry.ag.suitability.2040.2069.rcp8p5<- rast("Data/suitability.rasters/curtis.forestry.2040.2069.rcp8p5.ag.suitability.tif")
plot(forestry.ag.suitability.2040.2069.rcp8p5)
# get difference #
forestry.ag.suitability.2040.2069.rcp8p5.change<- forestry.ag.suitability.2040.2069.rcp8p5 - forestry.ag.suitability.current
plot(forestry.ag.suitability.2040.2069.rcp8p5.change)
writeRaster(forestry.ag.suitability.2040.2069.rcp8p5.change, "Data/change.rasters/curtis.forestry.2040.2069.rcp8p5.ag.suitability.change.tif", overwrite = T)

## 2070 - 2099 ##
forestry.ag.suitability.2070.2099.rcp8p5<- rast("Data/suitability.rasters/curtis.forestry.2070.2099.rcp8p5.ag.suitability.tif")
plot(forestry.ag.suitability.2070.2099.rcp8p5)
# get difference #
forestry.ag.suitability.2070.2099.rcp8p5.change<- forestry.ag.suitability.2070.2099.rcp8p5 - forestry.ag.suitability.current
plot(forestry.ag.suitability.2070.2099.rcp8p5.change)
writeRaster(forestry.ag.suitability.2070.2099.rcp8p5.change, "Data/change.rasters/curtis.forestry.2070.2099.rcp8p5.ag.suitability.change.tif", overwrite = T)


#### Convenience functions ####
##Fast mean 
#https://stackoverflow.com/questions/10397574/efficiently-compute-mean-and-standard-deviation-from-a-frequency-table

fastmean <- function(dat) {
  with(dat, sum(area*value)/sum(area)) 
}

##Fast SD 
#https://stackoverflow.com/questions/10397574/efficiently-compute-mean-and-standard-deviation-from-a-frequency-table
fastSD <- function(dat) {
  mu <- fastmean(dat)
  with(dat, sqrt(sum(area*(value-mu)^2)/(sum(area)-1) ) )
}


#### Global suitability area ####
## get list of layers
layers <- list.files("Data/classified.rasters/", full.names = TRUE)[1:14]
## storage
lyr.dat <- data.frame()

## loop through the layers pulling out variables and calculating the 
## expanse of each suitability.
for (i in 1:length(layers)) {
  lyr <- layers[i]
  r.lyr <- rast(lyr)
  cat(lyr, '\n', i, "out of", length(layers), '\n')
  
  if (grepl(pattern = "current", lyr) == TRUE){time <- "Current"}
  if (grepl(pattern = "2010", lyr) == TRUE){time <- "2010-2039"}
  if (grepl(pattern = "2040", lyr) == TRUE){time <- "2040-2069"}
  if (grepl(pattern = "2070", lyr) == TRUE){time <- "2070-2099"}
  if (grepl(pattern = "rcp8p5", lyr) == TRUE){rcp <- "rcp8.5"}
  if (grepl(pattern = "rcp2p6", lyr) == TRUE){rcp <- "rcp2.6"}
  if (grepl(pattern = "current", lyr) == TRUE){rcp <- NA}
  if (grepl(pattern = "no.forestry", lyr) == TRUE){land.cover <- "Non-forestry"} else {
    land.cover <- "Forestry"}
  
  tot <- r.lyr %>% expanse(unit = "ha")
  Unsuitable <- r.lyr %>% ifel(. != 1, NA, .) %>% expanse(unit = "ha")
  Marginal <- r.lyr %>% ifel(. != 2, NA, .) %>% expanse(unit = "ha")
  Moderately <- r.lyr %>% ifel(. != 3, NA, .) %>% expanse(unit = "ha")
  Highly <- r.lyr %>% ifel(. != 4, NA, .) %>% expanse(unit = "ha")
  
  lyr.add <- data.frame(suitability = c("Unsuitable", "Marginal", "Moderately", "Highly"),
                        area.ha = c(Unsuitable$area, Marginal$area, Moderately$area, Highly$area),
                        total.ha = tot$area, RCP = rcp, time = time, land.cover)
  
  lyr.dat <- rbind(lyr.dat, lyr.add)
  
}

#write.csv(lyr.dat, "Data/SuitableArea/SuitableArea.Summary.csv")
lyr.dat <- read.csv("Data/SuitableArea/SuitableArea.Summary.csv")


#### Global gain/loss/persist area ####
layers <- list.files("Data/CurtisLayers/good.land", full.names = TRUE)[1:12]
all.gain.loss.dat <- data.frame()
for (j in 1:length(layers)) {
  lyr <- layers[j]
  cat(lyr, '\n', j, "out of", length(layers), '\n')
  
  if (grepl(pattern = "2010", lyr) == TRUE){time <- "2010-2039"}
  if (grepl(pattern = "2040", lyr) == TRUE){time <- "2040-2069"}
  if (grepl(pattern = "2070", lyr) == TRUE){time <- "2070-2099"}
  if (grepl(pattern = "rcp8p5", lyr) == TRUE){rcp <- "rcp8.5"}
  if (grepl(pattern = "rcp2p6", lyr) == TRUE){rcp <- "rcp2.6"}
  if (grepl(pattern = "no.forestry", layers[j]) == TRUE){land.cover <- "no.forestry"} else {
    land.cover <- "forestry"}
  
  
  cat("Expansing", '\n')
  mask.lyr <- rast(lyr)
  gain.loss.ex <- expanse(mask.lyr, byValue = TRUE, unit = "ha")
  
  gain.loss.i.add <- gain.loss.ex %>% 
    mutate(rcp = rcp, time = time, land.cover = land.cover,
           suitability = case_when(value == 1 ~ "Loss", 
                                   value == 2 ~ "Persist",
                                   value == 3 ~ "Gain",
                                   value == "NO.FOREST" ~ "NO.FOREST"))
  
  all.gain.loss.dat <- rbind(all.gain.loss.dat, gain.loss.i.add)       
}
#write.csv(all.gain.loss.dat, "Data/CountryArea/global.gain.loss.raw.csv")

#### Travel time to population centers summary ####
tt.layers <- list.files("Data/TravelTime", full.names = TRUE)[c(3,4,7,8)]
tt.dat <- data.frame()
for (i in 1:length(tt.layers)) {
  lyr <- tt.layers[i]
  cat(lyr, '\n', i, "out of", length(tt.layers), '\n')
  if (grepl(pattern = "2010", lyr) == TRUE){time <- "2010-2039"}
  if (grepl(pattern = "2040", lyr) == TRUE){time <- "2040-2069"}
  if (grepl(pattern = "2070", lyr) == TRUE){time <- "2070-2099"}
  if (grepl(pattern = "rcp8p5", lyr) == TRUE){rcp <- "rcp8.5"}
  if (grepl(pattern = "rcp2p6", lyr) == TRUE){rcp <- "rcp2.6"}
  if (grepl(pattern = "non.forestry", lyr) == TRUE){land.cover <- "non.forestry"} else {
    land.cover <- "forestry"}
  tt.ex <- rast(lyr) %>% expanse(byValue = TRUE, unit = "ha") %>%
    mutate(rcp = rcp, time = time, land.cover = land.cover)
  tt.mean <- fastmean(tt.ex)
  tt.sd <- fastSD(tt.ex)
  
  ## quantile to find
  q5 <- 0.05
  q25 <- 0.25
  q50 <- 0.5
  q75 <- 0.75
  q95 <- 0.95
  tt.ex$cumfreq <- cumsum(tt.ex$area)/sum(tt.ex$area)
  q5.lyr <- tt.ex$value[tt.ex$cumfreq >= q5][1]
  q50.lyr <- tt.ex$value[tt.ex$cumfreq >= q50][1]
  q95.lyr <- tt.ex$value[tt.ex$cumfreq >= q95][1]
  q25.lyr <- tt.ex$value[tt.ex$cumfreq >= q25][1]
  q75.lyr <- tt.ex$value[tt.ex$cumfreq >= q75][1]
  
  tt.add <- data.frame(rcp = rcp, time = time, land.cover = land.cover,
                       mean = tt.mean, sd = tt.sd,
                       q5 = q5.lyr, q50 = q50.lyr, q95 = q95.lyr,
                       q25 = q25.lyr, q75 = q75.lyr)
  
  tt.dat <- rbind(tt.dat, tt.add)
  
  write.csv(tt.ex, paste0("Data/TravelTime/", land.cover,
                          ".", time, ".", rcp, "raw.csv"))
}

#write.csv(tt.dat, "Data/TravelTime/travel.time.sum.csv")

#### Distance to agriculture summary ####
d.layers <- list.files("Data/Distance", full.names = TRUE)[c(4,5,8, 9)]
d.dat <- data.frame()
for (i in 1:length(d.layers)) {
  lyr <- d.layers[i]
  cat(lyr, '\n', i, "out of", length(d.layers), '\n')
  if (grepl(pattern = "2010", lyr) == TRUE){time <- "2010-2039"}
  if (grepl(pattern = "2040", lyr) == TRUE){time <- "2040-2069"}
  if (grepl(pattern = "2070", lyr) == TRUE){time <- "2070-2099"}
  if (grepl(pattern = "rcp8p5", lyr) == TRUE){rcp <- "rcp8.5"}
  if (grepl(pattern = "rcp2p6", lyr) == TRUE){rcp <- "rcp2.6"}
  if (grepl(pattern = "non.forestry", lyr) == TRUE){land.cover <- "non.forestry"} else {
    land.cover <- "forestry"}
  d.ex <- rast(lyr) %>% expanse(byValue = TRUE, unit = "ha") %>%
    mutate(rcp = rcp, time = time, land.cover = land.cover)
  d.mean <- fastmean(d.ex)
  d.sd <- fastSD(d.ex)
  
  ## quantile to find
  q5 <- 0.05
  q50 <- 0.5
  q95 <- 0.95
  q75 <- 0.75
  q95 <- 0.95
  d.ex$cumfreq <- cumsum(d.ex$area)/sum(d.ex$area)
  q5.lyr <- d.ex$value[d.ex$cumfreq >= q5][1]
  q50.lyr <- d.ex$value[d.ex$cumfreq >= q50][1]
  q95.lyr <- d.ex$value[d.ex$cumfreq >= q95][1]
  q25.lyr <- d.ex$value[d.ex$cumfreq >= q25][1]
  q75.lyr <- d.ex$value[d.ex$cumfreq >= q75][1]
  
  d.add <- data.frame(rcp = rcp, time = time, land.cover = land.cover,
                      mean = d.mean, sd = d.sd, 
                      q5 = q5.lyr, q50 = q50.lyr, q95 = q95.lyr,
                      q25 = q25.lyr, q75 = q75.lyr)
  
  d.dat <- rbind(d.dat, d.add)
  
  write.csv(d.ex, paste0("Data/Distance/", land.cover,
                         ".", time, ".", rcp, "raw.csv"))
}
write.csv(d.dat, "Data/Distance/distance.sum.csv")

#### COUNTRY LEVEL ANALYSES ####
world <- sf::st_read("Data/GADMworld/gadm_410-levels.gpkg", layer = "ADM_0")

#### Country classified area ####
## Get country outlines
#world <- ne_download(scale = "large", returnclass = 'sf')
layers <- list.files("Data/classified.rasters/", full.names = TRUE)[1:14]
timber.countries <- c("USA", "Russia", "China", "Brazil", "Canada")
timber.codes <- c("USA", "RUS", "CHN", "BRA", "CAN")

c.lyr.dat <- data.frame()

for (i in 1:length(timber.countries)) {
  
  country.i <- timber.countries[i]
  code.i <- timber.codes[i]
  
  country.border <- vect(paste0("Data/GADM/GADM_",country.i,"/gadm41_", code.i, "_0.shp" ))
  cat("Working on ", country.i, ": ", i, "out of", length(timber.countries), '\n')
  
  for (j in 1:length(layers)) {
    lyr <- layers[j]
    cat(lyr, '\n', j, "out of", length(layers), '\n')
    
    if (grepl(pattern = "current", lyr) == TRUE){time <- "Current"}
    if (grepl(pattern = "2010", lyr) == TRUE){time <- "2010-2039"}
    if (grepl(pattern = "2040", lyr) == TRUE){time <- "2040-2069"}
    if (grepl(pattern = "2070", lyr) == TRUE){time <- "2070-2099"}
    if (grepl(pattern = "rcp8p5", lyr) == TRUE){rcp <- "rcp8.5"}
    if (grepl(pattern = "rcp2p6", lyr) == TRUE){rcp <- "rcp2.6"}
    if (grepl(pattern = "current", lyr) == TRUE){rcp <- NA}
    if (grepl(pattern = "no.forestry", lyr) == TRUE){land.cover <- "Non-forestry"} else {
      land.cover <- "Forestry"}
    
    c.extent <- terra::ext(country.border)
    mask.c.area <- rast(lyr) %>% crop(., c.extent) %>%mask(., country.border)  
    
    tot <- mask.c.area %>% expanse(unit = "ha")
    mask.ex <- mask.c.area %>% expanse(byValue = TRUE, unit = "ha")
    
    c.lyr.add <- mask.ex %>% mutate(total.ha = tot$area, RCP = rcp, time = time, land.cover, country = country.i)
    
    c.lyr.dat <- rbind(c.lyr.dat, c.lyr.add)                    
  }
  write.csv(c.lyr.dat, paste0("Data/CountryArea/temp.out.", i, ".csv"))
  
}
c.lyr.dat <- c.lyr.dat %>% mutate(suitability = case_when(value == 1 ~ "Unsuitable",
                                                          value == 2 ~ "Marginal",
                                                          value == 3 ~ "Moderately",
                                                          value == 4 ~ "Highly")) %>%
  rename("area.ha" = "area") %>% 
  select(-c(layer, value))

#write.csv(c.lyr.dat, "Data/CountryArea/Timber.top5.suitability.csv")
c.lyr.dat <- read.csv("Data/CountryArea/Timber.top5.suitability.csv")




#### All countries suitability ####
layers <- list.files("Data/classified.rasters/", full.names = TRUE)[1:14]

world <- sf::st_read("Data/GADMworld/gadm_410-levels.gpkg", layer = "ADM_0")

all.c.suit.dat <- data.frame()

i <- 6
for (i in 1:nrow(world)) {
  
  country.i.lyr <- world[i,]
  country.i.id <- country.i.lyr$COUNTRY
  country.i.gid <- country.i.lyr$GID_0
  
  #country.i.cont <- country.i.lyr$continent
  #country.i.reg <- country.i.lyr$region.wb
  
  cat("Working on ", country.i.id, ": ", i, "out of", nrow(world), '\n')
  
  c.extent <- terra::ext(vect(country.i.lyr$geom))
  country.border <- vect(country.i.lyr$geom)
  
  for (j in 1:length(layers)) {
    lyr <- layers[j]
    cat(lyr, '\n', j, "out of", length(layers), '\n')
    
    if (grepl(pattern = "current", lyr) == TRUE){time <- "Current"}
    if (grepl(pattern = "2010", lyr) == TRUE){time <- "2010-2039"}
    if (grepl(pattern = "2040", lyr) == TRUE){time <- "2040-2069"}
    if (grepl(pattern = "2070", lyr) == TRUE){time <- "2070-2099"}
    if (grepl(pattern = "rcp8p5", lyr) == TRUE){rcp <- "rcp8.5"}
    if (grepl(pattern = "rcp2p6", lyr) == TRUE){rcp <- "rcp2.6"}
    if (grepl(pattern = "current", lyr) == TRUE){rcp <- NA}
    
    
    cat("Masking", '\n')
    
    mask.c.area.curr <- rast(lyr) %>% 
      crop(., c.extent) %>%mask(., country.border)  
    
    cat("Expansing", '\n')
    
    c.suit.ex <- expanse(mask.c.area.curr, byValue = TRUE, unit = "ha")
    
    if(nrow(c.suit.ex) == 0) {c.suit.ex <- data.frame(layer = "NO.FOREST", value = "NO.FOREST", area = "NO.FOREST")}
    
    c.suit.i.add <- c.suit.ex %>% 
      mutate(country = country.i.id, country.gid = country.i.gid, rcp = rcp, time = time,
             suitability = case_when(value == 0 ~ "Unsuitable", value > 0 & value < 33 ~ "Marginal",
                                     value > 32 & value < 75 ~ "Moderately", value > 74 ~ "Highly",
                                     value == "NO.FOREST" ~ "NO.FOREST"))
    
    all.c.suit.dat <- rbind(all.c.suit.dat, c.suit.i.add)       
  }
}



write.csv(all.c.suit.dat, "Data/CountryArea/all.country.suitability.area.raw.csv")
#### All countries productive gain/loss/persist ####


all.crops.suitability.current.forest <-rast("Data/CurtisLayers/curtis.forestry.current.ag.suitability.classified 2.tif")
all.crops.suitability.current.noforest <-rast("Data/CurtisLayers/curtis.no.forestry.current.ag.suitability.classified 2.tif")

all.crops.good.current.forest<- clamp(all.crops.suitability.current.forest, lower = 3, values = F, 
                                      filename = "Data/CurtisLayers/good.land/steps/all.crops.good.land.forest.current.tif", overwrite = T)
all.crops.good.current.noforest<- clamp(all.crops.suitability.current.noforest, lower = 3, values = F, 
                                        filename = "Data/CurtisLayers/good.land/steps/all.crops.good.land.no.forestry.current.tif", overwrite = T)

layers <- list.files("Data/CurtisLayers/", full.names = TRUE)[c(1:6,8:13)]
for (i in 1:length(layers)) {
  lyr <- rast(layers[i])
  cat(layers[i], '\n', i, "out of", length(layers), '\n')
  
  if (grepl(pattern = "2010", layers[i]) == TRUE){time <- "2010.2039"}
  if (grepl(pattern = "2040", layers[i]) == TRUE){time <- "2040.2069"}
  if (grepl(pattern = "2070", layers[i]) == TRUE){time <- "2070.2099"}
  if (grepl(pattern = "rcp8p5", layers[i]) == TRUE){rcp <- "rcp8p5"}
  if (grepl(pattern = "rcp2p6", layers[i]) == TRUE){rcp <- "rcp2p6"}
  if (grepl(pattern = "no.forestry", layers[i]) == TRUE){land.cover <- "no.forestry"} else {
    land.cover <- "forestry"}
  
  if (land.cover == "forestry") {current <- all.crops.good.current.forest} else {
    current <- all.crops.good.current.noforest}
  
  lyr.good <- clamp(lyr, lower = 3, values = F,
                    filename = paste0("Data/CurtisLayers/good.land/steps/all.crops.good.land.",time, ".", rcp, ".", land.cover, ".tif"), overwrite = T)
  # get areas good.land in the future but not in the present
  lyr.newly.good <- mask(lyr.good, current, inverse = T,
                         filename = paste0("Data/CurtisLayers/good.land/steps/all.crops.newly.good.land.",time, ".", rcp, ".", land.cover, ".tif"), overwrite = T)
  # get areas good.land in the present but not the future
  lyr.not.good <-  mask(current,lyr.good, inverse = T, 
                        filename = paste0("Data/CurtisLayers/good.land/steps/all.crops.not.good.land.",time, ".", rcp, ".", land.cover, ".tif"), overwrite = T)
  # get areas that remain good.land in the future
  lyr.remain.good <- mask(lyr.good,current, 
                          filename = paste0("Data/CurtisLayers/good.land/steps/all.crops.remain.good.land.",time, ".", rcp, ".", land.cover, ".tif"), overwrite = T)
  
  ## change to values 1 = loss, 2 = persist, 3 = gain
  lyr.not.good.2 <- subst(lyr.not.good, c(1:4), 1)
  lyr.remain.good.2 <- subst(lyr.remain.good, c(1:4), 2)
  lyr.newly.good.2 <- subst(lyr.newly.good, c(1:4), 3)
  
  ## put together ##
  all.crops.land.changes <- max(lyr.not.good.2,lyr.remain.good.2,lyr.newly.good.2, na.rm = T) 
  
  writeRaster(all.crops.land.changes, 
              filename = paste0("Data/CurtisLayers/good.land/gain.loss.", time, ".", rcp, ".", land.cover, ".tif"), overwrite = T) 
}


world <- sf::st_read("Data/GADMworld/gadm_410-levels.gpkg", layer = "ADM_0")

all.c.gain.loss.dat <- data.frame()
layers <- list.files("Data/CurtisLayers/good.land", full.names = TRUE)[1:12]

for (i in 1:nrow(world)) {
  
  country.i.lyr <- world[i,]
  country.i.id <- country.i.lyr$COUNTRY
  country.i.gid <- country.i.lyr$GID_0
  
  cat("Working on ", country.i.id, ": ", i, "out of", nrow(world), '\n')
  
  c.extent <- terra::ext(vect(country.i.lyr$geom))
  country.border <- vect(country.i.lyr$geom)
  
  for (j in 1:length(layers)) {
    lyr <- layers[j]
    cat(lyr, '\n', j, "out of", length(layers), '\n')
    
    if (grepl(pattern = "2010", lyr) == TRUE){time <- "2010-2039"}
    if (grepl(pattern = "2040", lyr) == TRUE){time <- "2040-2069"}
    if (grepl(pattern = "2070", lyr) == TRUE){time <- "2070-2099"}
    if (grepl(pattern = "rcp8p5", lyr) == TRUE){rcp <- "rcp8.5"}
    if (grepl(pattern = "rcp2p6", lyr) == TRUE){rcp <- "rcp2.6"}
    if (grepl(pattern = "no.forestry", layers[j]) == TRUE){land.cover <- "no.forestry"} else {
      land.cover <- "forestry"}
    
    cat("Masking", '\n')
    
    mask.c.lyr<- rast(lyr) %>% 
      crop(., c.extent) %>% mask(., country.border) 
    
    cat("Expansing", '\n')
    
    c.gain.loss.ex <- expanse(mask.c.lyr, byValue = TRUE, unit = "ha")
    
    if(nrow(c.gain.loss.ex) == 0) {c.suit.ex <- data.frame(layer = "NO.FOREST", value = "NO.FOREST", area = "NO.FOREST")}
    
    c.gain.loss.i.add <- c.gain.loss.ex %>% 
      mutate(country = country.i.id, country.gid = country.i.gid, 
             rcp = rcp, time = time, land.cover = land.cover,
             suitability = case_when(value == 1 ~ "Loss", 
                                     value == 2 ~ "Persist",
                                     value == 3 ~ "Gain",
                                     value == "NO.FOREST" ~ "NO.FOREST"))
    
    all.c.gain.loss.dat <- rbind(all.c.gain.loss.dat, c.gain.loss.i.add)       
  }
}

write.csv(all.c.gain.loss.dat, "Data/CountryArea/all.country.gain.loss.raw.csv")


#### All countries change ####

world <- sf::st_read("Data/GADMworld/gadm_410-levels.gpkg", layer = "ADM_0")
world.l <- sf::st_layers("Data/GADMworld/gadm_410-levels.gpkg")

world <- rnaturalearth::ne_countries(returnclass = "sf")
c.world <- world[15,]

c.world$name
terra::ext(vect(c.world$geometry))


c.lyr.dat[1,]
plot(world)

all.c.dat <- data.frame()

i <- 6
for (i in 1:nrow(world)) {
  
  country.i.lyr <- world[i,]
  country.i.id <- country.i.lyr$COUNTRY
  country.i.gid <- country.i.lyr$GID_0
  
  #country.i.cont <- country.i.lyr$continent
  #country.i.reg <- country.i.lyr$region.wb
  
  cat("Working on ", country.i.id, ": ", i, "out of", nrow(world), '\n')
  
  c.extent <- terra::ext(vect(country.i.lyr$geom))
  country.border <- vect(country.i.lyr$geom)
  
  cat("Masking", '\n')
  
  mask.c.area.curr <- rast("Data/CurtisLayers/Raw_suitability_scores/curtis.forestry.current.ag.suitability.tif") %>% 
    crop(., c.extent) %>%mask(., country.border)  
  
  mask.c.area.26 <- rast("Data/CurtisLayers/Raw_suitability_scores/curtis.forestry.2070.2099.rcp2p6.ag.suitability.tif") %>% 
    crop(., c.extent) %>%mask(., country.border) 
  
  mask.c.area.85 <- rast("Data/CurtisLayers/Raw_suitability_scores/curtis.forestry.2070.2099.rcp8p5.ag.suitability.tif") %>% 
    crop(., c.extent) %>%mask(., country.border) 
  
  cat("Differencing", '\n')
  
  diff.85 <- mask.c.area.85 - mask.c.area.curr
  diff.26 <- mask.c.area.26 - mask.c.area.curr
  
  cat("Expansing", '\n')
  
  diff.85.ex <- expanse(diff.85, byValue = TRUE, unit = "ha")
  diff.26.ex <- expanse(diff.26, byValue = TRUE, unit = "ha")
  
  if(nrow(diff.26.ex) == 0) {diff.26.ex <- data.frame(layer = "NO.FOREST", value = "NO.FOREST", area = "NO.FOREST", rcp = "rcp2p6")}
  if(nrow(diff.85.ex) == 0) {diff.85.ex <- data.frame(layer = "NO.FOREST", value = "NO.FOREST", area = "NO.FOREST", rcp = "rcp8p5")}
  if(nrow(diff.85.ex) != 0) {diff.85.ex$rcp <- "rcp8p5"}
  if(nrow(diff.26.ex) != 0) {diff.26.ex$rcp <- "rcp2p6"}
  
  c.i.add <- rbind(diff.26.ex, diff.85.ex) %>% 
    mutate(country = country.i.id, country.gid = country.i.gid,
           change = case_when(value > 0 ~ "Increase", value < 0 ~ "Decrease", value == 0 ~ "No change",
                              value == "NO.FOREST" ~ "NO.FOREST"))
  
  all.c.dat <- rbind(all.c.dat, c.i.add)       
}

all.c.dat.sum <- all.c.dat %>% filter(area != "NO.FOREST") %>% 
  group_by(rcp, country, country.gid, change) %>% summarise(area = sum(as.numeric(as.character(area))))

write.csv(all.c.dat, "Data/CountryArea/all.country.change.area.2070.2099.raw.csv")
write.csv(all.c.dat.sum, "Data/CountryArea/all.country.change.area.2070.2099.sum.csv")

#### Country frontier area #####
front.layers <- list.files("Data/Agri.Frontiers/", full.names = TRUE)
timber.countries <- c("USA", "Russia", "China", "Brazil", "Canada")
timber.codes <- c("USA", "RUS", "CHN", "BRA", "CAN")

front.dat <- data.frame()

for (i in 1:length(timber.countries)) {
  
  country.i <- timber.countries[i]
  code.i <- timber.codes[i]
  
  country.border <- vect(paste0("Data/GADM/GADM_",country.i,"/gadm41_", code.i, "_0.shp" ))
  cat("Working on ", country.i, ": ", i, "out of", length(timber.countries), '\n')
  
  for (j in 1:length(front.layers)) {
    lyr <- front.layers[j]
    cat(lyr, '\n', j, "out of", length(front.layers), '\n')
    
    if (grepl(pattern = "2040", lyr) == TRUE){time <- "2040-2069"}
    if (grepl(pattern = "2070", lyr) == TRUE){time <- "2070-2099"}
    if (grepl(pattern = "rcp8p5", lyr) == TRUE){rcp <- "rcp8.5"}
    if (grepl(pattern = "rcp2p6", lyr) == TRUE){rcp <- "rcp2.6"}
    
    c.extent <- terra::ext(country.border)
    
    cat("Masking", '\n')
    
    front.c.area <- rast(lyr) %>% 
      crop(., c.extent) %>%mask(., country.border)  
    
    cat("Expansing", '\n')
    
    front.c.area.ex <- expanse(front.c.area, byValue = TRUE, unit = "ha")
    
    front.c.area.ex$rcp <- rcp
    front.c.area.ex$time <- time
    front.c.area.ex$country <- country.i
    
    front.dat <- rbind(front.dat, front.c.area.ex)       
    
  }
}
write.csv(front.dat, "Data/CountryArea/Timber.top5.frontier.area.raw.csv")

front.dat.sum <- front.dat %>% group_by(rcp,country, time) %>% summarise(area = sum(area))
write.csv(front.dat.sum, "Data/CountryArea/Timber.top5.frontier.area.sum.csv")
