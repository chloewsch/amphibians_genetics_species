## Determinants of genetic diversity and species richness of North American amphibians
# Schmidt, Munshi-South, Dray, Garroway, J. Biogeography (2022) 

## 1. Build data
library(tidyverse)
library(sf)
library(raster)

# Data ----
amph <- read.csv("Schmidt_amphibian_data_JBiogeo.csv", head = TRUE) # (complete dataset)
sites <- st_as_sf(amph, coords = c("lon", "lat"), crs=4326)

## Note: site level species diversity done in GIS: 
# spatial join 1:1 target = points, join = polygons, 10km buffer around points

## Environmental variables
# PET & AET: https://figshare.com/articles/dataset/Global_Aridity_Index_and_Potential_Evapotranspiration_ET0_Climate_Database_v2/7504448/4
# PET
PET <- raster("global-et0_annual.tif/et0_yr/et0_yr.tif")

# AET
AET <- raster("AET_YR/aet_yr/w001001.adf")

# Land cover
# 30m: http://www.cec.org/north-american-environmental-atlas/land-cover-30m-2015-landsat-and-rapideye/
landcov <- raster("north_america_2015/NA_NALCMS_2015_LC_30m_LAEA_mmu5pix_/NA_NALCMS_2015_LC_30m_LAEA_mmu5pix_.tif")

sites_laea <- st_transform(sites, crs = crs(landcov))

# Extract values per site within 10km, 25km, 40km, 80km ----
buffers <- list(10*1000, 25*1000, 40*1000, 80*1000)

## 1) AET
aetmulti <- lapply(buffers, function(x) raster::extract(AET, sites, fun=mean, buffer=x, na.rm=TRUE, df=TRUE))
aetdf <- as.data.frame(aetmulti)
aetdf <- aetdf[,-c(1,3,5)] # remove 'ID' columns; keep data only
names(aetdf) <- c("AET_10", "AET_25", "AET_40", "AET_80")

## 2) PET
petmulti <- lapply(buffers, function(x) raster::extract(PET, sites, fun=mean, buffer=x, na.rm=TRUE, df=TRUE))
petdf <- as.data.frame(petmulti)
petdf <- petdf[,-c(1,3,5)] # remove 'ID' columns; keep data only
names(petdf) <- c("PET_10", "PET_25", "PET_40", "PET_80")


## 3) Heterogeneity
## SIMPSON'S INDEX
simpson_diversity <- function(cat.vect){
  px  = table(cat.vect)/length(cat.vect)
  D1 = 1-sum(px^2)
  return(D1)
}

# One by one because the objects are large
lc_10 <- raster::extract(landcov, sites_laea, buffer = 10*1000) # this outputs a list, every item is the summary for 1 site
lc_10SI <- lapply(lc_10, function(x) simpson_diversity(as.factor(x)))
lc_10SI <- unlist(lc_10mSI)
rm(lc_10)

lc_25 <- raster::extract(landcov, sites_laea, buffer = 25*1000) 
lc_25SI <- lapply(lc_25, function(x) simpson_diversity(as.factor(x)))
lc_25SI <- unlist(lc_25SI)
rm(lc_25)

lc_40 <- raster::extract(landcov, sites_laea, buffer = 40*1000) 
lc_40SI <- lapply(lc_40, function(x) simpson_diversity(as.factor(x)))
lc_40SI <- unlist(lc_40SI)
rm(lc_40)

lc_80 <- raster::extract(landcov, sites_laea, buffer = 80*1000)
lc_80SI <- lapply(lc_80, function(x) simpson_diversity(as.factor(x)))
lc_80SI <- unlist(lc_80SI)
rm(lc_80)

hetdf <- do.call(cbind, list(lc_10SI, lc_25SI, lc_40SI, lc_80SI))
names(hetdf) <- c("het_10", "het_25", "het_40", "het_80")

sgdata.amph <- do.call(cbind.data.frame, list(amph, aetdf, petdf, hetdf))

write.csv(sgdata.amph, "Schmidt_amphibian_data_JBiogeo.csv", row.names = F, quote = F)