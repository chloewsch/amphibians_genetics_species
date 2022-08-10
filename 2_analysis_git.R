## Determinants of genetic diversity and species richness of North American amphibians
# Schmidt, Munshi-South, Dray, Garroway, J. Bioeography (2022) 

## 2. Analysis
library(tidyverse)

# Structural equation models
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(sf)
library(spdep)

# dbMEM
library(vegan)
library(adespatial)

# plots
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(extrafont)
library(viridis)

# Data prep ----
sgdata.amph.sem <- read.csv("Schmidt_amphibian_data_JBiogeo.csv")
sgdata.amph.sem$species <- as.factor(sgdata.amph.sem$species)

# Scale & center variables:
sgdata.amph.sem$genetic_diversity <- scale(sgdata.amph.sem$gene_diversity)
sgdata.amph.sem$species_richness <- scale(sgdata.amph.sem$poplvlsd)

# AET scales
sgdata.amph.sem$water_10 <- scale(sgdata.amph.sem$AET_10)
sgdata.amph.sem$water_25 <- scale(sgdata.amph.sem$AET_25)
sgdata.amph.sem$water_40 <- scale(sgdata.amph.sem$AET_40)
sgdata.amph.sem$water_80 <- scale(sgdata.amph.sem$AET_80)

# PET scales
sgdata.amph.sem$energy_10 <- scale(sgdata.amph.sem$PET_10)
sgdata.amph.sem$energy_25 <- scale(sgdata.amph.sem$PET_25)
sgdata.amph.sem$energy_40 <- scale(sgdata.amph.sem$PET_40)
sgdata.amph.sem$energy_80 <- scale(sgdata.amph.sem$PET_80)

# Heterogeneity scales
sgdata.amph.sem$heterogeneity_10 <- scale(sgdata.amph.sem$het_10)
sgdata.amph.sem$heterogeneity_25 <- scale(sgdata.amph.sem$het_25)
sgdata.amph.sem$heterogeneity_40 <- scale(sgdata.amph.sem$het_40)
sgdata.amph.sem$heterogeneity_80 <- scale(sgdata.amph.sem$het_80)

# Structural equation models ----
## 10km scale model ----
A0.10 <- psem(lmer(genetic_diversity ~ energy_10 + water_10 + heterogeneity_10 + 
                       (energy_10 + water_10 + heterogeneity_10||species), data = sgdata.amph.sem),
                lm(species_richness ~ energy_10 + water_10 + heterogeneity_10, data = sgdata.amph.sem)
) # singular fit with correlated slopes; don't estimate correlation (||)

summary(A0.10) # no further links suggested

# Check for spatially autocorrelated residuals:
xy <- dplyr::select(sgdata.amph.sem, lon, lat)
xysf <- st_as_sf(xy, coords = c("lon", "lat"), crs=4326)

# Generate neighborhood matrix
nb1 <- dnearneigh(xysf, 0, 529) # minimum distance (km) to connect all sites at least once

# Genetic diversity:
GDres <- residuals(lmer(genetic_diversity ~ energy_10 + water_10 + heterogeneity_10 + 
                          (energy_10 + water_10 + heterogeneity_10||species), data = sgdata.amph.sem))
moran.test(GDres, nb2listw(nb1)) # No

# Species richness:
SRres <- residuals(lm(species_richness ~ energy_10 + water_10 + heterogeneity_10, data = sgdata.amph.sem))
moran.test(SRres, nb2listw(nb1)) # Yes, but low: 0.05

## 25 km scale model ----
A0.25 <- psem(lmer(genetic_diversity ~ energy_25 + water_25 + heterogeneity_25 + 
                     (energy_25 + water_25 + heterogeneity_25||species), data = sgdata.amph.sem),
              lm(species_richness ~ energy_25 + water_25 + heterogeneity_25, data = sgdata.amph.sem)
)

summary(A0.25) # no further links suggested

# Check for spatially autocorrelated residuals:
# Genetic diversity:
GDres <- residuals(lmer(genetic_diversity ~ energy_25 + water_25 + heterogeneity_25 + 
                          (energy_25 + water_25 + heterogeneity_25||species), data = sgdata.amph.sem))
moran.test(GDres, nb2listw(nb1)) # No

# Species richness:
SRres <- residuals(lm(species_richness ~ energy_25 + water_25 + heterogeneity_25, data = sgdata.amph.sem))
moran.test(SRres, nb2listw(nb1)) # Yes, but low: 0.05

## 40 km scale model ----
A0.40 <- psem(lmer(genetic_diversity ~ energy_40 + water_40 + heterogeneity_40 + 
                     (energy_40 + water_40 + heterogeneity_40||species), data = sgdata.amph.sem),
              lm(species_richness ~ energy_40 + water_40 + heterogeneity_40, data = sgdata.amph.sem)
)

summary(A0.40) # no further links suggested

# Check for spatially autocorrelated residuals:
# Genetic diversity:
GDres <- residuals(lmer(genetic_diversity ~ energy_40 + water_40 + heterogeneity_40 + 
                          (energy_40 + water_40 + heterogeneity_40||species), data = sgdata.amph.sem))
moran.test(GDres, nb2listw(nb1)) # No

# Species richness:
SRres <- residuals(lm(species_richness ~ energy_40 + water_40 + heterogeneity_40, data = sgdata.amph.sem))
moran.test(SRres, nb2listw(nb1)) # Yes, but low: 0.06

## 80 km scale model ----
A0.80 <- psem(lmer(genetic_diversity ~ energy_80 + water_80 + heterogeneity_80 + 
                      (energy_80 + water_80 + heterogeneity_80||species), data = sgdata.amph.sem),
               lm(species_richness ~ energy_80 + water_80 + heterogeneity_80, data = sgdata.amph.sem)
)

summary(A0.80) # no further links suggested

# Check for spatially autocorrelated residuals:
# Genetic diversity:
GDres <- residuals(lmer(genetic_diversity ~ energy_80 + water_80 + heterogeneity_80 + 
                          (energy_80 + water_80 + heterogeneity_80||species), data = sgdata.amph.sem))
moran.test(GDres, nb2listw(nb1)) # No

# Species richness:
SRres <- residuals(lm(species_richness ~ energy_80 + water_80 + heterogeneity_80, data = sgdata.amph.sem))
moran.test(SRres, nb2listw(nb1)) # Yes, but low: 0.05


# Partialling out the spatial components of variation in genetic & species diversity ----

## Extract dbMEMs ##
# Check for linear gradient
amphxy <- select(sgdata.amph.sem, lon, lat)
anova(lm(sgdata.amph.sem$gene_diversity ~ ., data=amphxy)) #yes
anova(lm(sgdata.amph.sem$poplvlsd ~ ., data=amphxy)) #yes

# Detrend data
amph.det.gd  <-  resid(lm(sgdata.amph.sem$gene_diversity ~ .,   data=amphxy)) # gene diversity
amph.det.sp  <-  resid(lm(sgdata.amph.sem$poplvlsd ~ ., data=amphxy)) # species richness

# Construct the matrix of dbMEM variables
# dbMEM
amph.dbmem <- dbmem(amphxy, silent=FALSE)

## pull the spatial weights matrix out of the dbmem object
alistw <- attributes(amph.dbmem)$listw

# convert to dataframe
amph.dbmem <- as.data.frame(amph.dbmem)

# Select significant MEMs: Global significance test
### gene diversity
GDlm <- lm(amph.det.gd ~., amph.dbmem)
summary(GDlm)

(GDr2da <- RsquareAdj(GDlm)$adj.r.squared)
# forward selection
GDmemfwd <- forward.sel(amph.det.gd, as.matrix(amph.dbmem), 
                        adjR2thresh = GDr2da)
# sort & extract selected MEMs
(GDmems <- sort(GDmemfwd[,2]))
GDmem.red <- amph.dbmem[,GDmems]

### species richness
SPlm <- lm(amph.det.sp ~., amph.dbmem)
summary(SPlm)

(SPr2da <- RsquareAdj(SPlm)$adj.r.squared)
# forward selection
SPmemfwd <- forward.sel(amph.det.sp, as.matrix(amph.dbmem), 
                        adjR2thresh = SPr2da)
# sort & extract selected MEMs
(SPmems <- sort(SPmemfwd[,2]))
SPmem.red <- amph.dbmem[,SPmems]

shared.df <- GDmems[GDmems %in% SPmems]
shared <- amph.dbmem[,shared.df]

## Variation partitioning ##
gd_mem <- lm(scale(sgdata.amph.sem$gene_diversity) ~., data= GDmem.red)
sp_mem <- lm(scale(sgdata.amph.sem$poplvlsd) ~., data = SPmem.red)

shareGD <- lm(scale(sgdata.amph.sem$gene_diversity) ~., data= shared) ## same as gd_mem because all GD MEMs are shared
shareSR <- lm(scale(sgdata.amph.sem$poplvlsd) ~., data= shared)


## Plot ##
(gATS <- summary(gd_mem)$adj.r.squared) # total spatial
(gANS <- 1- gATS) # nonspatial 
(gASS <- summary(shareGD)$adj.r.squared) # shared spatial
(gANSS <- gATS - gASS) # non-shared spatial
gASS/gATS # % shared spatial out of total spatial

# SD: 
(ATS <- summary(sp_mem)$adj.r.squared) # total spatial
(ANS <- 1- ATS) # nonspatial
(ASS <- summary(shareSR)$adj.r.squared) # shared spatial
(ANSS <- ATS - ASS) # non-shared spatial
ASS/ATS # % shared spatial out of total spatial

variable <- c(rep("gene_diversity", 3), rep("species_diversity", 3))
variation <- factor(rep(c("shared_spatial", "non_shared_spatial", "var_non_spatial"), 2), 
                    levels = c("shared_spatial", "non_shared_spatial", "var_non_spatial"))
value <- c(gASS, gANSS, gANS, ASS, ANSS, ANS)
sdgdR2tab <- data.frame(variable, variation, value)

pal <- c("#023020", "#228B22", "#AFE1AF")
q <- ggplot(sdgdR2tab, aes(fill = variation, y=value, x=variable)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = pal,
                    name= "variation", labels = c("shared spatial", "non-shared spatial", "non-spatial")) +
  scale_x_discrete(labels=c("genetic diversity", "species richness")) +
  labs(title = "", x= "", y= "proportion of variation") +
  theme(text=element_text(size=16,  family="Lato Black"))

# Predicted maps ----
## GD and SR ----
## Moran's I for each MEM to assign cut-off for broad scale Moran's I (I > 0.25) ##
moranig <- moran.randtest(as.matrix(GDmem.red), alistw)
morani_gd <- moranig$obs[moranig$obs > 0.25] # Moran's I values

moranis <- moran.randtest(as.matrix(SPmem.red), alistw)
morani_sp <- moranis$obs[moranis$obs > 0.25]

gmem.broad <- amph.dbmem[,GDmems[c(1:length(morani_gd))]]
smem.broad <- amph.dbmem[,SPmems[c(1:length(morani_sp))]]


## Predicted values for maps
fit <- lm(scale(sgdata.amph.sem$gene_diversity)~gmem.broad) 
fittedGD <- predict(fit)

fitsp <- lm(scale(sgdata.amph.sem$poplvlsd)~., data = smem.broad)
fittedSP <- predict(fitsp)


## basemap
world <- ne_countries(scale="medium", returnclass = 'sf')
canadausa <- world[world$sovereignt=="Canada" | world$sovereignt=="United States of America" & world$type !="Dependency",]
rm(world)

# Remove Hawaii:
USA <- canadausa %>% 
  filter(sov_a3=="US1") %>% 
  st_cast("POLYGON") %>% 
  mutate(sub_id = as.factor(paste("piece", 1:nrow(.), sep = "")))

USA_nohi <- USA %>% 
  filter(!sub_id %in% c('piece1', 'piece2', 'piece3',
                        'piece4', 'piece5', 'piece6', 'piece7'))
USA_nohic <- USA_nohi %>% 
  st_combine() %>% 
  st_as_sf()

rm(USA, USA_nohi, USA_nohic)

# Replace USA in region data with new one:
st_geometry(canadausa[2,]) <- st_geometry(USA_nohic)

# Convert data to sf for plot:
plotdata <- sgdata.amph.sem %>% 
  select(lon, lat, gene_diversity, poplvlsd) %>% 
  mutate(fittedGD = fittedGD,
         fittedSP = fittedSP)
plotdata <- st_as_sf(plotdata, coords = c('lon', 'lat'), crs = 4326)

# Genetic diversity
gdpoint <- ggplot() + 
  geom_sf(data = canadausa, fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +  #remove grey background
  geom_sf(data = plotdata, aes(fill = fittedGD), shape=21, size=0.2) +
  scale_fill_viridis(option = "D", breaks=c(0,1),labels=c("low","high"),
                     limits=c(0,1)) +
  labs(fill="") +
  geom_sf(data = plotdata, aes(color = fittedGD), show.legend = FALSE, size=2) +
  scale_color_viridis(option = "D") +
  guides(fill = guide_colourbar(ticks = FALSE)) +
  #labs(title = "Genetic diversity") +
  coord_sf(crs = 'ESRI:102008') +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Species richness 
sppoint <- ggplot() + 
  geom_sf(data = canadausa, fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +  #remove grey background
  geom_sf(data = plotdata, aes(fill = fittedSP), shape=21, size=0.2) +
  scale_fill_viridis(option = "D", breaks=c(0,1),labels=c("low","high"),
                     limits=c(0,1)) +
  labs(fill="") +
  geom_sf(data = plotdata, aes(color = fittedSP), show.legend = FALSE, size=2) +
  scale_color_viridis(option = "D") +
  guides(fill = guide_colourbar(ticks = FALSE)) +
  #labs(title = "Species richness") +
  coord_sf(crs = 'ESRI:102008') +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## FST ----
amphfst <- amph %>%
  drop_na(site_fst)
xyfst <- amphfst %>% select(lon, lat)

amph.det.fst  <-  resid(lm(amphfst$site_fst ~ ., data=xyfst))

# Construct the matrix of dbMEM variables
# dbMEM
amph.dbmem <- dbmem(xyfst, silent=FALSE)

## pull the spatial weights matrix out of the dbmem object
alistw <- attributes(amph.dbmem)$listw

# convert to dataframe
amph.dbmem <- as.data.frame(amph.dbmem)

# Select significant MEMs: Global significance test
### gene diversity
fstlm <- lm(amph.det.fst ~., amph.dbmem)
summary(fstlm)

(fstr2da <- RsquareAdj(fstlm)$adj.r.squared)
# forward selection
fstmemfwd <- forward.sel(amph.det.fst, as.matrix(amph.dbmem), 
                         adjR2thresh = fstr2da)
# sort & extract selected MEMs
(fstmems <- sort(fstmemfwd[,2]))
fstmem.red <- amph.dbmem[,fstmems]


moranif <- moran.randtest(as.matrix(fstmem.red), alistw)
morani_fst <- moranif$obs[moranif$obs > 0.25] # Moran's I values

fmem.broad <- amph.dbmem[,fstmems[c(1:length(morani_fst))]]

fitf <- lm(scale(amphfst$site_fst)~fmem.broad) 
fittedFST <- predict(fitf)

plotdataF <- amphfst %>% 
  select(lon, lat, site_fst) %>% 
  mutate(fittedFST = fittedFST)

plotdataF <- st_as_sf(plotdataF, coords = c('lon', 'lat'), crs = 4326)
fstpoint <- ggplot() +
  geom_sf(data = canadausa, fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +  #remove grey background
  geom_sf(data = plotdataF, aes(fill = fittedFST), shape=21, size=0.2) +
  scale_fill_viridis(option = "D", breaks=c(0,1),labels=c("low","high"),
                     limits=c(0,1)) +
  labs(fill="") +
  geom_sf(data = plotdataF, aes(color = fittedFST), show.legend = FALSE, size=2) +
  scale_color_viridis(option = "D") +
  guides(fill = guide_colourbar(ticks = FALSE)) +
  #labs(title = "Population differentiation") +
  coord_sf(crs = 'ESRI:102008') +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Effect of heterogeneity on FST ----
## 10 fst
mod10 <- lm(scale(amphfst$het_10) ~ MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11, data = fstmem.red)
summary(mod10)# 2, 5, 7

FSTmod10 <- lmer(scale(amphfst$site_fst) ~ scale(amphfst$het_10) + MEM4 + MEM9 + MEM10 + MEM11 +
                     (scale(amphfst$het_10)|amphfst$species), data = fstmem.red)
summary(FSTmod10)

## 25 fst
mod25 <- lm(scale(amphfst$het_25) ~ MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11, data = fstmem.red)
summary(mod25) # 2, 5, 7

FSTmod25 <- lmer(scale(amphfst$site_fst) ~ scale(amphfst$het_25) + MEM4 + MEM9 + MEM10 + MEM11 +
                   (scale(amphfst$het_25)|amphfst$species), data = fstmem.red)
summary(FSTmod25)

## 40 fst
mod40 <- lm(scale(amphfst$het_40) ~ MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11, data = fstmem.red)
summary(mod40) # 2, 5, 7, 10

FSTmod40 <- lmer(scale(amphfst$site_fst) ~ scale(amphfst$het_40) + MEM4 + MEM9 + MEM11 +
                   (scale(amphfst$het_40)|amphfst$species), data = fstmem.red)
summary(FSTmod40)

## 80 fst
mod80 <- lm(scale(amphfst$het_80) ~ MEM2 + MEM4 + MEM5 + MEM7 + MEM9 + MEM10 + MEM11, data = fstmem.red)
summary(mod80) # 2, 4, 5, 7, 10

FSTmod80 <- lmer(scale(amphfst$site_fst) ~ scale(amphfst$het_80) + MEM9 + MEM11 +
                    (scale(amphfst$het_80)|amphfst$species), data = fstmem.red)
summary(FSTmod80)

# Raw data plots ####
## gene diversity
gdraw <- ggplot() + 
  geom_sf(data = canadausa, fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +  #remove grey background
  geom_sf(data = plotdata, aes(color = gene_diversity), size=2) +
  scale_color_viridis(option = "D") +
  labs(title = "Gene diversity",
       color = "") +
  coord_sf(crs = 'ESRI:102008') +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## species richness
spraw <- ggplot() + 
  geom_sf(data = canadausa, fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +  #remove grey background
  geom_sf(data = plotdata, aes(color = poplvlsd), size=2) +
  scale_color_viridis(option = "D") +
  labs(title = "Species richness",
       color = "") +
  coord_sf(crs = 'ESRI:102008') +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## FST
fstraw <- ggplot() + 
  geom_sf(data = canadausa, fill = "#b7c0c7", colour = "#b7c0c7", alpha = 1, size = 0.9) +
  scale_y_continuous(breaks = NULL) + #remove y-axis
  scale_x_continuous(breaks = NULL) + #remove x-axis
  xlab("") + #remove x-axis label
  ylab("") + #remove y-axis label
  theme(panel.background = element_blank()) +  #remove grey background
  geom_sf(data = plotdataF, aes(color = global_fst), size=2) +
  scale_color_viridis(option = "D") +
  labs(title = "Genetic differentiation",
       color = "") +
  coord_sf(crs = 'ESRI:102008') +
  theme(panel.background = element_rect(fill='transparent',colour=NA), 
        plot.background = element_rect(fill='transparent',colour=NA),
        text=element_text(size=14,  family="Lato Black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())