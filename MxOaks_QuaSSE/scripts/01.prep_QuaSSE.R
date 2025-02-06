#####################
##### Load Data #####
#####################

# Load packages
library(diversitree)
library(ape)
library(phytools)
library(picante)
library(treeio)
library(tidyverse)
library(ggtree)
library(tidytree)
library(sf)
library(terra)
sf_use_s2(FALSE)

# Set functions
source('PrepData/SCRIPTS/functions.R')
`%notin%` = Negate(`%in%`)

## PART 1: LOAD DATA -----------------------------------------------------------
tree = read.tree(file = 'PrepData/DATA/20240325_p04_finaltree.tre')

tree$tip.label = tree$tip.label %>%
  str_split(pattern = "\\|") %>% 
  lapply(., head, n = 1) %>% 
  lapply(., str_replace_all, pattern = "_", replacement = " ") %>% 
  unlist() %>% 
  str_replace_all("'","")

# Load states data 
gbif = st_read('PrepData/OUT/20240325_gbif_cleaned.shp')

## PART 2: PREP DATA -----------------------------------------------------------
# Mutually prune occurrence data and tree tips for species we have data for
sp = intersect(x = unique(tree$tip.label), y = unique(gbif$species)) # Species in common

# Prune GBIf Data to those that are in the tree
gbif = gbif %>% 
  filter(species %in% sp)
gbif$species = gsub(x = gbif$species, " ", "_")

## PART 3: EXTRACT QUANTITATIVE STATE INFORMATION ------------------------------
## Extract data for maximum, minimum and mean latitude
# Make columns for latitude and longitude
gbif = gbif %>% 
  mutate(long = st_coordinates(gbif)[,1],
         lat = st_coordinates(gbif)[,2])


## PART 4: EXTRACT TRI ---------------------------------------------------------
elev = terra::rast("PrepData/DATA/wc2.1_30s_elev.tif")
elev = terra::crop(x = elev, y = st_bbox(gbif))
tri = terra::terrain(x = elev, v = "TRI")

# Extract TRI / species
gbif$elev = extract(x = elev, y = gbif)$wc2.1_30s_elev
gbif$tri = extract(x = tri, y = gbif)$TRI

# Summarize data
terrain_data = gbif %>% 
  drop_na(tri, elev) %>% 
  st_drop_geometry() %>% 
  group_by(species) %>% 
  summarise_at(.vars = vars("tri", "elev"),
               .funs = list(mean = mean,
                            min = min,
                            max = max))

## PART 4: SAVE DATA -----------------------------------------------------------
write_csv(x = terrain_data, file = "MX_oaks_QUASSE/DATA/QuaSSE_terrain.csv")

