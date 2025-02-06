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

# Set functions
source('MxOaks_Prep/scripts/00.functions.R')
`%notin%` = Negate(`%in%`)

## PART 1: LOAD DATA -----------------------------------------------------------
tree = read.tree(file = 'MxOaks_Prep/data/20240325_p04_finaltree.tre') # As phylo object)
tree$tip.label = remove_lab(tree) # Remove MORTON number
tree$tip.label = str_replace_all(string = tree$tip.label, pattern = "_", replacement = " ")

# Load states data 
gbif = st_read('MxOaks_Prep/out/20240325_gbif_cleaned.shp') 

# Load shapefile
shapefile = st_read(dsn = 'MxOaks_DEC/data/DEC_7_bioregions.shp')
# Rename id column to "state"
shapefile = shapefile %>% 
  rename(state = id)
# Combine shapefiles into one 
shapefile = shapefile %>% 
  split(.$state) %>% 
  lapply(FUN = st_union) %>% 
  do.call(c,.) %>% 
  st_cast() %>% 
  st_as_sf() %>% 
  dplyr::mutate(state = 1:nrow(.))

# Relabel states
states = c("4", "2", "5", "6", "3", "1")
shapefile = shapefile %>% 
  st_as_sf() %>% 
  mutate(state = as.character(states)) %>% 
  rename(geometry = x) %>% 
  arrange(state)

## PART 2: SAVE DATA -----------------------------------------------------------
get_states(occurrences = gbif, shapefile = shapefile, threshold = 0.3,
           analysis = "DEC", out.dir = 'MxOaks_DEC/data/', out.name = "DEC_states")

# Save tree
tree$tip.label = gsub(x = tree$tip.label, " ", "_")
write.tree(phy = tree, file = "MxOaks_DEC/data/DEC_tree.tre")


# State codes as understood by RevBayes and the DEC model
# [
#   1000000 = Eastern U.S.
#   0100000 = Texas
#   0010000 = Sierra Madre Oriental
#   0001000 = Tropics
#   0000100 = Sierra Madre Occidental
#   0000010 = Arizona
#   0000001 = California
# ]
