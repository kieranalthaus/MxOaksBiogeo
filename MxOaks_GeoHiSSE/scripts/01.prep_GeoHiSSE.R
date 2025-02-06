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
`%notin%` = Negate(`%in%`)

# Load Tree
tree = read.tree(file = 'MxOaks_Dating/out/20240321_p04_finaltree.tre') # As phylo object

## Load states data
gbif = st_read('MxOaks_Prep/out/20240325_gbif_cleaned.shp') %>% 
  mutate(species = str_replace_all(species, " ", "_"))

# Edit the labels to just be species names
tree$tip.label = tree$tip.label %>%
  str_split(pattern = "\\|") %>% 
  lapply(., head, n = 1) %>% 
  lapply(., str_replace_all, pattern = "_", replacement = " ") %>% 
  unlist() %>% 
  str_replace_all(pattern = "'", "")

#####################
##### Prep Data #####
#####################

# Save data
get_states(occurrences = gbif, 
           shapefile = countries, 
           threshold = 0.3,
           analysis = "SSE",
           out.dir = "MxOaks_GeoHiSSE/data/",
           out.name = "GeoHiSSE_mntns_states")

# Save tree
write.tree(phy = tree_pruned, file = "MxOaks_GeoHiSSE/data/GeoHiSSE_tree.tre")











