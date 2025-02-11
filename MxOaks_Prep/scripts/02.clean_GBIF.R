#### Load Packages ####
library(rgbif)
library(tidyverse)
library(sf)
library(maps)
library(scrubr)
library(data.table)
source("PrepData/SCRIPTS/functions.R")

# Turn off spherical geometry
sf_use_s2(FALSE)

#### Clean Point Data ####
# Load Data
gbif = fread(file = 'PrepData/DATA/gbif/20240325_gbif_raw.csv')
checklist = read.csv(file = "PrepData/DATA/20230326_gbif_taxon_checklist.csv")
tree = ape::read.tree(file = "PrepData/DATA/20240325_p04_finaltree.tre")
tree$tip.label = remove_lab(tree) # remove MORTON labels

# Select only some subset of columns
gbif = gbif %>% 
  select(species,
         decimalLatitude,
         decimalLongitude,
         infraspecificEpithet, 
         coordinateUncertaintyInMeters, 
         verbatimScientificName, 
         taxonKey, 
         scientificName)

##########################  
## Clean species names ##
##########################
tree.sp = data.frame(treelab = tree$tip.label)
tree.sp$treelab_form = gsub(x = tree.sp$treelab, "_", " ") # Remove underscore

# Which species aren't represented in the GBIF species column?
tree.sp$treelab_form[tree.sp$treelab_form %notin% unique(gbif$species)]

#' Weird, 19 species that didn't match. Some of them are 100% spelling problems 
#' and problems with how GBIF organizes varieties and stuff. I think I'll have
#' to format this one-by-one to sort it out...

# Fix Quercus shumardii var. schneckii labels
gbif$species = ifelse(test = gbif$infraspecificEpithet == "schneckii",
       yes = "Quercus shumardii var. schneckii",
       no = gbif$species)

# Fix Quercus elliottii
gbif[str_detect(gbif$verbatimScientificName, "elliottii"),]$species = "Quercus elliottii"
# Fix Quercus gentryi
gbif[str_detect(gbif$verbatimScientificName, "gentryi"),]$species = "Quercus gentryi"
# Fix Quercus albocinta 
gbif[str_detect(gbif$verbatimScientificName, "albocincta"),]$species = "Quercus albocinta"
# Fix Quercus graciliformis
gbif[str_detect(gbif$verbatimScientificName, "graciliformis"),]$species = "Quercus graciliformis"
# Fix Quercus bolanyosensis 
gbif[str_detect(gbif$verbatimScientificName, "bolanyosensis"),]$species = "Quercus bolanyosensis"
# Fix Quercus crassifolia
gbif[str_detect(gbif$verbatimScientificName, "crassifolia"),]$species = "Quercus crassifolia s.l"
# Fix Quercus ocoteifolia
gbif[str_detect(gbif$verbatimScientificName, "ocoteifolia"),]$species = "Quercus ocoteifolia"
# Fix Quercus lowilliamsii
gbif[str_detect(gbif$verbatimScientificName, "lowilliamsii"),]$species = "Quercus lowilliamsii" # <----- THIS SPECIES MAY NEED TO BE REMOVED
# Fix Quercus bumelioides
gbif[str_detect(gbif$verbatimScientificName, "bumelioides"),]$species = "Quercus bumelioides"
# Fix Quercus eugeniifolia
gbif[gbif$species == "Quercus seemannii",]$species = "Quercus eugeniifolia"
# Fix Quercus cupreata
gbif[str_detect(gbif$verbatimScientificName, "cupreata"),]$species = "Quercus cupreata"
# Fix Quercus agrifolia var. oxyadenia
gbif[str_detect(gbif$infraspecificEpithet, "oxyadenia"),]$species = "Quercus agrifolia var. oxyadenia"
# Fix Quercus parvula var. parvula
gbif[str_detect(gbif$infraspecificEpithet, "parvula"),]$species = "Quercus parvula var. parvula"
# Fix Quercus parvula var. shrevei
gbif[str_detect(gbif$infraspecificEpithet, "shrevei"),]$species = "Quercus parvula var. shrevei"
# Fix Quercus aff. transmontana
gbif[str_detect(gbif$verbatimScientificName, "transmontana"),]$species = "Quercus aff. transmontana"
# Fix Quercus cf. striatula 
gbif[str_detect(gbif$verbatimScientificName, "striatula"),]$species = "Quercus cf. striatula"
# Fix Quercus cf. depressipes
gbif[str_detect(gbif$verbatimScientificName, "depressipes"),]$species = "Quercus cf. depressipes"
# Fix Quercus nudinervis
gbif[str_detect(gbif$verbatimScientificName, "nudinervis"),]$species = "Quercus nudinervis"

gbif[str_detect(gbif$verbatimScientificName, "sagraeana"),]

# Remove duplicates
gbif = gbif %>% 
  unite(species, decimalLatitude, decimalLongitude,
        col = "dupes",
        remove = FALSE) %>% 
  unique() %>% 
  select(-dupes)

# Remove problematic points
gbif = gbif %>% 
  coord_impossible() %>% 
  coord_unlikely() %>% 
  coord_incomplete()


##########################  
#### Spatial Cleaning ####
##########################

# Turn data into sf object
gbif_sf = st_as_sf(gbif, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# Load country shapefile
base_sf = st_read(dsn = 'PrepData/DATA/countries.shp')
base_sf = maps::map(database = "world", regions = "Cuba", fill = TRUE, plot = FALSE) %>% 
  st_as_sf() %>% 
  st_transform(crs = 4326) %>% 
  rename(geometry = geom) %>% 
  rbind(base_sf)

# Filter to points that occur in this spatial extent
gbif_sf = gbif_sf[base_sf,]

## Save data
write_csv(x = gbif_sf, file = "PrepData/OUT/20240325_gbif_cleaned.csv", append = FALSE)
st_write(obj = gbif_sf,dsn = "PrepData/OUT/20240325_gbif_cleaned.shp", append = FALSE)