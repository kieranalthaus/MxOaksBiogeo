# load packages
library(terra)
library(sf)
`%notin%` = Negate(`%in%`)

## 1. Load data ----------------------------------------------------------------
# load occurrence data
occs.sf = st_read("PrepData/OUT/20240325_gbif_cleaned.shp")

# path to climate rasters
env.list = list.files("",  full.names = TRUE)
# path to soil rasters
soil.list = list.files("", full.names = TRUE, pattern = ".nc4")
# path to elevation rasters
elev = ""

## Read in data
# Read in just some of the soil data
soil = rast(x = soil.list[c(1,2,19,20,21,22,23,24,25,27,28)])
env = rast(x = env.list)
elev = rast(x = elev)

# Crop rasters to extent of occurrence
env = crop(x = env, y = st_bbox(occs.sf))
soil = crop(x = soil, y = st_bbox(occs.sf))

# Extract ecological data
env.data = terra::extract(x = env, y = occs.sf) %>%
  as.data.frame() %>%
  dplyr::select(-ID)
soil.data = terra::extract(x = soil, y = occs.sf) %>%
  as.data.frame() %>%
  dplyr::select(-ID)
elev.data = terra::extract(x = elev, y = occs.sf) %>%
  as.data.frame() %>%
  dplyr::select(-ID)

# Merge data
eco.data = cbind(env.data, soil.data, elev.data)

# Species column
eco.data$species = occs.sf$species

# Create longitude and latitude columns
eco.data = eco.data %>%
  mutate(longitude = st_coordinates(occs.sf)[,1],
         latitude = st_coordinates(occs.sf)[,2])

# Remove NAs
eco.data = drop_na(eco.data)

# Save data
# write_csv(x = eco.data, 
#           file = "20240429_oak_envdata_complete.csv", 
#           append = FALSE)












