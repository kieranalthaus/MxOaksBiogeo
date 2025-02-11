#### Load Packages ####
library(rgbif)
library(tidyverse)
library(ape)

# ### Load taxa from tree
tree = read.tree(file = 'PrepData/DATA/20240325_p04_finaltree.tre')

# Reformat tip labels
tree$tip.label = tree$tip.label %>%
  str_split(pattern = "\\|") %>%
  lapply(., head, n = 1) %>%
  unlist()

### Check my taxonomy against GBIF ###
checklist = name_backbone_checklist(name_data = tree$tip.label) %>% 
  filter(verbatim_name != "Quercus_crassifolia_s.l") %>% 
  select(usageKey, scientificName, canonicalName, speciesKey, verbatim_name)

# For some reason Q. crassifolia kept being left behind...
crass = name_backbone_checklist(name_data = "Quercus_crassifolia") %>% 
  select(usageKey, scientificName, canonicalName, speciesKey, verbatim_name)

checklist = rbind(checklist, crass)

#### Create query ####
# Get genus key
usageKey = checklist$usageKey

# Format Query
user = ""
pwd = ""
email = ""

query = occ_download_prep(pred("hasCoordinate", TRUE),
                          pred("hasGeospatialIssue", FALSE),
                          pred_in("taxonKey", usageKey),
                          format = "SIMPLE_CSV",
                          user = user,
                          pwd = pwd,
                          email = email)

# Submit Query
occ_download(body = query$json_request,
             user = user,
             pwd = pwd,
             email = email)

### 2024-03-25 DOWNLOAD ####
#' DOI: 10.15468/dl.e7jn4e
#' DOWNLOAD KEY: 0033983-240321170329656
#' URL: https://www.gbif.org/occurrence/download/0033983-240321170329656

