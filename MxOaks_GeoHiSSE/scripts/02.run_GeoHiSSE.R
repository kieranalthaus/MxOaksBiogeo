GeoHiSSE_Caetano_models = function(phy, data, f, outdir, outname){
  library(hisse)
  library(ape)
  library(tidyverse)
  
  ##----------------------------- CLADOGENIC MODELS ------------------------------
  ## MODEL 1 ---------------------------------------------------------------------
  # CID - ORIGINAL GEOSSE
  turnover = c(1,1,1) 
  eps = c(1,1)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 0, 
                                     include.jumps = FALSE,
                                     separate.extirpation = FALSE)
  mod1 = GeoHiSSE(phy = phy, data = data, f = f, 
                  turnover = turnover, eps = eps, hidden.states = FALSE,
                  trans.rate = trans.rate, assume.cladogenetic = TRUE)
  
  saveRDS(object = mod1, file = paste0(outdir,outname,"_mod01.rds"))
  
  ## MODEL 2 ---------------------------------------------------------------------
  # CD - Original GeoSSE
  turnover = c(1,2,3)
  eps = c(1,2)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 0, 
                                     include.jumps = FALSE,
                                     separate.extirpation = FALSE)
  mod2 = GeoHiSSE(phy = phy, data = data, f = f,
                  turnover = turnover, eps = eps, hidden.states = FALSE,
                  trans.rate = trans.rate, assume.cladogenetic = TRUE)
  
  saveRDS(object = mod2, file = paste0(outdir,outname,"_mod02.rds"))
  
  ## MODEL 3 ---------------------------------------------------------------------
  # CID - GeoHiSSE, two hidden rate classes, null model
  turnover = c(1,1,1,2,2,2)
  eps = c(1,1,2,2)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 1, make.null = TRUE,
                                     include.jumps = FALSE,
                                     separate.extirpation = FALSE)
  mod3 = GeoHiSSE(phy = phy, data = data, f = f,
                  turnover = turnover, eps = eps, hidden.states = TRUE,
                  trans.rate = trans.rate, assume.cladogenetic = TRUE)
  
  saveRDS(object = mod3, file = paste0(outdir,outname,"_mod03.rds"))
  
  ## MODEL 4 ---------------------------------------------------------------------
  # Heterogeneous diversification, tied to range evolution
  # Assumes 6 distinct diversification rates
  turnover = c(1,2,3,4,5,6)
  eps = c(1,2,3,4)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 1,
                                     include.jumps = FALSE,
                                     separate.extirpation = FALSE)  
  mod4 = GeoHiSSE(phy = phy, data = data, f = f,
                  turnover = turnover, eps = eps, hidden.states = TRUE,
                  trans.rate = trans.rate, assume.cladogenetic = TRUE)
  
  saveRDS(object = mod4, file = paste0(outdir,outname,"_mod04.rds"))
  
  ## MODEL 5 ---------------------------------------------------------------------
  # Heterogeneous diversification, not tied to range evolution
  # Assumes 2 distinct diversification rates (2 hidden states)
  turnover = c(1,1,1,2,2,2)
  eps = c(1,1,2,2)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 1, 
                                     include.jumps = FALSE,
                                     separate.extirpation = FALSE)
  mod5 = GeoHiSSE(phy = phy, data = data, f = f,
                  turnover = turnover, eps = eps, hidden.states = TRUE,
                  trans.rate = trans.rate, assume.cladogenetic = TRUE)
  
  saveRDS(object = mod5, file = paste0(outdir,outname,"_mod05.rds"))
  
  #----------------------------- EXTRIPATION MODELS ------------------------------
  #' IN THIS BATCH OF MODELS, EXTRIPATION IS NOT LINKED TO RANGE REDUCTION
  # RANGE REDUCTION IS DIFFERENT FORM THE EXTRIPATION OF AN ENDEMIC LINEAGE
  
  ## MODEL 6 ---------------------------------------------------------------------
  # Disepersal paramters vary, no range-dependent diversification
  turnover = c(1,1,1)
  eps = c(1,1)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 0, 
                                     include.jumps = FALSE,
                                     separate.extirpation = TRUE)
  mod6 = GeoHiSSE(phy = phy, data = data, f = f,
                  turnover = turnover, eps = eps, hidden.states = FALSE,
                  trans.rate = trans.rate, assume.cladogenetic = TRUE)
  
  saveRDS(object = mod6, file = paste0(outdir,outname,"_mod06.rds"))
  
  ## MODEL 7 ---------------------------------------------------------------------
  #  GeoSSE model, with range effect on diversification
  turnover = c(1,2,3)
  eps = c(1,2)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 0, 
                                     include.jumps = FALSE,
                                     separate.extirpation = TRUE)
  mod7 = GeoHiSSE(phy = phy, data = data, f = f,
                  turnover = turnover, eps = eps, hidden.states = FALSE,
                  trans.rate = trans.rate, assume.cladogenetic = TRUE)
  
  saveRDS(object = mod7, file = paste0(outdir,outname,"_mod07.rds"))
  
  ## MODEL 8 ---------------------------------------------------------------------                
  #' Heterogeneous diversification, not tied to range evolution,
  #' Assumes two distinct diversification rates (2 hidden states)
  turnover = c(1,1,1,2,2,2)
  eps = c(1,1,2,2)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 1, make.null = TRUE,
                                     include.jumps = FALSE, separate.extirpation = TRUE)
  mod8 = GeoHiSSE(phy = phy, data = data, f = f, 
                  turnover = turnover, eps = eps, hidden.states = TRUE,
                  trans.rate = trans.rate, assume.cladogenetic = TRUE)
  
  saveRDS(object = mod8, file = paste0(outdir,outname,"_mod08.rds"))
  
  ## MODEL 9 --------------------------------------------------------------------
  #' Heterogeneous diversification, tied to range evolution
  #' Assumes 6 distinct diversification rates
  turnover = c(1,2,3,4,5,6)
  eps = c(1,2,3,4)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 1, include.jumps = FALSE,
                                     separate.extirpation = TRUE)
  mod9 = GeoHiSSE(phy = phy, data = data, f = f, 
                   turnover = turnover, eps = eps, hidden.states = TRUE,
                   trans.rate = trans.rate, assume.cladogenetic = TRUE)
  
  saveRDS(object = mod9, file = paste0(outdir,outname,"_mod09.rds"))
  
  ## MODEL 10 --------------------------------------------------------------------
  #' #' Heterogeneous diversification, not tied to range evolution
  #' #' Assumes two distinct diversification rates
  #' turnover = c(1,1,1,2,2,2)
  #' eps = c(1,1,2,2)
  #' trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 1, 
  #'                                    include.jumps = FALSE,
  #'                                    separate.extirpation = TRUE)
  #' mod10 = GeoHiSSE(phy = phy, data = data, f = f, 
  #'                  turnover = turnover, eps = eps, hidden.states = TRUE,
  #'                  trans.rate = trans.rate, assume.cladogenetic = TRUE)
  #' 
  #' saveRDS(object = mod10, file = paste0(outdir,outname,"_mod10.rds"))
  
  #----------------------------- ANAGENIC MODELS ---------------------------------
  #' ANAGENIC MODELS WHERE CHANGES ONLY HAPPEN ALONG BRANCHES, 
  #' AND THEREFORE RANGE CHANGE CAN'T LEAD TO SPECIATION
  
  ## MODEL 11 --------------------------------------------------------------------
  #' Transitions only. No character effect on diversification
  turnover = c(1,1,1)
  eps = c(1,1)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 0, include.jumps = FALSE,
                                     separate.extirpation = TRUE)
  mod11 = GeoHiSSE(phy = phy, data = data, f = f, 
                   turnover = turnover, eps = eps, hidden.states = FALSE,
                   trans.rate = trans.rate, assume.cladogenetic = FALSE)
  
  saveRDS(object = mod11, file = paste0(outdir,outname,"_mod11.rds"))
  
  ## MODEL 12 --------------------------------------------------------------------                
  #' Character effect on diversification
  turnover = c(1,2,3)
  eps = c(1,2)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 0,
                                     include.jumps = FALSE,
                                     separate.extirpation = TRUE)
  mod12 = GeoHiSSE(phy = phy, data = data, f = f,
                   turnover = turnover, eps = eps, hidden.states = FALSE,
                   trans.rate = trans.rate, assume.cladogenetic = FALSE)
  
  saveRDS(object = mod12, file = paste0(outdir,outname,"_mod12.rds"))
  
  ## MODEL 13 --------------------------------------------------------------------
  #' Character effect on diversification
  #' 2 hidden states
  turnover = c(1,2,3,4,5,6)
  eps = c(1,2,3,4)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 1, 
                                     include.jumps = FALSE,
                                     separate.extirpation = TRUE)
  mod13 = GeoHiSSE(phy = phy, data = data, f = f, 
                   hidden.states = TRUE, turnover = turnover, eps = eps,
                   trans.rate = trans.rate, assume.cladogenetic = FALSE)
  
  saveRDS(object = mod13, file = paste0(outdir,outname,"_mod13.rds"))
  
  ## MODEL 14 --------------------------------------------------------------------
  #' No character effect on diversification, multiple shits
  turnover = c(1,1,1,2,2,2)
  eps = c(1,1,2,2)
  trans.rate = TransMatMakerGeoHiSSE(hidden.traits = 1,
                                     include.jumps = FALSE,
                                     separate.extirpation = TRUE,
                                     make.null = TRUE)
  mod14 = GeoHiSSE(phy = phy, data = data, f = f,
                   hidden.states = TRUE, turnover = turnover, eps = eps,
                   trans.rate = trans.rate, assume.cladogenetic = FALSE)
  
  saveRDS(object = mod14, file = paste0(outdir,outname,"_mod14.rds"))
  
}

### RUN MODELS -----------------------------------------------------------------
#' widespread area 01 = 0
#' endemic area 00 = 1
#' endemic area 11 = 2

# Load tree
tree = read.tree(file = "MxOaks_Dating/out/20240325_p04_finaltree.tre")
tree$tip.label = tree$tip.label %>% 
  strsplit("\\|") %>% 
  lapply(function(x) x[1]) %>% 
  unlist()

# Load states
states = readRDS(file = "MxOaks_GeoHiSSE/data/GeoHiSSE_mntns_states.RData")
states = rbind(states, data.frame(species = "Quercus_sagraeana",
                                  range = 1))
states = states[states$species %in% tree$tip.label,]

## Run model
GeoHiSSE_Caetano_models(phy = tree,
                        data = states,
                        f = c(0.67, 0.26, 0.07), 
                        outdir = "MxOaks_GeoHiSSE/out/",
                        outname = "20240428_mntn_mods")