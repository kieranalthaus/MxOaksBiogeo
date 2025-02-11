#' @param occurrences an sf object of occurrence points
#' @param shapefile an sf shapefile. Must have a "state" column indicating state 0 or 1
#' @param threshold a value between 0 and 1 at which the occurrence cutoff will be applied
#' @param analysis either "DEC" or "SSE", which specifies the output format of the data

get_states = function(occurrences, shapefile, threshold, analysis, out.dir, out.name, type = "NULL"){
  sf_use_s2(FALSE)
  # Create output dataframe, with columns n+1 of input df
  out = data.frame(matrix(data = NA, ncol = ncol(occurrences)+1, 
                          dimnames = list(c(),c(colnames(occurrences), "state"))))

## SSE FORMATING ---------------------------------------------------------------
  if (analysis == "SSE") {
    suppressMessages(
    for(i in 1:nrow(shapefile)){ # Loop over individual areas
      # Subset regions into smaller area i
      tmp_shp = shapefile[i,]
      # Filter occurrences by those that occur in area i
      tmp_occ = occurrences[tmp_shp,]
      if (nrow(tmp_occ) == 0) next
      # Transfer state info from shapefile to occurrences
      tmp_occ$state = tmp_shp$state
      # Bind data to output
      out = rbind(out, tmp_occ)
    }) # END OF SSE FORLOOP

    # Format states
    states = out %>% 
      drop_na(state) %>% 
      mutate(species = gsub(x = species, " ", "_")) %>% 
      group_by(species) %>% 
      reshape2::dcast(species~state) %>% 
      drop_na(species) %>%
      column_to_rownames(var = "species") %>% 
      apply(MARGIN = 1, FUN = proportions) %>% 
      t() %>% 
      `>`(threshold) %>% 
      `*`(1)
    
    if(any(rowSums(states) == 3)) {
      print("ERROR: Maximum number of areas exeeded")
      break
    }
    
    # Format states 50-50
    states_df = df_50_50(states)
    
    if(analysis == "SSE" & type == "REV"){
      # Save data for use in RevBayes
      states.rev = lapply(X = states_df, FUN = function(x)
        x %>% 
          apply(MARGIN = 1, which.max) %>% 
          `-`(1) %>% 
          as.character() %>% 
          `names<-`(rownames(states))
      )
      
      # Write state data
      ape::write.nexus.data(x = states.rev$df1, file = paste0(out.dir,out.name,"_rev1",".nex"))
      ape::write.nexus.data(x = states.rev$df2, file = paste0(out.dir,out.name,"_rev2",".nex"))
      cat("WARNING: Before reading NEXUS file into RevBayes, replace the FORMAT row with: \n FORMAT DATATYPE=STANDARD symbols='012345' MISSING=? GAP=-;")
      
    } else {
      
      # Save data for use in R
      states_r = states_df %>%
        lapply(X = ., FUN = function(x)
          apply(X = x, MARGIN = 1, FUN = which.max))
      
      
      saveRDS(states_r$df1, file = paste0(out.dir,out.name,"_1",".RData"))
      saveRDS(states_r$df2, file = paste0(out.dir,out.name,"_2",".RData"))
      # saveRDS(states_df$df1, file = paste0(out.dir,out.name,"_1",".RData"))
      # saveRDS(states_r$df2, file = paste0(out.dir,out.name,"_2",".RData"))
      }
    
## DEC FORMATING ---------------------------------------------------------------
  } else{
    
    for(i in 1:nrow(shapefile)){ # Loop over individual areas
      # Subset regions into smaller area i
      tmp_shp = shapefile[i,]
      # Filter occurrences by those that occur in area i
      tmp_occ = occurrences[tmp_shp,]
      # Skip forloop iteration if no species occur in region i
      if (nrow(tmp_occ) == 0) next
      # Transfer state info from shapefile to occurrences
      tmp_occ$state = tmp_shp$state
      # Bind data to output
      out = rbind(out, tmp_occ)
    } # END OF SSE FORLOOP
  
    # Format species name data
    out$species = out$species %>% 
      str_replace_all(pattern = " ",
                      replacement = "_")
    
    # Format states
    state = out %>% 
      drop_na(state,species) %>% 
      group_by(species) %>% 
      reshape2::dcast(species~state) %>%
      column_to_rownames(var = "species") %>% 
      apply(MARGIN = 1, FUN = proportions) %>% 
      t() %>% 
      `>`(0.3) %>% 
      `*`(1)
    
    sp = rownames(state)
    
    state = state %>% 
      apply(MARGIN = 2, as.character) %>% 
      as.data.frame() %>% 
      unite(1:ncol(state), remove = TRUE, col = "state", sep = "") %>% 
      .[,1] %>% 
      `names<-`(sp)
  
    # Write state data
    ape::write.nexus.data(x = state, file = paste0(out.dir,out.name,".nex"))
    cat("WARNING: Before reading NEXUS file into RevBayes, replace the FORMAT row with: \n FORMAT DATATYPE=STANDARD MISSING=? GAP=- LABELS='01';")
   } #  END OF ELSE STATEMENT
} # FUNCTION END


#' @param df a dataframe with rownames and only numeric columns
#' @description
#' This function takes in a dataframe with multiple 1's per row, and splits
#' it into two dataframes, each df with a single 1 / row
df_50_50 <- function(df) {
  df1 <- df # Copy of the original dataframe for the first version
  df2 <- df # Copy for the second version
  
  for (i in 1:nrow(df)) {
    row <- df[i, ]
    if (sum(row == 1) > 1) { # More than one 1 in the row
      ones_indices <- which(row == 1) # Indices of 1's
      # For df1, drop the first 1 (by making it 0)
      df1[i, ones_indices[1]] <- 0
      # For df2, drop the second 1
      df2[i, ones_indices[2]] <- 0
    }
  }
  
  list(df1 = df1, df2 = df2)
}


remove_lab = function(tree){
  return(tree$tip.label %>%
    str_split(pattern = "\\|") %>% 
    lapply(., head, n = 1) %>% 
    unlist())
}


#' @param phy a phylogeny 
#' @param data a dataframe with three columns: 1) Sp name, 2) discrete character 3) continuous character
#' @param outdir the name of the output directory
#' @param outname the name of the output RData file
#' @description:
#' Models the correlated evolution of a discrete character, and it's affect on 
#' the rage of evolution of a continuous character. Included in this modeling 
#' framework are a host of character-dependent and character-independent (CID)
#' models of trait evolution

DIY_hOUwie_models = function(phy, data, outdir, outname, nSim){
  
  library(OUwie)
  library(data.table)
  library(corHMM)
  library(parallel)
  library(expm)
  library(phylolm)
  library(expm)
  library(doMC)
  library(ape)
  library(phytools)
  library(caper)
  source("hOUwie/scripts/functions/hOUwie.internal.R")
  source("hOUwie/scripts/functions/hOUwie.R")
  
  states = length(unique(data[,2]))
  
##
## CID Models ------------------------------------------------------------------
##

## Model 1: CID, OU model with a single optimum for all species
  try(expr = {mod1.cid = OUwie::hOUwie(
    phy = phy,
    data = data,
    rate.cat = 2,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OU1", states, 2, T),
    null.model = TRUE,
    nSim = nSim,
    adaptive_sampling = TRUE
  )})

## Model 2: CID, OU model with different state means and a single alpha and sigma-squared on all regimes
  try({mod2.cid = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 2,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUM", states, 2, T),
    null.model = TRUE,
    nSim = nSim,
    adaptive_sampling = TRUE
  )})

## Model 3: CID, OU model with different state means and multiple alphas
  try({mod3.cid = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 2,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUMA", states, 2, T),
    null.model = TRUE,
    nSim = nSim,
    adaptive_sampling = TRUE,
  )})

## Model 4: CID, OU model with different state means and multiple variances
  try({mod4.cid = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 2,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUMV", states, 2, T),
    null.model = TRUE,
    nSim = nSim,
    adaptive_sampling = TRUE,
  )})

## Model 5: CID, OU model with different state means, alphas and sigma2
  try({mod5.cid = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 2,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUMVA", states, 2, T),
    null.model = TRUE,
    nSim = nSim,
    adaptive_sampling = TRUE,
  )})

##
## CD Models, 0 hidden states --------------------------------------------------
##

## Model 6: CD, OU1 model with a single optimum for all species
  try({mod6.cd = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 1,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OU1", states, 1, F),
    null.model = FALSE,
    nSim = nSim,
    adaptive_sampling = TRUE
  )})

## Model 7: CD, OU model with different state means and a single alpha and sigma-squared on all regimes
  try({mod7.cd = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 1,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUM", states, 1, F),
    null.model = FALSE,
    nSim = nSim,
    adaptive_sampling = TRUE,
  )})

## Model 8: CID, OU model with different state means and multiple alphas
  try({mod8.cd = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 1,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUMA", states, 1, F),
    null.model = FALSE,
    nSim = nSim,
    adaptive_sampling = TRUE,
  )})

## Model 9: CID, OU model with different state means and multiple variances
  try({mod9.cd = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 1,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUMV", states, 1, F),
    null.model = FALSE,
    nSim = nSim,
    adaptive_sampling = TRUE,
  )})

## Model 10: CID, OU model with different state means, alphas and sigma2
  try({mod10.cd = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 1,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUMVA", states, 1, F),
    null.model = FALSE,
    nSim = nSim,
    adaptive_sampling = TRUE,
  )})

##
## CD Models, 2 hidden states --------------------------------------------------
##

## Model 11: CD, OU model with different state means and a single alpha and sigma-squared on all regimes
  try({mod11.h = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 2,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUM", states, 2, F),
    null.model = FALSE,
    nSim = nSim,
    adaptive_sampling = TRUE,
  )})

## Model 12: CID, OU model with different state means and multiple alphas
  try({mod12.h = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 2,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUMA", states, 2, F),
    null.model = FALSE,
    nSim = nSim,
    adaptive_sampling = TRUE,
  )})

## Model 13: CID, OU model with different state means and multiple variances
  try({mod13.h = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 2,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUMV", states, 2, F),
    null.model = FALSE,
    nSim = nSim,
    adaptive_sampling = TRUE,
  )})

## Model 14: CID, OU model with different state means, alphas and sigma2
  try({mod14.h = hOUwie(
    phy = phy,
    data = data,
    rate.cat = 2,
    discrete_model = "ARD",
    continuous_model = getOUParamStructure("OUMVA", states, 2, F),
    null.model = FALSE,
    nSim = nSim,
    adaptive_sampling = TRUE
  )})

# Get model list
  models.list = list(cid.ou1 = mod1.cid,
                     cid.oum = mod2.cid,
                     cid.ouma = mod3.cid,
                     cid.oumv = mod4.cid,
                     cid.oumva = mod5.cid,
                     cd.ou1 = mod6.cd,
                     cd.oum = mod7.cd,
                     cd.ouma = mod8.cd,
                     cd.oumv = mod9.cd,
                     cd.oumva = mod10.cd,
                     mod11.h = mod11.h,
                     mod12.h = mod12.h,
                     mod13.h = mod13.h,
                     mod14.h = mod14.h
  )

  # Save models
saveRDS(object = models.list, file = paste0(outdir,outname,".RData"))
}


# Function to get significant p-values between comparisons
pmat = function(dstat.mat, tree, section){
  
  # Create output matrix
  pval.mat = matrix(data = NA,
                    nrow = length(unique(dstat.mat$P1)),
                    ncol = length(unique(dstat.mat$P1)),
                    dimnames = list(
                      c(unique(dstat.mat$P1)),
                      c(unique(dstat.mat$P1))
                    ))
  pval.mat = pval.mat[,!is.na(colnames(pval.mat))]
  pval.mat = pval.mat[!is.na(rownames(pval.mat)),]
  
  i = 1
  j = 1
  # Run forloop for p-values
  for(i in 1:nrow(pval.mat)){
    for(j in 1:ncol(pval.mat)){
      
      namei = rownames(pval.mat)[i]
      namej = colnames(pval.mat)[j]
      
      # If it's a diagonal, make it
      if(namei == namej){
        pval.mat[i,j] = 0
      } else {
        
        # Filter dataframe to species P2i and P3j
        tmp = dstat.mat %>% 
          filter(P2 %in% c(namei, namej) & P3 %in% c(namei, namej))
        
        # Get the p-value
        pval = tmp[which.min(tmp$p.value),]$p.value
        
        if(length(pval) == 0){
          pval.mat[i,j] = 0
        } else{
          pval.mat[i,j] =  pval 
        }
      } # ELSE
    } # J
    print(i)
  } # I
  
  
  if(section == "white"){
    # Adding a new row with the name "Quercus_rubra"
    new_row <- matrix(0, nrow=1, ncol=74)  # New row with 74 columns initialized with 0
    rownames(new_row) <- "Quercus_rubra"
    
    # Adding a new column with the name "Quercus_rubra"
    new_col <- matrix(0, nrow=75, ncol=1)  # New column with 75 rows (including the new row) initialized with 0
    colnames(new_col) <- "Quercus_rubra"
    
    matrix_pval = rbind(pval.mat, new_row)
    matrix_pval = cbind(matrix_pval, new_col)
    
  } else { # If red oaks 
    # Adding a new row with the name "Quercus_lobata"
    new_row <- matrix(0, nrow=1, ncol=ncol(pval.mat))  # New row with 74 columns initialized with 0
    rownames(new_row) <- "Quercus_lobata"
    
    # Adding a new column with the name "Quercus_lobata"
    new_col <- matrix(0, nrow=nrow(pval.mat)+1, ncol=1)  # New column with 75 rows (including the new row) initialized with 0
    colnames(new_col) <- "Quercus_lobata"
    
    matrix_pval = rbind(pval.mat, new_row)
    matrix_pval = cbind(matrix_pval, new_col)
  }
  
  # Organize output by tree order
  tree = ladderize(tree)
  matrix_pval = matrix_pval[tree$tip.label,tree$tip.label]
  return(matrix_pval)
  
} # END OF FUNCTION

### D-STAT FUNCTION ------------------------------------------------------------
dstat = function(dstat.mat, tree, section){
  pairwise.mat = matrix(data = NA,
                        nrow = length(unique(dstat.mat$P1)),
                        ncol = length(unique(dstat.mat$P1)),
                        dimnames = list(
                          c(unique(dstat.mat$P1)),
                          c(unique(dstat.mat$P1))
                        ))
  
  pairwise.mat = pairwise.mat[,!is.na(colnames(pairwise.mat))]
  pairwise.mat = pairwise.mat[!is.na(rownames(pairwise.mat)),]
  
  i = 1
  j = 1
  for(i in 1:nrow(pairwise.mat)){
    for(j in 1:ncol(pairwise.mat)){
      
      namei = rownames(pairwise.mat)[i]
      namej = colnames(pairwise.mat)[j]
      
      # If it's a diagonal, make it
      if(namei == namej){
        pairwise.mat[i,j] = 0
      } else {
        
        # Filter dataframe to species P2i and P3j
        tmp = dstat.mat %>% 
          filter(P2 %in% c(namei, namej) & P3 %in% c(namei, namej))
        # Get most significant D-stat
        dstat = tmp[which.min(tmp$p.value),]$Dstatistic
        
        if(length(dstat) == 0){
          pairwise.mat[i,j] = 0
        } else{
          pairwise.mat[i,j] =  dstat 
        }
      } # ELSE
    } # J
    print(i)
  } # I
  
  if(section == "white"){
    # Adding a new row with the name "Quercus_rubra"
    new_row <- matrix(0, nrow=1, ncol=74)  # New row with 74 columns initialized with 0
    rownames(new_row) <- "Quercus_rubra"
    
    # Adding a new column with the name "Quercus_rubra"
    new_col <- matrix(0, nrow=75, ncol=1)  # New column with 75 rows (including the new row) initialized with 0
    colnames(new_col) <- "Quercus_rubra"
    
    new.mat = rbind(pairwise.mat, new_row)
    new.mat = cbind(new.mat, new_col)
    
  } else { # If red oaks 
    # Adding a new row with the name "Quercus_lobata"
    new_row <- matrix(0, nrow=1, ncol=ncol(pairwise.mat))  # New row with 74 columns initialized with 0
    rownames(new_row) <- "Quercus_lobata"
    
    # Adding a new column with the name "Quercus_lobata"
    new_col <- matrix(0, nrow=nrow(pairwise.mat)+1, ncol=1)  # New column with 75 rows (including the new row) initialized with 0
    colnames(new_col) <- "Quercus_lobata"
    
    new.mat = rbind(pairwise.mat, new_row)
    new.mat = cbind(new.mat, new_col)
  }
  
  # Organize output by tree order
  tree = ladderize(tree)
  new.mat = new.mat[tree$tip.label,tree$tip.label]
  return(new.mat)
  
  
} # END OF FUCNTION

# Function to label eac`h pair as being in the same or different group
same_dif_dstat = function(dstats, states){
  for(i in 1:nrow(dstats)){
    
    test = states[states$species %in% dstats[i,2:3],]$range |>
      unique() %>% 
      length() %>% 
      `==`(1)
    
    if(test){
      dstats$region[i] = "same"
    }
    else {
      dstats$region[i] = "dif"
    }
    # print(i)
  }
  return(dstats)
}





