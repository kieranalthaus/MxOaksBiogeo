library(diversitree)
`%notin%` = Negate(`%in%`)

### DEFINE FUNCTIONS FOR QUASSE ###
library(hisse)
library(diversitree)

#' @name make.oaks a helper function for creating QuaSSE models
#' @param lambda the function describing the relationship between speciation and 
#' the quantitative trait
#' @param mu the function describing the relationship between extinction and the 
#' quantitative trait
make.oaks = function(lambda, mu) {
  make.quasse(tree = tree, 
              states = states,
              states.sd = states.sd,
              lambda = lambda,
              mu = mu)}

#' @name nodrift a helper function for constraining QuaSSE models
#' @param f likelihood function describing the model
nodrift = function(f){
  constrain(f = f, formulae = drift ~ 0)
}

#' @name QuaSSE_Model carries out a set of 32 QuaSSE models with different 
#' speciation and extinction functions, and with and without the drift parameter
#' @param tree the phylogenetic tree
#' @param states a named numeric vector of the quantitative trait
#' @param outdir the path of the output directory
#' @param outname the name suffix of each set of saved models

QuaSSE_Model = function(tree, states, outdir, outname) {
  
  ### SET UP MODEL PARAMETERS---------------------------------------------------
  
  # Standard deviation of states
  states.sd = sd(states)
  
  # Create Starting point
  p = starting.point.quasse(tree = tree, 
                            states = states, 
                            states.sd = states.sd)
  
  # Define ranges
  xr = range(states) + c(-1,1) * 20 * p["diffusion"]
  
  # Define lienar function
  linear.x = make.linear.x(x0 = xr[1], # Lower limit for linear function
                           x1 = xr[2]) # Upper limit for linear function
  
  control = list(parscale = 0.1, reltol = 0.001)
  
  ### CONSTANT EXTINCTION MODELS -----------------------------------------------
  # Set up models with constant extinction rates, and varying speciation rates
  # both with and without drift parameters
  
  # Constant speciation
  f.c = make.oaks(lambda = constant.x, # NULL MODEL
                  mu = constant.x)
  
  # Linear speciation
  f.l = make.oaks(lambda = linear.x,
                  mu = constant.x)
  
  # Sigmoidal speciation
  f.s = make.oaks(lambda = sigmoid.x,
                  mu = constant.x)
  
  # Modal speciation
  f.m = make.oaks(lambda = noroptimal.x,
                  mu = constant.x)
  
  ## BEING RUNNING MODELS 
  print("Beginning constant models, no drift")
  mle.c = find.mle(func = nodrift(f.c),
                   x.init = p,
                   lower = 0,
                   verbose = 0,
                   control = control
  )
  
  # Set starting values for following models
  p.c = mle.c$par
  p.l = c(p.c[1], l.m = 0, p.c[2:3])
  p.s = p.m = c(p.c[1], p.c[1], mean(xr), 1, p.c[2:3])
  names(p.s) = argnames(nodrift(f.s))
  names(p.m) = argnames(nodrift(f.m))
  
  mle.l = find.mle(func = nodrift(f.l),
                   x.init = p.l,
                   lower = 0,
                   verbose = 0,
                   control = control
  )
  mle.s = find.mle(func = nodrift(f.s),
                   x.init = p.s,
                   lower = 0,
                   verbose = 0,
                   control = control
  )
  # p.h[1] = p.h[1]+1
  mle.m = find.mle(func = nodrift(f.m),
                   x.init = p.m,
                   lower = 0,
                   verbose = 0,
                   control = control 
  )
  
  print("Beginning constant models, with drift")
  
  # Drift models
  mle.d.l = find.mle(func = f.l, x.init = coef(mle.l, TRUE), control = control, verbose = 0)
  mle.d.s = find.mle(func = f.s, x.init = coef(mle.s, TRUE), control = control, verbose = 0)
  mle.d.m = find.mle(func = f.m, x.init = coef(mle.m, TRUE), control = control, verbose = 0)
  
  print("Saving constant models")
  models = list(constant = mle.c, linear = mle.l, sig = mle.s, modal = mle.m, lineard = mle.d.l, sigd = mle.d.s, modald = mle.d.m)
  saveRDS(object = models, file = paste0(outdir, outname, ".RData"))
  
}

######## TERRAIN RUGEDNESS INDEX MODELS ########
## PART 1: LOAD DATA -----------------------------------------------------------
# Load quantitative trait data
terrain_data = read_csv(file = "MxOaks_QuaSSE/data/QuaSSE_terrain.csv")
states = setNames(object = log(terrain_data$tri_mean), nm = terrain_data$species)
ststaes.sd = sd(states)

# Load ultrametric phylogeny
tree = read.tree('MxOaks_DEC/data/DEC_tree.tre')

## Run models
# QuaSSE_Models(tree = tree,
#               states = states, 
#               outdir = "OUT/20240326_out/",
#               outname = "20240326.models")

## Run simple models with fixed extinction rates
QuaSSE_Model(tree = tree, 
             states = states,
             outdir = "MxOaks_QuaSSE/out/",
             outname = "20240326_quasse_mods")

models = readRDS(file = "MxOaks_QuaSSE/out/20240326_quasse_TRI_mods.RData")

anova(models$constant,
      lin = models$linear,
      sig = models$sig,
      mod = models$modal,
      lin.d = models$lineard, # <----- The best fit model for TRI models 
      sig.d = models$sigd, 
      mod.d = models$modald)


######## ELEVATION MODELS ########
## PART 1: LOAD DATA -----------------------------------------------------------
# Load quantitative trait data
terrain_data = read_csv(file = "MxOaks_QuaSSE/out/QuaSSE_terrain.csv")
states = setNames(object = log(terrain_data$elev_mean), nm = terrain_data$species)
ststaes.sd = sd(states)

# Load ultrametric phylogeny
tree = read.tree('MxOaks_QuaSSE/data/DEC_tree.tre')

## Run models
# QuaSSE_Models(tree = tree, states = states, outdir = "OUT/20240326_out/",
#               outname = "20240326.models")

## Run simple models with fixed extinction rates
QuaSSE_Model(tree = tree, 
             states = states,
             outdir = "MxOaks_QuaSSE/out/",
             outname = "20240326_quasse_elev_mods")

models = readRDS(file = "MxOaks_QuaSSE/out/20240326_quasse_elev_mods.RData")

anova(models$constant,
      lin = models$linear,
      sig = models$sig,
      mod = models$modal,
      lin.d = models$lineard,
      sig.d = models$sigd, 
      mod.d = models$modald)



