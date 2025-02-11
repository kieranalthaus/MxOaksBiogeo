# Install packages
library(ade4)
library(RANN)
library(ecospat)
library(tidyverse)
library(factoextra)
library(phytools)
library(geoscale)
library(RPANDA)
library(geiger)
library(geometry)
`%notin%` = Negate(`%in%`)

## 1. Load Data ----------------------------------------------------------------
# Niche data
eco.data = read_csv(file = "MxOaks_NicheEvolution/data/20240429_oak_envdata_complete.csv")
# Phylogenetic Tree
tree = ape::read.tree("MxOaks_Prep/data/20240325_p04_finaltree.tre")
# Reformat tip names
tree$tip.label = tree$tip.label %>%
  strsplit(split = "\\|") %>%
  lapply(FUN = function(x)
    x[1]) %>%
  unlist()

## Reformat data
# Remove species with fewer than 5 occurrences
remove = eco.data %>% 
  group_by(species) %>% 
  tally() %>% 
  filter(n < 5) %>% 
  pull(species)
# Filter dataset
eco.data = eco.data %>% 
  filter(species %notin% remove)
# Get species names
sp = unique(eco.data$species)[order(unique(eco.data$species))]

## 2. Run PCA ------------------------------------------------------------------
## First, PCA of entire niche space
ecopca <- dudi.pca(eco.data[,1:30], 
                   center = TRUE,
                   scale = TRUE, 
                   scannf = FALSE,
                   nf = 3)
scores.clim <- ecopca$li # Store scores

## Visualize variable contributions
fig1 = fviz_contrib(X = ecopca,
             choice = "var",
             axes = 1,
             top = 10)
fig2 = fviz_contrib(X = ecopca,
             choice = "var",
             axes = 2,
             top = 10)
# Save plots
# ggsave(plot = cowplot::plot_grid(fig1, fig2),
#        filename = "20240517_pca_contributions.pdf")

## Visualize variable loadings on biplot
figures1 = fviz_pca_var(X = ecopca,
             col.var = "black",
             geom = c("arrow", "text"),
             col.circle = "transparent")
# Save plot
# ggsave("Figures/20240517_pca_biplot.pdf",
#        figures1)

## 3. Plot PCA -----------------------------------------------------------------
# define clades
crown.clades = c(
  'Lobatae' = 'agrifolia|emoryi',
  'Quercus' = 'lobata|arizonica',
  'Virentes' = 'fusiformis|geminata',
  'Protobalanus' = 'cedrosensis|chrysolepis'
) 

# get node numbers for each clade
crownNodes.clades = sapply(X = crown.clades, 
                           FUN = function(x)
                             findMRCA(tree, grep(x, tree$tip.label, value = T)))

# filter phylogenies for each section
section.tips = phangorn::Descendants(x = tree, 
                                     node = crownNodes.clades,
                                     type = "tip")
# protobalanus phylogeny
gold.remove = tree$tip.label[unlist(section.tips[1:3])]
gold.tree = drop.tip(phy = tree, tip = c(gold.remove, "Quercus_sadleriana"))

# lobatae phylogeny
red.remove = c(tree$tip.label[unlist(section.tips[c(2:4)])], "Quercus_sadleriana")
red.tree = drop.tip(phy = tree, tip = red.remove)

# quercus Phylogeny
white.remove = c(tree$tip.label[unlist(section.tips[c(1,3,4)])], "Quercus_sadleriana")
white.tree = drop.tip(phy = tree, tip = white.remove)

# virentes Phylogeny
live.remove = c(tree$tip.label[unlist(section.tips[c(1,2,4)])])
live.tree = drop.tip(phy = tree, tip = live.remove)

## Plot PCA
xlim = round(range(eco[,2]))
ylim = round(range(eco[,3]))

{ #pdf(file = "../Figures/20241104_env_pca_sections.pdf",
   #  height = 7,
    # width = 7)
  par(mar = c(1,3,0,1),
      mfrow = c(2,2),
      oma = c(2, 0, 0.2, 0.2))
  
  plot(x = eco[,2],
       y = eco[,3],
       col = alpha("grey", 0.7),
       pch = 20,
       cex = 3,
       bty = "n",
       axes = FALSE,
       xlab = "",
       ylab = "")
  phytools::phylomorphospace(tree = live.tree, 
                             X = eco[live.tree$tip.label,2:3], 
                             xlim = xlim,
                             ylim = ylim,
                             label = "off", 
                             axes = FALSE, 
                             xlab = "PC1", 
                             ylab = "PC2", 
                             pch = 20, 
                             lwd = 1,
                             node.size = 0,
                             bty = "n",
                             add = TRUE)
  points(x = eco[live.tree$tip.label,2],
         y = eco[live.tree$tip.label,3], 
         pch = 20,
         cex = 3,
         col = "forestgreen")
  points(x = eco["Quercus_sadleriana", 2],
         y = eco["Quercus_sadleriana", 3],
         pch = 20,
         cex = 3,
         col = "steelblue")
  abline(v = 0, 
         h = 0,
         lty = 2, 
         lwd = 1)
  axis(side = 2, 
       lwd = 2,
       las = 2)
  box(lwd = 2)
  
  plot(x = eco[,2],
       y = eco[,3],
       col = alpha("grey", 0.7),
       pch = 20,
       cex = 3,
       bty = "n",
       axes = FALSE,
       xlab = "",
       ylab = "")
  phytools::phylomorphospace(tree = white.tree, 
                             X = eco[white.tree$tip.label,2:3], 
                             xlim = xlim,
                             ylim = ylim,
                             label = "off", 
                             axes = FALSE, 
                             xlab = "PC1", 
                             ylab = "PC2", 
                             pch = 20, 
                             lwd = 1,
                             node.size = 0,
                             bty = "n",
                             add = TRUE)
  points(x = eco[white.tree$tip.label,2],
         y = eco[white.tree$tip.label,3], 
         pch = 20,
         cex = 3,
         col = "dodgerblue")
  abline(v = 0, 
         h = 0,
         lty = 2, 
         lwd = 1)
  box(lwd = 2)
  
  plot(x = eco[,2],
       y = eco[,3],
       col = alpha("grey", 0.7),
       pch = 20,
       cex = 3,
       bty = "n",
       axes = FALSE,
       xlab = "",
       ylab = "")
  phytools::phylomorphospace(tree = red.tree, 
                             X = eco[red.tree$tip.label,2:3], 
                             xlim = xlim,
                             ylim = ylim,
                             label = "off", 
                             axes = FALSE, 
                             xlab = "PC1", 
                             ylab = "PC2", 
                             pch = 20, 
                             lwd = 1,
                             node.size = 0,
                             bty = "n",
                             add = TRUE)
  points(x = eco[red.tree$tip.label,2],
         y = eco[red.tree$tip.label,3], 
         pch = 20,
         cex = 3,
         col = "red4")
  abline(v = 0, 
         h = 0,
         lty = 2, 
         lwd = 1)
  axis(1, lwd = 2)
  axis(2, lwd = 2, las = 2)
  box(lwd = 2)
  
  plot(x = eco[,2],
       y = eco[,3],
       col = alpha("grey", 0.7),
       pch = 20,
       cex = 3,
       bty = "n",
       axes = FALSE,
       xlab = "",
       ylab = "")
  phytools::phylomorphospace(tree = gold.tree, 
                             X = eco[gold.tree$tip.label,2:3], 
                             xlim = xlim,
                             ylim = ylim,
                             label = "off", 
                             axes = FALSE, 
                             xlab = "PC1", 
                             ylab = "PC2", 
                             pch = 20, 
                             lwd = 1,
                             node.size = 0,
                             bty = "n",
                             add = TRUE)
  points(x = eco[gold.tree$tip.label,2],
         y = eco[gold.tree$tip.label,3], 
         pch = 20,
         cex = 3,
         col = "gold")
  abline(v = 0, 
         h = 0,
         lty = 2, 
         lwd = 1)
  axis(side = 1, 
       lwd = 2)
  legend(x = -3,
         y = 13,
         legend = c(
           expression(paste("sect. ", italic("Virentes"))),
           expression(paste("sect. ", italic("Ponticae"))),
           expression(paste("sect. ", italic("Quercus"))),
           expression(paste("sect. ", italic("Lobatae"))),
           expression(paste("sect. ", italic("Protobalanus")))
         ),
         col = c("forestgreen",
                 "steelblue",
                 "dodgerblue",
                 "red4",
                 "gold"),
         bty = "o",
         box.lwd = 2,
         bg = "white",
         box.col = "transparent",
         pch = 20,
         pt.cex = 2
  )
  box(lwd = 2)
  dev.off()
}


## 4. Niche Breadth ------------------------------------------------------------
# Create output dataframe
out = data.frame(species = unique(eco.data$species),
                 PC1 = NA,
                 PC2 = NA,
                 PC1_b = NA,
                 PC2_b = NA,
                 breadth = NA)

# Calculate total niche breadth for each species
i = 1
for(i in 1:nrow(out)){
  # Extract climate data
  sp = eco.data %>% 
    dplyr::filter(species == out$species[i]) %>% 
    dplyr::select(-species, -latitude, -longitude, -wc2.1_30s_elev)
  
  # Get scores
  scores <- suprow(ecopca, sp)$li
  
  # Get PC 1 & 2 Axes
  out[i,c("PC1", "PC2")] = colMeans(scores[,1:3]) # For first three axes
  
  # Calculate niche breadth
  out[i,"breadth"] = prod(apply(scores[,1:2],2, var))
  out[i,"PC1_b"] = var(scores[,1])
  out[i,"PC2_b"] = var(scores[,2])

}

## Calculate elevation niche breadth for each species
elev.df = eco.data %>% 
  group_by(species) %>% 
  # Remove outlier elevation points
  filter(!wc2.1_30s_elev %in% boxplot.stats(wc2.1_30s_elev)$out) %>% 
  # Summarize elevational data
  summarise(
    elev.mean = mean(wc2.1_30s_elev),
    elev.breadth = var(wc2.1_30s_elev),
    elev.lower.25 = quantile(wc2.1_30s_elev)[2],
    elev.upper.25 = quantile(wc2.1_30s_elev)[4],
    lat.mean = mean(latitude))

## Merge absolute niche and elevational niche data
out.df = merge(x = out, y = elev.df)
# Edit species column
out.df$species = gsub(x=out$species, " ","_")
# Add rownmaes
rownames(out.df) = out.df$species

## Save niche data
# write.csv(x = out.df,
#           file = "20240816_oak_niche_pca.csv")
niche.data = out.df

## 5. Niche PHYLOLM ------------------------------------------------------------
# Read niche data
niche.data = read.csv(file = "niche_evolution/20240816_oak_niche_pca.csv")

# Remove non-overlapping taxa in data and phylogeny
names2remove = name.check(phy = tree, data = niche.data)
tree = drop.tip(phy = tree, tip = names2remove$tree_not_data)
niche.data = niche.data[!rownames(niche.data) %in% obj$data_not_tree,]

## Phylogenetic regression relating niche breadth to elevation and latitude
model = phylolm::phylolm(
  formula = log(breadth) ~ elev.mean + elev.mean*lat.mean, 
  data = niche.data,
  phy = tree,
  model = "lambda")

# Model summary
summary(model)


## 6. Niche Pagel's Lambda  ----------------------------------------------------
# Read niche data
niche.data = read.csv(file = "MxOaks_NicheEvolution/data/20240816_oak_niche_pca.csv")

# Calculate phylogenetic signal for pc1
lambda.pc1 = phylosig(tree = tree, 
                      x = setNames(niche.data$PC1, niche.data$species),
                      method = "lambda", 
                      test = TRUE,
                      niter = 100)

# Create 95% CI for PC1
lambda.seq = seq(from = 0, to = 1, by = 0.001)
pc1.llik = unlist(lapply(X = lambda.seq,
       FUN = lambda.pc1$lik))
pc1.llik = setNames(pc1.llik, lambda.seq)
pc1.llik[which(abs(pc1.llik - (lambda.pc1$logL-2)) == min(abs(pc1.llik - (lambda.pc1$logL-2))))]

# Calculate phylogenetic signal for pc2
lambda.pc2 = phylosig(tree = tree, 
                      x = setNames(niche.data$PC2, niche.data$species),
                      method = "lambda",
                      test = TRUE,
                      niter = 100)
pc2.llik = unlist(lapply(X = lambda.seq,
                         FUN = lambda.pc2$lik))
pc2.llik = setNames(pc2.llik, lambda.seq)
pc2.llik[which(abs(pc2.llik - (lambda.pc2$logL-2)) == min(abs(pc2.llik - (lambda.pc2$logL-2))))]

## 7. Mantel Test: Niche Overlap vs. Phylogenetic Distance  --------------------
# Create output matrix
out = matrix(data = NA,
             nrow = length(unique(eco.data$species)),
             ncol = length(unique(eco.data$species)),
             dimnames = list(
               c(unique(eco.data$species)),
               c(unique(eco.data$species))
             ))

i = 1
for(i in 1:nrow(out)){ # For each row...
  # Extract climate data
  spi = eco.data %>% 
    dplyr::filter(species == colnames(out)[i]) %>% 
    dplyr::select(-species, -latitude, -longitude, -wc2.1_30s_elev)
  
  # Get scores
  scoresi <- suprow(ecopca, spi)$li
  
  # Calculate hulls
  hulli = convhulln(scoresi)
  
  for(j in 1:nrow(out)){
    spj = eco.data %>% 
      dplyr::filter(species == colnames(out)[j]) %>% 
      dplyr::select(-species, -latitude, -longitude, -wc2.1_30s_elev)
    
    # Get scores
    scoresj = suprow(ecopca, spj)$li
    
    # Calculate hulls
    hullj = convhulln(scoresj)  
    
    # Overlap
    out[i,j] = tryCatch({intersectn(hulli, hullj, options = "Tv Qj")$ch$vol}, 
                        error = function(e) NA)
  }
  print(i)
}

# Change dimnames
rownames(out) = gsub(x = rownames(out), " ", "_")
colnames(out) = gsub(x = colnames(out), " ", "_")

# Subset for non-matching rows
out = out[,colnames(out) %in% tree$tip.label]
out = out[rownames(out) %in% tree$tip.label,]
tree = drop.tip(tree, tree$tip.label[tree$tip.label %notin% colnames(out)])

# Create distance matrix with phylogeny
tree.mat = cophenetic.phylo(tree)

# Run mantel test
mantel =  mantel(xdis = tree.mat,
                 ydis = out, 
                 permutations = 1000)

## 

## 8. Niche Disparity through time  --------------------------------------------
# DTT for pc1
dtt.pc1 = geiger::dtt(phy = tree,
                      data = setNames(niche.data$PC1, niche.data$species),
                      nsim = 1000, 
                      CI = 0.95, 
                      plot = FALSE,
                      calculateMDIp = TRUE)
# DTT for pc2
dtt.pc2 = geiger::dtt(phy = tree,
                      data = setNames(niche.data$PC2, niche.data$species),
                      nsim = 1000, 
                      CI = 0.95, 
                      plot = FALSE,
                      calculateMDIp = TRUE)

## 9. Models of Niche Evolution  -----------------------------------------------

rownames(niche.data) = niche.data$species
# Check names
names2remove = geiger::name.check(phy = tree, 
                                  data = niche.data)

# Filter names
model.tree = drop.tip(phy = tree,
                tip = names2remove$tree_not_data)
niche.data.sim = niche.data %>% 
  dplyr::filter(species %notin% names2remove$data_not_tree)

# Create two PC datasets
pc1 = setNames(object = niche.data.sim$PC1,
               nm = niche.data.sim$species)
pc2 = setNames(object = niche.data.sim$PC2,
               nm = niche.data.sim$species)

## Fit Models for PC1
# Brownian motion
fitBM = fitContinuous(phy = model.tree,
                      dat = pc1,
                      model = "BM")
# Single-rate OU
fitOU = fitContinuous(phy = model.tree,
                      dat = pc1,
                      model = "OU",
                      bounds = list(alpha = c(0,10)))
# Early-burst
fitEB = fitContinuous(phy = model.tree,
                      dat = pc1, 
                      model = "EB")
# Diversity dependent (linear)
fitDDlin = fit_t_comp(phylo = model.tree, 
                      data = pc1, 
                      model = "DDlin")
# Diversity dependent (exponential)
fitDDexp = fit_t_comp(phylo = model.tree, 
                      data = pc1, 
                      model = "DDexp")
# Bounded Brownian Motion
fitBBM = bounded_bm(tree = model.tree, 
                    x = pc1, 
                    lims = range(pc1), 
                    lik.func = "parallel",
                    root = "mle", 
                    levs = 200)
# Check AICw
aic.pc1 = setNames(object = c(AIC(fitBM, fitOU, fitEB,fitBBM)$AIC, fitDDlin$aic, fitDDexp$aic),
                   nm = c("BM", "OU", "EB", "BBM", "DDlin", "DDexp"))
aic.w(aic = aic.pc1[-2])

## Fit Models for PC2
# Brownian motion
fitBM = fitContinuous(phy = model.tree,
                      dat = pc2,
                      model = "BM")
# Single-rate OU
fitOU = fitContinuous(phy = model.tree,
                      dat = pc2,
                      model = "OU",
                      bounds = list(alpha = c(0,10)))
# Early-burst
fitEB = fitContinuous(phy = model.tree,
                      dat = pc2, 
                      model = "EB")
# Diversity dependent (linear)
fitDDlin = fit_t_comp(phylo = model.tree, 
                      data = pc2, 
                      model = "DDlin")
# Diversity dependent (exponential)
fitDDexp = fit_t_comp(phylo = model.tree, 
                      data = pc2, 
                      model = "DDexp")
# Bounded Brownian Motion
fitBBM = bounded_bm(tree = model.tree, 
                    x = pc2, 
                    lims = range(pc2), 
                    lik.func = "parallel",
                    root = "mle", 
                    levs = 200)
# Check AICw
aic.pc2 = setNames(object = c(AIC(fitBM, fitOU, fitEB,fitBBM)$AIC, fitDDlin$aic, fitDDexp$aic),
                   nm = c("BM", "OU", "EB", "BBM", "DDlin", "DDexp"))
aic.w(aic = aic.pc2[-2])

## 10. Simulate Disparity Through Time  -----------------------------------------
# Set parameters from best fit models
sig_bbm = 2.50
x0_bbm = -7.66
sig_exp = 4e-8
z0 = 0.63
r = 0.102

## Simulate 100 BBM models
sim_bbm = fastBM(tree = model.tree, 
                 a = x0_bbm,
                 sig2 = sig_bbm,
                 bounds = range(pc1), 
                 nsim = 100)

# Create empty output for DTT
bbm_dtt = matrix(data = NA,
                 nrow = length(tree$tip.label),
                 ncol = 100)

# Run DTT 100 times on each simulated dataset
for(i in 1:100){
  x = dtt(phy = model.tree, data = sim_bbm[,i], plot = FALSE)
  bbm_dtt[,i] = x$dtt
}

## Simulate 100 DDexp models
# Create empty output for DDexp  
out = matrix(data = NA,
             nrow = 177,
             ncol = 100,
             dimnames = list(
               c(tree$tip.label),
               c()
             ))

# Simulate 100 DDexp datasets
for(i in 1:100){
  out[,i] = sim_t_comp(phylo = tree, 
                       pars = c(
                         sig2 = sig_exp,
                         r = r
                       ),
                       root.value = z0,
                       model = "DDexp")
}

# Create empty output for DTT
ddexp_dtt = matrix(data = NA,
                   nrow = length(tree$tip.label),
                   ncol = 100)

# Run DTT 100 times on each simulated dataset
for(i in 1:100){
  x = dtt(phy = tree, data = out[,i], plot = FALSE)
  ddexp_dtt[,i] = x$dtt
}

## 11. Plot DTT & Biogeographic Timing -----------------------------------------
pdf(file = "../Figures/20250102_ltt_dtt_PCA_stacked.pdf",
    height = 4.85,
    width = 6.59
)
par(mfrow = c(2,2),
    oma = c(1,4,0,0),
    mar = c(0,0,1,2))
plot(x = dtt.pc1$times,
     y = dtt.pc1$dtt,
     type = 'l',
     axes = FALSE,
     ylim = c(-0.5,2),
     lwd = 0,
     xlab = "",
     ylab = "",
     cex.lab = 1
)
matplot(x = dtt.pc1$times,
        y = bbm_dtt,
        type = 'l',
        col = geiger:::.transparency("lightgray", 0.5),
        lty = 1,
        add = TRUE,
        axes = F)
segments(x0 = (th-t1[2:7, 'Start'][c(-6,-5)])/h, 
         y0 = -0.09, 
         y1 = length(tree$tip.label), 
         lty = 'dashed', 
         col = 'gray',
         lwd = 2)
rect(
  xleft = 1,
  xright = rev(c(50.99321, 48.24821, 30.58121, 19.68121, 0)/h),
  ytop = 0, 
  ybottom = -0.09,
  col = c("white", "gray40"),
  border = "black",
  lwd = 2)
lines(x = dtt.pc2$times,
      y = dtt.pc2$dtt,
      type = 'l',
      lwd = 2)
poly = geiger:::.dtt.polygon(dtt.pc1$sim, dtt.pc1$times, alpha = 1 - 0.85)
polygon(poly[, "x"], poly[, "y"], col = geiger:::.transparency("lightgray", 0.5), border = NA)
axis(2, 
     at = c(0,2,1),
     pos = 0, 
     las = 2,
     lwd = 2, 
     cex.axis = 2)

# For PC2
plot(x = dtt.pc2$times,
     y = dtt.pc2$dtt,
     type = 'l',
     axes = FALSE,
     ylim = c(-0.5,2),
     lwd = 0,
     xlab = "",
     ylab = "",
     cex.lab = 2
)
matplot(x = dtt.pc1$times,
        y = ddexp_dtt,
        type = 'l',
        col = geiger:::.transparency("lightgray", 0.5),
        lty = 1,
        add = TRUE,
        axes = F)
segments(x0 = (th-t1[2:7, 'Start'][c(-6,-5)])/h, 
         y0 = -0.09, 
         y1 = length(tree$tip.label), 
         lty = 'dashed', 
         col = 'gray',
         lwd = 2)
rect(
  xleft = 1,
  xright = rev(c(50.99321, 48.24821, 30.58121, 19.68121, 0)/h),
  ytop = 0, 
  ybottom = -0.09,
  col = c("white", "gray40"),
  border = "black",
  lwd = 2)
lines(x = dtt.pc1$times,
      y = dtt.pc1$dtt,
      type = 'l',
      lwd = 2)
poly = geiger:::.dtt.polygon(dtt.pc2$sim, dtt.pc2$times, alpha = 1 - 0.85)
polygon(poly[, "x"], poly[, "y"], col = geiger:::.transparency("lightgray", 0.5), border = NA)
axis(2, 
     at = c(0,2,1),
     pos = 0, 
     las = 2,
     lwd = 2, 
     cex.axis = 2)

biogeo_timing = matrix(
  data = c(
    "#E16462FF", 1.90, "SMOR", 13.57, 18.90, 19.68, 26.98,
    "#B12A90FF", 1.30, "TMVB", 11.38, 18.90, 17.37, 23.92,
    "#F0F921FF", 0.70, "SMOC", 7.47, 10.72, 10.52, 15.07,
    "#0D0887FF", 0.10, "TROP", 7.50, 11.00, 6.70, 11.89
  ),
  byrow = TRUE,
  ncol = 7,
  dimnames = list(
    c(),
    c("col", "y", "region", "red_min", "red_max", "white_min", "white_max")
  )
)

pch = 15
cex = 2
plot(x = dtt.pc1$times,
     y = dtt.pc1$dtt,
     type = 'l',
     axes = FALSE,
     ylim = c(-0.5,2),
     lwd = 0,
     xlab = "",
     ylab = "",
     cex.lab = 1,
     col = "white"
)
w <- cex*strwidth("m")/2
segments(x0 = (th-t1[2:7, 'Start'][c(-6,-5)])/h, 
         y0 = -0.09, 
         y1 = length(tree$tip.label), 
         lty = 'dashed', 
         col = 'gray',
         lwd = 2)
rect(
  xleft = 1,
  xright = rev(c(50.99321, 48.24821, 30.58121, 19.68121, 0)/h),
  ytop = 0, 
  ybottom = -0.09,
  col = c("white", "gray40"),
  border = "black",
  lwd = 2)
segments(
  x0 = smor[1]/h,
  y0 = as.numeric(biogeo_timing[1,2]),
  x1 = smor[2]/h,
  y1 = as.numeric(biogeo_timing[1,2]),
  lwd = 5,
  col = smor.col
)
segments(
  x0 = tmvb[1]/h,
  y0 = as.numeric(biogeo_timing[2,2]),
  x1 = tmvb[2]/h,
  y1 = as.numeric(biogeo_timing[2,2]),
  lwd = 5,
  col = tmvb.col
)
segments(
  x0 = smoc[1]/h,
  y0 = as.numeric(biogeo_timing[3,2]),
  x1 = smoc[2]/h,
  y1 = as.numeric(biogeo_timing[3,2]),
  lwd = 5,
  col = smoc.col
)
arrows(
  x0 = (th - as.numeric(biogeo_timing[,"white_min"]))/h,
  x1 = (th - as.numeric(biogeo_timing[,"white_max"]))/h,
  y0 = as.numeric(biogeo_timing[,"y"]),
  y1 = as.numeric(biogeo_timing[,"y"]),
  angle = 90, 
  length = 0.05,
  col = 'dodgerblue',
  code = 3,
  lwd = 2
)
arrows(
  x0 = (th - as.numeric(biogeo_timing[,"red_min"]))/h,
  x1 = (th - as.numeric(biogeo_timing[,"red_max"]))/h,
  y0 = as.numeric(biogeo_timing[,"y"]),
  y1 = as.numeric(biogeo_timing[,"y"]),
  angle = 90, 
  length = 0.05,
  col = 'red4',
  code = 3,
  lwd = 2
)
fossils = matrix(
  data = c(
    th-15, 1.3, "Puebla",
    th-5.3, 0.1, "Panama / Guatemala"
  ),
  byrow = TRUE,
  ncol = 3,
  dimnames = list(
    c(),
    c("x","y","region")
  )
)
points(
  x = as.numeric(fossils[,"x"])/h,
  y = as.numeric(fossils[,"y"]),
  pch = 4,
  cex = 1,
  lwd = 3
)
axis(
  side = 4,
  at = biogeo_timing[,2],
  labels = c("","","",""),
  lwd = 3
)
dev.off()
