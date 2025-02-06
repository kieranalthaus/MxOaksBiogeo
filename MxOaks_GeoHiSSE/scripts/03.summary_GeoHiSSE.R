# packages
library(hisse)
library(diversitree)

# Load Results
model.list = list.files(path = "MxOaks_GeoHiSSE/out/", pattern = "20240428_mntn_mods",full.names = T)
models = lapply(model.list, readRDS)

# Which model has the lowest AIC
best_model = lapply(models, function(x) x$AIC) %>% 
  unlist() %>% which.min()

models[[best_model]]

# Marginal reconstruction of model history
# recon = list()
# for(i in 1:length(models)){
#   
#   recon[[i]] = MarginReconGeoSSE(phy = models[[i]]$phy,
#                                  data = models[[i]]$data,
#                                  f = c(0.18, .67, 0.15), 
#                                  pars = models[[i]]$solution,
#                                  hidden.states = ncol(models[[i]]$trans.matrix)/3,
#                                  root.type = models[[i]]$root.type,
#                                  root.p = models[[i]]$root.p,
#                                  AIC = models[[i]]$AIC, 
#                                  n.cores = 7)
#   print(i)
#   
# }
# saveRDS(object = recon, file = "MX_oaks_GEOSSE/OUT/20240428_model_reconstruction.RData")
recon = readRDS("MxOaks_GeoHiSSE/out/20240428_model_reconstruction.RData")

# Store all reconstructions into a single list object
recon.models = list(mod1 = recon[[1]],
                    mod2 = recon[[2]],
                    mod3 = recon[[3]],
                    mod4 = recon[[4]],
                    mod5 = recon[[5]],
                    mod6 = recon[[6]],
                    mod7 = recon[[7]],
                    mod8 = recon[[8]],
                    mod9 = recon[[9]],
                    # mod10 = recon[[10]], # This model is a repeat of model 8 
                    mod11 = recon[[11]],
                    mod12 = recon[[12]],
                    mod13 = recon[[13]],
                    mod14 = recon[[14]]
)

# Get model AIC weights
aic.weights = GetAICWeights(hisse.results = recon.models, criterion = "AIC")

# Get model average rates of diversification 
model.ave.rates = GetModelAveRates(x = recon.models, type = "both")


## SUMMARIZE FINDINGS
# Explore AIC weights
aic.weights[order(aic.weights, decreasing = T)]

# Summarize tip rates
geosse_tiprate = apply(model.ave.rates$tips[,2:4], MARGIN = 1, which.max)
geosse_noderate = apply(model.ave.rates$nodes[,2:4], MARGIN = 1, which.max)

# Geometric mean of tip rates
by(data = model.ave.rates$tips$speciation, 
   INDICES = geosse_tiprate,
   FUN = function(x)
     exp(mean(log(x))))
by(data = model.ave.rates$tips$extinction, 
   INDICES = geosse_tiprate,
   FUN = function(x)
     exp(mean(log(x))))
by(data = model.ave.rates$tips$net.div,
   INDICES = geosse_tiprate,
   FUN = function(x)
     exp(mean(log(x))))

by(data = model.ave.rates$tips$net.div,
   INDICES = geosse_tiprate,
   FUN = quantile,
   probs = c(0.05,0.95))

# Geometric mean of all rates
# names(model.ave.rates$nodes)[1] = names(model.ave.rates$tips)[1]
# model.rates = do.call(rbind, model.ave.rates)
# groups = apply(model.rates[,2:4], MARGIN = 1, which.max)
# 
# by(data = model.rates$net.div,
#    INDICES = groups,
#    FUN = function(x)
#      exp(mean(log(x))))
# 
# by(data = model.rates$speciation,
#    INDICES = groups,
#    FUN = quantile,
#    probs = c(0.05,0.95))


## PLOT RESULTS ##
## Boxplots for tip rates
# {pdf(file = "MxOaks_GeoHiSSE/out/0240502_geohisse_rates.pdf",
#      height = 5.09,
#      width = 8.72)
par(mar = c(5, 5, 2, 2),
    mfrow = c(1,2))
# Plot speciation rates
plot(x = as.factor(geosse_tiprate), 
     y = model.ave.rates$tips$speciation, 
     col = "white",
     ylab = expression(paste("Speciation Rate (", lambda, ")")),
     xlab = "Region",
     axes = F,
     cex.lab = 1.5)
axis(side = 1, 
     at = c(1,2,3),
     label = c("Lowland", "Montane", "Widespread"), cex.axis = 1)
axis(side = 2,
     at = c(0, 0.2, 0.4, 0.6),
     las = 2)
# Plot extinction rates
# plot(x = as.factor(x), 
#      y = model.ave.rates$tips$extinction,
#      col = "white",
#      xlab = "Region",
#      ylab = expression(paste("Extinction Rate (", mu, ")")),
#      axes = F,
#      ylim = c(0,0.03),
#      cex.lab = 1.5)
# axis(side = 1, 
#      at = c(1,2,3),
#      label = c("Lowland", "Mountain", "Both"))
# axis(side = 2,
#      at = c(0,0.01,0.02, 0.03),
#      labels = c(0,0.01,0.02, 0.03))

# Plot net diversification rates
plot(x = as.factor(geosse_tiprate), 
     y = model.ave.rates$tips$net.div,
     col = "white",
     xlab = "Region",
     ylab = expression(paste("Net Diversification (", lambda, "+", mu, ")")),
     axes = FALSE,
     ylim = c(0, 0.8),
     cex.lab = 1.5)
axis(side = 1,
     at = c(1,2,3),
     label = c("Lowland", "Montane", "Widespread"))
axis(side = 2,
     at = c(0, 0.2, 0.4, 0.6, 0.8),
     las = 2)
dev.off()
# }

