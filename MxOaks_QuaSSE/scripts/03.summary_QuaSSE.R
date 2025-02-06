# packages
library(hisse)
library(diversitree)

# read in model results
elev = readRDS(file = "MxOaks_QuaSSE/out/20240326_quasse_elev_mods.RData") # elevation 
tri = readRDS(file = "MxOaks_QuaSSE/out/20240326_quasse_TRI_mods.RData") # TRI

## ANOVA -- ELEVATION
anova(elev$constant,
      linear = elev$linear,
      sig = elev$sig,
      modal = elev$modal,
      linear_d = elev$lineard,
      sig_d = elev$sigd,
      modal_d = elev$modald) # <- Lowest AIC & p-value

# Extract coefficients
ecc = coefficients(elev$constant)
ecld = coefficients(elev$lineard)
ecsd = coefficients(elev$sigd)
ecmd = coefficients(elev$modald)

# Plot results for 3 drift models and constant
x = seq(from = 2, 
        to = 10, 
        by = 0.01)

# Set linear function
linear.x = make.linear.x(x0 = 2, x1 = 10)

pdf(file = "Figures/20240507_quasse_tri_elev.pdf", height = 4, width = 8)
par(mar = c(5, 5, 2, 2),
    mfrow = c(1,2))
plot(x = x, y = constant.x(x, ecc[1]), type = "l", xlim = c(2,10), 
     ylim = c(0,0.15), 
     col = "grey40",
     axes = F,
     xlab = "log(Elevation)",
     ylab = expression(paste("Speciation Rate (", lambda, ")")),
     cex.lab = 1.4)
lines(x = x, y = linear.x(x = x, c = ecld[1], m = ecld[2]), col = "grey40")
lines(x = x, y = sigmoid.x(x = x, y0 = ecsd[1], y1 = ecsd[2], xmid = ecsd[3], r = ecsd[4]), col = "grey40")
lines(x = x, y = noroptimal.x(x = x, y0 = ecmd[1], y1 = ecmd[2], xmid = ecmd[3], s2 = ecmd[4]), col = "black", lwd = 3)
axis(1)
axis(2, las = 2)


## ANOVA -- TRI
anova(tri$constant,
      linear = tri$linear,
      sig = tri$sig,
      modal = tri$modal,
      linear_d = tri$lineard, # <- Lowest AIC & p-value
      sig_d = tri$sigd,
      modal_d = tri$modald)

# Extract coefficients
tcc = coefficients(tri$constant)
tcld = coefficients(tri$lineard)
tcsd = coefficients(tri$sigd)
tcmd = coefficients(tri$modald)

# Plot results for 3 drift models and constant
tri = seq(from = -2, 
        to = 10, 
        by = 0.01)

# Set linear function
linear.x = make.linear.x(x0 = -2, x1 = 10)

# Plot
# par(mar = c(5, 5, 2, 2))
plot(x = x, y = constant.x(x, tcc[1]), type = "l", xlim = c(2,8),
     ylim = c(0,0.4), 
     col = "grey40",
     axes = F,
     xlab = "log(Terrain Rugedness)",
     ylab = expression(paste("Speciation Rate (", lambda, ")")),
     cex.lab = 1.4)
lines(x = x, y = linear.x(x = x, c = tcld[1], m = tcld[2]), col = "black", lwd = 3)
lines(x = x, y = sigmoid.x(x = x, y0 = tcsd[1], y1 = tcsd[2], xmid = tcsd[3], r = tcsd[4]), col = "grey40")
lines(x = x, y = noroptimal.x(x = x, y0 = tcmd[1], y1 = tcmd[2], xmid = tcmd[3], s2 = tcmd[4]), col = "grey40", lwd = 1)
axis(1)
axis(2, las = 2)
dev.off()









height = 6.54
width = 3.33

pdf(file = "../Figures/20241004_quasse_geohisse_v2.pdf",
    height = height, 
    width = width)
#     height = 7.28,
#     width = 4.17,
# )
par(mar = c(5, 5, 2, 2),
    mfrow = c(2,1))
plot(x = tri, y = constant.x(tri, tcc[1]), type = "l", xlim = c(2,8),
     ylim = c(0,0.4), 
     col = "grey40",
     axes = F,
     xlab = "log(Terrain Rugedness)",
     ylab = expression(paste("Speciation Rate (", lambda, ")")),
     cex.lab = 1.4)
lines(x = tri, y = linear.x(x = tri, c = tcld[1], m = tcld[2]), col = "black", lwd = 3)
lines(x = tri, y = sigmoid.x(x = tri, y0 = tcsd[1], y1 = tcsd[2], xmid = tcsd[3], r = tcsd[4]), col = "grey40")
lines(x = tri, y = noroptimal.x(x = tri, y0 = tcmd[1], y1 = tcmd[2], xmid = tcmd[3], s2 = tcmd[4]), col = "grey40", lwd = 1)
axis(1, lwd = 2)
axis(2, las = 2, lwd = 2)
box(lwd = 2)
grid()
plot(x = as.factor(geosse_tiprate),
     y = model.ave.rates$tips$net.div,
     col = "white",
     xlab = "Region",
     ylab = expression(paste("Net Diversification (", lambda, "-", mu, ")")),
     axes = FALSE,
     ylim = c(0, 0.8),
     cex.lab = 1.5)
axis(side = 1,
     at = c(1,2,3),
     lwd = 2,
     labels = FALSE)
text(x = c(1,2,3),
     y = par("usr")[3]-0.075,
     label = c("Lowland", "Montane", "Widespread"),
     srt = 22.5,
     xpd = NA,
     adj = 0.9)
axis(side = 2,
     at = c(0, 0.2, 0.4, 0.6, 0.8),
     lwd = 2,
     las = 2)
box(lwd = 2)
grid()
dev.off()




