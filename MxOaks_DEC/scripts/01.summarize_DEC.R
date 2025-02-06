# Load packages
library(RevGadgets)
library(tidyverse)
library(geoscale)
data("timescales")
`%notin%` = Negate(`%in%`)
source("MX_oaks_DEC/SCRIPTS/plot_anc_range.util.R")

# Get time data 
t1 = subset(timescales$ICS2013, timescales$ICS2013[, 'Type'] == 'Epoch')

# Set input file paths
fp = "MxOaks_DEC/m0/"
label_fn = "MxOaks_DEC/data/range_labels.txt"

# Set output file paths
tree_fn = "MxOaks_DEC/m0/m0.ase.tre"

# Get state labels
state_labels = read.table(file = label_fn, sep = ",", header = TRUE) %>% 
  apply(MARGIN = 1, function(x)
    x[1]
  ) %>% 
  str_trim(side = "left") %>% 
  `names<-`(0:(length(.)-1))

# process the ancestral states
ase <- processAncStates(path = tree_fn,
                        state_labels = state_labels)


## Prepare for plotting
# set colors
range_colors = c("#FCA636FF","#E16462FF", "#B12A90FF", "#0D0887FF", "#F0F921FF", "#6A00A8FF")

# Create matrix for node values
range.data = data.frame(matrix(data = 0,
                             nrow = nrow(ase@data),
                             ncol = 6,
                             dimnames = list(
                               c(), # Rows
                               c(state_labels[2:7]) # Columns
                             )))

# extract range info from each tip and node
range_str = ase@data %>% 
  dplyr::select(end_state_1) %>% 
  dplyr::pull() %>% 
  strsplit(split = "\\+")

# Extract the MAP regions each species is in
for(i in 1:length(range_str)){
  # Which columns match
  x = which(colnames(range.data) %in% range_str[[i]])
  # Add pie values
  if(length(x) > 1){
    range.data[i,x] = 0.5 
  } else {
    range.data[i,x] = 1
  }  
}

# Add node labels
range.data$node = ase@data$node

# Get node numbers
nodes = phylobase::phylo4d(x = ase@phylo) %>% 
  as("data.frame") %>% 
  dplyr::filter(node.type != "tip") %>% 
  .$node

# Remove tip nodes
node.data = dplyr::filter(range.data, node %in% nodes)
node.data = node.data[order(node.data$node),]
tip.data = dplyr::filter(range.data, node %notin% nodes)
tip.data$node = as.numeric(tip.data$node)
tip.data = tip.data[order(tip.data$node),]

# Pie chart
node.pie = nodepie(node.data, cols = 1:6, color = range_colors, alpha = 1)
tip.pie = nodepie(tip.data, cols = 1:6, color = range_colors, alpha = 1)

# Get the size of the node that is proportional to it's MAP
node.size = ase@data %>% 
  dplyr::filter(node %in% names(node.pie)) %>% 
  .$end_state_1_pp %>% 
  as.numeric()

## Plot tree
h<-max(nodeHeights(ase@phylo))
plot(ladderize(ase@phylo), show.tip.label = F, edge.width = 1.5, y.lim = c(-15,179))
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
labs<-seq(0,h,by=5)
col = alpha(colour = "grey40", alpha = 0.1)
labs<-seq(0,-h-10,by=-10)
{
pdf("../Figures/20241231_DEC_ase2.pdf",
     height=6.87,
     width = 3.19)
par(mar = c(0,0,0,0))
  plot(ladderize(ase@phylo),
       show.tip.label = TRUE,
       edge.width = 1.5,
       y.lim = c(-15,179),
       cex = 0.2,
       label.offset = 0.5)
  for(i in 1:ase@phylo$Nnode){
    plotrix::floating.pie(pp$xx[i+Ntip(ase@phylo)],
                          pp$yy[i+Ntip(ase@phylo)],
                          radius = node.size[i]/1.2,
                          x=setNames(as.numeric(node.data[i,-7]), colnames(node.data)[-7]),
                          col=range_colors,
                          border="transparent")
  }
  for(i in 1:length(ase@phylo$tip.label)){
    plotrix::floating.pie(pp$xx[i],
                          pp$yy[i],
                          pch = 1,
                          radius=0.3,
                          x=setNames(as.numeric(tip.data[i,-7]), colnames(tip.data)[-7]),
                          col=range_colors,
                          border="transparent")
  }
  segments(x0 = h-t1[2:7, 'Start'][c(-5,-6)], 
           y0 = -4, 
           y1 = length(ase@phylo$tip.label) + 1, 
           lty = 'dashed', 
           col = 'gray')
  rect(xleft = c(53.58121),
       xright = rev(c(50.99321, 48.24821, 30.58121, 19.68121, 0)),
       ytop = -5, 
       ybottom = -3,
       col = c("white", "gray40"),
       border = "black",
       lwd = .75)
  text(x = h-t1[2:7, 'Start'][c(-6,-5)], 
       y = -6.5, 
       labels = round(t1[2:7, 'Start'], 1)[c(-5,-6)],
       cex = 0.4,
       srt = 0)
  text(x = c(51.5, 48.5, 37.41108, 21.12708, 7.64208, -7.40792),
       # x = h-t1[2:7, 'Midpoint']-2,
       y = -1,
       # labels = t1[2:7, 'Name'],
       labels = c("Ple.", "Pli.", "Miocene", "Oligocene", "Eocene"),
       srt = 0,
       adj = 0,
       xpd = TRUE,
       cex = 0.5)
  dev.off()
}

pdf(file = "../Figures/20241005_dec_legend.pdf",
    height = 6.93,
    width = 5.82)
par(mar = c(0,0,0,0))
{plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(
  "center",
  xjust = 0,
  legend=c(0.2,0.4,0.6,0.8),
  pch=21,
  pt.cex=c(0.4,0.6,0.8,1)*20,
  bty = "n",
  cex = 5,
  pt.lwd = 4,
  y.intersp = 1.5
)
}
dev.off()

range_colors = c("#FCA636FF","#E16462FF", "#B12A90FF", "#0D0887FF", "#F0F921FF", "#6A00A8FF")
pdf(file = "../Figures/20241007_dec_legend.pdf")
par(mar = c(0,0,0,0))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center",
       legend = c("Eastern North America",
                  "Sierra Madre Occidental",
                  "Trans-Mexican Volcanic Belt",
                  "Tropics",
                  "Sierra Madre Oriental",
                  "California"),
       bty = "n",
       fill = range_colors,
       border = "transparent",
       cex = 2
)
dev.off()