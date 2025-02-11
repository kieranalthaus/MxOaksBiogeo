library(magrittr)
library(dplyr)
library(phytools)
source("PrepData/SCRIPTS/functions.R")
`%notin%` = Negate(`%in%`)

 ## 1. Load data ---------------------------------------------------------------
#' A positive D-statistic (i.e. an excess of ABBA) points to introgression 
#' between P2 and P3, whereas a negative D-statistic (i.e. an excess of BABA) 
#' points to introgression between P1 and P3

# Load trees
white.tree = read.tree(file = "../DSUITE/data/dsuite_white_redout.nwk")
red.tree = read.tree(file = "../DSUITE/data/dsuite_red_whiteout.nwk")

# relabel tips
red.tree$tip.label = red.tree$tip.label %>%
  strsplit("\\|") %>% 
  lapply(function(x) x[1]) %>% 
  unlist()

white.tree$tip.label = white.tree$tip.label %>%
  strsplit("\\|") %>% 
  lapply(function(x) x[1]) %>% 
  unlist()

# ladderize
red.tree = ladderize(red.tree)
white.tree = ladderize(white.tree)

## Load D-stats
# For white oaks
white_dstat = read.table(file = "../DSUITE/out/dsuite_white_filter_redout_mx_white_redout_Dmin.txt", 
           sep = "\t",
           header = TRUE)

white_dstat[,c("P1", "P2", "P3")] = apply(X = white_dstat[,c("P1", "P2", "P3")],
      MARGIN = 2,
      FUN = function(x)
        x %>% 
        strsplit(split = "\\|") |>
        lapply(
          FUN = function(x) x[1]) |>
        unlist()
        )

# For red oaks
red_dstat = read.table(file = "../DSUITE/out/dsuite_red_filter_whiteout_mx_red_whiteout_Dmin.txt", 
                       sep = "\t",
                       header = TRUE)

red_dstat[,c("P1", "P2", "P3")] = apply(X = red_dstat[,c("P1", "P2", "P3")],
                                          MARGIN = 2,
                                          FUN = function(x)
                                            x %>% 
                                            strsplit(split = "\\|") |>
                                            lapply(
                                              FUN = function(x) x[1]) |>
                                            unlist()
)

## 2. Merge in bioregions ------------------------------------------------------
# Load states data
states = read.table("MX_oaks_DEC/DATA/DEC_states.txt", 
                    sep = ",",
                    header = TRUE)

# Get state data for red and white oaks 
white.states = states |>
  dplyr::filter(species %in% unique(white_dstat$P1))
red.states = states  |>
  dplyr::filter(species %in% unique(red_dstat$P1))

# Merge states data
white_w_states = same_dif_dstat(dstats = white_dstat, states = white.states)
red_w_states = same_dif_dstat(dstats = red_dstat, states = red.states)

# Remove NAs
white_w_states = white_w_states[!is.na(white_w_states$Dstatistic),]
red_w_states = red_w_states[!is.na(red_w_states$Dstatistic),]

## 3. Question 1 ---------------------------------------------------------------
##' (1) Is there more post-speciation gene flow among species in the same
#' geographic area?
# create new objects for this question
white_q1 = white_w_states
red_q1 = red_w_states

# p-value adjust
white_q1$p.value = p.adjust(p = white_q1$p.value, method = "bonferroni")
red_q1$p.value = p.adjust(p = red_q1$p.value, method = "bonferroni")

# just significant comparisons
white_q1 = filter(white_q1, p.value < 0.05)
red_q1 = filter(red_q1, p.value < 0.05)

# summarize for white oaks
white_q1 = white_q1 %>% 
  select(P2, P3, p.value, Dstatistic, region) %>% 
  group_by(P2,P3) %>% 
  mutate(MAX = max(Dstatistic)) %>% 
  filter(Dstatistic == MAX) %>% 
  ungroup() %>% 
  select(-MAX)

# summarize for red oaks
red_q1 = red_q1 %>% 
  select(P2, P3, p.value, Dstatistic, region) %>% 
  group_by(P2,P3) %>% 
  mutate(MAX = max(Dstatistic)) %>% 
  filter(Dstatistic == MAX) %>% 
  ungroup() %>% 
  select(-MAX)

table(white_q1$region) # 18-24 / diff-same
table(red_q1$region) # 353-59 / diff-same

## 4. D-stat permute geography -------------------------------------------------
## First, for white oaks
# number of white oaks
white_n = length(white.states$range)

# Randomly sample ranges with replacement 1000 times
random_states = replicate(n = 1000,
          expr = sample(x = white.states$range,
                        size = white_n, 
                        replace = FALSE))

#' Create null distribution of geographic areas, and calculate the number of 
#' significant d-statistics per region
white_null = list()
for(i in 1:ncol(random_states)){
  tmp.states = white.states
  tmp.states$range = random_states[,i]
  
  white_w_states = same_dif_dstat(dstats = white_dstat, states = tmp.states)
  
  # Remove NAs
  white_w_states = white_w_states[!is.na(white_w_states$Dstatistic),]
  
  white_q1 = white_w_states
  
  white_q1$p.value = p.adjust(p = white_q1$p.value, method = "bonferroni")
  
  white_q1 = filter(white_q1, p.value < 0.05)
  
  white_q1 = white_q1 %>% 
    select(P2, P3, p.value, Dstatistic, region) %>% 
    group_by(P2,P3) %>% 
    mutate(MAX = max(Dstatistic)) %>% 
    filter(Dstatistic == MAX) %>% 
    ungroup() %>% 
    select(-MAX)
  
  white_null[[i]] = table(white_q1$region)["same"]
  print(i/ncol(random_states))
}

# Save output
saveRDS(object = white_null, file = "../DSUITE/out/20241106_white_dmin_null.RData")

## For red oaks
# number of red oaks
red_n = length(red.states$range)

# Randomly sample ranges with replacement 1000 times
random_states = replicate(n = 1000,
                          expr = sample(x = red.states$range,
                                        size = red_n, 
                                        replace = FALSE))

#' Create null distribution of geographic areas, and calculate the number of 
#' significant d-statistics per region
red_null = list()

for(i in 1:ncol(random_states)){
  tmp.states = red.states
  tmp.states$range = random_states[,i]
  
  red_w_states = same_dif_dstat(dstats = red_dstat, states = tmp.states)
  
  # Remove NAs
  red_w_states = red_w_states[!is.na(red_w_states$Dstatistic),]
  
  red_q1 = red_w_states
  
  red_q1$p.value = p.adjust(p = red_q1$p.value, method = "bonferroni")
  
  red_q1 = filter(red_q1, p.value < 0.05)
  
  red_q1 = red_q1 %>% 
    select(P2, P3, p.value, Dstatistic, region) %>% 
    group_by(P2,P3) %>% 
    mutate(MAX = max(Dstatistic)) %>% 
    filter(Dstatistic == MAX) %>% 
    ungroup() %>% 
    select(-MAX)
  
  red_null[[i]] = table(red_q1$region)["same"]
}

# save data
saveRDS(object = red_null, file = "../DSUITE/out/20241106_red_dmin_null.RData")

# Load data
white_null = readRDS(file = "../DSUITE/out/20241106_white_dmin_null.RData")
red_null = readRDS(file = "../DSUITE/out/20241106_red_dmin_null.RData")

# Conduct t.tests
white.ttest = t.test(mu = 24, x = unlist(white_null), alternative = "two.sided")
red.ttest = t.test(mu = 59, x = unlist(red_null), alternative = "two.sided")

# Quantiles for white oaks
proportions(table(unlist(white_null) >= 24))
mean(white_null >= 24)

# Quantiles for red oaks
proportions(table(red_null >= 59))
mean(red_null >= 59)


# Plot histograms 
pdf(file = "../Figures/20241107_dstat_null_histograms.pdf",
    height = 7.35,
    width = 4.99)
par(mfrow = c(2,1))
red_range = range(unlist(red_null))
hist(x = unlist(red_null), breaks = 50, xlim = range(unlist(red_null)), axes = F, xlab = "Significant Within-region D-statistics", main = "")
axis(1, at = seq(from = red_range[1], to = red_range[2], by = 10))
axis(2, las = 1)
abline(v = 59, col = "red", lwd = 3)
abline(v = mean(unlist(red_null)), lty = 2, col = 'black', lwd = 3)

white_range = range(unlist(white_null))
hist(x = unlist(white_null), breaks = 50, xlim = c(white_range[1], 24), axes = F, xlab = "Significant Within-region D-statistics", main = "")
axis(1, at = seq(from = white_range[1], to = 24, by = 2))
axis(2, las = 1)
abline(v = 24, col = "red", lwd = 3)
abline(v = mean(unlist(white_null)), lty = 2, col = 'black', lwd = 3)
legend(x = 10,
       y = 150,
       legend = c("Empirical Value",
                  "Null Mean"),
       lty = c(1,2),
       lwd = 2,
       cex = 1,
       col = c("red", "black"),
       bty = "n")
dev.off()

## 5. Question 2 ---------------------------------------------------------------
##' (2) Is there more post-speciation gene flow among close relatives?

## Red oaks first
# Create matrix
red.tree.mat = cophenetic(drop.tip(red.tree, "Quercus_lobata"))

# adjust p-value
red_dstat$p.value.adj = p.adjust(p = red_dstat$p.value, method = "bonferroni")

# Format data
red_dstat.mat = red_dstat %>% 
  select(P2, P3, p.value, Dstatistic, p.value.adj) %>% 
  mutate(sig = case_when(
    p.value.adj <= 0.05 ~ 1,
    p.value.adj > 0.05 ~ 0,
  )) %>% 
  group_by(P2, P3) %>% 
  mutate(dupe = duplicated(P2,P3)) %>% 
  filter(dupe == F) %>%
  ungroup() %>% 
  select(P2, P3, sig, Dstatistic) %>% 
  reshape2::dcast(
    formula = P2~P3,
    value.var = "Dstatistic"
  ) %>% 
  drop_na(P2) %>% 
  select(-`NA`) %>% 
  column_to_rownames(var = "P2")

# replcae NA w/ 0
red_dstat.mat[is.na(red_dstat.mat)] = 0

# run mantel test 
mantel.test(m1 = red.tree.mat, 
            m2 = red_dstat.mat, 
            nperm = 1000, 
            graph = TRUE,
            alternative = "less")

# run mantel test w/ vegan package
notin = colnames(red_dstat.mat)[colnames(red_dstat.mat) %notin% rownames(red_dstat.mat)]
notindex = which(rownames(red.tree.mat) == notin)
red.tree.mat = red.tree.mat[-notindex,-notindex]
vegan::mantel(xdis = red.tree.mat,
              ydis = red_dstat.mat,
              permutations = 10000)

## white oaks
# Create matrix
white.tree.mat = cophenetic(drop.tip(white.tree, "Quercus_rubra"))

# adjust p-values
white_dstat$p.value.adj = p.adjust(p = white_dstat$p.value, method = "bonferroni")

# Format data
white_dstat.mat = white_dstat %>% 
  select(P2, P3, p.value, Dstatistic, p.value.adj) %>% 
  mutate(sig = case_when(
    p.value.adj <= 0.05 ~ 1,
    p.value.adj > 0.05 ~ 0,
  )) %>% 
  group_by(P2, P3) %>% 
  mutate(dupe = duplicated(P2,P3)) %>% 
  filter(dupe == F) %>%
  ungroup() %>% 
  select(P2, P3, sig, Dstatistic) %>% 
  reshape2::dcast(
    formula = P2~P3,
    value.var = "Dstatistic"
  ) %>% 
  drop_na(P2) %>% 
  select(-`NA`) %>% 
  column_to_rownames(var = "P2")

# replace NAs with 0
white_dstat.mat[is.na(white_dstat.mat)] = 0

# Two functions for mantel tests
vegan::mantel(xdis = white.tree.mat,
              ydis = white_dstat.mat,
              permutations = 10000)

mantel.test(m1 = white.tree.mat, 
            m2 = white_dstat.mat, 
            nperm = 1000, 
            graph = TRUE,
            alternative = "less")