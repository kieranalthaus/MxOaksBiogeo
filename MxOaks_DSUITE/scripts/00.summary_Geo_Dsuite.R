library(magrittr)
library(dplyr)
library(phytools)

## 1. Load Data ----------------------------------------------------------------
#' A positive D-statistic (i.e. an excess of ABBA) points to introgression 
#' between P2 and P3, whereas a negative D-statistic (i.e. an excess of BABA) 
#' points to introgression between P1 and P3
#' The D_min values only looks at introgression between P2 and P3

#' First, load the results for white and red oaks into two lists, organized 
#' internally by region (CA, ENA, etc.)
# Load D-stats for white oaks
white.path = list.files(path = "../DSUITE/out/20240913_dsuite_output/",
                     pattern = "^_SETS_white.*_Dmin\\.txt$", 
                     full.names = TRUE)
white.D = lapply(X = white.path,
                 FUN = read.table,
                 sep = "\t",
                 header = TRUE)
names(white.D) = c("CA", "ENA", "SMOC", "SMOR", "TMVB", "TROP")

# Load D-stats for red oaks
red.path = list.files(path = "../DSUITE/out/20240913_dsuite_output/",
                        pattern = "^_SETS_red.*_Dmin\\.txt$",
                        full.names = TRUE)
red.D = lapply(X = red.path,
                 FUN = read.table,
                 sep = "\t",
                 header = TRUE)
names(red.D) = c("CA", "ENA", "SMOC", "SMOR", "TMVB", "TROP")

## 2. D-stats by region --------------------------------------------------------
# For white oaks
# For each region...
lapply(X = white.D, 
       FUN = function(x){
         # Record the rows in the subset
         nrow = nrow(x)
         x %>% 
           # Filter to significant tests
           filter(p.value < 0.05) %>% 
           # How many significant tests are there?
           nrow() %>%
           # What proportion of the total set is significant?
           `/`(nrow) 
       }
)
# For white oaks
# For each region...
lapply(X = red.D,
       FUN = function(x){
         # Record the rows in the subset
         nrow = nrow(x)
         x %>% 
           # Filter to significant tests
           filter(p.value < 0.05) %>% 
           # How many significant tests are there?
           nrow() %>%
           # What proportion of the total set is significant?
           `/`(nrow) 
       }
)
