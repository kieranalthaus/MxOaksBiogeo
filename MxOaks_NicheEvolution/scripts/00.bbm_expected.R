library(phytools)
library(ape)

threshold = 2.7

# Load data and prepare tree
h<-max(nodeHeights(tree))
mt<-make.era.map(tree,seq(0,h,length.out=100))
st<-map.to.singleton(mt)

# Function to find first passage time
first_passage_time <- function(trait_values, times, threshold) {
  exceed <- which(trait_values >= threshold)[1]
  if (is.na(exceed)) return(NA)
  return(times[exceed])
}

out = c()
for(j in 1:n_sims){
# Simulate traits
sim_traits <- fastBM(tree = st, 
                     sig2 = fitBBM$sigsq, 
                     a = fitBBM$x0, 
                     internal = TRUE,
                     bounds = fitBBM$bounds)

# Get node heights
node_heights <- nodeHeights(st)

# Calculate first passage times for each branch
n_nodes <- nrow(st$edge)
passage_times <- numeric(n_nodes)

for (i in 1:n_nodes) {
  parent <- st$edge[i, 1]
  child <- st$edge[i, 2]
  start_time <- node_heights[i, 1]
  end_time <- node_heights[i, 2]
  branch_length <- end_time - start_time
  
  # Interpolate trait values along the branch
  trait_values <- seq(sim_traits[parent], sim_traits[child], length.out=100)
  times <- seq(start_time, end_time, length.out=100)
  
  passage_times[i] <- first_passage_time(trait_values, times, threshold = threshold)
}
  out[j] = min(passage_times, na.rm = TRUE)
  print(j)

}



# # Print results
# for (i in 1:n_nodes) {
#   if (!is.na(passage_times[i])) {
#     cat("Branch", i, "(Node", st$edge[i,1], "to", st$edge[i,2], ") reached 10 at time:", passage_times[i], "\n")
#   }
# }


# # Visualize
par(mar=c(5.1,4.1,2.1,2.1))
phenogram(st, sim_traits, spread.cost=c(1,0), ftype='off', lwd = 2, col="grey", fsize=0.7)
phenogram(st, sim_traits, spread.cost=c(1,0), ftype='off',lwd = 0.5, col="white", fsize=0.7, add = T)
clip(0, h, min(x), max(x))
grid()

# Add points for first passage times
for (i in 1:n_nodes) {
  if (!is.na(passage_times[i])) {
    points(passage_times[i], threshold, col="red", pch=19)
  }
}
abline(h=threshold, col="blue", lty=2)

