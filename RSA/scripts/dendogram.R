# ************ Dendogram
# AK 28.08.2023

library(dendextend)
library(circlize)

my_data = test_df
# Select the columns for clustering
data_for_clustering <- my_data[, c("var1", "var2", "mean_all")]

# Perform hierarchical clustering
dist_matrix <- dist(data_for_clustering)  # Calculate the distance matrix
hc <- hclust(dist_matrix)  # Perform hierarchical clustering

# Plot the dendrogram
#plot(hc, main = "Dendrogram of Clustering")

# Convert hierarchical clustering object to a dendrogram
dendro <- as.dendrogram(hc)

dendro <- dendro %>% 
  color_branches(k=4) %>% 
  color_labels

plot(dendro)

# Create a circular dendrogram
# plot the radial plot
par(mar = rep(0,4))
# circlize_dendrogram(dend, dend_track_height = 0.8) 
circlize_dendrogram(dendro, labels_track_height = NA, dend_track_height = .4) 
