##CLUSTERING WITH MCLUST

clustered_data <- Mclust (tsnedf)
clustered_data <- Mclust (tsnedf) ## clustering analysis, storing in variable clustered_data

## visualisation of clustering analysis
par(mfrow = c(1, 2))
plot(clustered_data, what = 'BIC', xlab = 'Number of Components')
plot(clusters, what = "uncertainty", xlab = 'clusters', ylab = 'clustering uncertainty')
png("mclust_plot.png", width = 800, height = 600)
## plot results 
par(mfrow = c(1, 2))
# generating a plot of the Bayesian information criterion (BIC) values
## what = 'BIC', specifies that the BIC values should be plotted as done in the paper
## xlab = specifies the label for the x-axis of the plot.
plot(clustered_data, what = 'BIC', xlab = 'Number of Components')
#  plot where each data point is colored according to its assigned class
plot(clusters, what = "uncertainty", xlab = 'clusters', ylab = 'clustering uncertainty')
## what = "BIC" argument tells plot to display the BIC values, which are a measure of model fit for the clustering model, for different numbers of clusters.
## This creates a plot with the BIC value on the left and the classification results on the right.
# Close the PNG file and save
dev.off()