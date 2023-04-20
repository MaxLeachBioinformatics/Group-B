##if packages are not already installed in R then install, installation shown in comments


#devtools::install_version('flexmix', '2.3-13') this version of flexmix is needed to be downloaded before SCDE version 1.99.1
#devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)
library(scde)
library(parallel)
#install_version('mclust', version = '5.4.10', repos = 'http://lib.stat.cmu.edu/R/CRAN/')
library(mclust)
##install.packages("Rtsne")
library(Rtsne)
##install.packages("igraph")
library(igraph)
##install.packages("FactoMineR")
library(FactoMineR)
##install.packages("parallelly")
library(parallelly)


data <- read.table('GSE67835_CPM_named_clean.tsv', header = TRUE, row.names = 1, sep ='\t')
###Convert all non gene columns to numeric
data[,1:466] <- sapply(data[,1:466], as.integer)

##clean.counts function from the scde package to remove low-count genes and low-quality cells
clean_data <- clean.counts(data, min.lib.size=466, min.reads = 1, min.detected = 1)

##fit error model to each cell
##  fit an error model to each cell and filter out cells that did not show positive correlation with the expected expression magnitudes.


model <- scde.error.models(counts = clean_data, n.cores = 4, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- model$corr.a > 0 ##only want these
table(valid.cells) ##should print out 466 into terminal

model <- model[valid.cells, ]
write.csv(as.data.frame(model), 'model.csv')  ##creating a file to have a copy in case any crashing occurs

p.self.fail <- scde.failure.probability(models = model, counts = clean_data)
n.simulations <- 500; k <- 0.9;
cell.names <- colnames(clean_data); names(cell.names) <- cell.names;
dl <- mclapply(1:n.simulations,function(i) {
  scd1 <- do.call(cbind,lapply(cell.names,function(nam) {
    x <- clean_data[,nam];
    # replace predicted drop outs with NAs
    x[!as.logical(rbinom(length(x),1,1-p.self.fail[,nam]*k))] <- NA;
    x;
  }))
  rownames(scd1) <- rownames(clean_data); 
  # calculate correlation on the complete observation pairs
  cor(log10(scd1+1),use="pairwise.complete.obs");
}, mc.cores <- parallelly::availableCores()/2 ##parallel::detectCores() / 2  ##mc.cores = parallely:AvailableCores()/2
)
##pairwisedistancematrix
# calculate average distance across sampling rounds
direct.dist <- as.dist(1-Reduce("+",dl)/length(dl))
distances <- as.data.frame(as.matrix(direct.dist))
distances[is.na(distances)] <- 0 ## replacings nans with 0
write.csv(as.data.frame(distances), 'distancematrix.csv')
distancematrix <- as.matrix(direct.dist)


##data needs to be normalised for clustering
normalised_data <- normalize_input(distancematrix)
## Dimensionality reduction using t-SNE:
tsne_results <- Rtsne(normalised_data, perplexity = 10, theta = 0.0)       ##tsne_results <- Rtsne(normalised_data, theta = 0.0)
tsnedf <- as.data.frame(tsne_results$Y) ##converting the Y matrix into a data frame, which can be more easily manipulated and plotted


##CLUSTERING WITH MCLUST

clustered_data <- Mclust (tsnedf)

par(mfrow = c(1, 2))
plot(clustered_data, what = 'BIC', xlab = 'Number of Components')
plot(clustered_data, what = 'classification')

# Set up the PNG file, can alter height and width if need be
png("mclust_plot.png", width = 800, height = 600)
## plot results 
par(mfrow = c(1, 2))
plot(clustered_data, what = 'BIC', xlab = 'Number of Components')
plot(clustered_data, what = 'classification')
## what = "BIC" argument tells plot to display the BIC values, which are a measure of model fit for the clustering model, for different numbers of clusters.
## This creates a plot with the BIC value on the left and the classification results on the right.
# Close the PNG file and save
dev.off()


##MINIMUM SPANNING TREES

#convert the distance matrix into a object (graph) that can be used by igraph using `graph.adjacency`
graph <- graph.adjacency(distancematrix, mode = "undirected", weighted = TRUE)

mst <- minimum.spanning.tree(graph)

##visualising them 

plot(mst, edge.label = E(mst)$weight)

##saving plot as a png

png("mst.png")
plot(mst, edge.label = E(mst)$weight)
dev.off()


## PCA

## create a PCA object 
pca_result <- PCA(normalised_data, graph = TRUE)

png('pca_result.png')
pca_result <- PCA(normalised_data, choix = "ind")  ##  graph = TRUE argument may produce a 3D plot, which cannot be saved as a 2D PNG file
dev.off()


