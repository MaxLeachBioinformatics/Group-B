# Authors:Friday
# Modified Date:25-04-2023 BST 00:00
# Version: 1.1

##in terminal
#sudo apt install cmake
##sudo add-apt-repository ppa:dns/gnu
##sudo apt install libgsl-dev

## in R
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.16")

#BiocManager::install(c(tsne", "tsne", "Rtsne"))
#BiocManager::install(c("factoextra", "MASS", "ggplot2", "mclust"))

#BiocManager::install(c("Rfast", "RcppGSL", "RcppZiggurat"),force = TRUE)
#devtools::install_github("jlmelville/smallvis", subdir = "smallvis")

#BiocManager::install(c("FactoMineR"))

#################################################################################
#################################################################################
require(scde)
require(tsne)
require(Rtsne)
library(smallvis)

require(factoextra)
require(MASS)
require(ggplot2)
require(mclust)

library(FactoMineR)

library(igraph)
library(fields)
require(RColorBrewer)

# Set working directory
setwd("/home/yyl23/Documents/Lectures/BS7120_SteeredResearchProject/Script")

# load the data
load("reciprocal.distance_unbiased.RData")
load("direct.dropout.pairwise.biased_500n.RData")
df <- read.csv("GSE67835_transformedRenamed_noLog_rounded.tsv", sep = "\t", header = TRUE)

# Set Gene_Name column as row names and remove original column from data frame
rownames(df) <- df$Gene_Name
df$Gene_Name <- NULL


########### Extra#############################################################
## Calculating optimal number of principle component for the initial dimensional reduction in the default t-SNE setting; 
PC <- prcomp(t(log10(df + 1)), center=TRUE, scale=FALSE)
expl_var <- PC$sdev^2/sum(PC$sdev^2)
barplot(expl_var[1:50], ylab="EXPLAINED VARIANCE",names.arg=paste0("PC",seq(1:50)), col="darkgreen", main="VARIANCE EXPLAINED BY PRINCIPAL COMPONENTS")

# Bootstrapping
N_perm <- 5
expl_var_perm <- matrix(NA, ncol = length(PC$sdev), nrow = N_perm)
for(k in 1:N_perm)
{
  df_perm <- apply(df,2,sample)
  PC_perm <- prcomp(t(log10(df_perm+1)), center=TRUE, scale=FALSE)
  expl_var_perm[k,] <- PC_perm$sdev^2/sum(PC_perm$sdev^2)
}
plot(expl_var[1:50]~seq(1:50), ylab="EXPLAINED VARIANCE",
     col="darkgreen", type='o', xlab="PRINCIPAL COMPONENTS")
lines(colMeans(expl_var_perm)[1:50]~seq(1:50),col="red")
legend("topright", c("Explained by PCS", "Explained by chance"),
       fill=c("darkgreen","red"), inset=0.02)

pval <- apply(t(expl_var_perm) >= expl_var,1,sum) / N_perm
plot(pval[1:50]~seq(1:50),col="darkred",type='o', xlab="PRINCIPAL COMPONENTS",ylab="PVALUE")
optPC<-head(which(pval>=0.05),1)-1
mtext(paste0("OPTIMAL NUMBER OF PRINCIPAL COMPONENTS = ",optPC))
#############################################################################


## T-Distributed Stochastic Neighbor Embeddin (t-SNE)
# Factorise cell by the prefix
cell.labels <- substr(colnames(df), 0, 3)
groups <- factor(cell.labels, levels = c("AST", "NEU", "MIC", "OLI", "OPC", "END", "HYB", "FQU", "FRE"))

# Assign colours to each group
sgCol <- rainbow(length(levels(groups)))[groups]

# Define the modified ecb function
ecb <- function(x, y) {
  plot(x, t = 'n')
  text(x, labels = groups, col = sgCol)
  return(NULL)
}


# Compute standard t-SNE with early exaggeration only for correlation matrix using bias model, optimal perplexity, K~N^(1/2)
set.seed(42)  # tsne has some stochastic steps (gradient descent) so need to set seed 
tsne_res_bias <- tsne(direct.dist, epoch_callback = ecb, k = 2, initial_dims = 30, perplexity=20, epoch=1000)

# Compute t-SNE for unbiased correlations; force global topolgical clustering by increasing perplexity
set.seed(42)
tsne_res <- tsne(reciprocal.dist, epoch_callback = ecb, k = 2, initial_dims = 30, perplexity=80, epoch=1000)

# Create the plot and save it to a PNG file
pdf("baised_tsne_named_clustering.pdf")
set.seed(42) 
tsne_res_bias <- tsne(direct.dist, epoch_callback = ecb, k = 2, initial_dims = 30, perplexity=20, epoch=1000)
dev.off()
pdf("unbaised_tsne_named_clustering.pdf")
set.seed(42)
tsne_res <- tsne(reciprocal.dist, epoch_callback = ecb, k = 2, initial_dims = 30, perplexity=80, epoch=1000)
dev.off()

# saving tsne data
write.csv(tsne_res, file="unbaised_comparing-dens.csv")
save(tsne_res,file="unbaised_comparing-dens.RData")

# save black and white tsne
png("unbaised_tsne_bw.png", width = 600, height = 600)
plot(tsne_res, main="t-SNE, unbiased error model with reciprocal dropout weighting \n(stimu=500;perplex=50;epoch.iter=1000, error=1.16)")
dev.off()

# Calculate geometric density on the tsne result and save to file
geo_plot <- ggplot(as.data.frame(tsne_res),aes(x=V1,y=V2)) + geom_point()+geom_density2d(adjust=.4)
ggsave("geom_density.png", plot = geo_plot)

# Create a density perspective plot and save it
png("density_perspective.png", width = 600, height = 600)
dens <- densityMclust(tsne_res) # Generating a Gaussian plot
plot(dens,what = 'density',type='persp') # Density perspective plot
dev.off()

## Barnes-Hut approximation of t-SNE to compare trading off accuracy with speed 
set.seed(42)
tsne_out <- Rtsne(reciprocal.dist, dims = 2, is_distance=TRUE, perplexity=30, max_iter=2000, check_duplicates = TRUE, theta = 0.5, stop_lying_iter = 500, exaggeration_factor = 12, verbose = TRUE)
plot(tsne_out$Y, col=sgCol, pch=16, main='Barnes-Hut approximation of t-SNE with early exaggeration\nunbiased error model with reciprocal dropout weighting\n(stimu=500;perplex=30;max_iter=2000;stop_lying_iter=500;exaggeration_factor=12;error=1.37)')


## Force Clustering experiment
# Default (early)
set.seed(42)
tsne_out <- smallvis(reciprocal.dist, perplexity = 50, eta = 100, exaggeration_factor = 12, stop_lying_iter = 100, epoch_callback = ecb)
# Late only
tsne_out <- smallvis(reciprocal.dist, perplexity = 50, eta = 100, late_exaggeration_factor = 12, start_late_lying_iter = 100, epoch_callback = ecb)
# Standard and late to force clustering
set.seed(42)
tsne_out <- smallvis(reciprocal.dist, perplexity = 50, eta = 100, exaggeration_factor = 12, stop_lying_iter = 100, late_exaggeration_factor = 1.5, start_late_lying_iter = 100, epoch_callback = ecb)



## Model-based-clustering
mc <- Mclust(tsne_res) 
#plot(mc) # classification plot to show which has aligned with axis; uncertainty plot showing which part of the dataaset is poorly aligned 
summary(mc) # best model Mclust VII (spherical, varying volume) model with 6 components:  BIC -4715.816

# Computes the BIC (Bayesian Information Criterion) for parameterized mixture models given the
loglikelihood, the dimension of the data, and number of mixture components in the model
# Save the BIC graph by mclust
png("mclust_bic.png", width = 600, height = 600)
plot(mc, what = 'BIC', xlab = 'Number of Components')
dev.off()

# Save the BIC graph by fviz_mclust
png("fviz_mclust_bic.png", width = 600, height = 600)
fviz_mclust_bic(mc)
dev.off()

png("mclust_plot.png", width = 1000, height = 600)
par(mfrow = c(1, 2))
plot(mc, what = 'classification') # classification plot to show which cluster has aligned with axis
plot(mc, what = 'uncertainty') # uncertainty plot to show which part of the data is poorly aligned
dev.off()



## Principal component analysis (PCA) 
# A linear dimensionality reduction method to summarize the data into 2 dimensions and then visually identify obvious clusters.

# load df and factorise group
df <- read.csv("GSE67835_transformedRenamed_nolog_notrounded.tsv", sep = "\t")
rownames(df) <- df$Gene_Name
df$Gene_Name <- NULL
cell.labels <- substr(colnames(df),0,3)
groups <- factor(cell.labels,levels=c("AST","NEU","MIC","OLI","OPC", "END", "HYB", "FQU","FRE"))  # Note that principal component analysis generally separates our subpopulations based on their expected group labels. 

# Perform multidimensional scaling (MDS) and convert the matrix of Euclidean distances into a centered Gram matrix, which can be directly used to perform PCA via eigendecomposition
cg_matrix <- cmdscale(reciprocal.dist)

# Run PCA
set.seed(42)
pca_result <- PCA(cg_matrix, scale.unit = TRUE, ncp=5, graph =T)

# pca graph of variables are so clouded with each cells layered on-top of each other which is difficult to interpret
fviz_pca_ind(pca_result, geom.ind = "point", col.ind = groups, palette = c("#FF0000", "#FFAA00", "#AAFF00", "#00FF00", "#00FFAA", "#00AAFF", "#0000FF", "#AA00FF", "#FF00AA"), addEllipses = TRUE, legend.title = "Cell_types")
# Each axis tells you roughly how much variation of the dataset is represented by each of the dimensions; astrocyte group is separated out especially in dimension 2 from the other group as well as fetal quiescent from the adult cell types, indicating an unique gene expression patterns compared to the others.
# the centroid point for the group did not overlap for AST in PCA space; while over gene expression patterns were relatively similarly

png("MDS.png", width = 800, height = 800)
fviz_pca_ind(pca_result, geom.ind = "point", col.ind = groups, palette = c("#FF0000", "#FFAA00", "#AAFF00", "#00FF00", "#00FFAA", "#00AAFF", "#0000FF", "#AA00FF", "#FF00AA"), addEllipses = TRUE, legend.title = "Cell_types")
dev.off()


# Run actual PCA using raw gene expression data
# Set global parameters:
options(ggrepel.max.overlaps = Inf)

# Read dataframe
df3 <- read.csv("GSE67835_transformedRenamed_nolog_notrounded.tsv", sep = "\t")

# filter out genes that are not seen in a sufficient number of cells (colSums(df3>0)>1.8e3 for  low-gene cells often empty wells e.g.>10)
df3 <- df3[rowSums(df3>0)>5, ]

# transpose df
df4 <- data.frame(t(df3[-1]))
colnames(df4) <- df3[, 1]

# transform to make more data normal
df5 <- log10(as.matrix(df4)+1)

# Factorise group
cell.labels_raw <- substr(rownames(df5),0,3)
groups_raw <- factor(cell.labels_raw,levels=c("AST","NEU","MIC","OLI","OPC", "END", "HYB", "FQU","FRE"))

# standardized PCA by scaling to unit variance
pca_result_raw <- PCA(df5, scale.unit = TRUE, ncp = 18, graph =T)

# View the tables of eigenvalues and the percentages of inertia associated with each axis, each cells and distance from the cloud's centre of gravity, the result of the each component including its coordinate value, its quality of representation on the axis given by the squared cosine 
summary(pca_result_raw)

# Description of the dimensions
#dimdesc(pca_result_raw)

#plot(pca_result_raw, cex=0.6, invisible = "none", select="cos2 0.5", title="Individual PCA graph")
#plot(pca_result_raw, cex=0.6, select="contrib 10", title="Individual PCA graph")

# show the grap, named and coloured 
fviz_pca_ind(pca_result_raw, geom = c("point", "text"), label = "all", invisible = "none", labelsize = 2, pointsize = 0.5,, xlim = c(-60, 100), ylim = c(-25, 43))
fviz_pca_ind(pca_result_raw, geom = c("point"), label = "all", invisible = "none", pointsize = 1, col.ind = groups_raw, palette = c("#FF0000", "#FFAA00", "#AAFF00", "#00FF00", "#00FFAA", "#00AAFF", "#0000FF", "#AA00FF", "#FF00AA"), addEllipses = TRUE, ellipse.level = 0.95, legend.title = "Cell_types", xlim = c(-80, 100), ylim = c(-30, 45))

# Save the graphs
png("PCA_named.png", width = 800, height = 800)
fviz_pca_ind(pca_result_raw, geom = c("point", "text"), label = "all", invisible = "none", labelsize = 2, pointsize = 0.5,, xlim = c(-60, 100), ylim = c(-25, 43))
dev.off()
png("PCA_grouped.png", width = 800, height = 800)
fviz_pca_ind(pca_result_raw, geom = c("point"), label = "all", invisible = "none", pointsize = 1, col.ind = groups_raw, palette = c("#FF0000", "#FFAA00", "#AAFF00", "#00FF00", "#00FFAA", "#00AAFF", "#0000FF", "#AA00FF", "#FF00AA"), addEllipses = TRUE, ellipse.level = 0.95, legend.title = "Cell_types", xlim = c(-80, 100), ylim = c(-30, 45))
dev.off()



## MINIMUM SPANNING TREES for all fetal cells, community identification and computation of longest paths, and relative expression for MHCI genes
# Read the data
df <- read.csv("GSE67835_transformedRenamed_noLog_rounded.tsv", sep = "\t", header = TRUE)

# Set Gene_Name column as row names and remove original column from data frame
rownames(df) <- df$Gene_Name
df$Gene_Name <- NULL

# subset data to contain only columns of fetal cells using grepl to match a regular expression to a target and returns TRUE if a match is found and FALSE otherwise.
# Vectorised the function and pass a vector of strings to match and return with a vector of boolean values
df_fetal <- df[ , grepl("FQU|FRE*", colnames(df))]
fetal.cell.labels <- substr(colnames(df_fetal),0,3)
fetal.groups <- factor(fetal.cell.labels,levels=c("FQU","FRE"))

# Map the colours
fetal_color <- fetal.groups[as.numeric(fetal.groups)]
#table(fetal_color)

# read the fetal pairwise distance csv file
reciprocal.dist.fetal <- read.csv("reciprocal_dist_fetal_group.csv", sep = ",")
rownames(reciprocal.dist.fetal) <- reciprocal.dist.fetal$X
reciprocal.dist.fetal$X <- NULL

# Convert correlation matrix in to pairwise distance for pairs defined by the upper triangle of the distance matrix
xy <- t(combn(colnames(reciprocal.dist.fetal), 2))
fetal.distance <- data.frame(xy, dist=reciprocal.dist.fetal[xy])

# Create a graph adjacency based on correlation distances between genes in pairwise fashion.
graph_obj <- graph.adjacency(as.matrix(reciprocal.dist.fetal), mode = "undirected", weighted = TRUE, diag=FALSE)
#E(graph_obj)$weight

# Convert the graph adjacency object into a minimum spanning tree
mst <- mst(graph_obj)

# Plot the tree object
par(mfrow=c(1,2), bg="white", mar = c(0, 0, 0, 0))
set.seed(42)
plot.igraph(mst, edge.color="grey70", vertex.size=3,
            vertex.label=NA, vertex.color=fetal_color)
legend("left",bty = "n",legend=levels(fetal.groups),fill= c("#FFAA00", "#00AAFF"), border=NA,  cex = 0.8)
title("Minimal Spanning Tree", line = -6, cex.main = 0.8)


# identify communities in the tree object based on 'edge between-ness'
mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1
set.seed(42)
plot(
  mst.clustering, mst,
  edge.color="grey70", vertex.size=3,
  vertex.label=NA, vertex.color=fetal_color)
title("All fetal communities", line = -6, cex.main = 0.8)
dev.off()


# Find the longest path in the MST
Simple <- all_simple_paths(mst, "FQU_SRR1974938")
# specify max distance and colour the path in red
LP <- which.max(sapply(Simple, length))
EL1 <- rep(Simple[[LP]], each=2)[-1]
EL1 <- EL1[-length(EL1)]
ECol <- rep("gray", ecount(mst))
ECol[get.edge.ids(mst, EL1)] <- "red"
set.seed(42)
plot(mst, edge.color=ECol, vertex.size=3,
     vertex.label=NA, vertex.color=fetal_color, main="Longest path")
legend("left",bty = "n",legend=levels(fetal.groups),fill= c("#FFAA00", "#00AAFF"), border=NA,  cex = 0.8)


# Find FAT3 expression level in all fetal cells 
df_fetal.FAT3 <- df_fetal[grepl("FAT3", rownames(df_fetal)), ]
df_fetal.FAT3.inverse <- data.frame(t(df_fetal.FAT3[-1]))
# transform to make more data normal
df_fetal.FAT3.inverse <- log10(as.matrix(df_fetal.FAT3.inverse)+1)

# set the gradient colour palette and its content
my_palette <- colorRampPalette(c("white", "#F0E442", "red"))(n = 1000)
node.colors <- (df_fetal.FAT3.inverse) / (max(abs(df_fetal.FAT3.inverse))) * 1000
par(mar=c(0.5, 0.5, 0.5, 6))
set.seed(42)
plot(mst, vertex.color=my_palette[node.colors], vertex.size=3,vertex.label=NA)
title("FAT3", line=-2, cex = 0.5)

# Modified legend code to create a smooth gradient legend
par(new = TRUE)
image.plot(axes = FALSE, legend.only = TRUE, col = my_palette, zlim = c(0, 1), legend.shrink = 0.2, legend.width = 1, 
           axis.args = list(at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), cex.axis = 0.7), 
           legend.args = list(text = "Relative\nExpression:\n", cex=0.8))
dev.off()


# Find EGR1 expression level in all fetal cells; filter out NEGR1 & NEGR1.IT1
df_fetal.EGR1 <- df_fetal[grepl("^EGR1$", rownames(df_fetal)), ]
df_fetal.EGR1.inverse <- data.frame(t(df_fetal.EGR1[-1]))
# transform to make more data normal
df_fetal.EGR1.inverse <- log10(as.matrix(df_fetal.EGR1.inverse)+1)

# set the gradient colour palette and its content
my_palette <- colorRampPalette(c("white", "#F0E442", "red"))(n = 1000)
node.colors <- (df_fetal.EGR1.inverse) / (max(abs(df_fetal.EGR1.inverse))) * 1000
# set right-sided to be more spacious for legend
par(mar=c(0.5, 0.5, 0.5, 6))
set.seed(42)
plot(mst, vertex.color=my_palette[node.colors], vertex.size=3,vertex.label=NA)
title("EGR1", line=-2, cex = 0.5)

par(new = TRUE)
image.plot(axes = FALSE, legend.only = TRUE, col = my_palette, zlim = c(0, 1), legend.shrink = 0.2, legend.width = 1, 
           axis.args = list(at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), cex.axis = 0.7), 
           legend.args = list(text = "Relative\nExpression:\n", cex=0.8))
dev.off()


## Looking into fetal cell populations' lack of expression of HLA-A, -B, and -C, and MHCI-associated genes
# Find all HLA-A expression level to generate a relative expression to all cell types
df.HLA.A <- df[grepl("^HLA-A$", rownames(df)), ]

# subset fetal cell types only
df_fetal.HLA.A <- df_fetal[grepl("^HLA-A$", rownames(df_fetal)), ]
df_fetal.HLA.A.inverse <- as.matrix(data.frame(t(df_fetal.HLA.A[-1])))

# Set the gradient colour palette and its content
my_cb_palette <- colorRampPalette(c("white","#0072B2"))(n=1000)
# Calculate the relative expression of HLA-A of fetal cells to max of all cell types
node.colors <- (df_fetal.HLA.A.inverse) / (max(abs(df.HLA.A)))*1000

# Plot layer 1 graph
par(mar=c(0.5, 0.5, 0.5, 6))
set.seed(42)
plot(mst, vertex.color=my_cb_palette[node.colors], vertex.size=3,vertex.label=NA)
title("HLA-A", line=-2, cex = 0.5)

# plot legend
par(new = TRUE)
image.plot(axes = FALSE, legend.only = TRUE, col = my_cb_palette, zlim = c(0, 1), legend.shrink = 0.3, legend.width = 1, 
           axis.args = list(at = c(0, 0.5, 1), labels = c("0", (max(abs(df.HLA.A)))/2, (max(abs(df.HLA.A)))), cex.axis = 0.7), 
           legend.args = list(text = "Absolute\nCount:\n", cex=0.8))
dev.off()


# HLA-B expression level
df.HLA.B <- df[grepl("^HLA-B$", rownames(df)), ]
df_fetal.HLA.B <- df_fetal[grepl("^HLA-B$", rownames(df_fetal)), ]
df_fetal.HLA.B.inverse <- as.matrix(data.frame(t(df_fetal.HLA.B[-1])))

my_cb_palette <- colorRampPalette(c("white","#0072B2"))(n=1000)
node.colors <- (df_fetal.HLA.B.inverse) / (max(abs(df.HLA.B)))*1000

par(mar=c(0.5, 0.5, 0.5, 6))
set.seed(42)
plot(mst, vertex.color=my_cb_palette[node.colors], vertex.size=3,vertex.label=NA)
title("HLA-B", line=-2, cex = 0.5)

par(new = TRUE)
image.plot(axes = FALSE, legend.only = TRUE, col = my_cb_palette, zlim = c(0, 1), legend.shrink = 0.3, legend.width = 1, 
           axis.args = list(at = c(0, 0.5, 1), labels = c("0", (max(abs(df.HLA.B)))/2, (max(abs(df.HLA.B)))), cex.axis = 0.7), 
           legend.args = list(text = "Absolute\nCount:\n", cex=0.8))
dev.off()


# HLA-C expression level
df.HLA.C <- df[grepl("^HLA-C$", rownames(df)), ]
df_fetal.HLA.C <- df_fetal[grepl("^HLA-C$", rownames(df_fetal)), ]
df_fetal.HLA.C.inverse <- as.matrix(data.frame(t(df_fetal.HLA.C[-1])))

my_cb_palette <- colorRampPalette(c("white","#0072B2"))(n=1000)
node.colors <- (df_fetal.HLA.C.inverse) / (max(abs(df.HLA.C)))*1000

par(mar=c(0.5, 0.5, 0.5, 6))
set.seed(42)
plot(mst, vertex.color=my_cb_palette[node.colors], vertex.size=3,vertex.label=NA)
title("HLA-C", line=-2, cex = 0.5)

par(new = TRUE)
image.plot(axes = FALSE, legend.only = TRUE, col = my_cb_palette, zlim = c(0, 1), legend.shrink = 0.3, legend.width = 1, 
           axis.args = list(at = c(0, 0.5, 1), labels = c("0", (max(abs(df.HLA.C)))/2, (max(abs(df.HLA.C)))), cex.axis = 0.7), 
           legend.args = list(text = "Absolute\nCount:\n", cex=0.8))
dev.off()


# MHCI-associated genes, TAPBP expression level
df.TAPBP <- df[grepl("^TAPBP$", rownames(df)), ]
df_fetal.TAPBP <- df_fetal[grepl("^TAPBP$", rownames(df_fetal)), ]
df_fetal.TAPBP.inverse <- as.matrix(data.frame(t(df_fetal.TAPBP[-1])))

my_cb_palette <- colorRampPalette(c("white","#CC79A7"))(n=1000)
node.colors <- (df_fetal.TAPBP.inverse) / (max(abs(df.TAPBP)))*1000

par(mar=c(0.5, 0.5, 0.5, 6))
set.seed(42)
plot(mst, vertex.color=my_cb_palette[node.colors], vertex.size=3,vertex.label=NA)
title("TAPBP", line=-2, cex = 0.5)

par(new = TRUE)
image.plot(axes = FALSE, legend.only = TRUE, col = my_cb_palette, zlim = c(0, 1), legend.shrink = 0.3, legend.width = 1, 
           axis.args = list(at = c(0, 0.5, 1), labels = c("0", (max(abs(df.TAPBP)))/2, (max(abs(df.TAPBP)))), cex.axis = 0.7), 
           legend.args = list(text = "Absolute\nCount:\n", cex=0.8))
dev.off()


# ERAP1 expression level
df.ERAP1 <- df[grepl("^ERAP1$", rownames(df)), ]
df_fetal.ERAP1 <- df_fetal[grepl("^ERAP1$", rownames(df_fetal)), ]
df_fetal.ERAP1.inverse <- as.matrix(data.frame(t(df_fetal.ERAP1[-1])))

my_cb_palette_2 <- colorRampPalette(c("white","#009E73"))(n=1000)
node.colors <- (df_fetal.ERAP1.inverse) / (max(abs(df.ERAP1)))*1000

par(mar=c(0.5, 0.5, 0.5, 6))
set.seed(42)
plot(mst, vertex.color=my_cb_palette_2[node.colors], vertex.size=3,vertex.label=NA)
title("ERAP1", line=-2, cex = 0.5)

par(new = TRUE)
image.plot(axes = FALSE, legend.only = TRUE, col = my_cb_palette_2, zlim = c(0, 1), legend.shrink = 0.3, legend.width = 1, 
           axis.args = list(at = c(0, 0.5, 1), labels = c("0", (max(abs(df.ERAP1)))/2, (max(abs(df.ERAP1)))), cex.axis = 0.7), 
           legend.args = list(text = "Absolute\nCount:\n", cex=0.8))
dev.off()
