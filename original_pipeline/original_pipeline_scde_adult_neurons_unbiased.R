# Authors:Friday
# Modified Date:24-04-2023 BST 00:00
# Version: 1.0
#################################################################################
#################################################################################
# load Single-Cell Differential Expression Analysis
require(scde)
require(parallel)
# load boot package for the weighted correlation implementation
require(boot)

# Similar procedure to biased model
#setwd("/home/y/yyl23/R_analysis/BS7120")
setwd("/home/yyl23/Documents/Lectures/BS7120_SteeredResearchProject/Script")

# Since we fit models separately for each cell in SCDE, ideally we would pass in non-normalized counts.
df <- read.csv("GSE67835_transformedRenamed_noLog_rounded.tsv", sep = "\t", header = TRUE)

# Set Gene_Name column as row names and remove original column from data frame
rownames(df) <- df$Gene_Name
df$Gene_Name <- NULL

# filtering low count cells and genes without any expression
counts <- df[rowSums(df)>0,]

# subsetting adult neurons
df_neu <- counts[ , grepl("NEU*", colnames(counts))]

# Generate error model WITHOUT supplying groups, therefore all cells are considered as a whole
scde.fitted.model.neu.unbiased <- scde.error.models(counts = df_neu, n.cores = parallelly::availableCores() / 2, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- scde.fitted.model.neu.unbiased$corr.a > 0
table(valid.cells)
scde.fitted.model.neu.unbiased <- scde.fitted.model.neu.unbiased[valid.cells, ]

# Save data
save(scde.fitted.model.neu.unbiased, file="scde_neurons_unbiased.RData")
write.csv(as.data.frame(as.matrix(scde.fitted.model.neu.unbiased)), file="scde_neurons_unbiased.csv")

# get expression magnitude estimates
o.fpm <- scde.expression.magnitude(scde.fitted.model.neu.unbiased, counts = df_neu)

# Reciprocal weighting to generate correlation matrix
k <- 0.95;
reciprocal.dist <- as.dist(1 - do.call(rbind, mclapply(cell.names, function(nam1) {
  unlist(lapply(cell.names, function(nam2) {
    # reciprocal probabilities
    f1 <- scde.failure.probability(models = scde.fitted.model.neu.unbiased[nam1,,drop = FALSE], magnitudes = o.fpm[, nam2])
    f2 <- scde.failure.probability(models = scde.fitted.model.neu.unbiased[nam2,,drop = FALSE], magnitudes = o.fpm[, nam1])
    # weight factor
    pnf <- sqrt((1-f1)*(1-f2))*k +(1-k); 
    boot::corr(log10(cbind(df_neu[, nam1], df_neu[, nam2])+1), w = pnf)
  }))
},mc.cores = parallelly::availableCores()/2)), upper = FALSE)

# saving output
save(reciprocal.dist,file="reciprocal.distance_neurons.RData")
write.csv(as.data.frame(as.matrix(reciprocal.dist)), file="reciprocal.distance_neurons.csv")
