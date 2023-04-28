# Authors:Friday
# Modified Date:01-04-2023 BST 00:00
# Version: 1.0
#################################################################################
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.16")
# adding .libPaths() path[0] to the argument
#BiocManager::install(c("Cairo","edgeR","pcaMethods"), lib = "/alice-home/1/y/yyl23/R/x86_64-pc-linux-gnu-library/4.2")

## May encounter Error in checkSlotAssignment(object, name, value) if not revert flexmix to
# the version prior to the update on April 28, 2017
#require(devtools)
#install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")

# in bash, reinstall scde in order for a few C programs that depend on flexmix to recompile properly.
# download the 1.99.2 stable release from https://github.com/hms-dbmi/scde/releases
# R CMD INSTALL 1.99.2.tar.gz 

#install.packages(c("scde"))
#################################################################################
#################################################################################
# load Single-Cell Differential Expression Analysis
require(scde)
require(parallel)
# load boot package for the weighted correlation implementation
require(boot)

# Similar procedure to biased model
setwd("/home/y/yyl23/R_analysis/BS7120")

# Since we fit models separately for each cell in SCDE, ideally we would pass in non-normalized counts.
df <- read.csv("GSE67835_transformedRenamed_noLog_rounded.tsv", sep = "\t", header = TRUE)

# Set Gene_Name column as row names and remove original column from data frame
rownames(df) <- df$Gene_Name
df$Gene_Name <- NULL

# filtering low count cells and genes without any expression
counts <- df[rowSums(df)>0,]
nGenes <- length(counts[,1])

# Generate error model WITHOUT supplying groups, therefore all cells are considered as a whole
scde.fitted.model.unbiased <- scde.error.models(counts = counts, n.cores = parallelly::availableCores() / 2, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- scde.fitted.model.unbiased$corr.a > 0
table(valid.cells)
scde.fitted.model.unbiased <- scde.fitted.model.unbiased[valid.cells, ]

# Save data
save(scde.fitted.model.unbiased, file="scde_fit_unbiased.RData")
write.csv(as.data.frame(as.matrix(scde.fitted.model.unbiased)), file="scde_fit_unbiased.csv")


# get expression magnitude estimates
o.fpm <- scde.expression.magnitude(scde.fitted.model.unbiased, counts = counts)

# detectCores() may return a missing value/1 and it does not give the number of “allowed” cores
# This construct is guaranteed to always return at least one core.
# ncores <- max(1L, detectCores(), na.rm = TRUE)
# Alternative use availableCores() instead, which guarantees to always return at least one core.

# Use the reciprocal weighting of the Pearson correlation to give increased weight 
# to pairs of observations where a gene expressed (on average) at a level x1 observed 
# in a cell c1 would not be likely to fail in a cell c2
# Using direct weighting without modifying "nam" would lead a mismatch error
n.simulations <- 500; k <- 0.95;
cell.names <- colnames(counts); names(cell.names) <- cell.names;

reciprocal.dist <- as.dist(1 - do.call(rbind, mclapply(cell.names, function(nam1) {
  unlist(lapply(cell.names, function(nam2) {
    # reciprocal probabilities
    f1 <- scde.failure.probability(models = scde.fitted.model.unbiased[nam1,,drop = FALSE], magnitudes = o.fpm[, nam2])
    f2 <- scde.failure.probability(models = scde.fitted.model.unbiased[nam2,,drop = FALSE], magnitudes = o.fpm[, nam1])
    # weight factor
    pnf <- sqrt((1-f1)*(1-f2))*k +(1-k); 
    boot::corr(log10(cbind(counts[, nam1], counts[, nam2])+1), w = pnf)
  }))
},mc.cores = parallelly::availableCores()/2)), upper = FALSE)

# saving output
save(reciprocal.dist,file="reciprocal.distance_unbiased.RData")
write.csv(as.data.frame(as.matrix(reciprocal.dist)), file="reciprocal_dist_unbiased.csv")
