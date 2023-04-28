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

# Set working directory
setwd("/home/y/yyl23/R_analysis/BS7120")

# Since we fit models separately for each cell in SCDE, ideally we would pass in non-normalized counts.
df <- read.csv("GSE67835_transformedRenamed_noLog_rounded.tsv", sep = "\t", header = TRUE)

# Set Gene_Name column as row names and remove original column from data frame
rownames(df) <- df$Gene_Name
df$Gene_Name <- NULL

# A factor object can be defined to index the different cell types.
# We run both the error modelling with and without the reference cell types
# so we can see how different methods perform in recapitulating these labels
cell.labels <- substr(colnames(df),0,3)
groups <- factor(cell.labels,levels=c("AST","NEU","MIC","OLI","OPC", "END", "HYB", "FQU","FRE"))
names(groups) <- colnames(df)
table(groups)
#AST NEU MIC OLI OPC END HYB FQU FRE 
# 62 131  16  38  18  20  46 110  25 

# filtering low count cells and genes without any expression
counts <- df[rowSums(df)>0,]
nGenes <- length(counts[,1])

# Generate error model
scde.fitted.model <- scde.error.models(counts = counts, groups = groups, n.cores = parallelly::availableCores() / 2, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
# supply the groups argument, so that the error models for the two cell types are fit independently (using nine different sets of "robust" genes). 
# The fitting process relies on a subset of robust genes that are detected in multiple cross-cell comparisons. 
# Due to the number of cross-cell comparisons, this step is fairly computationally intensive.

# # Here, corr.a and corr.b are slope and intercept of the correlated component fit, conc.* refer to the concomitant fit, corr.theta is the NB over-dispersion, and fail.r is the background Poisson rate (fixed).
head(scde.fitted.model)

#Particularly poor cells may result in abnormal fits, most commonly showing negative corr.a, and should be removed.
# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- scde.fitted.model$corr.a > 0
table(valid.cells)
scde.fitted.model <- scde.fitted.model[valid.cells, ]

# Save data
save(scde.fitted.model, file="scde_fit.RData")
write.csv(as.data.frame(as.matrix(scde.fitted.model)), file="scde.fitted.bias_model.csv")

# genes with a non-zero expression are not detected in some cells due to failure to amplify the RNA.
p.self.fail <- scde.failure.probability(models = scde.fitted.model, counts = counts)

# Direct weighting down weights the contribution of a given gene to the cell-to-cell distance based on the probability that the given measurement is a drop-out event (i.e. belongs to the drop-out component)
# Stimulate drop-outs using 500 (or more) sampling rounds and generate cell-to-cell distances in terms of Pearson correlation matrix
n.simulations <- 500; k <- 0.9;
cell.names <- colnames(counts); names(cell.names) <- cell.names;

dl <- mclapply(1:n.simulations,function(i) {
  scd1 <- do.call(cbind,lapply(cell.names,function(nam) {
    x <- counts[,nam];
    # replace predicted drop outs with NA values
    x[!as.logical(rbinom(length(x),1,1-p.self.fail[,nam]*k))] <- NA;
    x;
  }))
  rownames(scd1) <- rownames(counts); 
  # calculate correlation on the complete observation pairs using the remaining points
  cor(log10(scd1+1),use="pairwise.complete.obs");
}, mc.cores = parallelly::availableCores() / 2)
# mclapply is a parallelized version of lapply, provided mc.cores > 1: for mc.cores == 1 it simply calls lapply. 
# We use N/2 with one level of parallelism in outer loop to keep overhead happening once.
# availableCores() handles the above case automatically, and it guarantees to always return at least one core

# calculate average distance across sampling rounds
direct.dist <- as.dist(1-Reduce("+",dl)/length(dl))

# saving output
save(direct.dist,file="direct.dropout.pairwise.biased_500n.RData")
write.csv(as.data.frame(as.matrix(direct.dist)), file="direct_dropout_adj_pair_dist_biased_500n.csv")