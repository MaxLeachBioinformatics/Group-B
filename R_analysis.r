# Load SCDE package
lirary(scde)
# load data and read table
sq_data <- read.table("GSE67835_CPM_named_clean.tsv", header=TRUE, row.names=1)
# Convert all non gene nam columns to numeric  
sq_data[, 1:466] <- sapply(sq_data[, 1:466], as.integer)
# clean up the data removig low-count genes and low-quality cells
cleaned <- clean.counts(sq_data, min.lib.size=1000, min.reads = 1, min.detected = 1)
# Fit error models for the different cell types
model <- scde.error.models(counts = cleaned, n.cores = 4, threshold.segmentaticross-fitting cells.sfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
# filter out cells that don't show positive correlation with the expected expression magnitudes (very poor fits)
valid.cells <- model$corr.a > 0
# show table of valid cell out of the 467 cells.... (464 cells are valid)
table(valid.cells)
model <- model[valid.cells, ]
# write the model to a dataframe and export at csv file
write.csv(as.data.frame (model), ("fit_model.csv"))

# get self-fail probabilities (at a given observed count)
p.self.fail <- scde.failure.probability(models = model, counts = cleaned)
# simulate drop-outs
n.simulations <- 200; k <- 0.9;
cell.names <- colnames(cleaned); names(cell.names) <- cell.names;
dl <- mclapply(1:n.simulations,function(i) {
  scd1 <- do.call(cbind,lapply(cell.names,function(nam) {
    x <- cleaned[,nam]
    # replace predicted drop outs with NAs
    x[!as.logical(rbinom(length(x),1,1-p.self.fail[,nam]*k))] <- NA;
    x;
  }))
   rownames(scd1) <- rownames(cleaned); 
  # calculate correlation on the complete observation pairs
  cor(log10(scd1+1),use="pairwise.complete.obs");
}, mc.cores = 1)
# calculate average distance across sampling rounds
direct.dist <- as.dist(1-Reduce("+",dl)/length(dl))
# write output to a dataframe 
scde_output <- as.data.frame(as.matrix(direct.dist))
# export as csv file
write.csv(scde_output, "PairwiseD_matrix.csv") 
