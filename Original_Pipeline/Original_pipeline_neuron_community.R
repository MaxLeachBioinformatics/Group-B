# Authors:Friday
# Modified Date:25-04-2023 BST 00:00
# Version: 1.0

#BiocManager::install(c("M3C"))

#################################################################################
#################################################################################
## Heterogeneity analysis through clustering of neurons
library(ggplot2)
require(M3C) # wrapper for Rtsne and ggplot2
require(gridExtra)
library(igraph)

setwd("/home/yyl23/Documents/Lectures/BS7120_SteeredResearchProject/Script")

load("reciprocal.distance_unbiased.RData")
df <- read.csv("GSE67835_transformedRenamed_noLog_rounded.tsv", sep = "\t", header = TRUE)

# Set Gene_Name column as row names and remove original column from data frame
rownames(df) <- df$Gene_Name
df$Gene_Name <- NULL

# filtering low count cells and genes without any expression
df <- df[rowSums(df)>0,]

# Factorise cell by the prefix
cell.labels <- substr(colnames(df), 0, 3)
groups <- factor(cell.labels, levels = c("AST", "NEU", "MIC", "OLI", "OPC", "END", "HYB", "FQU", "FRE"))

tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Cell types", labels = groups)+ ggtitle("t-SNE") + theme(plot.title = element_text(hjust = 0.5))
tsne(reciprocal.dist, seed = 42, perplex = 12, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df[row.names(df)=='PVALB',])),
     controlscale = TRUE, scale = 2) + ggtitle("parvalbumin (PVALB)") + theme(plot.title = element_text(hjust = 0.5))

load("reciprocal.distance_neurons.RData")

# subsetting adult neurons
df_neu <- df[ , grepl("NEU*", colnames(df))]
#write.csv(as.data.frame(as.matrix(df_neu)), file="neu_count_matrix.csv")

# Rt-SNE, 131^(1/2) = 11.4
tsne(reciprocal.dist, seed = 42, perplex = 12, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Cell types") + ggtitle("t-SNE of neurons") + theme(plot.title = element_text(hjust = 0.5))

# GABAergic neurons - PVALB+, CRHBP+, TAC1+, LHX6+, SOX6+ 
tsne(reciprocal.dist, seed = 42, perplex = 12, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='PVALB',])),
     controlscale = TRUE, scale = 2) + ggtitle("parvalbumin (PVALB)") + theme(plot.title = element_text(hjust = 0.5))
tsne(reciprocal.dist, seed = 42, perplex = 12, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='CRHBP',])),
     controlscale = TRUE, scale = 2) + ggtitle("corticotropin releasing factor binding protein (CRHBP)") + theme(plot.title = element_text(hjust = 0.5))
tsne(reciprocal.dist, seed = 42, perplex = 12, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='TAC1',])),
     controlscale = TRUE, scale = 2) + ggtitle("tachykinin precursor 1 (TAC1)") + theme(plot.title = element_text(hjust = 0.5))
tsne(reciprocal.dist, seed = 42, perplex = 12, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='LHX6',])),
     controlscale = TRUE, scale = 2) + ggtitle("LIM homeobox 6 (LHX6)") + theme(plot.title = element_text(hjust = 0.5))
tsne(reciprocal.dist, seed = 42, perplex = 12, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='SOX6',])),
     controlscale = TRUE, scale = 2) + ggtitle("SRY (sex determining region Y)-box 6 (SOX6)") + theme(plot.title = element_text(hjust = 0.5))

# Excitatory communities - SLC17A7+, NEUROD6+, SATB2+
tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='SLC17A7',])),
     controlscale = TRUE, scale = 2, high = "#0072B2") + ggtitle("solute carrier family 17 member 7 (SLC17A7)") + theme(plot.title = element_text(hjust = 0.5))
tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='NEUROD6',])),
     controlscale = TRUE, scale = 2, high = "#0072B2") + ggtitle("neuronal differentiation 6 (NEUROD6)") + theme(plot.title = element_text(hjust = 0.5))
tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='SATB2',])),
     controlscale = TRUE, scale = 2, high = "#0072B2") + ggtitle("SATB homeobox 2 (SATB2)") + theme(plot.title = element_text(hjust = 0.5))

# Inhibitory/interneurons -  GAD1+, CCK+, CALB2+
tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='GAD1',])),
     controlscale = TRUE, scale = 2, high = "#009E73") + ggtitle("glutamate decarboxylase 1 (GAD1)") + theme(plot.title = element_text(hjust = 0.5))
tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='CCK',])),
     controlscale = TRUE, scale = 2, high = "#009E73") + ggtitle("cholecystokinin (CCK)") + theme(plot.title = element_text(hjust = 0.5))
tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='CALB2',])),
     controlscale = TRUE, scale = 2, high = "#009E73") + ggtitle("calbindin 2 (CALB2)") + theme(plot.title = element_text(hjust = 0.5))


# Layered enriched populations
# PAX6, RELN coexpression (layer 1)
plot1 <- tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='PAX6',])),
              controlscale = TRUE, scale = 2, high = "#56B4E9") + ggtitle("Layer 1: paired box 6 (PAX6)") + theme(plot.title = element_text(hjust = 0.5))
plot2 <- tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='RELN',])),
              controlscale = TRUE, scale = 2, high = "#56B4E9") + ggtitle("Layer 1: reelin (RELN)") + theme(plot.title = element_text(hjust = 0.5))
grid.arrange(plot1, plot2, ncol=2)
ggsave("Layered_1_enriched_neuron_populations.pdf", arrangeGrob(plot1, plot2))

#VIP, TAC3 (in neurons enriched in layers 2-4)
plot3 <- tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='VIP',])),
              controlscale = TRUE, scale = 2, high = "#D55E00") + ggtitle("Layer 2-4: vasoactive intestinal peptide(VIP)") + theme(plot.title = element_text(hjust = 0.5))
plot4 <- tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='TAC3',])),
              controlscale = TRUE, scale = 2, high = "#D55E00") + ggtitle("tachykinin precursor 3(TAC3)") + theme(plot.title = element_text(hjust = 0.5))
grid.arrange(plot3, plot4, ncol=2)
ggsave("Layered_2-4_enriched_neuron_populations.pdf", arrangeGrob(plot3, plot4))

# Layer 3 - coexpressing SPARC, SV2C
plot5 <- tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='SPARC',])),
              controlscale = TRUE, scale = 2, high = "#CC79A7") + ggtitle("Layer 3: secreted protein acidic and cysteine rich(SPARC)") + theme(plot.title = element_text(hjust = 0.5))
plot6 <- tsne(reciprocal.dist, seed = 42, perplex = 30, axistextsize = 10, legendtextsize = 10, dotsize = 1.5, legendtitle = "Expression Level", labels = scale(as.numeric(df_neu[row.names(df_neu)=='SV2C',])),
              controlscale = TRUE, scale = 2, high = "#CC79A7") + ggtitle("synaptic vesicle glycoprotein 2C(SV2C)") + theme(plot.title = element_text(hjust = 0.5))
grid.arrange(plot5, plot6, ncol=2)
ggsave("Layered_3_enriched_neuron_populations.pdf", arrangeGrob(plot5, plot6))



# Create a graph adjacency based on correlation distances between genes in pairwise fashion.
graph_obj <- graph.adjacency(as.matrix(reciprocal.dist), mode = "undirected", weighted = TRUE, diag=FALSE)

# Convert the graph adjacency object into a minimum spanning tree
mst <- mst(graph_obj)

# identify communities in the mst object based on the Walktrap community finding algorithm, see Pascal Pons, Matthieu Latapy: Computing communities in large networks using random walks, http://arxiv.org/abs/physics/0512106
mst.communities <- walktrap.community(mst, steps = 4)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1
set.seed(42)
plot(
  mst.clustering, mst,
  edge.color="grey70", vertex.size=3,
  vertex.label=NA)
title("Minimal Spanning Tree Neuronal Community", line = 1, cex.main = 0.8)

# saving communities
capture.output(mst.communities$names, file = "mst.communities_names.txt")
capture.output(mst.communities$membership, file = "mst.communities_membership.txt")

#GABAergic(inhibtory) neurons: PVALB+, CRHBP+, TAC1+, LHX6+, SOX6+ 
#Excitatory: SLC17A7+, NEUROD6+, SATB2+
#Inhibitory/interneurons: GAD1+, CCK+, CALB2+

#### Subsequent analysis in jypter_notebook ####