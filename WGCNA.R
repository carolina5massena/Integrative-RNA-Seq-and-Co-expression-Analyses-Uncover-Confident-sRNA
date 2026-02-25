# ======================================================================================
# USER INPUTS (edit only this section)
# ======================================================================================

# ---- File paths ----
expr_counts_tsv   <- "salmon.merged.gene_counts.tsv"   # expression matrix with gene_name column
traits_csv        <- "table_treatment.csv"             # sample trait table

# ---- Required column names in your files ----
gene_id_col       <- "gene_name"   # column in expr_counts_tsv that contains gene names
trait_sample_col  <- "rownames"    # column in traits_csv that matches sample IDs (rownames(datExpr))
treatment_col     <- "treatment"   # categorical column to dummy-encode

# ---- Output directory and filenames ----
out_dir <- "."  # e.g., "wgcna_outputs"

out_sample_dendro_png        <- "sample_dendrogram_outlier_check.png"
out_sample_dendro_traits_png <- "sample_dendrogram_with_traits.png"
out_softpower_png            <- "soft_thresholding_power_selection.png"
out_gene_dendro_png          <- "gene_dendrogram_with_modules.png"
out_module_trait_pdf         <- "Module_Trait_Heatmap.pdf"
out_cyt_edges                <- "CytoscapeInput-edges.txt"
out_cyt_nodes                <- "CytoscapeInput-nodes.txt"
out_gene_module_csv          <- "gene_module_membership.csv"
out_top_hubs_csv             <- "top_hub_genes_per_module.csv"

# ---- Compute / threading ----
# Set NULL to let WGCNA decide; or set e.g. 8
n_threads <- NULL

# ---- Outlier visualization cut (for the red line only) ----
sample_tree_cut_height <- 140

# ---- Soft-threshold selection ----
powers <- c(1:10, seq(from = 12, to = 20, by = 2))
softPower <- 10  # IMPORTANT: user sets this after looking at the plot

# ---- Network + correlation settings ----
networkType <- "unsigned"     # "unsigned" or "signed"
adjacency_type <- "unsigned"  # usually same as networkType
corFnc <- "bicor"             # "bicor" (robust) or "cor"
tomType <- "unsigned"         # usually matches networkType

# ---- Module detection + merging ----
minModuleSize <- 30
deepSplit <- 2
pamRespectsDendro <- FALSE
MEDissThres <- 0.25  # merge if eigengene corr > 1 - MEDissThres (i.e., > 0.75)

# ---- Trait dummy encoding ----
remove_first_dummy <- FALSE   # keep all levels
remove_selected_columns <- TRUE

# ---- Cytoscape export ----
cytoscape_edge_threshold <- 0.02  # export edges with weight > threshold

# ======================================================================================
# END OF USER INPUTS
# ======================================================================================


# Create output directory if needed
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Helper to build output paths
op <- function(x) file.path(out_dir, x)

# Threading
allowWGCNAThreads(nThreads = n_threads)

# -----------------------------
# Expression loading
# -----------------------------
SAbatch <- read.table(expr_counts_tsv,
                      header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE) %>%
  group_by(.data[[gene_id_col]]) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  column_to_rownames(var = gene_id_col)

SAdatExpr0 <- as.data.frame(t(SAbatch))

# -----------------------------
# Sample dendrogram plot
# -----------------------------
sampleTree <- hclust(dist(SAdatExpr0), method = "average")

png(op(out_sample_dendro_png), width = 12, height = 9, units = "in", res = 300)
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = sample_tree_cut_height, col = "red")
dev.off()

datExpr <- SAdatExpr0

# -----------------------------
# Trait loading + matching
# -----------------------------
traitData <- read.csv(traits_csv, stringsAsFactors = FALSE)

Samples <- rownames(datExpr)
traitRows <- match(Samples, traitData[[trait_sample_col]])
datTraits <- traitData[traitRows, , drop = FALSE]
rownames(datTraits) <- datTraits[[trait_sample_col]]

datTraits_dummy <- fastDummies::dummy_cols(
  datTraits,
  select_columns = treatment_col,
  remove_first_dummy = remove_first_dummy,
  remove_selected_columns = remove_selected_columns
)

datTraits_final <- datTraits_dummy[, !(names(datTraits_dummy) %in% trait_sample_col), drop = FALSE]
rownames(datTraits_final) <- rownames(datTraits)

# Plot dendrogram + traits
sampleTree2 <- hclust(dist(datExpr), method = "average")
traitColors <- numbers2colors(datTraits_final, signed = FALSE)

png(op(out_sample_dendro_traits_png), width = 12, height = 9, units = "in", res = 300)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits_final),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# -----------------------------
# Soft thresholding
# -----------------------------
sft <- pickSoftThreshold(datExpr,
                         networkType = networkType,
                         powerVector = powers,
                         verbose = 5)

png(op(out_softpower_png), width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.80, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# User sets softPower in config

adjacency <- adjacency(datExpr, power = softPower, type = adjacency_type, corFnc = corFnc)
TOM <- TOMsimilarity(adjacency, TOMType = tomType)
dissTOM <- 1 - TOM

geneTree <- flashClust(as.dist(dissTOM), method = "average")

dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = deepSplit,
  pamRespectsDendro = pamRespectsDendro,
  minClusterSize = minModuleSize
)
dynamicColors <- labels2colors(dynamicMods)

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
moduleColors <- merge$colors
MEs <- merge$newMEs

png(op(out_gene_dendro_png), width = 12, height = 9, units = "in", res = 300)
plotDendroAndColors(geneTree, cbind(dynamicColors, moduleColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# -----------------------------
# Module-trait correlations
# -----------------------------
nSamples <- nrow(datExpr)
moduleTraitCor <- cor(MEs, datTraits_final, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

pdf(file = op(out_module_trait_pdf), width = 10, height = 10)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits_final),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

# -----------------------------
# Cytoscape + exports
# -----------------------------
probes <- names(datExpr)
inModule <- is.finite(match(moduleColors, unique(moduleColors[moduleColors != "grey"])))
modProbes <- probes[inModule]
modTOM <- dissTOM[inModule, inModule]
diag(modTOM) <- 0

exportNetworkToCytoscape(
  1 - modTOM,
  edgeFile = op(out_cyt_edges),
  nodeFile = op(out_cyt_nodes),
  weighted = TRUE,
  threshold = cytoscape_edge_threshold,
  nodeNames = modProbes,
  nodeAttr = moduleColors[inModule]
)

geneModuleMembership <- data.frame(Gene = probes, Module = moduleColors)
write.csv(geneModuleMembership, op(out_gene_module_csv), row.names = FALSE)

topHubs <- chooseTopHubInEachModule(datExpr, colorh = moduleColors, omitColors = "grey", power = softPower)
write.csv(topHubs, op(out_top_hubs_csv), row.names = FALSE)

message("WGCNA analysis complete. Check outputs in: ", normalizePath(out_dir))

