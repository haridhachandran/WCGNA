# ================================
# WGCNA Step-by-Step Network Analysis
# ================================

# 1. Setup -----------------------------------------------------------------

# Display current directory and set working directory
cat("Current working directory:", getwd(), "\n")
workingDir <- "."  # Change if needed
setwd(workingDir)

# Load required package
library(WGCNA)

# Important options
options(stringsAsFactors = FALSE)

# Enable multi-threading (skip if using RStudio or IDE with own threading)
enableWGCNAThreads()

# Load prepared data (modify filename if needed)
lnames <- load(file = "FemaleLiver-01-dataInput.RData")
cat("Loaded variables:", lnames, "\n")

# 2. Soft-threshold power selection ----------------------------------------

powers <- c(1:10, seq(12, 20, 2))

# Network topology analysis
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot scale-free fit index and mean connectivity
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))

cex1 <- 0.9

# Scale-free fit plot
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], labels = powers, col = "red", cex = cex1)
abline(h = 0.90, col = "red")

# Mean connectivity plot
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red", cex = cex1)

# 3. Construct adjacency matrix --------------------------------------------

softPower <- 6  # Choose based on plot or prior knowledge
adjacency <- adjacency(datExpr, power = softPower)

# 4. Compute Topological Overlap Matrix (TOM) ------------------------------

TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# 5. Hierarchical clustering ------------------------------------------------

geneTree <- hclust(as.dist(dissTOM), method = "average")

sizeGrWindow(12, 9)
plot(geneTree, main = "Gene clustering on TOM-based dissimilarity",
     xlab = "", sub = "", labels = FALSE, hang = 0.04)

# 6. Module identification using Dynamic Tree Cut --------------------------

minModuleSize <- 30

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
cat("Module sizes:\n")
print(table(dynamicMods))

# 7. Assign colors to modules -----------------------------------------------

dynamicColors <- labels2colors(dynamicMods)
cat("Module colors and sizes:\n")
print(table(dynamicColors))

sizeGrWindow(8, 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# 8. Calculate module eigengenes --------------------------------------------

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

# 9. Merge modules with close eigengenes ------------------------------------

MEDissThres <- 0.25
abline(h = MEDissThres, col = "red")

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# 10. Plot dendrogram with merged module colors -----------------------------

sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# 11. Save results ----------------------------------------------------------

moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree,
     file = "FemaleLiver-02-networkConstruction-stepByStep.RData")

cat("Network construction and module detection completed and saved.\n")
