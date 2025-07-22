# WCGNA

Weighted Gene Co-expression Network Analysis (WGCNA) pipeline for gene expression data

This scripted workflow performs weighted gene co-expression network analysis on expression data (datExpr), aiming to identify modules of co-expressed genes and their relationships.


**Key steps:**

•	Setting working directory and loading essential packages

•	Selection of soft-thresholding power to approximate scale-free topology

•	Construction of weighted adjacency and topological overlap matrices (TOM)

•	Hierarchical clustering and dynamic tree cutting to define gene modules

•	Assigning module colors for easy visualization

•	Calculating module eigengenes and clustering them

•	Merging similar modules based on eigengene similarity

•	Visualization of dendrograms and module colors at every crucial step

•	Saving results for downstream analyses such as module-trait correlation


**Typical applications:**

•	Identification of gene modules linked to clinical traits or experimental conditions

•	Discovery of hub genes within biologically relevant modules

•	Integration with other omics data or phenotypic traits in systems biology studies

