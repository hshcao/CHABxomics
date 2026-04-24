xdf=read.delim("mG19_genus_abundance.txt",sep = '\t')
xdf=xdf[c(-1:-2, (-nrow(xdf)+1),-nrow(xdf)),]
colnames(xdf)=gsub("_([679])\\.","_0\\1\\.", colnames(xdf))
#write.table(xdf,"mG19_genus_abundance_allSsamples.tsv",sep = "\t",quote = FALSE)
df=xdf[,-1]
rownames(df)=xdf[,1]

df=df[,order(colnames(df))]
library(heatmap3)
pdf("mG19_genus_abundance_allSamples.pdf",8,8)
heatmap3(as.matrix(df), Colv = NA,labRow = "",margins = c(7,1))
dev.off()

############################################################################
################### WGCNA analysis of metabolome with metadata #############
rm(list=ls())
library(WGCNA)
allowWGCNAThreads()
ALLOW_WGCNA_THREADS=14
options(stringsAsFactors = FALSE)
collectGarbage()

otu=as.data.frame(t(df))

# hellinger transformation
sotu=sqrt(otu)
gsg = goodSamplesGenes(otu, verbose = 3);
gsg$allOK

### if gsg$allOK is FALSE (as opposed to TRUE), do this block of code and test at the end.
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(otu)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(otu)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  otu = otu[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(otu, verbose = 3);
gsg$allOK

# check for outliers
sampleTree = hclust(dist(otu), method = "average");
sizeGrWindow(12,9)
#pdf("Pi_sample_clusters.pdf", width = 8, height = 8);
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "NPi_samples", sub="", xlab="", 
     cex.lab = 1.8, cex.axis = 1.5, cex.main = 1.2, cex=1.5
)
#dev.off()

# Plot a line to show the cut
abline(h = 30, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 30, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = otu[keepSamples, ] # datExpr = sqrt(otu[keepSamples, ])
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##### 1.c Loading enviromental factor data
wdf=read.delim("water_chemistry.csv",row.names = 1, header = T, sep = ",")
datTraits=wdf[!rownames(wdf)=="mG19_07.23.2",]
collectGarbage()

# Before we continue with network construction and module detection, 
# we visualize how the clinical traits relate to the sample dendrogram
# Re-cluster samples
#pdf("water_condition_metaGenome_clusters.pdf",8,8)
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",
                    cex.colorLabels = 0.85, cex.dendroLabels = 1, 
                    cex.rowText = 0.8,las=2
)
#dev.off()

## 2.a Automatic network construction and module detection

# 2.a.1 Choosing the soft-thresholding power: analysis of network topology

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# 2.a.2 One-step network construction and module detection
net = blockwiseModules(datExpr, power = 10,
                       TOMType = "unsigned", minModuleSize = 30,#unsigned
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "diel-Msamples",
                       verbose = 3);
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05);

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Taihu_2019mG-networkConstruction-auto.RData");

# 3.a Quantifying module-trait associations
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(4, 7, 2, 1)+0.1);
# Display the correlation values within a heatmap plot
pdf("Taihu_2019_mG-water_chemistry_heatmap.pdf",9,6)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50), #blueWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("mG_genus water-chemistry correlation"));
library(heatmap3)
library(mixOmics)
color.edge <- color.GreenRed(200)
heatmap3(moduleTraitCor, col = color.edge, main = "mG_genus water-chemistry correlation")
dev.off()

# 3.b Gene relationship to trait and important modules: Gene Significance and Module
# Membership
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$NO2); # ====>>> here to SELECT TRAITs for regression
names(weight) = "NO2" # ====>>> here to SELECT TRAITs for regression
# names (colors) of the modules to pick the module
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

# 3.c Intramodular analysis: identifying genes with high GS and MM
pdf("mG19_network-module_correlation_with_NO2.pdf",6,6)
for (i in 1:length(modNames)){
  module = modNames[i] ## ====>>>>> pick the MODULE here
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  #sizeGrWindow(7, 7);
  #par(mfrow = c(1,1));
  plott<-verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                            abs(geneTraitSignificance[moduleGenes, 1]),
                            xlab = paste("Module Membership in", module, "module"),
                            ylab = "Genus significance for NO2",
                            main = paste("Module membership vs. genus significance\n"),
                            abline = T, abline.color = module,cex.axis = 1.2, 
                            cex = 1.2, pch=16, cex.main = 0.9, cex.lab = 1.2,
                            displayAsZero=0.5, col = module)
}
dev.off()

# 6.b Exporting to Cytoscape

# Recalculate topological overlap if needed
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

TOM = TOMsimilarityFromExpr(datExpr, power = 10);

for (i in 1:length(modNames)){
  modules = modNames[i]
  # Select module probes
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  #modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("mG19_network_genus_list_edge_", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("mG19_network_genus_list_node_", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,# 
                                 nodeNames = modProbes,
                                 includeColNames = TRUE,
                                 #altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]
  )
}
##################################### WGCNA Analysis ############################################
#################################################################################################
############# pathway abundance and correlation with water conditions ################
#ydf=read.delim("mG19_pathID-Abundance.txt",header = T, row.names = 1, sep = "\t")
ydf=read.delim("mG19_MARonly-pathID_Abundance.txt",header = T, row.names = 1, sep = "\t")
colnames(ydf)=gsub("mG","mG19_", colnames(ydf))
ydf=ydf[,!colnames(ydf)=="mG19_11.3.1"]
#df=ydf[,order(colnames(ydf))]
library(heatmap3)
pdf("mG19_MARonly_pathway_abundance_allSamples.pdf",8,8)
heatmap3(as.matrix(ydf)), Colv = NA, labRow = "",margins = c(7,1))
dev.off()

rm(list=ls())
library(WGCNA)
allowWGCNAThreads()
ALLOW_WGCNA_THREADS=14
options(stringsAsFactors = FALSE)
collectGarbage()

otu=as.data.frame(t(ydf))

# hellinger transformation
sotu=sqrt(otu)
gsg = goodSamplesGenes(otu, verbose = 3);
gsg$allOK

### if gsg$allOK is FALSE (as opposed to TRUE), do this block of code and test at the end.
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(otu)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(otu)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  otu = otu[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(otu, verbose = 3);
gsg$allOK

# check for outliers
sampleTree = hclust(dist(otu), method = "average");
sizeGrWindow(12,9)
#pdf("Pi_sample_clusters.pdf", width = 8, height = 8);
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "NPi_samples", sub="", xlab="", 
     cex.lab = 1.8, cex.axis = 1.5, cex.main = 1.2, cex=1.5
)
#dev.off()

# Plot a line to show the cut
abline(h = 30, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 30, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = otu[keepSamples, ] # datExpr = sqrt(otu[keepSamples, ])
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##### 1.c Loading enviromental factor data
wdf=read.delim("water_chemistry.csv",row.names = 1, header = T, sep = ",")
datTraits=wdf[!rownames(wdf)=="mG19_11.3.1",]
collectGarbage()

# Before we continue with network construction and module detection, 
# we visualize how the clinical traits relate to the sample dendrogram
# Re-cluster samples
#pdf("water_condition_pathway_clusters.pdf",8,8)
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",
                    cex.colorLabels = 0.85, cex.dendroLabels = 1, 
                    cex.rowText = 0.8,las=2
)
#dev.off()

## 2.a Automatic network construction and module detection

# 2.a.1 Choosing the soft-thresholding power: analysis of network topology

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# 2.a.2 One-step network construction and module detection
net = blockwiseModules(datExpr, power = 8,
                       TOMType = "unsigned", minModuleSize = 30,#unsigned
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "diel-Msamples",
                       verbose = 3);
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05);

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Taihu_2019mG-pathway_networkConstruction-auto.RData");

# 3.a Quantifying module-trait associations
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(4, 7, 2, 1)+0.1);
# Display the correlation values within a heatmap plot
pdf("Taihu_2019mG-water_chemistry_pathway_heatmap.pdf",9,6)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50), #blueWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("mG_genus water-chemistry correlation"));
library(heatmap3)
library(mixOmics)
color.edge <- color.GreenRed(200)
heatmap3(moduleTraitCor, col = color.edge, main = "mG_water-chemistry_pathway_correlation")
dev.off()

# 3.b Gene relationship to trait and important modules: Gene Significance and Module
# Membership
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Temp); # ====>>> here to SELECT TRAITs for regression
names(weight) = "Temp" # ====>>> here to SELECT TRAITs for regression
# names (colors) of the modules to pick the module
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

# 3.c Intramodular analysis: identifying genes with high GS and MM
pdf("mG19_pathway_network-module_correlation_with_Temp.pdf",6,6)
for (i in 1:length(modNames)){
  module = modNames[i] ## ====>>>>> pick the MODULE here
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  #sizeGrWindow(7, 7);
  #par(mfrow = c(1,1));
  plott<-verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                            abs(geneTraitSignificance[moduleGenes, 1]),
                            xlab = paste("Module Membership in", module, "module"),
                            ylab = "Pathway significance for Temp",
                            main = paste("Module membership vs. pathway significance\n"),
                            abline = T, abline.color = module,cex.axis = 1.2, 
                            cex = 1.2, pch=16, cex.main = 0.9, cex.lab = 1.2,
                            displayAsZero=0.5, col = module)
}
dev.off()

# 6.b Exporting to Cytoscape

# Recalculate topological overlap if needed
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

TOM = TOMsimilarityFromExpr(datExpr, power = 10);

for (i in 1:length(modNames)){
  modules = modNames[i]
  # Select module probes
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  #modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("mG19_network_genus_list_edge_", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("mG19_network_genus_list_node_", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,# 
                                 nodeNames = modProbes,
                                 includeColNames = TRUE,
                                 #altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]
  )
}
