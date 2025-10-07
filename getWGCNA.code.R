###############################################################################################
##################################### WGCNA Analysis ##########################################
rm(list=ls())
#install.packages("WGCNA")
#BiocManager::install("impute")
#BiocManager::install("preprocessCore")
library(WGCNA)
allowWGCNAThreads()
#ALLOW_WGCNA_THREADS=8
options(stringsAsFactors = FALSE)
collectGarbage()


# try data transformed with 
## rows are samples and columns are OTUs
otu0=read.csv("feature-water-pro.csv", header=T)
#otu0=read.csv("54-OTU_table_newID.txt", header=T, sep="\t",row.names = 1)
#absotu=as.data.frame(t(otu0));
otu01=otu0[,-1]
rownames(otu01)=otu0[,1]
otu1=t(apply(otu01, 1, function(x) x/sum(x)))
otu1=otu1[order(row.names(otu1)),]
otu=sqrt(otu1[, colSums(otu1>0)>=10])
#otu=sqrt(otu1[c(4:6,19:24,16:18),colSums(otu1>0)>=10]) # for the Pi samples
dim(otu)
otu=otu[order(row.names(otu)),]

# hellinger transformation
#sotu=sqrt(otu)
gsg = goodSamplesGenes(otu, verbose = 3);
gsg$allOK

# ### if gsg$allOK is FALSE (as opposed to TRUE), do this block of code and test at the end.
# if (!gsg$allOK)
# {
#   # Optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes)>0)
#     printFlush(paste("Removing genes:", paste(names(otu)[!gsg$goodGenes], collapse = ", ")));
#   if (sum(!gsg$goodSamples)>0)
#     printFlush(paste("Removing samples:", paste(rownames(otu)[!gsg$goodSamples], collapse = ", ")));
#   # Remove the offending genes and samples from the data:
#   otu = otu[gsg$goodSamples, gsg$goodGenes]
# }
# gsg = goodSamplesGenes(otu, verbose = 3);
# gsg$allOK

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
abline(h = 0.09, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 0.09, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = otu[keepSamples, ] # datExpr = sqrt(otu[keepSamples, ])
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##### 1.c Loading clinical trait data
geo=read.csv("Env_Data.csv",header = T)
#geo=read.csv("geochemical_factors.csv",header = T, row.names = 1)
#colnames(geo)[3:7]=c("NH4","SiO4","Nitrate","SRP","Nitrite")
#geo=geo[grep("_S[456]",row.names(geo),perl = T),]
#geo=geo[grep("_S[123]",row.names(geo),perl = T),]
geo1=geo[,c(-1, -10,-11)]
rownames(geo1)=geo[,1]
colnames(geo1)=c("rain","humidity","airT","Wind","Solar","Salinity","waterT","Chl")
geo2=geo1[order(row.names(geo1)),]
dim(geo2)
#datTraits=geo[-1,] # remove the outlier sample
datTraits=geo2
#rownames(datExpr)=rownames(datTraits);
collectGarbage();

# Before we continue with network construction and module detection, 
# we visualize how the clinical traits relate to the sample dendrogram
# Re-cluster samples
#pdf("NPi_sample_geo_clusters.pdf",8,8)
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
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,#unsigned
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "arctiaPiTOM",
                       verbose = 3);
write.table(table(net$colors),file = "Plastics_taxa-in-modules.txt",sep = "\t")

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
     file = "Plastics-networkConstruction-auto.RData");

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
par(mar = c(4, 7, 0, 0)+0.1);
# Display the correlation values within a heatmap plot
#pdf("Pi_modules_geochem_corr_heatmap.pdf",10,10)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50), #blueWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.05,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"));
#dev.off()
# 3.b Gene relationship to trait and important modules: Gene Significance and Module
# Membership
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$rain); # ====>>> here to SELECT TRAITs for regression
names(weight) = "rain" # ====>>> here to SELECT TRAITs for regression
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
module = "brown" ## ====>>>>> pick the MODULE here
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
#pdf("Pi_pink-module_geochem_corr.pdf",8,8)
plott<-verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                          abs(geneTraitSignificance[moduleGenes, 1]),
                          xlab = paste("Module Membership in", module, "module"),
                          ylab = "OTU significance for SRP",
                          main = paste("Module membership vs. OTU significance\n"),
                          abline = T, abline.color = module,
                          cex = 1.2, pch=16, cex.main = 0.9,
                          cex.lab = 1.2, cex.axis = 1.2, 
                          col = module)
#dev.off()

# ## get the total RA of signficant modules in each samples
# #picked.modules=c("red","green","cyan","lightyellow", "pink")
# picked.modules=c("brown","turquoise","red","pink", "grey")
# 
# modtRA=matrix(0,24,10)
# rownames(modtRA)=rownames(otu1)
# colnames(modtRA)=c(paste("Pi",picked.modules,sep = "_"),
#                         paste("NPi",1:5,sep = "_"))
# colnames(modtRA)[6:10]=paste("NPi",picked.modules,sep = "_")
# for (ii in 1:5){
#   #ii=1
#   module = picked.modules[ii]
#   modOTUs=colnames(otu)[moduleColors==module]
#   modOTUs.ra=otu1[,colnames(otu1) %in% modOTUs]
#   modtRA[,ii+5]=rowSums(modOTUs.ra)
# }
# 
# write.table(modtRA,"PiNpi_modules_OTUs_RA.txt",sep="\t",quote = F)
# ## To export the OTUs of of some module of interest
# # get module names
# modNames
# # pick the module you want and get its position index (position number MINUS one) in the vector
# modOTUs=colnames(otu)[moduleGenes]
# modOTUs.coeff=cbind("modcoeff"=as.numeric(geneModuleMembership[moduleGenes, column]),
#              "traitcoeff"=as.numeric(geneTraitSignificance[moduleGenes, 1])
#              )
# rownames(modOTUs.coeff)=modOTUs
# modOTUs.coeff=t(modOTUs.coeff)
# modOTUs.ra=otu1[,colnames(otu1) %in% modOTUs]
# #modOTUs.ra=modOTUs.ra[c(1:6,10:15,7:9,19:24,16:18),]
# modOTUs.ra=modOTUs.ra[c(1:3,10:15,7:9,4:6, 19:24,16:18),]
# modOTUs.merge<-merge(modOTUs.coeff,modOTUs.ra, all.x=T, all.y=T)
# rownames(modOTUs.merge)=c(rownames(modOTUs.ra),rownames(modOTUs.coeff))
# #write.table(modOTUs.merge,"Pi_mod-pink_coeffs.txt",sep = "\t",quote = F)
# piTaxa=read.table("56-OTU_taxanomy.txt",header = F,sep = "\t")
# pinkRA=read.table("Pi_mod-pink_coeffs.txt",header = T,sep = "\t")



# 3.d Summary output of network analysis results
# names(datExpr)
# names(datExpr)[moduleColors=="black"] #
# annot = read.csv(file = "GeneAnnotation.csv");
# dim(annot)
# names(annot)
# probes = names(datExpr)
# probes2annot = match(probes, annot$substanceBXH)
# # The following is the number or probes without annotation:
# sum(is.na(probes2annot))
# # Create the starting data frame
# geneInfo0 = data.frame(substanceBXH = probes,
#                        geneSymbol = annot$gene_symbol[probes2annot],
#                        LocusLinkID = annot$LocusLinkID[probes2annot],
#                        moduleColor = moduleColors,
#                        geneTraitSignificance,
#                        GSPvalue)
# # Order modules by their significance for weight
# modOrder = order(-abs(cor(MEs, weight, use = "p")));
# # Add module membership information in the chosen order
# for (mod in 1:ncol(geneModuleMembership))
# {
#   oldNames = names(geneInfo0)
#   geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
#                          MMPvalue[, modOrder[mod]]);
#   names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
#                        paste("p.MM.", modNames[modOrder[mod]], sep=""))
# }
# # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
# geneInfo = geneInfo0[geneOrder, ]
# write.csv(geneInfo, file = "geneInfo.csv")

# 6.b Exporting to Cytoscape

# Recalculate topological overlap if needed
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

TOM = TOMsimilarityFromExpr(datExpr, power = 6);

# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
#allmodules=c("red","green","purple","cyan","lightyellow","midnightblue",
#             "pink","magenta","brown","turquoise")
# allmodules=c("red","green","yellow","brown","grey",
#              "pink","blue","turquoise")
# for (i in 1:length(allmodules) ){
# #i=1
#   modules=allmodules[i]
    
modules = "brown"

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
                               edgeFile = paste("Plastic_brown_cytoscape-edges_", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("Plastic_brown_cytoscape-nodes_", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,# 
                               nodeNames = modProbes,
                               includeColNames = TRUE,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]
                               )

##################################### WGCNA Analysis ############################################
#################################################################################################
