################################################################################
# For:
#    eukaryotic algae: Chlorella and Thalassiosira
#    cyanobacteria: M and Dolichospermum
#    epibiotic bacteria: TA02/07/08, TF02/F05, TW01/08
#DO:
#  WGCNA
#  pheatmap
#  EBSeq
#  topGO
################################################################################
##### Chlorella #####
library(pheatmap)
library(RColorBrewer)
library(topGO)
## Chlorella
rm(list = ls())
df1=read.csv("chlorella_20_gene_functions.tsv.txt", sep = "\t",header = T)
df=read.csv("rsem_tpm_chlorella.tsv", sep = "\t",header = T, row.names = 1)
df=df[rowSums(df)>0,]
colnames(df)=gsub(".genes","",colnames(df))
rownames(df)=gsub("gene-ACK3TF_\\d+_","", rownames(df), perl=T)
gene_list=rownames(df)
rownames(df)=df1$geneID[match(rownames(df), df1$LocusID)]

pdf("Chlorella_TPM_by-date-site.pdf",5,5)
pheatmap(log10(df+1), fontsize = 7,cluster_rows = TRUE, cluster_cols = FALSE,
         color = colorRampPalette(brewer.pal(n = 4, name = "Oranges"))(50))
dev.off()

geneID2GO <- readMappings(file = "chlorella_protID_geneID_GOterm_paired2.tsv")
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% gene_list))
names(geneList) <- geneNames
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
      orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)
fullTerms <- Term(GOTERM[allRes$GO.ID])
allRes$Full_Term <- fullTerms
write.csv(allRes, file = "chlorella_topGO_results.csv", row.names = FALSE)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 20, useInfo = 'all')
################################################################################
## Thalassiosira ##
library(pheatmap)
library(RColorBrewer)
library(topGO)
library(DESeq2)

rm(list=ls())
df=read.csv("rsem_tpm_thalassiosira.tsv", sep = "\t",header = T, row.names = 1)
df=df[rowSums(df)>0,]
colnames(df)=gsub(".genes","",colnames(df))
head(rownames(df))
rownames(df)=gsub("gene-", "", rownames(df))
parts <- do.call(rbind, lapply(strsplit(rownames(df), "_"), function(x) {
  part1 <- paste(x[1:2], collapse = "_")
  part2 <- paste(x[-c(1:2)], collapse = "_")
  c(part1, part2)
}))
colnames(parts) <- c("locusID", "geneName")
df_genenames <- as.data.frame(parts, stringsAsFactors = FALSE)
rownames(df) = df_genenames$geneName
#pdf("Thalassiosira_TPM_by-date-site.pdf",5,10)
pheatmap(log10(df+1), fontsize = 6, cluster_rows = TRUE, cluster_cols = FALSE,
         color = colorRampPalette(brewer.pal(n = 4, name = "Oranges"))(50))
dev.off()

deg_list=list(
  mT6v7=c("THAPSDRAFT_37880", "THAPSDRAFT_37959", "THAPSDRAFT_bd1112",
  "THAPSDRAFT_bd1258", "THAPSDRAFT_bd598", "THAPSDRAFT_bd1426",
  "THAPSDRAFT_bd1039", "THAPSDRAFT_26246", "THAPSDRAFT_268516", 
  "THAPSDRAFT_41829", "THAPSDRAFT_bd1023", "THAPSDRAFT_261702","THAPSDRAFT_4710",
  "THAPS_270312","THAPS_6796","THAPSDRAFT_23653","THAPSDRAFT_22117",
  "THAPSDRAFT_27836","THAPSDRAFT_bd467","THAPSDRAFT_260392","THAPSDRAFT_21175",
  "THAPSDRAFT_26573","THAPSDRAFT_23988","THAPSDRAFT_bd1342","THAPSDRAFT_bd1387",
  "THAPSDRAFT_bd1313"))
deg_list[['mT6v11']]=c("THAPSDRAFT_37880","THAPS_270304","THAPSDRAFT_260392",
"THAPSDRAFT_21175","THAPSDRAFT_10048","THAPS_41392","THAPSDRAFT_268304",
"THAPSDRAFT_bd1863","THAPSDRAFT_29842","THAPS_23576","THAPSDRAFT_bd1658",
"THAPSDRAFT_12193","THAPSDRAFT_270340","THAPSDRAFT_24535","THAPSDRAFT_5009",
"THAPSDRAFT_bd388","THAPSDRAFT_10521","THAPSDRAFT_5026","THAPSDRAFT_bd1339",
"THAPSDRAFT_270233","THAPSDRAFT_33270","THAPSDRAFT_21292","THAPSDRAFT_270210",
"THAPSDRAFT_bd1387","THAPSDRAFT_bd1313","THAPSDRAFT_24694","THAPSDRAFT_24112",
"THAPSDRAFT_24325","THAPSDRAFT_2587","THAPSDRAFT_34487","THAPSDRAFT_34541",
"THAPSDRAFT_268862","THAPSDRAFT_8748","THAPSDRAFT_10331","THAPSDRAFT_24403",
"THAPSDRAFT_12191","THAPSDRAFT_bd1771","THAPSDRAFT_21336","THAPSDRAFT_25074",
"THAPSDRAFT_29228","THAPS_35726","THAPSDRAFT_12695","THAPSDRAFT_9394",
"THAPSDRAFT_21328","THAPSDRAFT_23801","THAPSDRAFT_269900")
deg_list['mT11']=c("THAPSDRAFT_bd1258","THAPSDRAFT_bd598","THAPSDRAFT_bd1426",
                   "THAPSDRAFT_bd1039","THAPSDRAFT_26246","THAPSDRAFT_268516",
                   "THAPSDRAFT_41829","THAPSDRAFT_bd467","THAPSDRAFT_260392",
                   "THAPSDRAFT_21175","THAPSDRAFT_24112","THAPSDRAFT_24325")

# count data for DESeq2
rm(list=ls())
df=read.csv("rsem_readCount_thalassiosira.tsv", sep = "\t",header = T, row.names = 1)
df=df[rowSums(df)>0,]
colnames(df)=gsub(".genes","",colnames(df))
rownames(df)=gsub("gene-","", rownames(df))
head(rownames(df))
parts <- do.call(rbind, lapply(strsplit(rownames(df), "_"), function(x) {
  part1 <- paste(x[1:2], collapse = "_")
  part2 <- paste(x[-c(1:2)], collapse = "_")
  c(part1, part2)
}))
colnames(parts) <- c("locusID", "geneName")
df_genenames <- as.data.frame(parts, stringsAsFactors = FALSE)
rownames(df) = df_genenames$locusID

cts=round(as.matrix(df))
#cts=cts[, c(1:3, 6:8)]
cts=cts[, c(6:8,13:15)]
coldata=data.frame(sample=colnames(cts),condition=rep(c("mT07","mT11"), each=3))
rownames(coldata)=coldata[,1]
all(rownames(coldata)==colnames(cts))

# running DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, name = "condition_mT11_vs_mT07")
resultsNames(dds)
resLFC <- lfcShrink(dds, coef = "condition_mT11_vs_mT07", type = "apeglm")
resLFC
resOrdered <- res[order(res$pvalue),]
summary(res)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")
resSig <- subset(resOrdered, padj < 0.1)
degs=list(mt6v7=rownames(resSig))

#gene_list=deg_list[["mT6v11"]]
gene_list=deg_list[['mT11']]
geneID2GO <- readMappings(file = "./geneIDGO/thalassiosira_geneID_GOterms_paired2.txt")
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% gene_list))
names(geneList) <- geneNames
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
allRes <- GenTable(GOdata, classicFisher = resultFisher,
                orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 50)
fullTerms <- Term(GOTERM[allRes$GO.ID])
allRes$Full_Term <- fullTerms
write.csv(allRes, file = "thalassiosira_mT11_topGO_results.csv", row.names = FALSE)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo='all')
################################################################################
######## DESeq2 & WGCNA ########
##_____ 1. DESeq2 + topGO _____
##--Microcystis flos aquae--
library(pheatmap)
library(RColorBrewer)
library(topGO)
library(DESeq2)
library(openxlsx)
# count data for DESeq2
rm(list=ls())
df=read.csv("rsem_readCount_Mflosaquae.tsv", sep = "\t",header = T, row.names = 1)
df=df[rowSums(df==0) <7,]
colnames(df)=gsub(".genes$","",colnames(df))
colnames(df)=gsub("\\.","_",colnames(df))
parts <- do.call(rbind, lapply(strsplit(rownames(df), "_"), function(x) {
  part1 <- paste(x[1:2], collapse = "_")
  part2 <- paste(x[-c(1:2)], collapse = "_")
  c(part1, part2)
}))
colnames(parts) <- c("locusID", "geneName")
parts=as.data.frame(parts)
rownames(df)=parts$locusID
dim(df)

list_up=list_down=list()
list_names=c('mt7v6', 'mt9v6', 'mt11v6', 'mt12v6')
cts0=round(as.matrix(df))
for (ll in 1:length(list_names)){
  if(ll==1){
    cts=cts0[, 1:6]
  } else if(ll==2){
    cts=cts0[, c(1:3,7:9)]
  } else if(ll==3){
    cts=cts0[,c(1:3, 13:15)]
  }else{
    cts=cts0[,c(1:3, 16:18)]
  }
  coldata=data.frame(sample=colnames(cts),condition=gsub("_(1|2|3)","",colnames(cts), perl = TRUE))
  rownames(coldata)=coldata[,1]
  coldata$condition <- factor(coldata$condition)
  
  all(rownames(coldata)==colnames(cts))
  # running DESeq2
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata,
                              design = ~ condition)
  keep <- rowSums(counts(dds)) >= 1
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds)
  resultsNames(dds)
  condName=resultsNames(dds)[2]
  res <- results(dds, name = condName)
  resOrdered <- res[order(res$pvalue),]
  res05 <- results(dds, alpha=0.05)
  res05_order = res05[order(res05$padj), ]
  summary(res05)
  sum(res05$padj < 0.05, na.rm=TRUE)
  resSig <- subset(resOrdered, padj < 0.05)
  list_up[[list_names[ll]]]=rownames(resSig)[resSig$log2FoldChange>0 & resSig$padj<0.05]
  list_down[[list_names[ll]]]=rownames(resSig)[resSig$log2FoldChange<0 & resSig$padj<0.05]
}

topGO_list=list()
goafile="./geneIDGO/Mflosaque_protID_geneID_GOterm_paired2.tsv"
geneID2GO <- readMappings(file = goafile)
geneNames <- names(geneID2GO)
for (i in 1:length(up_list)){
  #i=1
  genesDN=list_down[[list_names[i]]]
  genesUP=list_up[[list_names[i]]]
  for (j in 1:2){
    if (j==1) {
      gene_list=genesUP
      tmpl_name=paste0('up_',list_names[i])
    } else {
      gene_list=genesDN
      tmpl_name=paste0('dn_',list_names[i])
    }
    geneList <- factor(as.integer(geneNames %in% gene_list))
    names(geneList) <- geneNames
    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    #resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    allRes <- GenTable(GOdata, classicFisher = resultFisher,
                       orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 100)
    fullTerms <- Term(GOTERM[allRes$GO.ID])
    allRes$Full_Term <- fullTerms
    topGO_list[[tmpl_name]]=allRes
    #write.csv(allRes, file = "thalassiosira_mT11_topGO_results.csv", row.names = FALSE)
    showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo='all')
  }
}

wb <- createWorkbook()
# Add each data frame to a sheet
for (sheet_name in names(topGO_list)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, topGO_list[[sheet_name]])
}
# Save workbook
saveWorkbook(wb, "Mflosaquae_DEGs_topGO_4pairs.xlsx", overwrite = TRUE)

######## WGCNA and module function enrichment in GO terms and query pathways ####
detach("package:topGO", unload = TRUE)
detach("package:WGCNA", unload = TRUE)
rm(list=ls())
library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
collectGarbage()
list_mod_trait=list() # Ana-b15-b16 Cylb15-b16
list_mod_topGO=list()
##### 1. load get GOA files of all three genomes
goafile="./geneIDGO/Mflosaque_protID_geneID_GOterm_paired2.tsv"
### 4.1 Load water chemistry metadata and combined metadata with mG
meta=read.csv("../water_chemistry.csv",header = T)
rownames(meta)=gsub("\\.","_",meta[,1])
rownames(meta)=gsub("mG19_", "mT", rownames(meta))
datTraits=meta[,-1]
colnames(datTraits)
datTraits=datTraits[,c(1,8,9,7,11,4,3,2,5,6,10,12)]
### 4.2 WGCNA analysis for two genome 
genefiles=file.path(".",list.files(pattern="rsem_tpm_"))
g=3
## read in gene abundance from mT reads
geneAbd_mt=read.delim(genefiles[g],header = T,sep = "\t", row.names = 1)
colnames(geneAbd_mt)=gsub("\\.","_",colnames(geneAbd_mt))
colnames(geneAbd_mt)=gsub("_genes", "", colnames(geneAbd_mt))
parts <- do.call(rbind, lapply(strsplit(rownames(geneAbd_mt), "_"), function(x) {
  part1 <- paste(x[1:2], collapse = "_")
  part2 <- paste(x[-c(1:2)], collapse = "_")
  c(part1, part2)
}))
rownames(geneAbd_mt)=parts[,1]

otu=as.data.frame(t(geneAbd_mt[rowSums(geneAbd_mt>0)>10,]))
logotu=log10(t(as.matrix(otu)+1))
all(rownames(otu)==rownames(datTraits))

library(pheatmap)
pdf("mT19_Mflosaquaed_log10tpm_heatmap.pdf",6,6)
pheatmap(logotu, cluster_cols=FALSE, show_rownames=FALSE, main = "log10TPM_mT19_Mfa")
barplot(datTraits$T, names.arg = rownames(datTraits),
        las = 2, col = "lightblue", ylab = "Temperature (°C)", cex.names = 0.8)
dev.off()
## WGCNA analysis starts
gsg = goodSamplesGenes(otu, verbose = 3);
gsg$allOK
### if gsg$allOK is FALSE (as opposed to TRUE), do this block of code and test at the end.
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

## check for outliers
sampleTree = hclust(dist(otu), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,5,2,0))
pdf("mT19_Mflosaquae_wgcna_all-samples_cluster.pdf", width = 5, height = 4);
pt<-plot(sampleTree, main = "mT19_Mflosaquae_tpm", sub="", xlab="Sample", 
         cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8, cex=0.8, ylim=c(0,0.5)
)
dev.off()

# Plot a line to show the cut
abline(h = 6e+05, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 6e+05, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = otu[keepSamples, ] # datExpr = sqrt(otu[keepSamples, ])
datTraits=datTraits[keepSamples,]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Before we continue with network construction and module detection, 
# we visualize how the meta data relate to the sample dendrogram
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf("Mflosaquae_wgcna_tpm_meta_combined.pdf",5,5)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",
                    cex.colorLabels = 0.85, cex.dendroLabels = 1, 
                    cex.rowText = 0.8,las=2
)
dev.off()
## 2 Automatic network construction and module detection
# 2.1 Choosing the soft-thresholding power: analysis of network topology
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

## One-step network construction and module detection
collectGarbage()
#-: do not import topGO before runing this function below.
power=8
net = blockwiseModules(datExpr, power = power,
                       TOMType = "signed", #minModuleSize = 15,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3);
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05
);
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
# Quantifying module-trait associations
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#rownames(moduleTraitCor)=gsub("ME",paste(sites[4+f], "_",sep = ""), rownames(moduleTraitCor))
listname= "Mflosaquae" #paste(sites[4+f],species[g],sep = "_")
# add module-trait association to list
#list_mod_trait[[listname]]=moduleTraitCor
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "|",
                   signif(moduleTraitPvalue, 1), sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
pdf("mT19_Mflosaquae_wgcna_modules_trait_corr_heatmap_4x5.pdf",4,5)
library(colorspace)
par(mar = c(5, 5, 0, 0)+0.1);
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = diverging_hcl(30, palette='Blue-Red3'), #blueWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.25, cex.lab.x = 0.4, cex.lab.y = 0.5,
               zlim = c(-1,1),
               xLabelsAngle = 90,
               main = paste("mT19_Mflosaquae_metadata correlation"));
pheatmap(
  mat = moduleTraitCor,
  color = diverging_hcl(30, palette='Blue-Red3'),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 4,
  fontsize_col = 5,
  main = "Mfa",
  angle_col = 90,
  border_color = "white"
)
dev.off()

# pick modules and do topGO analysis
library(topGO)
library(dplyr)
goafile="./geneIDGO/Mflosaque_protID_geneID_GOterm_paired2.tsv"
geneID2GO<-readMappings(file = goafile)
allGenes <- unique(names(geneID2GO))

# modules associated with temperature
#pdf("wgcanMT_Mflosaque_3module-T_positiveCorr_heatmap.pdf",8,8)
#outfile="wgcanMT_Mflosaque_3module-T_positiveCorr_topGO.tsv"
#picked_module = c("midnightblue", "darkmagenta","darkturquoise")
# pdf("wgcanMT_Mflosaque_4module-T_negativeCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mflosaque_4module-T_negativeCorr_topGO.tsv"
# picked_module =c("paleturquoise","turquoise","cyan", "royalblue")

# # modules associated with PC
# pdf("wgcanMT_Mflosaque_2module-PC_postiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mflosaque_2module-PC_postiveCorr_topGO.tsv"
# picked_module = c("ivory","steelblue")

# # modules associated with NH4-N
# pdf("wgcanMT_Mflosaque_2module-NH4-N_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mflosaque_1module-NH4-N_positiveCorr_topGO.tsv"
# picked_module = c("plum1","turquoise")

# # modules associated with NO3-N
# pdf("wgcanMT_Mflosaque_3module-NO3-N_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mflosaque_3module-NO3-N_positiveCorr_topGO.tsv"
# picked_module = c("salmon", "royalblue", "blue")
# pdf("wgcanMT_Mflosaque_2module-NO3-N_negativeCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mflosaque_2module-NO3-N_negativeCorr_topGO.tsv"
# picked_module = c("midnightblue", "darkolivegreen")

#modules associated with NO2-N
# pdf("wgcanMT_Mflosaque_3module-NO2-N_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mflosaque_3module-NO2-N_positiveCorr_topGO.tsv"
# picked_module = c("lightyellow","sienna3","darkturquoise")

# modules associated with Pi
pdf("wgcanMT_Mflosaque_2module-Pi_positiveCorr_heatmap.pdf",8,8)
outfile="wgcanMT_Mflosaque_2module-Pi_positiveCorr_topGO.tsv"
picked_module = c("orangered4","purple","lightsteelblue1")

module_dfList <- list()
term_dfList <- list()
for (i in 1:length(picked_module)){
  module_genes = names(datExpr)[moduleColors == picked_module[i]]
  module_df=datExpr[,colnames(datExpr) %in% module_genes]
  module_dfList[[i]]=module_df
  old_names <- colnames(module_df)
  name_map <- setNames(parts[,2], parts[,1])
  new_names <- ifelse(old_names %in% names(name_map), name_map[old_names], old_names)
  colnames(module_df) <- new_names
  pheatmap(log10(t(module_df+1)), cluster_cols=FALSE, fontsize_row = 4, fontsize_col = 8,
           show_rownames = TRUE, main=paste(picked_module[i],"log10value",sep=" "))
  # topGO enrichment
  moduleGenes_factor <- as.integer(allGenes %in% module_genes)
  #moduleGenes_factor <- factor(as.integer(allGenes %in% module_genes))
  names(moduleGenes_factor)=allGenes
  moduleGenes_factor <- factor(moduleGenes_factor)
  table(moduleGenes_factor)
  GOdata <- new("topGOdata", ontology = "BP", allGenes = moduleGenes_factor,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata, classicFisher = resultFisher,
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 50)
  fullTerms <- Term(GOTERM[allRes$GO.ID])
  allRes$Full_Term <- fullTerms
  allRes=allRes[allRes$classicFisher <= 0.09,]
  if(nrow(allRes)>0){
    term_dfList[[i]]=allRes
    showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo='all', useFullNames=TRUE)
    mtext(picked_module[i], side = 3, line = -0.2, adj = 0.1, cex = 0.8, font = 1)
  }
}
is_empty <- function(x) length(x) == 0
if (!all(sapply(term_dfList, is_empty))) {
  term_dfList <- Filter(function(x) length(x) > 0, term_dfList)
  combined_term_df <- do.call(rbind,term_dfList)
  write.table(combined_term_df, outfile, sep = "\t", row.names = F, quote = F)
  message("list Empty elements removed and written to file.")
}else{
  message("list empty")
}
combined_module_df <- do.call(cbind, module_dfList)
old_names <- colnames(combined_module_df)
name_map <- setNames(parts[,2], parts[,1])
new_names <- ifelse(old_names %in% names(name_map), name_map[old_names], old_names)
colnames(combined_module_df) <- new_names
pheatmap(t(log10(combined_module_df + 1)), main = outfile, cluster_cols=FALSE, 
         show_rownames = TRUE, fontsize_row = 4, fontsize_col = 8)
dev.off()

################################################################################
######## DESeq2 & WGCNA ########
##_____ 1. DESeq2 + topGO _____
##--Microcystis panniformis
library(pheatmap)
library(RColorBrewer)
library(topGO)
library(DESeq2)
library(openxlsx)
# count data for DESeq2
rm(list=ls())
df=read.csv("rsem_readCount_Mpanniformis.tsv", sep = "\t",header = T, row.names = 1)
df=df[rowSums(df==0) <7,]
colnames(df)=gsub(".genes$","",colnames(df))
colnames(df)=gsub("\\.","_",colnames(df))
parts <- do.call(rbind, lapply(strsplit(rownames(df), "_"), function(x) {
  part1 <- paste(x[1:2], collapse = "_")
  part2 <- paste(x[-c(1:2)], collapse = "_")
  c(part1, part2)
}))
colnames(parts) <- c("locusID", "geneName")
parts=as.data.frame(parts)
rownames(df)=parts$locusID
dim(df)

list_up=list_down=list()
list_names=c('mt7v6', 'mt9v6', 'mt11v6', 'mt12v6')
cts0=round(as.matrix(df))
for (ll in 1:length(list_names)){
  if(ll==1){
    cts=cts0[, 1:6]
  } else if(ll==2){
    cts=cts0[, c(1:3,7:9)]
  } else if(ll==3){
    cts=cts0[,c(1:3, 13:15)]
  }else{
    cts=cts0[,c(1:3, 16:18)]
  }
  coldata=data.frame(sample=colnames(cts),condition=gsub("_(1|2|3)","",colnames(cts), perl = TRUE))
  rownames(coldata)=coldata[,1]
  coldata$condition <- factor(coldata$condition)
  
  all(rownames(coldata)==colnames(cts))
  # running DESeq2
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata,
                                design = ~ condition)
  keep <- rowSums(counts(dds)) >= 1
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds)
  resultsNames(dds)
  condName=resultsNames(dds)[2]
  res <- results(dds, name = condName)
  resOrdered <- res[order(res$pvalue),]
  res05 <- results(dds, alpha=0.05)
  res05_order = res05[order(res05$padj), ]
  summary(res05)
  sum(res05$padj < 0.05, na.rm=TRUE)
  resSig <- subset(resOrdered, padj < 0.05)
  list_up[[list_names[ll]]]=rownames(resSig)[resSig$log2FoldChange>0 & resSig$padj<0.05]
  list_down[[list_names[ll]]]=rownames(resSig)[resSig$log2FoldChange<0 & resSig$padj<0.05]
}

topGO_list=list()
goafile="./geneIDGO/Mpanniformis_protID_geneID_GOterm_paired2.tsv"
geneID2GO <- readMappings(file = goafile)
geneNames <- names(geneID2GO)
for (i in 1:length(list_up)){
  #i=1
  genesDN=list_down[[list_names[i]]]
  genesUP=list_up[[list_names[i]]]
  for (j in 1:2){
    if (j==1) {
      gene_list=genesUP
      tmpl_name=paste0('up_',list_names[i])
    } else {
      gene_list=genesDN
      tmpl_name=paste0('dn_',list_names[i])
    }
    geneList <- factor(as.integer(geneNames %in% gene_list))
    names(geneList) <- geneNames
    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    #resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    allRes <- GenTable(GOdata, classicFisher = resultFisher,
                       orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 100)
    fullTerms <- Term(GOTERM[allRes$GO.ID])
    allRes$Full_Term <- fullTerms
    topGO_list[[tmpl_name]]=allRes
    #write.csv(allRes, file = "thalassiosira_mT11_topGO_results.csv", row.names = FALSE)
    showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo='all')
  }
}

wb <- createWorkbook()
# Add each data frame to a sheet
for (sheet_name in names(topGO_list)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, topGO_list[[sheet_name]])
}
# Save workbook
saveWorkbook(wb, "Mpanniformis_DEGs_topGO_4pairs.xlsx", overwrite = TRUE)

######## WGCNA and module function enrichment in GO terms and query pathways ####
detach("package:topGO", unload = TRUE)
detach("package:WGCNA", unload = TRUE)
rm(list=ls())
library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
collectGarbage()
list_mod_trait=list() # Ana-b15-b16 Cylb15-b16
list_mod_topGO=list()
##### 1. load get GOA files of all three genomes
goafile="./geneIDGO/Mpanniformis_protID_geneID_GOterm_paired2.tsv"
### 4.1 Load water chemistry metadata and combined metadata with mG
meta=read.csv("../water_chemistry.csv",header = T)
rownames(meta)=gsub("\\.","_",meta[,1])
rownames(meta)=gsub("mG19_", "mT", rownames(meta))
datTraits=meta[,-1]
colnames(datTraits)
datTraits=datTraits[,c(1,8,9,7,11,4,3,2,5,6,10,12)]

### 4.2 WGCNA analysis for two genome 
genefiles=file.path(".",list.files(pattern="rsem_tpm_"))
g=4
## read in gene abundance from mT reads
geneAbd_mt=read.delim(genefiles[g],header = T,sep = "\t", row.names = 1)
colnames(geneAbd_mt)=gsub("\\.","_",colnames(geneAbd_mt))
colnames(geneAbd_mt)=gsub("_genes", "", colnames(geneAbd_mt))
parts <- do.call(rbind, lapply(strsplit(rownames(geneAbd_mt), "_"), function(x) {
  part1 <- paste(x[1:2], collapse = "_")
  part2 <- paste(x[-c(1:2)], collapse = "_")
  c(part1, part2)
}))
rownames(geneAbd_mt)=parts[,1]

otu=as.data.frame(t(geneAbd_mt[rowSums(geneAbd_mt>0)>10,]))
logotu=log10(t(as.matrix(otu)+1))
all(rownames(otu)==rownames(datTraits))

library(pheatmap)
pdf("mT19_Mpanniformis_log10tpm_heatmap.pdf",6,6)
pheatmap(logotu, cluster_cols=FALSE, show_rownames=FALSE, main = "log10TPM_mT19_MAR")
barplot(datTraits$T, names.arg = rownames(datTraits),
        las = 2, col = "lightblue", ylab = "Temperature (°C)", cex.names = 0.8)
dev.off()
## WGCNA analysis starts
gsg = goodSamplesGenes(otu, verbose = 3);
gsg$allOK
### if gsg$allOK is FALSE (as opposed to TRUE), do this block of code and test at the end.
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

## check for outliers
sampleTree = hclust(dist(otu), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,5,2,0))
pdf("mT19_Mpanniformis_wgcna_all-samples_cluster.pdf", width = 5, height = 4);
pt<-plot(sampleTree, main = "mT19_Mflosaquae_tpm", sub="", xlab="Sample", 
         cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8, cex=0.8, ylim=c(0,0.5)
)
dev.off()

# Plot a line to show the cut
abline(h = 6e+05, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 6e+05, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = otu[keepSamples, ] # datExpr = sqrt(otu[keepSamples, ])
datTraits=datTraits[keepSamples,]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Before we continue with network construction and module detection, 
# we visualize how the meta data relate to the sample dendrogram
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf("Mflosaquae_wgcna_tpm_meta_combined.pdf",5,5)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",
                    cex.colorLabels = 0.85, cex.dendroLabels = 1, 
                    cex.rowText = 0.8,las=2
)
dev.off()
## 2 Automatic network construction and module detection
# 2.1 Choosing the soft-thresholding power: analysis of network topology
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

## One-step network construction and module detection
collectGarbage()
#-: do not import topGO before runing this function below.
power=8
net = blockwiseModules(datExpr, power = power,
                       TOMType = "signed", #minModuleSize = 15,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3);
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05
);
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
# Quantifying module-trait associations
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#rownames(moduleTraitCor)=gsub("ME",paste(sites[4+f], "_",sep = ""), rownames(moduleTraitCor))

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "|",
                   signif(moduleTraitPvalue, 1), sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
pdf("mT19_Mpanniformis_wgcna_modules_trait_corr_heatmap_4x5.pdf",4,5)
library(colorspace)
par(mar = c(5, 5, 0, 0)+0.1);
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = diverging_hcl(30, palette='Blue-Red3'), #blueWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.25, cex.lab.x = 0.4, cex.lab.y = 0.5,
               zlim = c(-1,1),
               xLabelsAngle = 90,
               main = paste("mT19_Mpanniformis_metadata correlation"));
pheatmap(
  mat = moduleTraitCor,
  color = diverging_hcl(30, palette='Blue-Red3'),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 4,
  fontsize_col = 5,
  main = "Mpa",
  angle_col = 90,
  border_color = "white"
)
dev.off()

# pick modules and do topGO analysis
library(topGO)
library(dplyr)
goafile="./geneIDGO/Mpanniformis_protID_geneID_GOterm_paired2.tsv"
geneID2GO<-readMappings(file = goafile)
allGenes <- unique(names(geneID2GO))

# modules associated with temperature
# pdf("wgcanMT_Mpanniformis_3module-T_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mpanniformis_3module-T_positiveCorr_topGO.tsv"
# picked_module = c("royalblue", "saddlebrown","violet")
# pdf("wgcanMT_Mpanniformis_5module-T_negativeCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mpanniformis_5module-T_negativeCorr_topGO.tsv"
# picked_module =c("green","purple","midnightblue", "black","salmon")

# # modules associated with PC
# pdf("wgcanMT_Mpanniformis_2module-PC_postiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mpanniformis_2module-PC_postiveCorr_topGO.tsv"
# picked_module = c("darkorange2","white")
# 
# # modules associated with NH4-N
# pdf("wgcanMT_Mpanniformis_3module-NH4-N_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mpanniformis_3module-NH4-N_positiveCorr_topGO.tsv"
# picked_module = c("black","salmon","green")

# # modules associated with NO3-N
# pdf("wgcanMT_Mpanniformis_3module-NO3-N_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mpanniformis_3module-NO3-N_positiveCorr_topGO.tsv"
# picked_module = c("grey60", "salmon", "purple")
# pdf("wgcanMT_Mpanniformis_2module-NO3-N_negativeCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mpanniformis_2module-NO3-N_negativeCorr_topGO.tsv"
# picked_module = c("midnightblue", "darkolivegreen")

#modules associated with NO2-N
# pdf("wgcanMT_Mpanniformis_4module-NO2-N_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mpanniformis_4module-NO2-N_positiveCorr_topGO.tsv"
# picked_module = c("ivory","bisque4","darkgrey","darkred")

# # modules associated with Pi
# pdf("wgcanMT_Mpanniformis_3module-Pi_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mpanniformis_3module-Pi_positiveCorr_topGO.tsv"
# picked_module = c("yellowgreen","mediumpurple3","yellow")

module_dfList <- list()
term_dfList <- list()
for (i in 1:length(picked_module)){
  module_genes = names(datExpr)[moduleColors == picked_module[i]]
  module_df=datExpr[,colnames(datExpr) %in% module_genes]
  module_dfList[[i]]=module_df
  old_names <- colnames(module_df)
  name_map <- setNames(parts[,2], parts[,1])
  new_names <- ifelse(old_names %in% names(name_map), name_map[old_names], old_names)
  colnames(module_df) <- new_names
  pheatmap(log10(t(module_df+1)), cluster_cols=FALSE, fontsize_row = 4, fontsize_col = 8,
           show_rownames = TRUE, main=paste(picked_module[i],"log10value",sep=" "))
  # topGO enrichment
  moduleGenes_factor <- as.integer(allGenes %in% module_genes)
  #moduleGenes_factor <- factor(as.integer(allGenes %in% module_genes))
  names(moduleGenes_factor)=allGenes
  moduleGenes_factor <- factor(moduleGenes_factor)
  table(moduleGenes_factor)
  GOdata <- new("topGOdata", ontology = "BP", allGenes = moduleGenes_factor,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata, classicFisher = resultFisher,
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 50)
  fullTerms <- Term(GOTERM[allRes$GO.ID])
  allRes$Full_Term <- fullTerms
  allRes=allRes[allRes$classicFisher <= 0.09,]
  if(nrow(allRes)>0){
    term_dfList[[i]]=allRes
    showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo='all', useFullNames=TRUE)
    mtext(picked_module[i], side = 3, line = -0.2, adj = 0.1, cex = 0.8, font = 1)
  }
}
is_empty <- function(x) length(x) == 0
if (!all(sapply(term_dfList, is_empty))) {
  term_dfList <- Filter(function(x) length(x) > 0, term_dfList)
  combined_term_df <- do.call(rbind,term_dfList)
  write.table(combined_term_df, outfile, sep = "\t", row.names = F, quote = F)
  message("list Empty elements removed and written to file.")
}else{
  message("list empty")
}
combined_module_df <- do.call(cbind, module_dfList)
old_names <- colnames(combined_module_df)
name_map <- setNames(parts[,2], parts[,1])
new_names <- ifelse(old_names %in% names(name_map), name_map[old_names], old_names)
colnames(combined_module_df) <- new_names
pheatmap(t(log10(combined_module_df + 1)), main = outfile, cluster_cols=FALSE, 
         show_rownames = TRUE, fontsize_row = 4, fontsize_col = 8)
dev.off()

################################################################################
######## DESeq2 & WGCNA ########
##_____ 1. DESeq2 + topGO _____
##--Microcystis wesenbergii--
library(pheatmap)
library(RColorBrewer)
library(topGO)
library(DESeq2)
library(openxlsx)
# count data for DESeq2
rm(list=ls())
df=read.csv("rsem_readCount_Mwesenbergii.tsv", sep = "\t",header = T, row.names = 1)
df=df[rowSums(df==0) <7,]
colnames(df)=gsub(".genes$","",colnames(df))
colnames(df)=gsub("\\.","_",colnames(df))
parts <- do.call(rbind, lapply(strsplit(rownames(df), "_"), function(x) {
  part1 <- paste(x[1:2], collapse = "_")
  part2 <- paste(x[-c(1:2)], collapse = "_")
  c(part1, part2)
}))
colnames(parts) <- c("locusID", "geneName")
parts=as.data.frame(parts)
rownames(df)=parts$locusID
dim(df)

list_up=list_down=list()
list_names=c('mt7v6', 'mt9v6', 'mt11v6', 'mt12v6')
cts0=round(as.matrix(df))
for (ll in 1:length(list_names)){
  if(ll==1){
    cts=cts0[, 1:6]
  } else if(ll==2){
    cts=cts0[, c(1:3,7:9)]
  } else if(ll==3){
    cts=cts0[,c(1:3, 13:15)]
  }else{
    cts=cts0[,c(1:3, 16:18)]
  }
  coldata=data.frame(sample=colnames(cts),condition=gsub("_(1|2|3)","",colnames(cts), perl = TRUE))
  rownames(coldata)=coldata[,1]
  coldata$condition <- factor(coldata$condition)
  
  all(rownames(coldata)==colnames(cts))
  # running DESeq2
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata,
                                design = ~ condition)
  keep <- rowSums(counts(dds)) >= 1
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds)
  resultsNames(dds)
  condName=resultsNames(dds)[2]
  res <- results(dds, name = condName)
  resOrdered <- res[order(res$pvalue),]
  res05 <- results(dds, alpha=0.05)
  res05_order = res05[order(res05$padj), ]
  summary(res05)
  sum(res05$padj < 0.05, na.rm=TRUE)
  resSig <- subset(resOrdered, padj < 0.05)
  list_up[[list_names[ll]]]=rownames(resSig)[resSig$log2FoldChange>0 & resSig$padj<0.05]
  list_down[[list_names[ll]]]=rownames(resSig)[resSig$log2FoldChange<0 & resSig$padj<0.05]
}

topGO_list=list()
goafile="./geneIDGO/Mwesenbergii_protID_geneID_GOterm_paired2.tsv"
geneID2GO <- readMappings(file = goafile)
geneNames <- names(geneID2GO)
for (i in 1:length(list_up)){
  #i=1
  genesDN=list_down[[list_names[i]]]
  genesUP=list_up[[list_names[i]]]
  for (j in 1:2){
    if (j==1) {
      gene_list=genesUP
      tmpl_name=paste0('up_',list_names[i])
    } else {
      gene_list=genesDN
      tmpl_name=paste0('dn_',list_names[i])
    }
    geneList <- factor(as.integer(geneNames %in% gene_list))
    names(geneList) <- geneNames
    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                  annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    #resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    allRes <- GenTable(GOdata, classicFisher = resultFisher,
                       orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 100)
    fullTerms <- Term(GOTERM[allRes$GO.ID])
    allRes$Full_Term <- fullTerms
    topGO_list[[tmpl_name]]=allRes
    #write.csv(allRes, file = "thalassiosira_mT11_topGO_results.csv", row.names = FALSE)
    showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo='all')
  }
}

wb <- createWorkbook()
# Add each data frame to a sheet
for (sheet_name in names(topGO_list)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, topGO_list[[sheet_name]])
}
# Save workbook
saveWorkbook(wb, "Mwesenbergii_DEGs_topGO_4pairs.xlsx", overwrite = TRUE)

######## WGCNA and module function enrichment in GO terms and query pathways ####
detach("package:topGO", unload = TRUE)
detach("package:WGCNA", unload = TRUE)
rm(list=ls())
library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
collectGarbage()
list_mod_trait=list() # Ana-b15-b16 Cylb15-b16
list_mod_topGO=list()
##### 1. load get GOA files of all three genomes
goafile="./geneIDGO/Mwesenbergii_protID_geneID_GOterm_paired2.tsv"
### 4.1 Load water chemistry metadata and combined metadata with mG
meta=read.csv("../water_chemistry.csv",header = T)
rownames(meta)=gsub("\\.","_",meta[,1])
rownames(meta)=gsub("mG19_", "mT", rownames(meta))
datTraits=meta[,-1]
colnames(datTraits)
datTraits=datTraits[,c(1,8,9,7,11,4,3,2,5,6,10,12)]

### 4.2 WGCNA analysis for two genome 
genefiles=file.path(".",list.files(pattern="rsem_tpm_"))
g=5
## read in gene abundance from mT reads
geneAbd_mt=read.delim(genefiles[g],header = T,sep = "\t", row.names = 1)
colnames(geneAbd_mt)=gsub("\\.","_",colnames(geneAbd_mt))
colnames(geneAbd_mt)=gsub("_genes", "", colnames(geneAbd_mt))
parts <- do.call(rbind, lapply(strsplit(rownames(geneAbd_mt), "_"), function(x) {
  part1 <- paste(x[1:2], collapse = "_")
  part2 <- paste(x[-c(1:2)], collapse = "_")
  c(part1, part2)
}))
rownames(geneAbd_mt)=parts[,1]

otu=as.data.frame(t(geneAbd_mt[rowSums(geneAbd_mt>0)>10,]))
logotu=log10(t(as.matrix(otu)+1))
all(rownames(otu)==rownames(datTraits))

library(pheatmap)
pdf("mT19_Mwesenbergii_log10tpm_heatmap.pdf",6,6)
pheatmap(logotu, cluster_cols=FALSE, show_rownames=FALSE, main = "log10TPM_mT19_MAR")
barplot(datTraits$T, names.arg = rownames(datTraits),
        las = 2, col = "lightblue", ylab = "Temperature (°C)", cex.names = 0.8)
dev.off()
## WGCNA analysis starts
gsg = goodSamplesGenes(otu, verbose = 3);
gsg$allOK
### if gsg$allOK is FALSE (as opposed to TRUE), do this block of code and test at the end.
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

## check for outliers
sampleTree = hclust(dist(otu), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,5,2,0))
pdf("mT19_Mwesenbergii_wgcna_all-samples_cluster.pdf", width = 5, height = 4);
pt<-plot(sampleTree, main = "mT19_Mflosaquae_tpm", sub="", xlab="Sample", 
         cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8, cex=0.8, ylim=c(0,0.5)
)
dev.off()

# Plot a line to show the cut
abline(h = 6e+05, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 6e+05, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = otu[keepSamples, ] # datExpr = sqrt(otu[keepSamples, ])
datTraits=datTraits[keepSamples,]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Before we continue with network construction and module detection, 
# we visualize how the meta data relate to the sample dendrogram
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf("Mwesenbergii_wgcna_tpm_meta_combined.pdf",5,5)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",
                    cex.colorLabels = 0.85, cex.dendroLabels = 1, 
                    cex.rowText = 0.8,las=2
)
dev.off()
## 2 Automatic network construction and module detection
# 2.1 Choosing the soft-thresholding power: analysis of network topology
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

## One-step network construction and module detection
collectGarbage()
#-: do not import topGO before runing this function below.
power=8
net = blockwiseModules(datExpr, power = power,
                       TOMType = "signed", #minModuleSize = 15,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3);
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05
);
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
# Quantifying module-trait associations
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#rownames(moduleTraitCor)=gsub("ME",paste(sites[4+f], "_",sep = ""), rownames(moduleTraitCor))

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "|",
                   signif(moduleTraitPvalue, 1), sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
pdf("mT19_Mwesenbergii_wgcna_modules_trait_corr_heatmap_4x5.pdf",4,5)
library(colorspace)
par(mar = c(5, 5, 0, 0)+0.1);
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = diverging_hcl(30, palette='Blue-Red3'), #blueWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.25, cex.lab.x = 0.4, cex.lab.y = 0.5,
               zlim = c(-1,1),
               xLabelsAngle = 90,
               main = paste("mT19_Mflosaquae_metadata correlation"));
pheatmap(
  mat = moduleTraitCor,
  color = diverging_hcl(30, palette='Blue-Red3'),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 4,
  fontsize_col = 5,
  main = "Mfa",
  angle_col = 90,
  border_color = "white"
)
dev.off()

# pick modules and do topGO analysis
library(topGO)
library(dplyr)
goafile="./geneIDGO/Mwesenbergii_protID_geneID_GOterm_paired2.tsv"
geneID2GO<-readMappings(file = goafile)
allGenes <- unique(names(geneID2GO))

# modules associated with temperature
# pdf("wgcanMT_Mwesenbergii_4module-T_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mwesenbergii_4module-T_positiveCorr_topGO.tsv"
# picked_module = c("paleturquoise", "orangered4","tan","red")
# pdf("wgcanMT_Mwesenbergii_5module-T_negativeCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mwesenbergii_5module-T_negativeCorr_topGO.tsv"
# picked_module =c("green","yellowgreen","midnightblue", "lightcyan1","skyblue")

# # modules associated with PC
# pdf("wgcanMT_Mwesenbergii_1module-PC_postiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mwesenbergii_1module-PC_postiveCorr_topGO.tsv"
# picked_module = c("darkgreen")
# 
# # modules associated with NH4-N
# pdf("wgcanMT_Mwesenbergii_4module-NH4-N_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mwesenbergii_4module-NH4-N_positiveCorr_topGO.tsv"
# picked_module = c("midnightblue","purple","green","yellowgreen")

# # modules associated with NO3-N
# pdf("wgcanMT_Mwesenbergii_4module-NO3-N_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mwesenbergii_4module-NO3-N_positiveCorr_topGO.tsv"
# picked_module = c("skyblue", "lightgreen", "orange","pink")
# pdf("wgcanMT_Mwesenbergii_1module-NO3-N_negativeCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mwesenbergii_1module-NO3-N_negativeCorr_topGO.tsv"
# picked_module = c("paleturquoise")

#modules associated with NO2-N
# pdf("wgcanMT_Mwesenbergii_4module-NO2-N_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mwesenbergii_4module-NO2-N_positiveCorr_topGO.tsv"
# picked_module = c("darkturquoise","red","royalblue","tan")

# # modules associated with Pi
# pdf("wgcanMT_Mwesenbergii_2module-Pi_positiveCorr_heatmap.pdf",8,8)
# outfile="wgcanMT_Mwesenbergii_2module-Pi_positiveCorr_topGO.tsv"
# picked_module = c("darkred","yellow")

module_dfList <- list()
term_dfList <- list()
for (i in 1:length(picked_module)){
  module_genes = names(datExpr)[moduleColors == picked_module[i]]
  module_df=datExpr[,colnames(datExpr) %in% module_genes]
  module_dfList[[i]]=module_df
  old_names <- colnames(module_df)
  name_map <- setNames(parts[,2], parts[,1])
  new_names <- ifelse(old_names %in% names(name_map), name_map[old_names], old_names)
  colnames(module_df) <- new_names
  pheatmap(log10(t(module_df+1)), cluster_cols=FALSE, fontsize_row = 4, fontsize_col = 8,
           show_rownames = TRUE, main=paste(picked_module[i],"log10value",sep=" "))
  # topGO enrichment
  moduleGenes_factor <- as.integer(allGenes %in% module_genes)
  #moduleGenes_factor <- factor(as.integer(allGenes %in% module_genes))
  names(moduleGenes_factor)=allGenes
  moduleGenes_factor <- factor(moduleGenes_factor)
  table(moduleGenes_factor)
  GOdata <- new("topGOdata", ontology = "BP", allGenes = moduleGenes_factor,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata, classicFisher = resultFisher,
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 50)
  fullTerms <- Term(GOTERM[allRes$GO.ID])
  allRes$Full_Term <- fullTerms
  allRes=allRes[allRes$classicFisher <= 0.09,]
  if(nrow(allRes)>0){
    term_dfList[[i]]=allRes
    showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo='all', useFullNames=TRUE)
    mtext(picked_module[i], side = 3, line = -0.2, adj = 0.1, cex = 0.8, font = 1)
  }
}
is_empty <- function(x) length(x) == 0
if (!all(sapply(term_dfList, is_empty))) {
  term_dfList <- Filter(function(x) length(x) > 0, term_dfList)
  combined_term_df <- do.call(rbind,term_dfList)
  write.table(combined_term_df, outfile, sep = "\t", row.names = F, quote = F)
  message("list Empty elements removed and written to file.")
}else{
  message("list empty")
}
combined_module_df <- do.call(cbind, module_dfList)
old_names <- colnames(combined_module_df)
name_map <- setNames(parts[,2], parts[,1])
new_names <- ifelse(old_names %in% names(name_map), name_map[old_names], old_names)
colnames(combined_module_df) <- new_names
pheatmap(t(log10(combined_module_df + 1)), main = outfile, cluster_cols=FALSE, 
         show_rownames = TRUE, fontsize_row = 4, fontsize_col = 8)
dev.off()

################################################################################
######## DESeq2 & WGCNA ########
##_____ 1. DESeq2 + topGO _____
##--TA02--
library(pheatmap)
library(RColorBrewer)
library(topGO)
rm(list = ls())

epi_mTList=list()
genefiles=file.path(".",list.files(pattern="rsem_tpm_"))
genefiles
g=8
## read in gene abundance from mT reads
df=read.delim(genefiles[g],header = T,sep = "\t", row.names = 1)
rownames(df)=gsub("_\\d+$","",rownames(df))
#rownames(df)=gsub("fig\\|34072.190.peg.", "TA02_", rownames(df))
rownames(df)=gsub("fig\\|32033.31.peg.", "TA07_", rownames(df))
#rownames(df)=gsub("fig\\|41275.94.peg.", "TA08_",rownames(df))
df=df[rowSums(df==0) <10,]
colnames(df)=gsub(".genes$","",colnames(df))
colnames(df)=gsub("\\.","_",colnames(df))

dim(df)
dff=df

goa_files=list.files("./geneIDGO/", pattern = "geneID_GOterm_paired\\.tsv")
goa_files
n=7
path <- file.path("geneIDGO", goa_files[n])
ndf=read.delim(path, header = F, sep = "\t")
#ndf[,1]=gsub("fig\\|34072.190.peg.", "TA02_", ndf[,1])
ndf[,1]=gsub("fig\\|32033.31.peg.", "TA07_", ndf[,1])
ndf[,1]=gsub("fig\\|41275.94.peg.", "TA08_", ndf[,1])

ndf$gname=gsub(" .*$", "", ndf$V3,)
ndf$gname2=ifelse(grepl("^[A-Z]",ndf$gname), ndf$V1, paste(ndf$V1, ndf$gname, sep = ":"))
ndf$fullname=paste(ndf$gname2, ndf$V4,sep = "|")
ndf[grepl("tuf|groL|rpsJ", ndf$fullname),]
#ndf[360,7]='tuf:Elongation factor Tu (EF-Tu) (EC 3.6.5.3)'
#ndf[450,7]='tuf2:Elongation factor Tu (EF-Tu) (EC 3.6.5.3)'

matched_idx <- match(rownames(dff), ndf$V1)
new_names <- ifelse(!is.na(matched_idx), ndf$fullname[matched_idx], rownames(dff))
rownames(dff) <- new_names

epi_mTList[['TA02']]=dff

logdff=log10(as.matrix(dff)+1)

pdf("mT19_TA07_log10tpm_heatmap.pdf",8,12)
pheatmap(logdff, cluster_cols=FALSE, show_rownames=TRUE, main = "log10TPM_mT19",
         fontsize_row=5, fontsize_col = 6)
dev.off()

geneID2GO <- readMappings(file = "./geneIDGO/ta08_geneID_GOterm_paired2.tsv")
geneID2GO <- readMappings(file = "./geneIDGO/ta02_geneID_GOterm_paired2.tsv")
geneID2GO <- readMappings(file = "./geneIDGO/ta07_geneID_GOterm_paired2.tsv")
geneNames <- names(geneID2GO)

gene_list=as.character(rownames(df))
geneList <- factor(as.integer(geneNames %in% gene_list))
names(geneList) <- geneNames
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
            annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
allRes <- GenTable(GOdata, classicFisher = resultFisher,
         orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 50)
fullTerms <- Term(GOTERM[allRes$GO.ID])
allRes$Full_Term <- fullTerms
write.csv(allRes, file = "ta08_mT19_topGO_results.csv", row.names = FALSE)
write.csv(allRes, file = "ta02_mT19_topGO_results.csv", row.names = FALSE)
write.csv(allRes, file = "ta07_mT19_topGO_results.csv", row.names = FALSE)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo='all')

################################################################################
######## DESeq2 & WGCNA ########
##_____ 1. DESeq2 + topGO _____
##--TF02-TF05--
library(pheatmap)
library(RColorBrewer)
library(topGO)
rm(list = ls())
genefiles=file.path(".",list.files(pattern="rsem_tpm_"))
genefiles
g=10
## read in gene abundance from mT reads
df=read.delim(genefiles[g],header = T,sep = "\t", row.names = 1)
rownames(df)=gsub("_\\d+$","",rownames(df))
df=df[rowSums(df==0) <10,]
colnames(df)=gsub(".genes$","",colnames(df))
colnames(df)=gsub("\\.","_",colnames(df))
rownames(df)=gsub("fig\\|80864.49.peg.", "TF02_",rownames(df))

dim(df)
dff=df

goa_files=list.files("./geneIDGO/", pattern = "geneID_GOterm_paired\\.tsv")
goa_files
n=9
path <- file.path("geneIDGO", goa_files[n])
ndf=read.delim(path, header = F, sep = "\t")
ndf[,1]=gsub("fig\\|80864.49.peg.", "TF02_", ndf[,1])
ndf$gname=gsub(" .*$", "", ndf$V3,)
ndf$fullname=paste(ndf$gname, ndf$V4, sep = ":")

ndf[grepl("tuf", ndf$fullname),]
ndf[2057,7]='tuf:Elongation factor Tu (EF-Tu) (EC 3.6.5.3)'
ndf[2576,7]='tuf2:Elongation factor Tu (EF-Tu) (EC 3.6.5.3)'

matched_idx <- match(rownames(dff), ndf$V1)
new_names <- ifelse(!is.na(matched_idx), ndf$fullname[matched_idx], rownames(dff))
rownames(dff) <- new_names

logdff=log10(as.matrix(dff)+1)

#pdf("mT19_TF02_log10tpm_heatmap.pdf",7,4)
pheatmap(logdff, cluster_cols=FALSE, show_rownames=TRUE, main = "log10TPM_mT19",
         fontsize_row=5, fontsize_col = 8)
dev.off()

geneID2GO <- readMappings(file = "./geneIDGO/tf02_geneID_GOterm_paired2.tsv")
geneNames <- names(geneID2GO)

gene_list=as.character(rownames(df))
geneList <- factor(as.integer(geneNames %in% gene_list))
names(geneList) <- geneNames
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 50)
fullTerms <- Term(GOTERM[allRes$GO.ID])
allRes$Full_Term <- fullTerms
write.csv(allRes, file = "tf02_mT19_topGO_results.csv", row.names = FALSE)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo='all')
dev.off()
################################################################################
######## DESeq2 & WGCNA ########
##_____ 1. DESeq2 + topGO _____
##--TW01-TW08--
library(pheatmap)
library(RColorBrewer)
library(topGO)
rm(list = ls())
genefiles=file.path(".",list.files(pattern="rsem_tpm_"))
genefiles
g=13
## read in gene abundance from mT reads
df=read.delim(genefiles[g],header = T,sep = "\t", row.names = 1)
rownames(df)=gsub("_\\d+$","",rownames(df))
df=df[rowSums(df==0) <10,]
colnames(df)=gsub(".genes$","",colnames(df))
colnames(df)=gsub("\\.","_",colnames(df))
#rownames(df)=gsub("fig\\|80864.49.peg.", "TF02_",rownames(df))
rownames(df)=gsub("fig\\|237.436.peg.", "TW01_",rownames(df))
dim(df)
dff=df

goa_files=list.files("./geneIDGO/", pattern = "geneID_GOterm_paired\\.tsv")
goa_files
n=13
path <- file.path("geneIDGO", goa_files[n])
ndf=read.delim(path, header = F, sep = "\t")
ndf[,1]=gsub("fig\\|237.436.peg.", "TW01_", ndf[,1])
ndf$gname=gsub(" .*$", "", ndf$V3,)
ndf$fullname=paste(ndf$gname, ndf$V4, sep = ":")

ndf[grepl("nqrF|glmS", ndf$fullname),]
aa="glmS:Glutamine−−fructose−6−phosphate aminotransferase [isomerizing] (EC 2.6.1.16)"
bb="nqrF:Na(+)−translocating NADH−quinone reductase subunit F (Na(+)−NQR subunit F)"
ndf[437,7]=aa
ndf[1278,7]=bb

matched_idx <- match(rownames(dff), ndf$V1)
new_names <- ifelse(!is.na(matched_idx), ndf$fullname[matched_idx], rownames(dff))
rownames(dff) <- new_names

logdff=log10(as.matrix(dff)+1)

#pdf("mT19_TW01_log10tpm_heatmap.pdf",10,12)
pheatmap(logdff, cluster_cols=FALSE, show_rownames=TRUE, main = "log10TPM_mT19",
         fontsize_row=5, fontsize_col = 8)
dev.off()

geneID2GO <- readMappings(file = "./geneIDGO/tw01_geneID_GOterm_paired2.tsv")
geneNames <- names(geneID2GO)

gene_list=as.character(rownames(df))
geneList <- factor(as.integer(geneNames %in% gene_list))
names(geneList) <- geneNames
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 50)
fullTerms <- Term(GOTERM[allRes$GO.ID])
allRes$Full_Term <- fullTerms
write.csv(allRes, file = "tw01_mT19_topGO_results.csv", row.names = FALSE)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo='all')

################################################################################
######## heatmap for all epibiotic bacteria ########
rm(list = ls())
# Load libraries
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)

# Get file list
files <- list.files(pattern = "t.*_mT19_topGO_results.csv")
names(files) <- tools::file_path_sans_ext(files)

# Read and filter each file
term_lists <- lapply(files, function(f) {
  df <- read.csv(f)
  sig_terms <- subset(df, classicFisher < 0.06)$Full_Term
  unique(sig_terms)
})

# Remove entries with no significant terms (like ta08)
term_lists <- term_lists[lengths(term_lists) > 0]

# Proceed only if at least 2 non-empty files remain
if (length(term_lists) < 2) stop("Not enough non-empty files with significant terms.")

# All unique terms
all_terms <- unique(unlist(term_lists))

# Binary matrix of presence/absence
binary_matrix <- sapply(term_lists, function(x) as.numeric(all_terms %in% x))
rownames(binary_matrix) <- all_terms
colnames(binary_matrix) <- names(term_lists)

# Nature-style blue color palette
nature_colors <- colorRampPalette(c("#FFFFFF", "#1B3B6F"))(2)

# ---- HEATMAP (PDF) ----
pdf("heatmap_topGO_terms.pdf", width = 4, height = 4)
pheatmap(binary_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 5,
         fontsize_row = 4.5,
         fontsize_col = 5,
         legend = FALSE,
         color = nature_colors,
         border_color = "grey90",
         main = "TopGO Terms (p < 0.06)")
dev.off()

# ---- VENN DIAGRAM (PDF) ----
pdf("venn_topGO_terms.pdf", width = 3, height = 3)
venn.plot <- venn.diagram(
  x = term_lists,
  category.names = names(term_lists),
  filename = NULL,
  output = TRUE,
  cat.cex = 0.4,
  cex = 0.4,
  margin = 0.1,
  fill = brewer.pal(min(5, length(term_lists)), "Set2"),
  col = "grey30",
  lwd = 0.3  # Half-width border lines
)
grid::grid.draw(venn.plot)
dev.off()

#####################################################################
# make bubble plot for GO terms among TA TF TW
rm(list = ls())
library(tidyverse)

#---- Settings ----
font_size <- 6
point_size_range <- c(3, 10)  # bubble size scaling

#---- Load selected files ----
files <- c("tw01_mT19_topGO_results.csv",
           "tf02_mT19_topGO_results.csv",
           "ta02_mT19_topGO_results.csv",
           "ta08_mT19_topGO_results.csv")

# Short labels
file_labels <- c("TW01","TF02","TA02","TA08")
names(file_labels) <- files

# Read and combine
df <- files %>%
  map_df(~read.csv(.x) %>%
           filter(as.numeric(classicFisher) < 0.04) %>%   # stricter filter
           mutate(Source = file_labels[.x],
                  logFisher = -log10(as.numeric(classicFisher))),
         .id = "file_id")

#---- Bubble Plot ----
p <- ggplot(df, aes(x = Source, y = Full_Term)) +
  geom_point(aes(size = Significant, color = logFisher)) +
  scale_size(range = point_size_range) +
  scale_color_gradient(low = "steelblue", high = "coral") +
  theme_bw(base_size = font_size) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color="black"),
        axis.text.y = element_text(color="black"),
        axis.title = element_text(color="black")) +
  labs(x = "Sample", y = "GO term",
       size = "Significant genes",
       color = expression(-log[10]~Fisher))

print(p)
## pheatmap
rm(list = ls())
library(tidyverse)
library(pheatmap)

#---- Settings ----
font_size <- 6   # font size for heatmap labels
cutoff <- 0.04    # Fisher threshold

#---- Load selected files ----
files <- c("tw01_mT19_topGO_results.csv",
           "tf02_mT19_topGO_results.csv",
           "ta02_mT19_topGO_results.csv",
           "ta08_mT19_topGO_results.csv")

# Short labels
file_labels <- c("TW01","TF02","TA02","TA08")
names(file_labels) <- files

# Read and combine
df <- files %>%
  map_df(~read.csv(.x) %>%
           filter(as.numeric(classicFisher) < cutoff) %>%
           mutate(Source = file_labels[.x],
                  logFisher = -log10(as.numeric(classicFisher))),
         .id = "file_id")

#---- Reshape into matrix ----
mat <- df %>%
  select(Full_Term, Source, logFisher) %>%
  pivot_wider(names_from = Source, values_from = logFisher, values_fill = 0) %>%
  column_to_rownames("Full_Term") %>%
  as.matrix()

#---- Heatmap ----
pdf("epibiotic_TA02_TA08_TF02_TW01_GOterms.pdf",5,5)
pheatmap(mat,
         color = colorRampPalette(c("white", "coral"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize = font_size,
         fontsize_row = font_size*0.7,
         fontsize_col = font_size,
         border_color = NA,
         main = "-log10(Fisher) enrichment (GO terms)")
dev.off()

################################################################################
# ========================================================
# FINAL UPDATED R Script: Month-specific top 200 genes + heatmaps
# 
# Changes applied (exactly as requested):
#   - Gene name trimming now handles duplicates correctly:
#       • First try the short name you requested (last part for 3 parts, last two joined for 4 parts)
#       • If the short name is NOT unique (duplicate row.names error), fall back to the FULL original gene_id
#       • This guarantees no duplicate row names while keeping names as short as possible
#   - All other features unchanged:
#       • Top 200 (or all significant) genes per month via DESeq2
#       • Combined heatmap for 07-23 + 09-10
#       • All heatmaps saved as PDF (8 × 12 inches)
# ========================================================
rm(list = ls())
library(DESeq2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# 2. Read the raw count file
counts_raw <- read.table("rsem_readCount_taihu98.tsv", 
                         header = TRUE, sep = "\t", 
                         row.names = 1, check.names = FALSE)

# 3. GENE NAME TRIMMING with duplicate handling
original_ids <- rownames(counts_raw)

gene_labels <- sapply(original_ids, function(id) {
  parts <- strsplit(id, "_", fixed = TRUE)[[1]]
  n <- length(parts)
  
  if (n == 3) {
    # 3 parts → take the last part only
    parts[3]
  } else if (n == 4) {
    # 4 parts → take the last two parts joined by "_"
    paste(parts[3], parts[4], sep = "_")
  } else {
    # Fallback
    paste(tail(parts, 2), collapse = "_")
  }
})

# Fix duplicates: if a trimmed name appears more than once, use the FULL original gene_id
dup_names <- names(table(gene_labels))[table(gene_labels) > 1]
for (dup in dup_names) {
  idx <- which(gene_labels == dup)
  gene_labels[idx] <- original_ids[idx]   # ← use full original ID for duplicates
}

rownames(counts_raw) <- gene_labels

# 4. Clean column names (remove "mT19_" prefix)
colnames(counts_raw) <- sub("^mT19_", "", colnames(counts_raw))

# 5. Parse metadata
meta <- data.frame(sample = colnames(counts_raw), stringsAsFactors = FALSE) %>%
  separate(sample, into = c("month", "day", "bio", "tech"), 
           sep = "-", remove = FALSE, convert = TRUE) %>%
  mutate(
    month = as.character(month),
    timepoint = case_when(
      month == "6"  ~ "06-15",
      month == "7"  ~ "07-23",
      month == "9" & day == "10" ~ "09-10",
      month == "9" & day == "29" ~ "09-29",
      month == "11" ~ "11-3",
      month == "12" ~ "12-4"
    ),
    bio_sample = paste(timepoint, bio, sep = "_")
  )

# Average technical replicates per biological replicate
count_avg_long <- counts_raw %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(meta %>% select(sample, bio_sample), by = "sample") %>%
  group_by(gene, bio_sample) %>%
  summarise(count = mean(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = bio_sample, values_from = count) %>%
  column_to_rownames("gene")

count_matrix <- as.matrix(count_avg_long)

# colData for DESeq2
coldata <- meta %>%
  distinct(bio_sample, timepoint) %>%
  arrange(timepoint) %>%
  column_to_rownames("bio_sample") %>%
  mutate(timepoint = factor(timepoint, 
                            levels = c("06-15", "07-23", "09-10", 
                                       "09-29", "11-3", "12-4")))

# Run DESeq2 once for normalization + VST
dds_full <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                   colData = coldata,
                                   design = ~ timepoint)
dds_full <- DESeq(dds_full)
vsd <- vst(dds_full, blind = TRUE)

# 8. DESeq2 per timepoint → TSV files
timepoints <- levels(coldata$timepoint)
top_results <- list()

for (tp in timepoints) {
  cat("Processing", tp, "...\n")
  
  coldata$group <- ifelse(coldata$timepoint == tp, "this", "other")
  coldata$group <- factor(coldata$group, levels = c("other", "this"))
  
  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                colData = coldata,
                                design = ~ group)
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("group", "this", "other"), 
                 alpha = 0.05, lfcThreshold = 0)
  
  res_ordered <- res[order(res$log2FoldChange, decreasing = TRUE), ]
  sig <- res_ordered[!is.na(res_ordered$padj) & 
                       res_ordered$padj < 0.05 & 
                       res_ordered$log2FoldChange > 0, ]
  
  n_genes <- min(200, nrow(sig))
  if (n_genes == 0) next
  
  top_genes <- rownames(sig)[1:n_genes]
  
  norm_counts <- counts(dds, normalized = TRUE)
  mean_this <- rowMeans(norm_counts[top_genes, coldata$group == "this", drop = FALSE])
  
  out_df <- data.frame(
    gene = top_genes,
    as.data.frame(sig[1:n_genes, ]),
    mean_normalized_expr_this_month = mean_this[1:n_genes]
  )
  
  write.table(out_df, 
              file = paste0("top200_genes_", tp, ".tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  top_results[[tp]] <- list(genes = top_genes, df = out_df)
}

cat("\n=== All TSV files saved! ===\n")

# 9. HEATMAPS — saved as PDF (8 × 12 inches)
pdf_width  <- 6
pdf_height <- 12

# Single-month heatmaps (skip 07-23 and 09-10)
for (tp in timepoints) {
  if (tp %in% c("07-23", "09-10")) next
  if (is.null(top_results[[tp]])) next
  
  genes <- top_results[[tp]]$genes
  mat <- assay(vsd)[genes, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  
  col_order <- order(coldata$timepoint)
  mat_z <- mat_z[, col_order, drop = FALSE]
  
  ann_col <- data.frame(Timepoint = coldata$timepoint[col_order])
  rownames(ann_col) <- colnames(mat_z)
  
  title <- paste("Top genes highly expressed in", tp, 
                 "(row z-score; DESeq2 this vs others)")
  
  pdf(paste0("heatmap_top_genes_", tp, ".pdf"), width = pdf_width, height = pdf_height)
  pheatmap(mat_z,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           annotation_col = ann_col,
           main = title,
           fontsize = 7,
           fontsize_row = 6,
           color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
           border_color = NA)
  dev.off()
  
  cat("PDF heatmap saved for", tp, "\n")
}

# Combined heatmap for 07-23 + 09-10
if (!is.null(top_results[["07-23"]]) && !is.null(top_results[["09-10"]])) {
  genes07 <- top_results[["07-23"]]$genes
  genes09 <- top_results[["09-10"]]$genes
  combined_genes <- unique(c(genes07, genes09))
  
  mat <- assay(vsd)[combined_genes, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  
  col_order <- order(coldata$timepoint)
  mat_z <- mat_z[, col_order, drop = FALSE]
  
  ann_col <- data.frame(Timepoint = coldata$timepoint[col_order])
  rownames(ann_col) <- colnames(mat_z)
  
  title <- "Top genes highly expressed in 07-23 or 09-10 (row z-score; DESeq2 this vs others)"
  
  pdf("heatmap_top_genes_combined_07-23_09-10.pdf", width = pdf_width, height = pdf_height)
  pheatmap(mat_z,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           annotation_col = ann_col,
           main = title,
           fontsize = 7,
           fontsize_row = 6,
           color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
           border_color = NA)
  dev.off()
  
  cat("PDF heatmap saved for combined 07-23 + 09-10\n")
}


## @@ >> top 100 instead of 200 above

# 2. Read the raw count file
counts_raw <- read.table("rsem_readCount_taihu98.tsv", 
                         header = TRUE, sep = "\t", 
                         row.names = 1, check.names = FALSE)

# 3. GENE NAME TRIMMING with duplicate handling
original_ids <- rownames(counts_raw)

gene_labels <- sapply(original_ids, function(id) {
  parts <- strsplit(id, "_", fixed = TRUE)[[1]]
  n <- length(parts)
  
  if (n == 3) {
    parts[3]
  } else if (n == 4) {
    paste(parts[3], parts[4], sep = "_")
  } else {
    paste(tail(parts, 2), collapse = "_")
  }
})

# Fix duplicates: use full original ID when short name would collide
dup_names <- names(table(gene_labels))[table(gene_labels) > 1]
for (dup in dup_names) {
  idx <- which(gene_labels == dup)
  gene_labels[idx] <- original_ids[idx]
}

rownames(counts_raw) <- gene_labels

# 4. Clean column names
colnames(counts_raw) <- sub("^mT19_", "", colnames(counts_raw))

# 5. Parse metadata
meta <- data.frame(sample = colnames(counts_raw), stringsAsFactors = FALSE) %>%
  separate(sample, into = c("month", "day", "bio", "tech"), 
           sep = "-", remove = FALSE, convert = TRUE) %>%
  mutate(
    month = as.character(month),
    timepoint = case_when(
      month == "6"  ~ "06-15",
      month == "7"  ~ "07-23",
      month == "9" & day == "10" ~ "09-10",
      month == "9" & day == "29" ~ "09-29",
      month == "11" ~ "11-3",
      month == "12" ~ "12-4"
    ),
    bio_sample = paste(timepoint, bio, sep = "_")
  )

# Average technical replicates
count_avg_long <- counts_raw %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(meta %>% select(sample, bio_sample), by = "sample") %>%
  group_by(gene, bio_sample) %>%
  summarise(count = mean(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = bio_sample, values_from = count) %>%
  column_to_rownames("gene")

count_matrix <- as.matrix(count_avg_long)

# colData for DESeq2
coldata <- meta %>%
  distinct(bio_sample, timepoint) %>%
  arrange(timepoint) %>%
  column_to_rownames("bio_sample") %>%
  mutate(timepoint = factor(timepoint, 
                            levels = c("06-15", "07-23", "09-10", 
                                       "09-29", "11-3", "12-4")))

# Run DESeq2 once for VST
dds_full <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                   colData = coldata,
                                   design = ~ timepoint)
dds_full <- DESeq(dds_full)
vsd <- vst(dds_full, blind = TRUE)

# 8. DESeq2 per timepoint → top 100 genes (or all if fewer)
timepoints <- levels(coldata$timepoint)
top_results <- list()

for (tp in timepoints) {
  cat("Processing", tp, "...\n")
  
  coldata$group <- ifelse(coldata$timepoint == tp, "this", "other")
  coldata$group <- factor(coldata$group, levels = c("other", "this"))
  
  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                colData = coldata,
                                design = ~ group)
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("group", "this", "other"), 
                 alpha = 0.05, lfcThreshold = 0)
  
  res_ordered <- res[order(res$log2FoldChange, decreasing = TRUE), ]
  sig <- res_ordered[!is.na(res_ordered$padj) & 
                       res_ordered$padj < 0.05 & 
                       res_ordered$log2FoldChange > 0, ]
  
  n_genes <- min(100, nrow(sig))          # ← CHANGED TO 100
  if (n_genes == 0) next
  
  top_genes <- rownames(sig)[1:n_genes]
  
  norm_counts <- counts(dds, normalized = TRUE)
  mean_this <- rowMeans(norm_counts[top_genes, coldata$group == "this", drop = FALSE])
  
  out_df <- data.frame(
    gene = top_genes,
    as.data.frame(sig[1:n_genes, ]),
    mean_normalized_expr_this_month = mean_this[1:n_genes]
  )
  
  write.table(out_df, 
              file = paste0("top100_genes_", tp, ".tsv"),   # ← updated filename
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  top_results[[tp]] <- list(genes = top_genes, df = out_df)
}

cat("\n=== All TSV files saved! ===\n")

# 9. HEATMAPS — saved as PDF (8 × 12 inches)
pdf_width  <- 6
pdf_height <- 8

# Single-month heatmaps (skip 07-23 and 09-10)
for (tp in timepoints) {
  if (tp %in% c("07-23", "09-10")) next
  if (is.null(top_results[[tp]])) next
  
  genes <- top_results[[tp]]$genes
  mat <- assay(vsd)[genes, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  
  col_order <- order(coldata$timepoint)
  mat_z <- mat_z[, col_order, drop = FALSE]
  
  ann_col <- data.frame(Timepoint = coldata$timepoint[col_order])
  rownames(ann_col) <- colnames(mat_z)
  
  title <- paste("Top 100 genes highly expressed in", tp, 
                 "(row z-score; DESeq2 this vs others)")
  
  pdf(paste0("heatmap_top_genes_", tp, "top100_.pdf"), width = pdf_width, height = pdf_height)
  pheatmap(mat_z,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           annotation_col = ann_col,
           main = title,
           fontsize = 7,
           fontsize_row = 6,
           color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
           border_color = NA)
  dev.off()
  
  cat("PDF heatmap saved for", tp, "\n")
}

# Combined heatmap for 07-23 + 09-10
if (!is.null(top_results[["07-23"]]) && !is.null(top_results[["09-10"]])) {
  genes07 <- top_results[["07-23"]]$genes
  genes09 <- top_results[["09-10"]]$genes
  combined_genes <- unique(c(genes07, genes09))
  
  mat <- assay(vsd)[combined_genes, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  
  col_order <- order(coldata$timepoint)
  mat_z <- mat_z[, col_order, drop = FALSE]
  
  ann_col <- data.frame(Timepoint = coldata$timepoint[col_order])
  rownames(ann_col) <- colnames(mat_z)
  
  title <- "Top genes highly expressed in 07-23 or 09-10 (row z-score; DESeq2 this vs others)"
  
  pdf("heatmap_top_genes_combined_07-23_09-10_top100.pdf", width = pdf_width, height = pdf_height)
  pheatmap(mat_z,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           annotation_col = ann_col,
           main = title,
           fontsize = 7,
           fontsize_row = 6,
           color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
           border_color = NA)
  dev.off()
  
  cat("PDF heatmap saved for combined 07-23 + 09-10\n")
}


################################################################################
# ========================================================
# FINAL R Script: Bottom 100 & Bottom 200 genes (least expressed)
#                per timepoint + heatmaps (PDF)
# 
# Exactly as requested:
#   - Bottom 100 genes (most downregulated in that timepoint) → PDF 6×8 inches
#   - Bottom 200 genes (most downregulated in that timepoint) → PDF 6×12 inches
#   - TSV files: bottom100_genes_*.tsv and bottom200_genes_*.tsv
#   - Combined heatmap for 07-23 + 09-10 for both 100 and 200
#   - Gene names: short when unique, full original ID when duplicates
# ========================================================

library(DESeq2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# 2. Read the raw count file
counts_raw <- read.table("rsem_readCount_taihu98.tsv", 
                         header = TRUE, sep = "\t", 
                         row.names = 1, check.names = FALSE)

# 3. GENE NAME TRIMMING with duplicate handling
original_ids <- rownames(counts_raw)

gene_labels <- sapply(original_ids, function(id) {
  parts <- strsplit(id, "_", fixed = TRUE)[[1]]
  n <- length(parts)
  if (n == 3) {
    parts[3]
  } else if (n == 4) {
    paste(parts[3], parts[4], sep = "_")
  } else {
    paste(tail(parts, 2), collapse = "_")
  }
})

# Fix duplicates: use full original ID when short name collides
dup_names <- names(table(gene_labels))[table(gene_labels) > 1]
for (dup in dup_names) {
  idx <- which(gene_labels == dup)
  gene_labels[idx] <- original_ids[idx]
}

rownames(counts_raw) <- gene_labels

# 4. Clean column names
colnames(counts_raw) <- sub("^mT19_", "", colnames(counts_raw))

# 5. Parse metadata
meta <- data.frame(sample = colnames(counts_raw), stringsAsFactors = FALSE) %>%
  separate(sample, into = c("month", "day", "bio", "tech"), 
           sep = "-", remove = FALSE, convert = TRUE) %>%
  mutate(
    month = as.character(month),
    timepoint = case_when(
      month == "6"  ~ "06-15",
      month == "7"  ~ "07-23",
      month == "9" & day == "10" ~ "09-10",
      month == "9" & day == "29" ~ "09-29",
      month == "11" ~ "11-3",
      month == "12" ~ "12-4"
    ),
    bio_sample = paste(timepoint, bio, sep = "_")
  )

# Average technical replicates
count_avg_long <- counts_raw %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(meta %>% select(sample, bio_sample), by = "sample") %>%
  group_by(gene, bio_sample) %>%
  summarise(count = mean(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = bio_sample, values_from = count) %>%
  column_to_rownames("gene")

count_matrix <- as.matrix(count_avg_long)

# colData for DESeq2
coldata <- meta %>%
  distinct(bio_sample, timepoint) %>%
  arrange(timepoint) %>%
  column_to_rownames("bio_sample") %>%
  mutate(timepoint = factor(timepoint, 
                            levels = c("06-15", "07-23", "09-10", 
                                       "09-29", "11-3", "12-4")))

# Run DESeq2 once for VST
dds_full <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                   colData = coldata,
                                   design = ~ timepoint)
dds_full <- DESeq(dds_full)
vsd <- vst(dds_full, blind = TRUE)

# 8. Compute BOTTOM results once (most downregulated = most negative LFC)
timepoints <- levels(coldata$timepoint)
bottom_results <- list()

for (tp in timepoints) {
  cat("Computing bottom genes for", tp, "...\n")
  
  coldata$group <- ifelse(coldata$timepoint == tp, "this", "other")
  coldata$group <- factor(coldata$group, levels = c("other", "this"))
  
  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                colData = coldata,
                                design = ~ group)
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("group", "this", "other"), 
                 alpha = 0.05, lfcThreshold = 0)
  
  # Order by log2FoldChange ASCENDING → most negative (most downregulated) first
  res_ordered <- res[order(res$log2FoldChange, decreasing = FALSE), ]
  
  sig <- res_ordered[!is.na(res_ordered$padj) & 
                       res_ordered$padj < 0.05 & 
                       res_ordered$log2FoldChange < 0, ]
  
  bottom_results[[tp]] <- sig
}

cat("\n=== Bottom gene results computed! ===\n")

# 9. Generate heatmaps & TSV for BOTH bottom 100 and bottom 200
for (N in c(100, 200)) {
  cat("\n=== Processing BOTTOM", N, "genes ===\n")
  
  bottom_results_N <- list()
  pdf_width  <- 5
  pdf_height <- ifelse(N == 100, 8, 15)
  
  for (tp in timepoints) {
    sig <- bottom_results[[tp]]
    n_genes <- min(N, nrow(sig))
    if (n_genes == 0) next
    
    top_genes <- rownames(sig)[1:n_genes]   # already ordered most downregulated
    
    # Mean normalized expression in THIS month
    norm_counts <- counts(dds, normalized = TRUE)   # reuse last dds (same for all)
    mean_this <- rowMeans(norm_counts[top_genes, coldata$group == "this", drop = FALSE])
    
    out_df <- data.frame(
      gene = top_genes,
      as.data.frame(sig[1:n_genes, ]),
      mean_normalized_expr_this_month = mean_this[1:n_genes]
    )
    
    write.table(out_df, 
                file = paste0("bottom", N, "_genes_", tp, ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    bottom_results_N[[tp]] <- list(genes = top_genes, df = out_df)
  }
  
  # Single-month heatmaps (skip 07-23 and 09-10)
  for (tp in timepoints) {
    if (tp %in% c("07-23", "09-10")) next
    if (is.null(bottom_results_N[[tp]])) next
    
    genes <- bottom_results_N[[tp]]$genes
    mat <- assay(vsd)[genes, , drop = FALSE]
    mat_z <- t(scale(t(mat)))
    
    col_order <- order(coldata$timepoint)
    mat_z <- mat_z[, col_order, drop = FALSE]
    
    ann_col <- data.frame(Timepoint = coldata$timepoint[col_order])
    rownames(ann_col) <- colnames(mat_z)
    
    title <- paste("Bottom", N, "genes lowly expressed in", tp, 
                   "(row z-score; DESeq2 this vs others)")
    
    pdf(paste0("heatmap_bottom", N, "_genes_", tp, ".pdf"), 
        width = pdf_width, height = pdf_height)
    pheatmap(mat_z,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             annotation_col = ann_col,
             main = title,
             fontsize = 7,
             fontsize_row = 6,
             color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
             border_color = NA)
    dev.off()
    
    cat("PDF (", pdf_width, "×", pdf_height, ") saved for bottom", N, "in", tp, "\n")
  }
  
  # Combined heatmap for 07-23 + 09-10
  if (!is.null(bottom_results_N[["07-23"]]) && !is.null(bottom_results_N[["09-10"]])) {
    genes07 <- bottom_results_N[["07-23"]]$genes
    genes09 <- bottom_results_N[["09-10"]]$genes
    combined_genes <- unique(c(genes07, genes09))
    
    mat <- assay(vsd)[combined_genes, , drop = FALSE]
    mat_z <- t(scale(t(mat)))
    
    col_order <- order(coldata$timepoint)
    mat_z <- mat_z[, col_order, drop = FALSE]
    
    ann_col <- data.frame(Timepoint = coldata$timepoint[col_order])
    rownames(ann_col) <- colnames(mat_z)
    
    title <- paste("Bottom", N, "genes lowly expressed in 07-23 or 09-10",
                   "(row z-score; DESeq2 this vs others)")
    
    pdf(paste0("heatmap_bottom", N, "_genes_combined_07-23_09-10.pdf"), 
        width = pdf_width, height = pdf_height)
    pheatmap(mat_z,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             annotation_col = ann_col,
             main = title,
             fontsize = 7,
             fontsize_row = 6,
             color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
             border_color = NA)
    dev.off()
    
    cat("PDF (", pdf_width, "×", pdf_height, ") saved for bottom", N, "combined 07-23 + 09-10\n")
  }
}

################################################################################
# Acir310F SCRIPT: TOP + BOTTOM genes per timepoint + heatmaps
################################################################################
library(DESeq2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

rm(list = ls())

# -----------------------------
# 1. Read data
# -----------------------------
counts_raw <- read.table("rsem_readCount_Acir310F.tsv",
                         header = TRUE,
                         sep = "\t",
                         row.names = 1,
                         check.names = FALSE)

# -----------------------------
# 2. Rename genes
# -----------------------------
original_ids <- rownames(counts_raw)

make_gene_label <- function(id) {
  parts <- strsplit(id, "_", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(id)
  
  rs_part <- parts[2]
  suffix  <- paste(parts[3:length(parts)], collapse = "_")
  
  # case: ACIR310F_RS0108385_ACIR310F_RS0108385 -> RS0108385
  if (suffix == paste(parts[1], parts[2], sep = "_") || suffix == rs_part) {
    return(rs_part)
  }
  
  # case: ACIR310F_RS0116065_mutT -> mutT
  return(suffix)
}

gene_labels <- sapply(original_ids, make_gene_label, USE.NAMES = FALSE)

# duplicates -> RSxxxx_suffix
dup_labels <- names(table(gene_labels))[table(gene_labels) > 1]
for (lab in dup_labels) {
  idx <- which(gene_labels == lab)
  for (i in idx) {
    parts <- strsplit(original_ids[i], "_", fixed = TRUE)[[1]]
    if (length(parts) >= 3) {
      rs_part <- parts[2]
      suffix  <- paste(parts[3:length(parts)], collapse = "_")
      if (suffix == paste(parts[1], parts[2], sep = "_") || suffix == rs_part) {
        gene_labels[i] <- rs_part
      } else {
        gene_labels[i] <- paste(rs_part, suffix, sep = "_")
      }
    } else {
      gene_labels[i] <- original_ids[i]
    }
  }
}

# still duplicated -> original IDs
dup2 <- names(table(gene_labels))[table(gene_labels) > 1]
for (lab in dup2) {
  idx <- which(gene_labels == lab)
  gene_labels[idx] <- original_ids[idx]
}

rownames(counts_raw) <- gene_labels

# -----------------------------
# 3. Clean sample names
# -----------------------------
colnames(counts_raw) <- gsub("^Acir\\.mT19_", "", colnames(counts_raw))

# -----------------------------
# 4. Metadata
# -----------------------------
desired_tp_order <- c("06-15", "07-23", "09-10", "09-29", "11-03", "12-04")

meta <- data.frame(sample = colnames(counts_raw), stringsAsFactors = FALSE) %>%
  separate(sample,
           into = c("month", "day", "bio", "tech"),
           sep = "-",
           remove = FALSE,
           convert = TRUE) %>%
  mutate(
    month_num  = sprintf("%02d", month),
    day_num    = sprintf("%02d", day),
    bio        = as.integer(bio),
    tech       = as.integer(tech),
    timepoint  = paste(month_num, day_num, sep = "-"),
    bio_sample = paste(timepoint, bio, sep = "_")
  ) %>%
  mutate(
    timepoint = factor(timepoint, levels = desired_tp_order)
  ) %>%
  arrange(timepoint, bio, tech)

# -----------------------------
# 5. Average technical replicates
# -----------------------------
count_avg <- counts_raw %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  left_join(meta %>% select(sample, bio_sample), by = "sample") %>%
  group_by(gene, bio_sample) %>%
  summarise(count = mean(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = bio_sample, values_from = count) %>%
  column_to_rownames("gene")

count_matrix <- as.matrix(count_avg)
mode(count_matrix) <- "numeric"

# -----------------------------
# 6. colData
# -----------------------------
coldata <- meta %>%
  distinct(bio_sample, timepoint) %>%
  arrange(timepoint, bio_sample) %>%
  column_to_rownames("bio_sample")

# force exact column order match
sample_order_18 <- c(
  "06-15_1", "06-15_2", "06-15_3",
  "07-23_1", "07-23_2", "07-23_3",
  "09-10_1", "09-10_2", "09-10_3",
  "09-29_1", "09-29_2", "09-29_3",
  "11-03_1", "11-03_2", "11-03_3",
  "12-04_1", "12-04_2", "12-04_3"
)

count_matrix <- count_matrix[, sample_order_18, drop = FALSE]
coldata <- coldata[sample_order_18, , drop = FALSE]

stopifnot(identical(colnames(count_matrix), rownames(coldata)))

# -----------------------------
# 7. DESeq2 + VST
# -----------------------------
dds_full <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData   = coldata,
  design    = ~ timepoint
)

dds_full <- DESeq(dds_full)
vsd <- varianceStabilizingTransformation(dds_full, blind = FALSE)

# -----------------------------
# Heatmap for ALL genes
# same x-axis order and labels
# no gene names
# -----------------------------

plot_all_genes_heatmap <- function() {
  mat <- assay(vsd)[, sample_order_18, drop = FALSE]
  
  # row z-score
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  
  # displayed x labels only
  display_labels <- sample_order_18
  display_labels <- gsub("^11-03_", "11-3_", display_labels)
  display_labels <- gsub("^12-04_", "12-4_", display_labels)
  
  # keep only timepoint annotation bar
  ann_col <- data.frame(
    Timepoint = factor(
      gsub("_\\d+$", "", sample_order_18),
      levels = desired_tp_order
    )
  )
  rownames(ann_col) <- sample_order_18
  
  pdf("Acir_all_genes_heatmap.pdf", width = 5, height = 10)
  pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    labels_col = display_labels,
    annotation_col = ann_col,
    fontsize = 8,
    fontsize_col = 7,
    angle_col = 90,
    border_color = NA,
    main = "Acir all genes",
    color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)
  )
  dev.off()
}

# run it
plot_all_genes_heatmap()


# -----------------------------
# 8. DE results: each timepoint vs all others
# -----------------------------
timepoints <- desired_tp_order
result_list <- list()

for (tp in timepoints) {
  message("Computing DE for ", tp)
  
  coldata_tp <- as.data.frame(coldata)
  coldata_tp$group <- ifelse(coldata_tp$timepoint == tp, "this", "other")
  coldata_tp$group <- factor(coldata_tp$group, levels = c("other", "this"))
  
  dds_tp <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData   = coldata_tp,
    design    = ~ group
  )
  
  dds_tp <- DESeq(dds_tp)
  
  res <- results(dds_tp, contrast = c("group", "this", "other"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  norm_counts <- counts(dds_tp, normalized = TRUE)
  focal_samples <- rownames(coldata_tp)[coldata_tp$group == "this"]
  res_df$mean_expr_this_timepoint <- rowMeans(
    norm_counts[rownames(res_df), focal_samples, drop = FALSE]
  )
  
  result_list[[tp]] <- res_df
}

# -----------------------------
# 9. Heatmap function
# x-axis = 18 biological replicate columns
# -----------------------------
plot_gene_heatmap <- function(genes, prefix, tp, n_show, pdf_height) {
  genes <- intersect(genes, rownames(vsd))
  if (length(genes) == 0) {
    message("No genes to plot for ", prefix, " ", tp)
    return(NULL)
  }
  
  mat <- assay(vsd)[genes, sample_order_18, drop = FALSE]
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  
  # display labels only: remove leading zero from day in 11-03 and 12-04
  display_labels <- sample_order_18
  display_labels <- gsub("^11-03_", "11-3_", display_labels)
  display_labels <- gsub("^12-04_", "12-4_", display_labels)
  
  # keep only timepoint annotation bar
  ann_col <- data.frame(
    Timepoint = factor(
      gsub("_\\d+$", "", sample_order_18),
      levels = desired_tp_order
    )
  )
  rownames(ann_col) <- sample_order_18
  
  out_pdf <- paste0("Acir_", prefix, "_", n_show, "_genes_", tp, "_heatmap.pdf")
  
  pdf(out_pdf, width = 8, height = pdf_height)
  pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    labels_col = display_labels,
    annotation_col = ann_col,
    fontsize = 8,
    fontsize_row = ifelse(n_show == 100, 6, 4.5),
    fontsize_col = 7,
    angle_col = 90,
    border_color = NA,
    main = paste0("Acir ", prefix, " ", n_show, " genes: ", tp),
    color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)
  )
  dev.off()
}
# -----------------------------
# 10. Export top/bottom genes for all 6 timepoints
# -----------------------------
for (tp in timepoints) {
  res_df <- result_list[[tp]]
  
  res_use <- res_df %>%
    filter(!is.na(log2FoldChange), !is.na(padj))
  
  top_df <- res_use %>%
    arrange(desc(log2FoldChange), padj)
  
  bottom_df <- res_use %>%
    arrange(log2FoldChange, padj)
  
  # fallback if too few padj-defined genes
  if (nrow(res_use) < 50) {
    res_use <- res_df %>%
      filter(!is.na(log2FoldChange))
    top_df <- res_use %>% arrange(desc(log2FoldChange))
    bottom_df <- res_use %>% arrange(log2FoldChange)
  }
  
  for (N in c(100, 200)) {
    # top
    top_n <- top_df %>% slice_head(n = min(N, nrow(top_df)))
    write.table(
      top_n,
      file = paste0("Acir_top", N, "_genes_", tp, ".tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    plot_gene_heatmap(
      genes = top_n$gene,
      prefix = "top",
      tp = tp,
      n_show = N,
      pdf_height = ifelse(N == 100, 8, 15)
    )
    
    # bottom
    bottom_n <- bottom_df %>% slice_head(n = min(N, nrow(bottom_df)))
    write.table(
      bottom_n,
      file = paste0("Acir_bottom", N, "_genes_", tp, ".tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    plot_gene_heatmap(
      genes = bottom_n$gene,
      prefix = "bottom",
      tp = tp,
      n_show = N,
      pdf_height = ifelse(N == 100, 8, 15)
    )
  }
}

message("All Acir outputs completed.")


################################################################################
## plot taihu98 TPM values across 6 timepoints
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

rm(list = ls())

# -----------------------------
# 1. Read Taihu98 TPM file
# -----------------------------
expr_raw <- read.table("rsem_tpm_taihu98.tsv",
                       header = TRUE,
                       sep = "\t",
                       row.names = 1,
                       check.names = FALSE)

# -----------------------------
# 2. Keep Taihu98 gene IDs as-is
# -----------------------------
rownames(expr_raw) <- rownames(expr_raw)

# -----------------------------
# 3. Clean sample names
# -----------------------------
colnames(expr_raw) <- gsub("^.*mT19_", "", colnames(expr_raw))

# -----------------------------
# 4. Metadata
# -----------------------------
desired_tp_order <- c("06-15", "07-23", "09-10", "09-29", "11-03", "12-04")

meta <- data.frame(sample = colnames(expr_raw), stringsAsFactors = FALSE) %>%
  separate(sample,
           into = c("month", "day", "bio", "tech"),
           sep = "-",
           remove = FALSE,
           convert = TRUE) %>%
  mutate(
    month = as.integer(month),
    day   = as.integer(day),
    bio   = as.integer(bio),
    tech  = as.integer(tech),
    month_num  = sprintf("%02d", month),
    day_num    = sprintf("%02d", day),
    timepoint  = paste(month_num, day_num, sep = "-"),
    bio_sample = paste(timepoint, bio, sep = "_")
  ) %>%
  mutate(timepoint = factor(timepoint, levels = desired_tp_order)) %>%
  arrange(timepoint, bio, tech)

# -----------------------------
# 5. Average technical replicates
# -----------------------------
expr_avg <- expr_raw %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expr") %>%
  left_join(meta %>% select(sample, bio_sample), by = "sample") %>%
  group_by(gene, bio_sample) %>%
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = bio_sample, values_from = expr) %>%
  column_to_rownames("gene")

expr_matrix <- as.matrix(expr_avg)
mode(expr_matrix) <- "numeric"

# -----------------------------
# 6. Force the same 18-column x-axis order as Acir
# -----------------------------
sample_order_18 <- c(
  "06-15_1", "06-15_2", "06-15_3",
  "07-23_1", "07-23_2", "07-23_3",
  "09-10_1", "09-10_2", "09-10_3",
  "09-29_1", "09-29_2", "09-29_3",
  "11-03_1", "11-03_2", "11-03_3",
  "12-04_1", "12-04_2", "12-04_3"
)

expr_matrix <- expr_matrix[, sample_order_18, drop = FALSE]

# -----------------------------
# 7. Strong filtering:
# expression filter + variability filter
# -----------------------------
# expression filter: keep genes with TPM >= 1 in at least 2 samples
keep_expr <- rowSums(expr_matrix >= 1, na.rm = TRUE) >= 2
expr_matrix_filt <- expr_matrix[keep_expr, , drop = FALSE]

cat("Genes before filtering:           ", nrow(expr_matrix), "\n")
cat("After expression filter:          ", nrow(expr_matrix_filt), "\n")

# log2 transform
mat_log <- log2(expr_matrix_filt + 1)

# variability filter: keep genes with SD >= 0.5 after log2 transform
keep_var <- apply(mat_log, 1, sd, na.rm = TRUE) >= 0.5
mat_log_filt <- mat_log[keep_var, , drop = FALSE]

cat("After variability filter:         ", nrow(mat_log_filt), "\n")

# -----------------------------
# 8. Row z-score
# -----------------------------
mat_z <- t(scale(t(mat_log_filt)))
mat_z[is.na(mat_z)] <- 0

# optional clipping for cleaner color contrast
mat_z[mat_z >  2] <-  2
mat_z[mat_z < -2] <- -2

# -----------------------------
# 9. Display labels
# -----------------------------
display_labels <- sample_order_18
display_labels <- gsub("^11-03_", "11-3_", display_labels)
display_labels <- gsub("^12-04_", "12-4_", display_labels)

ann_col <- data.frame(
  Timepoint = factor(
    gsub("_\\d+$", "", sample_order_18),
    levels = desired_tp_order
  )
)
rownames(ann_col) <- sample_order_18

# -----------------------------
# 10. Heatmap
# same/similar color style as Acir
# -----------------------------
output_pdf <- "taihu98_all_genes_heatmap_filtered.pdf"
plot_title <- "taihu98 all genes (filtered, log2 TPM + z-score)"

# -----------------------------
# Heatmap with ultra-thin dendrogram (1/200 width)
# -----------------------------
nature_palette <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)

pdf(output_pdf, width = 5, height = 10)

p <- pheatmap(
  mat_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  labels_col = display_labels,
  annotation_col = ann_col,
  fontsize = 8,
  fontsize_col = 7,
  angle_col = 90,
  border_color = NA,
  main = plot_title,
  color = nature_palette,
  breaks = seq(-2, 2, length.out = 101),
  silent = TRUE
)

# -----------------------------
# Reduce dendrogram line width to 1/200 (~0.005)
# -----------------------------
library(grid)

for (i in seq_along(p$gtable$grobs)) {
  grob <- p$gtable$grobs[[i]]
  
  if (inherits(grob, "polyline")) {
    grob$gp$lwd <- 0.005
  }
  
  if ("grobs" %in% names(grob)) {
    for (j in seq_along(grob$grobs)) {
      if (inherits(grob$grobs[[j]], "polyline")) {
        grob$grobs[[j]]$gp$lwd <- 0.005
      }
    }
  }
  
  p$gtable$grobs[[i]] <- grob
}

grid.newpage()
grid.draw(p$gtable)

dev.off()

################################################################################
# Plot Acir310F all gene expression at 6 timepoints
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

rm(list = ls())

# -----------------------------
# 1. Read Acir310F TPM file
# -----------------------------
expr_raw <- read.table("rsem_tpm_Acir310F.tsv",
                       header = TRUE,
                       sep = "\t",
                       row.names = 1,
                       check.names = FALSE)

# -----------------------------
# 2. Gene renaming
# Rule:
# ACIR310F_RS0116065_mutT -> mutT preferred
# if duplicate, use RS0116065_mutT
# ACIR310F_RS0108385_ACIR310F_RS0108385 -> RS0108385
# -----------------------------
original_ids <- rownames(expr_raw)

make_gene_label <- function(id) {
  parts <- strsplit(id, "_", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(id)
  
  rs_part <- parts[2]
  suffix  <- paste(parts[3:length(parts)], collapse = "_")
  
  if (suffix == paste(parts[1], parts[2], sep = "_") || suffix == rs_part) {
    return(rs_part)
  }
  
  return(suffix)
}

gene_labels <- sapply(original_ids, make_gene_label, USE.NAMES = FALSE)

# duplicates -> RSxxxx_suffix
dup_labels <- names(table(gene_labels))[table(gene_labels) > 1]
for (lab in dup_labels) {
  idx <- which(gene_labels == lab)
  for (i in idx) {
    parts <- strsplit(original_ids[i], "_", fixed = TRUE)[[1]]
    if (length(parts) >= 3) {
      rs_part <- parts[2]
      suffix  <- paste(parts[3:length(parts)], collapse = "_")
      if (suffix == paste(parts[1], parts[2], sep = "_") || suffix == rs_part) {
        gene_labels[i] <- rs_part
      } else {
        gene_labels[i] <- paste(rs_part, suffix, sep = "_")
      }
    } else {
      gene_labels[i] <- original_ids[i]
    }
  }
}

# still duplicated -> original IDs
dup2 <- names(table(gene_labels))[table(gene_labels) > 1]
for (lab in dup2) {
  idx <- which(gene_labels == lab)
  gene_labels[idx] <- original_ids[idx]
}

rownames(expr_raw) <- gene_labels

# -----------------------------
# 3. Clean sample names
# handles Acir.mT19_11-3-1-1 or mT19_11-3-1-1
# -----------------------------
colnames(expr_raw) <- gsub("^.*mT19_", "", colnames(expr_raw))

# -----------------------------
# 4. Metadata
# -----------------------------
desired_tp_order <- c("06-15", "07-23", "09-10", "09-29", "11-03", "12-04")

meta <- data.frame(sample = colnames(expr_raw), stringsAsFactors = FALSE) %>%
  separate(sample,
           into = c("month", "day", "bio", "tech"),
           sep = "-",
           remove = FALSE,
           convert = TRUE) %>%
  mutate(
    month = as.integer(month),
    day   = as.integer(day),
    bio   = as.integer(bio),
    tech  = as.integer(tech),
    month_num  = sprintf("%02d", month),
    day_num    = sprintf("%02d", day),
    timepoint  = paste(month_num, day_num, sep = "-"),
    bio_sample = paste(timepoint, bio, sep = "_")
  ) %>%
  mutate(timepoint = factor(timepoint, levels = desired_tp_order)) %>%
  arrange(timepoint, bio, tech)

# -----------------------------
# 5. Average technical replicates
# -----------------------------
expr_avg <- expr_raw %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expr") %>%
  left_join(meta %>% select(sample, bio_sample), by = "sample") %>%
  group_by(gene, bio_sample) %>%
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = bio_sample, values_from = expr) %>%
  column_to_rownames("gene")

expr_matrix <- as.matrix(expr_avg)
mode(expr_matrix) <- "numeric"

# -----------------------------
# 6. Force 18-column x-axis order
# -----------------------------
sample_order_18 <- c(
  "06-15_1", "06-15_2", "06-15_3",
  "07-23_1", "07-23_2", "07-23_3",
  "09-10_1", "09-10_2", "09-10_3",
  "09-29_1", "09-29_2", "09-29_3",
  "11-03_1", "11-03_2", "11-03_3",
  "12-04_1", "12-04_2", "12-04_3"
)

expr_matrix <- expr_matrix[, sample_order_18, drop = FALSE]

# -----------------------------
# 7. log2 transform before plotting
# -----------------------------
keep_expr <- rowSums(expr_matrix >= 1, na.rm = TRUE) >= 2
expr_matrix_filt <- expr_matrix[keep_expr, , drop = FALSE]

mat_log <- log2(expr_matrix_filt + 1)

keep_var <- apply(mat_log, 1, sd, na.rm = TRUE) >= 0.5
mat_log <- mat_log[keep_var, , drop = FALSE]

mat_z <- t(scale(t(mat_log)))
mat_z[is.na(mat_z)] <- 0
# -----------------------------
# 8. Display labels
# vertical x labels
# show 11-3 and 12-4 instead of 11-03 and 12-04
# -----------------------------
display_labels <- sample_order_18
display_labels <- gsub("^11-03_", "11-3_", display_labels)
display_labels <- gsub("^12-04_", "12-4_", display_labels)

ann_col <- data.frame(
  Timepoint = factor(
    gsub("_\\d+$", "", sample_order_18),
    levels = desired_tp_order
  )
)
rownames(ann_col) <- sample_order_18

# -----------------------------
# Heatmap with ultra-thin dendrogram (1/200 width)
# -----------------------------
output_pdf <- "Acir_all_genes_heatmap_filtered.pdf"
plot_title <- ""

nature_palette <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)

pdf(output_pdf, width = 5, height = 10)

p <- pheatmap(
  mat_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  labels_col = display_labels,
  annotation_col = ann_col,
  fontsize = 8,
  fontsize_col = 7,
  angle_col = 90,
  border_color = NA,
  main = plot_title,
  color = nature_palette,
  breaks = seq(-2, 2, length.out = 101),
  silent = TRUE
)

# -----------------------------
# Reduce dendrogram line width to 1/200 (~0.005)
# -----------------------------
library(grid)

for (i in seq_along(p$gtable$grobs)) {
  grob <- p$gtable$grobs[[i]]
  
  if (inherits(grob, "polyline")) {
    grob$gp$lwd <- 0.005
  }
  
  if ("grobs" %in% names(grob)) {
    for (j in seq_along(grob$grobs)) {
      if (inherits(grob$grobs[[j]], "polyline")) {
        grob$grobs[[j]]$gp$lwd <- 0.005
      }
    }
  }
  
  p$gtable$grobs[[i]] <- grob
}

grid.newpage()
grid.draw(p$gtable)

dev.off()


################################
##############################
# new code for Acir Max or Min at Each time points with PDF
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(gridExtra)
  library(grid)
})
rm(list = ls())
counts_file <- "rsem_readCount_Acir310F.tsv"
tpm_file    <- "rsem_tpm_Acir310F.tsv"
outdir      <- "Acir310F_6tp_DESeq_mean_plots"
dir.create(outdir, showWarnings = FALSE)
desired_tp_order <- c("06-15", "07-23", "09-10", "09-29", "11-03", "12-04")

make_gene_label <- function(id) {
  parts <- strsplit(id, "_", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(id)
  rs_part <- parts[2]
  suffix  <- paste(parts[3:length(parts)], collapse = "_")
  if (suffix == paste(parts[1], parts[2], sep = "_") || suffix == rs_part) return(rs_part)
  suffix
}

clean_gene_labels <- function(ids) {
  original_ids <- ids
  gene_labels <- vapply(original_ids, make_gene_label, character(1))
  dup_labels <- names(table(gene_labels))[table(gene_labels) > 1]
  for (lab in dup_labels) {
    idx <- which(gene_labels == lab)
    for (i in idx) {
      parts <- strsplit(original_ids[i], "_", fixed = TRUE)[[1]]
      if (length(parts) >= 3) {
        rs_part <- parts[2]
        suffix  <- paste(parts[3:length(parts)], collapse = "_")
        if (suffix == paste(parts[1], parts[2], sep = "_") || suffix == rs_part) {
          gene_labels[i] <- rs_part
        } else {
          gene_labels[i] <- paste(rs_part, suffix, sep = "_")
        }
      } else gene_labels[i] <- original_ids[i]
    }
  }
  dup2 <- names(table(gene_labels))[table(gene_labels) > 1]
  for (lab in dup2) gene_labels[gene_labels == lab] <- original_ids[gene_labels == lab]
  gene_labels
}

clean_sample_name <- function(x) {
  x <- gsub("^Acir\\.mT19_", "", x)
  
  out <- vapply(x, function(s) {
    parts <- strsplit(s, "-", fixed = TRUE)[[1]]
    
    if (length(parts) != 4) {
      warning(paste("Unexpected sample format:", s))
      return(NA_character_)
    }
    
    # safely convert
    nums <- suppressWarnings(as.integer(parts))
    
    if (any(is.na(nums))) {
      warning(paste("Non-numeric sample parts:", s))
      return(NA_character_)
    }
    
    sprintf("%02d-%02d-%d-%d", nums[1], nums[2], nums[3], nums[4])
  }, character(1))
  
  return(out)
}
make_meta <- function(cols) {
  clean <- vapply(cols, clean_sample_name, character(1))
  tibble(sample_raw = clean) %>%
    separate(sample_raw, into = c("month", "day", "bio", "tech"), sep = "-", convert = TRUE) %>%
    mutate(timepoint = sprintf("%02d-%02d", month, day), bio_sample = paste0(timepoint, "_", bio), timepoint = factor(timepoint, levels = desired_tp_order)) %>%
    arrange(timepoint, bio, tech)
}

collapse_counts <- function(df, meta) {
  df %>% rownames_to_column("geneID") %>%
    pivot_longer(-geneID, names_to = "sample_raw", values_to = "count") %>%
    left_join(meta %>% select(sample_raw, bio_sample), by = "sample_raw") %>%
    group_by(geneID, bio_sample) %>% summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = bio_sample, values_from = count) %>%
    column_to_rownames("geneID") %>% as.matrix()
}

average_tpm <- function(df, meta) {
  df %>% rownames_to_column("geneID") %>%
    pivot_longer(-geneID, names_to = "sample_raw", values_to = "TPM") %>%
    left_join(meta %>% select(sample_raw, bio_sample), by = "sample_raw") %>%
    group_by(geneID, bio_sample) %>% summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = bio_sample, values_from = TPM) %>%
    column_to_rownames("geneID") %>% as.matrix()
}

counts_raw <- read.table(counts_file, header = TRUE, sep = "	", row.names = 1, check.names = FALSE)
tpm_raw    <- read.table(tpm_file, header = TRUE, sep = "	", row.names = 1, check.names = FALSE)
stopifnot(setequal(rownames(counts_raw), rownames(tpm_raw)))
tpm_raw <- tpm_raw[rownames(counts_raw), , drop = FALSE]

gene_name <- clean_gene_labels(rownames(counts_raw))
colnames(counts_raw) <- vapply(colnames(counts_raw), clean_sample_name, character(1))
colnames(tpm_raw)    <- vapply(colnames(tpm_raw), clean_sample_name, character(1))

library(dplyr)
library(tidyr)
library(tibble)

# metadata from already-cleaned column names
meta <- tibble(sample = colnames(counts_raw)) %>%
  separate(sample, into = c("month", "day", "bio", "tech"),
           sep = "-", remove = FALSE, convert = TRUE) %>%
  mutate(
    month = sprintf("%02d", month),
    day   = sprintf("%02d", day),
    timepoint = paste(month, day, sep = "-"),
    bio_sample = paste(timepoint, bio, sep = "_")
  )

collapse_counts <- function(count_df, meta) {
  count_df %>%
    rownames_to_column("geneID") %>%
    pivot_longer(
      cols = -geneID,
      names_to = "sample",
      values_to = "count"
    ) %>%
    left_join(meta %>% select(sample, bio_sample), by = "sample") %>%
    group_by(geneID, bio_sample) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    pivot_wider(names_from = bio_sample, values_from = count) %>%
    column_to_rownames("geneID") %>%
    as.matrix()
}

bio_order <- c(
  "06-15_1", "06-15_2", "06-15_3",
  "07-23_1", "07-23_2", "07-23_3",
  "09-10_1", "09-10_2", "09-10_3",
  "09-29_1", "09-29_2", "09-29_3",
  "11-03_1", "11-03_2", "11-03_3",
  "12-04_1", "12-04_2", "12-04_3"
)

counts_collapsed <- collapse_counts(as.data.frame(counts_raw), meta)[, bio_order]

average_tpm <- function(tpm_df, meta) {
  tpm_df %>%
    rownames_to_column("geneID") %>%
    pivot_longer(
      cols = -geneID,
      names_to = "sample",
      values_to = "TPM"
    ) %>%
    left_join(meta %>% select(sample, bio_sample), by = "sample") %>%
    group_by(geneID, bio_sample) %>%
    summarise(TPM = mean(TPM), .groups = "drop") %>%
    pivot_wider(names_from = bio_sample, values_from = TPM) %>%
    column_to_rownames("geneID") %>%
    as.matrix()
}

tpm_avg <- average_tpm(as.data.frame(tpm_raw), meta)[, bio_order]

write.table(cbind(geneID = rownames(counts_collapsed), gene_name = gene_name, as.data.frame(counts_collapsed)), file.path(outdir, "Acir310F_counts_collapsed_for_DESeq2.tsv"), sep = "	", quote = FALSE, row.names = FALSE)
write.table(cbind(geneID = rownames(tpm_avg), gene_name = gene_name, as.data.frame(tpm_avg)), file.path(outdir, "Acir310F_TPM_averaged_technical_replicates.tsv"), sep = "	", quote = FALSE, row.names = FALSE)

coldata <- meta %>% distinct(bio_sample, timepoint) %>% as.data.frame()
rownames(coldata) <- coldata$bio_sample
coldata <- coldata[bio_order, , drop = FALSE]


counts_collapsed <- round(counts_collapsed)
storage.mode(counts_collapsed) <- "integer"
dds <- DESeqDataSetFromMatrix(countData = counts_collapsed, colData = coldata, design = ~ timepoint)
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)
write.table(cbind(geneID = rownames(norm_counts), gene_name = gene_name, as.data.frame(norm_counts)), file.path(outdir, "Acir310F_normalized_counts_for_ranking.tsv"), sep = "	", quote = FALSE, row.names = FALSE)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

selection_list <- list()
res_list <- list()

for (tp in desired_tp_order) {
  coldata_tp <- coldata
  coldata_tp$group <- factor(
    ifelse(coldata_tp$timepoint == tp, "this", "other"),
    levels = c("other", "this")
  )
  
  dds_tp <- DESeqDataSetFromMatrix(
    countData = counts_collapsed,
    colData = coldata_tp,
    design = ~ group
  )
  dds_tp <- DESeq(dds_tp)
  
  res <- as.data.frame(results(dds_tp, contrast = c("group", "this", "other"))) %>%
    rownames_to_column("geneID") %>%
    left_join(id_map, by = "geneID") %>%
    mutate(
      timepoint = tp,
      sig_class = case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange > 0 ~ "Up",
        !is.na(padj) & padj < 0.05 & log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  res_list[[tp]] <- res
  
  tp_cols <- grep(paste0("^", tp, "_"), colnames(norm_counts), value = TRUE)
  
  rank_df <- tibble(
    geneID = rownames(norm_counts),
    mean_normalized_count_at_timepoint = rowMeans(norm_counts[, tp_cols, drop = FALSE])
  ) %>%
    left_join(id_map, by = "geneID")
  
  top_df <- rank_df %>%
    arrange(desc(mean_normalized_count_at_timepoint), geneID) %>%
    slice_head(n = 100) %>%
    mutate(block = "Top")
  
  bottom_df <- rank_df %>%
    arrange(mean_normalized_count_at_timepoint, geneID) %>%
    slice_head(n = 100) %>%
    mutate(block = "Bottom")
  
  sel <- bind_rows(top_df, bottom_df) %>%
    left_join(
      res %>% select(geneID, log2FoldChange, pvalue, padj, sig_class),
      by = "geneID"
    ) %>%
    mutate(timepoint = tp)
  
  selection_list[[tp]] <- sel
  
  write.table(
    sel,
    file.path(outdir, paste0("Acir310F_", tp, "_top_bottom_100_selected_genes.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}


selection_df <- bind_rows(selection_list)
write.table(selection_df, file.path(outdir, "Acir310F_all_timepoints_top_bottom_100_selected_genes.tsv"), sep = "	", quote = FALSE, row.names = FALSE)
write.table(bind_rows(res_list), file.path(outdir, "Acir310F_DESeq2_significance_all_genes.tsv"), sep = "	", quote = FALSE, row.names = FALSE)

z_tpm <- t(scale(t(tpm_avg))); z_tpm[is.na(z_tpm)] <- 0; z_tpm[z_tpm > 2.5] <- 2.5; z_tpm[z_tpm < -2.5] <- -2.5
ann_col <- data.frame(
  Timepoint = factor(sub("_[0-9]+$", "", bio_order), levels = desired_tp_order)
)
rownames(ann_col) <- bio_order
ann_colors <- list(Timepoint = c("06-15"="#6C8EBF","07-23"="#8FB3D9","09-10"="#B8CCE4","09-29"="#F4B183","11-03"="#E69138","12-04"="#C55A11"), Block = c("Top"="#D55E00","Bottom"="#0072B2"), Significance = c("Up"="#B2182B","NS"="#D9D9D9","Down"="#2166AC"))
my_breaks <- seq(-2.5, 2.5, length.out = 101); my_colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

build_panel <- function(tp, show_rownames = FALSE) {
  sel <- selection_list[[tp]]
  top_genes <- sel %>% filter(block == "Top") %>% arrange(desc(mean_normalized_count_at_timepoint)) %>% pull(geneID)
  bottom_genes <- sel %>% filter(block == "Bottom") %>% arrange(mean_normalized_count_at_timepoint) %>% pull(geneID)
  top_ord <- top_genes[hclust(dist(z_tpm[top_genes, bio_order, drop = FALSE]))$order]
  bottom_ord <- bottom_genes[hclust(dist(z_tpm[bottom_genes, bio_order, drop = FALSE]))$order]
  genes <- c(top_ord, bottom_ord)
  sig_map <- res_list[[tp]] %>% select(geneID, sig_class) %>% tibble::deframe()
  ann_row <- data.frame(Block = factor(c(rep("Top", length(top_ord)), rep("Bottom", length(bottom_ord))), levels = c("Top", "Bottom")), Significance = factor(sig_map[genes], levels = c("Up", "NS", "Down"))); rownames(ann_row) <- genes
  pheatmap(z_tpm[genes, bio_order, drop = FALSE], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = show_rownames, show_colnames = TRUE, fontsize = 8, fontsize_col = 7, fontsize_row = ifelse(show_rownames, 4.5, 0.1), border_color = NA, annotation_col = ann_col, annotation_row = ann_row, annotation_colors = ann_colors, breaks = my_breaks, color = my_colors, silent = TRUE, main = paste0(tp, " top and bottom 100 genes"))
}
for (tp in desired_tp_order) { ph <- build_panel(tp, show_rownames = TRUE); pdf(file.path(outdir, paste0("Acir310F_", tp, "_detailed_top_bottom100_TPM_heatmap.pdf")), width = 9, height = 13, useDingbats = FALSE); grid.newpage(); grid.draw(ph$gtable); dev.off() }
panel_list <- lapply(desired_tp_order, build_panel, show_rownames = FALSE)
pdf(file.path(outdir, "Acir310F_multipanel_top_bottom100_TPM_heatmaps.pdf"), width = 14, height = 16, useDingbats = FALSE)
grid.arrange(grobs = lapply(panel_list, `[[`, "gtable"), ncol = 2, top = textGrob("Acir310F top and bottom 100 expressed genes per timepoint", gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()
#####
##@@>> top100 and bottom100 plots for each time point
# -----------------------------
# TPM z-score matrix for plotting
# -----------------------------
z_tpm <- t(scale(t(tpm_avg)))
z_tpm[is.na(z_tpm)] <- 0
z_tpm[z_tpm > 2.5] <- 2.5
z_tpm[z_tpm < -2.5] <- -2.5

# column annotation
ann_col <- data.frame(
  Timepoint = factor(sub("_[0-9]+$", "", bio_order), levels = desired_tp_order)
)
rownames(ann_col) <- bio_order

# x-axis labels
display_labels <- bio_order
display_labels <- gsub("^11-03_", "11-3_", display_labels)
display_labels <- gsub("^12-04_", "12-4_", display_labels)

# -----------------------------
# function to plot one heatmap
# -----------------------------
plot_one_heatmap <- function(tp, which_block, selection_list, z_tpm, tpm_avg, ann_col, display_labels, outdir) {
  
  sel_tp <- selection_list[[tp]] %>%
    dplyr::filter(block == which_block)
  
  if (nrow(sel_tp) == 0) {
    message("No genes found for ", tp, " ", which_block)
    return(NULL)
  }
  
  gene_ids <- sel_tp$geneID
  gene_ids <- intersect(gene_ids, rownames(z_tpm))
  
  # filter low-TPM genes
  keep <- apply(tpm_avg[gene_ids, , drop = FALSE], 1, function(x) max(x) > 1)
  gene_ids <- gene_ids[keep]
  
  if (length(gene_ids) == 0) {
    message("All genes filtered out for ", tp, " ", which_block)
    return(NULL)
  }
  
  # preserve ranking
  gene_ids <- sel_tp$geneID[sel_tp$geneID %in% gene_ids]
  
  mat <- z_tpm[gene_ids, bio_order, drop = FALSE]
  
  gene_labels <- sel_tp$gene_name[match(rownames(mat), sel_tp$geneID)]
  gene_labels[is.na(gene_labels) | gene_labels == ""] <- rownames(mat)
  rownames(mat) <- gene_labels
  
  outfile <- file.path(outdir, paste0("Acir310F_", tp, "_", tolower(which_block), "100_TPM_heatmap.pdf"))
  
  pdf(outfile, width = 6, height = max(8, 0.12 * nrow(mat)))
  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    labels_col = display_labels,
    annotation_col = ann_col,
    fontsize = 8,
    fontsize_row = 6,
    fontsize_col = 7,
    angle_col = 90,
    border_color = NA,
    color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
    breaks = seq(-2.5, 2.5, length.out = 101),
    main = paste0("Acir310F ", tp, " ", which_block, " genes (TPM > 1)")
  )
  dev.off()
  
  message("Saved: ", outfile)
}

# -----------------------------
# generate 12 PDFs
# -----------------------------
for (tp in desired_tp_order) {
  plot_one_heatmap(
    tp = tp,
    which_block = "Top",
    selection_list = selection_list,
    z_tpm = z_tpm,
    tpm_avg = tpm_avg,
    ann_col = ann_col,
    display_labels = display_labels,
    outdir = outdir
  )
  
  plot_one_heatmap(
    tp = tp,
    which_block = "Bottom",
    selection_list = selection_list,
    z_tpm = z_tpm,
    tpm_avg = tpm_avg,
    ann_col = ann_col,
    display_labels = display_labels,
    outdir = outdir
  )
}


###############################
##############################
# top 75 bttom 75 genes
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

rm(list = ls())

# =============================
# 1. Input / output
# =============================
counts_file <- "rsem_readCount_Acir310F.tsv"
tpm_file    <- "rsem_tpm_Acir310F.tsv"
outdir      <- "Acir310F_top75_bottom75_heatmaps"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

desired_tp_order <- c("06-15", "07-23", "09-10", "09-29", "11-03", "12-04")
bio_order <- unlist(lapply(desired_tp_order, function(tp) paste0(tp, "_", 1:3)))

top_n <- 75
bottom_n <- 75
tpm_filter_threshold <- 1

# =============================
# 2. Helper functions
# =============================

make_gene_label <- function(id) {
  parts <- strsplit(id, "_", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(id)
  
  rs_part <- parts[2]
  suffix  <- paste(parts[3:length(parts)], collapse = "_")
  
  if (suffix == paste(parts[1], parts[2], sep = "_") || suffix == rs_part) {
    return(rs_part)
  }
  
  return(suffix)
}

clean_gene_labels <- function(ids) {
  original_ids <- ids
  gene_labels <- vapply(original_ids, make_gene_label, character(1))
  
  dup_labels <- names(table(gene_labels))[table(gene_labels) > 1]
  for (lab in dup_labels) {
    idx <- which(gene_labels == lab)
    for (i in idx) {
      parts <- strsplit(original_ids[i], "_", fixed = TRUE)[[1]]
      if (length(parts) >= 3) {
        rs_part <- parts[2]
        suffix  <- paste(parts[3:length(parts)], collapse = "_")
        if (suffix == paste(parts[1], parts[2], sep = "_") || suffix == rs_part) {
          gene_labels[i] <- rs_part
        } else {
          gene_labels[i] <- paste(rs_part, suffix, sep = "_")
        }
      } else {
        gene_labels[i] <- original_ids[i]
      }
    }
  }
  
  dup2 <- names(table(gene_labels))[table(gene_labels) > 1]
  for (lab in dup2) {
    idx <- which(gene_labels == lab)
    gene_labels[idx] <- original_ids[idx]
  }
  
  gene_labels
}

clean_sample_name <- function(x) {
  x <- gsub("^Acir\\.mT19_", "", x)
  
  out <- vapply(x, function(s) {
    parts <- strsplit(s, "-", fixed = TRUE)[[1]]
    
    if (length(parts) != 4) {
      warning(paste("Unexpected sample format:", s))
      return(NA_character_)
    }
    
    nums <- suppressWarnings(as.integer(parts))
    if (any(is.na(nums))) {
      warning(paste("Non-numeric sample parts:", s))
      return(NA_character_)
    }
    
    sprintf("%02d-%02d-%d-%d", nums[1], nums[2], nums[3], nums[4])
  }, character(1))
  
  return(out)
}

make_meta <- function(cols) {
  tibble(sample = cols) %>%
    separate(sample,
             into = c("month", "day", "bio", "tech"),
             sep = "-",
             remove = FALSE,
             convert = TRUE) %>%
    mutate(
      month = sprintf("%02d", month),
      day = sprintf("%02d", day),
      timepoint = paste(month, day, sep = "-"),
      bio_sample = paste(timepoint, bio, sep = "_")
    ) %>%
    mutate(
      timepoint = factor(timepoint, levels = desired_tp_order)
    ) %>%
    arrange(timepoint, bio, tech)
}

collapse_counts <- function(count_df, meta) {
  count_df %>%
    rownames_to_column("geneID") %>%
    pivot_longer(
      cols = -geneID,
      names_to = "sample",
      values_to = "count"
    ) %>%
    left_join(meta %>% select(sample, bio_sample), by = "sample") %>%
    group_by(geneID, bio_sample) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = bio_sample, values_from = count) %>%
    column_to_rownames("geneID") %>%
    as.matrix()
}

average_tpm <- function(tpm_df, meta) {
  tpm_df %>%
    rownames_to_column("geneID") %>%
    pivot_longer(
      cols = -geneID,
      names_to = "sample",
      values_to = "TPM"
    ) %>%
    left_join(meta %>% select(sample, bio_sample), by = "sample") %>%
    group_by(geneID, bio_sample) %>%
    summarise(TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = bio_sample, values_from = TPM) %>%
    column_to_rownames("geneID") %>%
    as.matrix()
}

plot_one_heatmap <- function(tp, which_block, selection_list, z_tpm, tpm_avg,
                             ann_col, display_labels, outdir,
                             desired_tp_order, bio_order,
                             tpm_filter_threshold = 1) {
  
  sel_tp <- selection_list[[tp]] %>%
    dplyr::filter(block == which_block)
  
  if (nrow(sel_tp) == 0) {
    message("No genes found for ", tp, " ", which_block)
    return(NULL)
  }
  
  gene_ids <- sel_tp$geneID
  gene_ids <- intersect(gene_ids, rownames(z_tpm))
  
  # Filter low TPM genes using focal timepoint only
  tp_sample_cols <- grep(paste0("^", tp, "_"), colnames(tpm_avg), value = TRUE)
  keep <- apply(tpm_avg[gene_ids, tp_sample_cols, drop = FALSE], 1, function(x) mean(x) > tpm_filter_threshold)
  gene_ids <- gene_ids[keep]
  
  if (length(gene_ids) == 0) {
    message("All genes filtered out for ", tp, " ", which_block)
    return(NULL)
  }
  
  # Preserve selection/ranking order
  gene_ids <- sel_tp$geneID[sel_tp$geneID %in% gene_ids]
  
  mat <- z_tpm[gene_ids, bio_order, drop = FALSE]
  
  gene_labels <- sel_tp$gene_name[match(rownames(mat), sel_tp$geneID)]
  gene_labels[is.na(gene_labels) | gene_labels == ""] <- rownames(mat)[is.na(gene_labels) | gene_labels == ""]
  rownames(mat) <- gene_labels
  
  outfile <- file.path(
    outdir,
    paste0("Acir310F_", tp, "_", tolower(which_block), "75_TPM_heatmap.pdf")
  )
  
  pdf_height <- max(7, 0.11 * nrow(mat))
  
  pdf(outfile, width = 6, height = pdf_height)
  pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    labels_col = display_labels,
    annotation_col = ann_col,
    fontsize = 8,
    fontsize_row = 6,
    fontsize_col = 7,
    angle_col = 90,
    border_color = NA,
    color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
    breaks = seq(-2.5, 2.5, length.out = 101),
    main = paste0("Acir310F ", tp, " ", which_block, " 75 genes")
  )
  dev.off()
  
  message("Saved: ", outfile)
}

# =============================
# 3. Read data
# =============================
counts_raw <- read.table(
  counts_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

tpm_raw <- read.table(
  tpm_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

stopifnot(setequal(rownames(counts_raw), rownames(tpm_raw)))
tpm_raw <- tpm_raw[rownames(counts_raw), , drop = FALSE]

# =============================
# 4. Clean gene and sample names
# =============================
gene_name <- clean_gene_labels(rownames(counts_raw))
id_map <- tibble(geneID = rownames(counts_raw), gene_name = gene_name)

colnames(counts_raw) <- clean_sample_name(colnames(counts_raw))
colnames(tpm_raw)    <- clean_sample_name(colnames(tpm_raw))

meta <- make_meta(colnames(counts_raw))

# =============================
# 5. Collapse technical replicates
# =============================
counts_collapsed <- collapse_counts(as.data.frame(counts_raw), meta)[, bio_order, drop = FALSE]
counts_collapsed <- round(counts_collapsed)
storage.mode(counts_collapsed) <- "integer"

tpm_avg <- average_tpm(as.data.frame(tpm_raw), meta)[, bio_order, drop = FALSE]

# Make sure rows align
common_genes <- intersect(rownames(counts_collapsed), rownames(tpm_avg))
counts_collapsed <- counts_collapsed[common_genes, , drop = FALSE]
tpm_avg <- tpm_avg[common_genes, , drop = FALSE]
id_map <- id_map %>% filter(geneID %in% common_genes)

# Export collapsed matrices
write.table(
  cbind(geneID = rownames(counts_collapsed),
        gene_name = id_map$gene_name[match(rownames(counts_collapsed), id_map$geneID)],
        as.data.frame(counts_collapsed)),
  file.path(outdir, "Acir310F_counts_collapsed_for_DESeq2.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  cbind(geneID = rownames(tpm_avg),
        gene_name = id_map$gene_name[match(rownames(tpm_avg), id_map$geneID)],
        as.data.frame(tpm_avg)),
  file.path(outdir, "Acir310F_TPM_averaged_technical_replicates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# =============================
# 6. DESeq2
# =============================
coldata <- meta %>%
  distinct(bio_sample, timepoint) %>%
  as.data.frame()

rownames(coldata) <- coldata$bio_sample
coldata <- coldata[bio_order, "timepoint", drop = FALSE]
coldata$timepoint <- factor(coldata$timepoint, levels = desired_tp_order)

stopifnot(all(colnames(counts_collapsed) == rownames(coldata)))

dds <- DESeqDataSetFromMatrix(
  countData = counts_collapsed,
  colData = coldata,
  design = ~ timepoint
)

dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE)
write.table(
  cbind(geneID = rownames(norm_counts),
        gene_name = id_map$gene_name[match(rownames(norm_counts), id_map$geneID)],
        as.data.frame(norm_counts)),
  file.path(outdir, "Acir310F_normalized_counts_for_ranking.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Optional transform; not required for the heatmaps below, but useful to have
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# =============================
# 7. Select top/bottom 75 per timepoint
# =============================
selection_list <- list()
res_list <- list()

for (tp in desired_tp_order) {
  coldata_tp <- coldata
  coldata_tp$group <- factor(
    ifelse(coldata_tp$timepoint == tp, "this", "other"),
    levels = c("other", "this")
  )
  
  dds_tp <- DESeqDataSetFromMatrix(
    countData = counts_collapsed,
    colData = coldata_tp,
    design = ~ group
  )
  dds_tp <- DESeq(dds_tp)
  
  res <- as.data.frame(results(dds_tp, contrast = c("group", "this", "other"))) %>%
    rownames_to_column("geneID") %>%
    left_join(id_map, by = "geneID") %>%
    mutate(
      timepoint = tp,
      sig_class = case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange > 0 ~ "Up",
        !is.na(padj) & padj < 0.05 & log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  res_list[[tp]] <- res
  
  tp_cols <- grep(paste0("^", tp, "_"), colnames(norm_counts), value = TRUE)
  
  rank_df <- tibble(
    geneID = rownames(norm_counts),
    mean_normalized_count_at_timepoint = rowMeans(norm_counts[, tp_cols, drop = FALSE])
  ) %>%
    left_join(id_map, by = "geneID")
  
  top_df <- rank_df %>%
    arrange(desc(mean_normalized_count_at_timepoint), geneID) %>%
    slice_head(n = top_n) %>%
    mutate(block = "Top")
  
  bottom_df <- rank_df %>%
    arrange(mean_normalized_count_at_timepoint, geneID) %>%
    slice_head(n = bottom_n) %>%
    mutate(block = "Bottom")
  
  sel <- bind_rows(top_df, bottom_df) %>%
    left_join(
      res %>% select(geneID, log2FoldChange, pvalue, padj, sig_class),
      by = "geneID"
    ) %>%
    mutate(timepoint = tp)
  
  selection_list[[tp]] <- sel
  
  write.table(
    sel,
    file.path(outdir, paste0("Acir310F_", tp, "_top_bottom_75_selected_genes.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

selection_df <- bind_rows(selection_list)
write.table(
  selection_df,
  file.path(outdir, "Acir310F_all_timepoints_top_bottom_75_selected_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  bind_rows(res_list),
  file.path(outdir, "Acir310F_DESeq2_significance_all_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# =============================
# 8. TPM z-score matrix for plotting
# =============================
z_tpm <- t(scale(t(tpm_avg)))
z_tpm[is.na(z_tpm)] <- 0
z_tpm[z_tpm > 2.5] <- 2.5
z_tpm[z_tpm < -2.5] <- -2.5

ann_col <- data.frame(
  Timepoint = factor(sub("_[0-9]+$", "", bio_order), levels = desired_tp_order)
)
rownames(ann_col) <- bio_order

display_labels <- bio_order
display_labels <- gsub("^11-03_", "11-3_", display_labels)
display_labels <- gsub("^12-04_", "12-4_", display_labels)

# =============================
# 9. Generate 12 PDFs
# =============================
for (tp in desired_tp_order) {
  plot_one_heatmap(
    tp = tp,
    which_block = "Top",
    selection_list = selection_list,
    z_tpm = z_tpm,
    tpm_avg = tpm_avg,
    ann_col = ann_col,
    display_labels = display_labels,
    outdir = outdir,
    desired_tp_order = desired_tp_order,
    bio_order = bio_order,
    tpm_filter_threshold = tpm_filter_threshold
  )
  
  plot_one_heatmap(
    tp = tp,
    which_block = "Bottom",
    selection_list = selection_list,
    z_tpm = z_tpm,
    tpm_avg = tpm_avg,
    ann_col = ann_col,
    display_labels = display_labels,
    outdir = outdir,
    desired_tp_order = desired_tp_order,
    bio_order = bio_order,
    tpm_filter_threshold = tpm_filter_threshold
  )
}

message("All outputs completed.")
