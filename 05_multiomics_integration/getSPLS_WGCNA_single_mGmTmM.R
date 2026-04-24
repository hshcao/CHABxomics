# prepare metagenome data
mg0=read.csv("mG19_merged_gene_families_TPM.tsv", sep = "\t", header = T)
colnames(mg0)
mg0=mg0[rowSums(is.na(mg0))<10,]
mg0=mg0[!grepl('\\||unclassified|UNMAPPED|unknown', mg0$Gene, perl=TRUE),]
colnames(mg0)=gsub("_([679])","_0\\1",colnames(mg0), perl = TRUE)
mg0=mg0[,sort(colnames(mg0))]

# prepare metagenome data
mt0=read.csv("Table5_mGmT_geneAb_across-TimeSite.tsv", sep = "\t", header = T)
colnames(mt0)
keep_cols=grepl("_TPM",colnames(mt0))
keep_cols[1:2]=c(TRUE, TRUE)
mt0=mt0[,keep_cols]
mt0=mt0[rowSums(is.na(mt0))<13,]

mergedmGT=merge(mg0, mt0, by.x = "Gene", by.y = "UniRef90_ID", all.x = TRUE)
mergedmGT=mergedmGT[,-19]
colnames(mergedmGT)=gsub("X([679])","X_0\\1",colnames(mergedmGT), perl = TRUE)
colnames(mergedmGT)=gsub("X_","mT19_", colnames(mergedmGT))
colnames(mergedmGT)=gsub("X","mT19", colnames(mergedmGT))
colnames(mergedmGT)=gsub("_TPM|19","", colnames(mergedmGT))
colnames(mergedmGT)=gsub("_","", colnames(mergedmGT))
colnames(mergedmGT)=gsub("\\.","_", colnames(mergedmGT))
colnames(mergedmGT)=gsub("_([123]$)","\\.\\1", colnames(mergedmGT))
write.table(mergedmGT, "Table5-2_mGmT_tpm_bytimeDate.tsv", sep = "\t", quote = F, row.names = F)
################################################################################
############## Correlation b/w species abundance and T #########################
## read in metadata
rm(list=ls())
mdf = read.csv('water_chemistry.csv', header = T, row.names = 1, check.names = F)
dim(mdf) # confirm the dimension of data
rownames(mdf) = gsub('mG19_', 'm', rownames(mdf))
df=mdf[1:12,c(1,4)]
df$npc=mdf$PC[1:12]
df$npc[7:9]=(as.numeric(df$npc[7:9])+as.numeric(df$npc[10:12]))/2
cor.test(df$T, df$npc)
range(mdf$T[13:18])
cor.test(mdf$T, mdf$PC)
plot(df$T,df$npc)
abline(df$T~df$npc, color="red")
library(ggplot2)
ggplot(df, aes(x = T, y = npc)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal()
## read in mG species data
sdf=read.delim("kaiju_equalR4m_taxonAbund_genus_cutoffR10.tsv",sep = '\t')
# Split the column by ";" and extract second-to-last
last_elements = sapply(strsplit(xdf$taxon_name, ";"), function(x) {
  x=x[x != ""]
  tail(x, 1)
})
rownames(xdf) <- last_elements
xdf=xdf[,-1]

X = as.data.frame(t(xdf))


############### end of do spls for mG genus #################
################################################################################
## draw the sPCA of metaData ##
rm(list=ls())
library(mixOmics)
X = read.csv('water_chemistry.csv', header = T, row.names = 1, check.names = F)
dim(X) # confirm the dimension of data
# run preliminary model
rownames(X) = gsub('mG19_', '', rownames(X))
x.pca = pca(X, center = T, scale = T)
plot(x.pca)

pdf("spls_PCA_water_metaData.pdf", 5, 5)
# plot final sPCA samples for first two components
plotIndiv(x.pca, comp = c(1, 2), ind.names = TRUE, 
          group = rownames(X),  # use class to colour each sample
          legend = FALSE, title = 'metaData PCA comp 1 - 2')
dev.off()
# plot variables against the sPCA components
plotVar(x.pca, comp = c(1, 2), var.names = TRUE,  
        title = 'metaData PCA comp 1 - 2')
# plot samples and variables against the sPCA components
pdf("spls_PCA_water_metaData_biplot.pdf", 5, 5)
biplot(x.pca, cex = 1,  pch.size = 1.5, ind.names.size = 3, arrow=FALSE,
       xlabs = paste(ncol(X), 1:nrow(X)), #simplify names
       var.arrow.size = 0.2, var.names.size = 2.8, legend = FALSE,
       group = rownames(X), #multidrug$cell.line$Class,  # colour by sample class
       title = 'metaData PCA comp 1 - 2')
dev.off()

set.seed(8589) # for reproducibility with this case study, remove otherwise
test.keepX <- 1:ncol(X) # set the number of variable values to be tested
tune.spca.res <- tune.spca(X, ncomp = 3, # generate the first 3 components
                           nrepeat = 5, # repeat the cross-validation 5 times
                           folds = 3, # use 3 folds for the cross-validation
                           test.keepX = test.keepX)
plot(tune.spca.res) # plot the optimisation output
tune.spca.res$choice.keepX # how many variables per component is optimal
final.spca <- spca(X, ncomp = 2, # based off figure 1, three components is best
                   keepX =c(14,4)) # tune.spca.res$choice.keepX)
final.spca$variates$X[,2] <- -final.spca$variates$X[,2]
# plot final sPCA samples for first two components
plotIndiv(final.spca, comp = c(1, 2), ind.names = TRUE, 
          group = rownames(X),  # use class to colour each sample
          legend = FALSE, title = 'Multidrug transporter, sPCA comp 1 - 2')
# plot variables against the sPCA components
plotVar(final.spca, comp = c(1, 2), var.names = TRUE,  
        title = 'Multidrug transporter, sPCA comp 1 - 2')
# plot samples and variables against the sPCA components
pdf("spls_sPCA_water_metaData_biplot.pdf", 3, 3)
biplot(final.spca, cex =0.7,  pch.size = 1.2, ind.names.size = 2,
       xlabs = paste(ncol(X), 1:nrow(X)), #simplify names
       var.arrow.size = 0.1, var.names.size = 1.5, legend = FALSE,
       group = rownames(X), #multidrug$cell.line$Class,  # colour by sample class
       title = 'Multidrug transporter, sPCA comp 1 - 2')
dev.off()
## the end of metaData PCA and sPAC plots ##
################################################################################
############### do spls for mG genus #################
## prepare  mG taxon abundance data
rm(list=ls())
xdf=read.delim("kaiju_equalR4m_taxonAbund_genus_cutoffR10.tsv",sep = '\t')
# Split the column by ";" and extract second-to-last
last_elements = sapply(strsplit(xdf$taxon_name, ";"), function(x) {
  x=x[x != ""]
  tail(x, 1)
})
rownames(xdf) <- last_elements
xdf=xdf[,-1]

X = as.data.frame(t(xdf))

library(mixOmics)

x.pca = pca(X, center = T, scale = T)
x.pca$variates$X[,1] <- -x.pca$variates$X[,1]
#x.pca$variates$X[,2] <- -x.pca$variates$X[,2]
plot(x.pca)

#pdf("spls_PCA_mGgenus_abundance.pdf", 5,5)
# plot final sPCA samples for first two components
plotIndiv(x.pca, comp = c(1, 2), ind.names = TRUE, 
          group = rownames(X),  # use class to colour each sample
          legend = FALSE, title = "")
#dev.off()
# plot variables against the sPCA components
plotVar(x.pca, comp = c(1, 2), var.names = TRUE,  
        title = 'metaData PCA comp 1 - 2')
# plot samples and variables against the sPCA components
pdf("spls_PCA_mGgenus_abudance_biplot.pdf", 8, 6)
biplot(x.pca, cex = 1,  pch.size = 1, ind.names.size = 2, 
       xlabs = paste(ncol(X), 1:nrow(X)), #simplify names
       var.arrow.col = 'grey90',var.names.angle = FALSE,
       var.arrow.size = 0.01, var.names.size = 0.5, legend = FALSE,
       group = rownames(X), #multidrug$cell.line$Class,  # colour by sample class
       title = 'mG PCA comp 1 - 2', arrow=FALSE,var.arrow.length = 0.1)

dev.off()


test.keepX <- 1:nrow(X) # set the number of variable values to be tested
tune.spca.res <- tune.spca(X, ncomp = 3, # generate the first 3 components
                           nrepeat = 5, # repeat the cross-validation 5 times
                           folds = 3, # use 3 folds for the cross-validation
                           test.keepX = test.keepX)
plot(tune.spca.res) # plot the optimisation output
tune.spca.res$choice.keepX # how many variables per component is optimal
final.spca <- spca(X, ncomp = 2, # based off figure 1, three components is best
                   keepX =c(4,1)) # tune.spca.res$choice.keepX)
# plot final sPCA samples for first two components
plotIndiv(final.spca, comp = c(1, 2), ind.names = TRUE, 
          group = rownames(X),  # use class to colour each sample
          legend = TRUE, title = 'Multidrug transporter, sPCA comp 1 - 2')
# plot variables against the sPCA components
plotVar(final.spca, comp = c(1, 2), var.names = TRUE,  
        title = 'Multidrug transporter, sPCA comp 1 - 2')
# plot samples and variables against the sPCA components
cim(final.spca)
pdf("sPCA_mGgenus.pdf", 5, 5)
biplot(final.spca, cex = 1,  pch.size = 1, ind.names.size = 2.5,
       xlabs = paste(ncol(X), 1:nrow(X)), #simplify names
       var.arrow.size = 0.2, var.names.size = 2.8, legend = FALSE,
       group = rownames(X), #multidrug$cell.line$Class,  # colour by sample class
       title = 'Multidrug transporter, sPCA comp 1 - 2')
dev.off()
############### draw the PCA & sPCS of mG genus #################
## prepare  mG taxon abundance data
rm(list=ls())
xdf=read.delim("mG19_genus_abundance.txt",sep = '\t')
xdf=xdf[c(-1:-2, (-nrow(xdf)+1),-nrow(xdf)),]
colnames(xdf)=gsub("_([679])\\.","_0\\1\\.", colnames(xdf))
colnames(xdf)=gsub("mG19_","mG", colnames(xdf))

df=xdf[,-1]
rownames(df)=xdf[,1]
df=df[,order(colnames(df))]
df = as.data.frame(t(df))
colnames(df) <- sapply(strsplit(colnames(df), ";"), function(x) tail(x, 1))

X = df

library(mixOmics)

x.pca = pca(X, center = T, scale = T)
x.pca$variates$X[,1] <- -x.pca$variates$X[,1]
x.pca$variates$X[,2] <- -x.pca$variates$X[,2]
plot(x.pca)

pdf("spls_PCA_mGgenus_abundance.pdf", 5,5)
# plot final sPCA samples for first two components
plotIndiv(x.pca, comp = c(1, 2), ind.names = TRUE, 
          group = rownames(X),  # use class to colour each sample
          legend = FALSE, title = "")
dev.off()
# plot variables against the sPCA components
plotVar(x.pca, comp = c(1, 2), var.names = TRUE,  
        title = 'metaData PCA comp 1 - 2')
# plot samples and variables against the sPCA components
pdf("spls_PCA_mGgenus_abudance_biplot.pdf", 8, 6)
biplot(x.pca, cex = 1,  pch.size = 1, ind.names.size = 2, 
       xlabs = paste(ncol(X), 1:nrow(X)), #simplify names
       var.arrow.col = 'grey90',var.names.angle = FALSE,
       var.arrow.size = 0.01, var.names.size = 0.5, legend = FALSE,
       group = rownames(X), #multidrug$cell.line$Class,  # colour by sample class
       title = 'mG PCA comp 1 - 2', arrow=FALSE,var.arrow.length = 0.1)
dev.off()


test.keepX <- 1:nrow(X) # set the number of variable values to be tested
tune.spca.res <- tune.spca(X, ncomp = 3, # generate the first 3 components
                           nrepeat = 5, # repeat the cross-validation 5 times
                           folds = 3, # use 3 folds for the cross-validation
                           test.keepX = test.keepX)
plot(tune.spca.res) # plot the optimisation output
tune.spca.res$choice.keepX # how many variables per component is optimal
final.spca <- spca(X, ncomp = 2, # based off figure 1, three components is best
                   keepX =c(4,1)) # tune.spca.res$choice.keepX)
# plot final sPCA samples for first two components
plotIndiv(final.spca, comp = c(1, 2), ind.names = TRUE, 
          group = rownames(X),  # use class to colour each sample
          legend = TRUE, title = 'Multidrug transporter, sPCA comp 1 - 2')
# plot variables against the sPCA components
plotVar(final.spca, comp = c(1, 2), var.names = TRUE,  
        title = 'Multidrug transporter, sPCA comp 1 - 2')
# plot samples and variables against the sPCA components
cim(final.spca)
pdf("sPCA_mGgenus.pdf", 5, 5)
biplot(final.spca, cex = 1,  pch.size = 1, ind.names.size = 2.5,
       xlabs = paste(ncol(X), 1:nrow(X)), #simplify names
       var.arrow.size = 0.2, var.names.size = 2.8, legend = FALSE,
       group = rownames(X), #multidrug$cell.line$Class,  # colour by sample class
       title = 'Multidrug transporter, sPCA comp 1 - 2')
dev.off()
##### the end of mG genus abundance #####
################################################################################
############### draw the PCA & sPCS of mG genes #################
## prepare  mG gene abundance data
rm(list=ls())
xdf=read.delim("summary_table_5_across-time.tsv",sep = '\t')
rownames(xdf) =xdf[[1]]
xdf=xdf[, seq(3,71,4)]
colnames(xdf)=gsub("X([679])","X0\\1", colnames(xdf))
colnames(xdf)=gsub("^X","mG", colnames(xdf))
colnames(xdf)=gsub("_abundance","", colnames(xdf))
# remove NAs from columns
xdf <- xdf[, colSums(is.na(xdf)) < nrow(xdf)]
xdf[is.na(xdf)] <- 0
df = as.data.frame(t(xdf))
df = df[order(rownames(df)),]
X = df
x.pca = pca(X, center = T, scale = T)
x.pca$variates$X[,1] <- -x.pca$variates$X[,1]
x.pca$variates$X[,2] <- -x.pca$variates$X[,2]
plot(x.pca)

pdf("spls_PCA_mGgene_abundance.pdf", 6,5)
# plot final sPCA samples for first two components
plotIndiv(x.pca, comp = c(1, 2), ind.names = TRUE, 
          group = rownames(X),  # use class to colour each sample
          legend = FALSE, title = "")
dev.off()
# plot variables against the sPCA components
plotVar(x.pca, comp = c(1, 2), var.names = TRUE,  
        title = 'metaData PCA comp 1 - 2')
# plot samples and variables against the sPCA components
pdf("spls_sPCA_mGgene_abudance_biplot.pdf", 8, 6)
biplot(x.pca, cex = 1,  pch.size = 1, ind.names.size = 2, cutoff = 0.15,
       xlabs = paste(ncol(X), 1:nrow(X)), #simplify names
       var.arrow.size = 0.01, var.names.size = 0.5, 
       var.arrow.col = 'grey90', var.arrow.length = 0.1,
       legend = FALSE, alpha=0.9,
       group = rownames(X), #multidrug$cell.line$Class,  # colour by sample class
       title = 'mG PCA comp 1 - 2')
dev.off()
##### the end of mG gene abundance #####
################################################################################
############### draw the PCA & sPCA of mT genes #################
## prepare  mT gene abundance data
rm(list=ls())
library(mixOmics)
xdf=read.delim("summary_table_5_across-time.tsv",sep = '\t')
rownames(xdf) =xdf[[1]]
xdf=xdf[, seq(5,73,4)]
colnames(xdf)=gsub("X([679])","X_0\\1", colnames(xdf))
colnames(xdf)=gsub("^X","mTgene", colnames(xdf))
colnames(xdf)=gsub("_TPM","", colnames(xdf))
# remove NAs from columns
xdf <- xdf[, colSums(is.na(xdf)) < nrow(xdf)]
xdf[is.na(xdf)] <- 0
df = as.data.frame(t(xdf))
X = df
X = X[,-c(17,143,222,412,446,461,487,765,929,1115,1146,1228,1395,2129,2166,2409,2589,2794,2843,3132,3149,3251,3273,3383,3389,3450,3460,3501,3594,3624,3632,3756,3856,3859,4028,4132,4534,4538,4546,4737,4859,4984,5000,5056,5189,5207,5243,5329,5518,5581,5582,5600,5721,5811,5844,5887,5943,5955,6034,6059,6252,6339,6425,6429,6456,6471,6578,6716,6751,6901,6947,7049,7262,7302,7655,7810,7849,7911,8189,8735,8822,9056,9398,9600,9737,9826,9886,9979,9987,10048,10230,10455,10494,10533,10534,10809,11045,11073,11101,11932,11933,12144,12292,12571,12667,13026,13032,13058,13195,13262,13337,13369,13702,13802,13847,13905,14190,14390,14670,14793,15206,15620,15745,15857,16104,16223,16287,16295,16339,16511,16607,16626,16680,16686,16971,17130,17291,17292,17469,17601,17615,17677,17749)]

x.pca = pca(X, center = T, scale = T)
x.pca$variates$X[,1] <- -x.pca$variates$X[,1]
#x.pca$variates$X[,2] <- -x.pca$variates$X[,2]
plot(x.pca)

pdf("spls_PCA_mTgene_abundance.pdf", 5,5)
# plot final sPCA samples for first two components
plotIndiv(x.pca, comp = c(1, 2), ind.names = TRUE, 
          group = rownames(X),  # use class to colour each sample
          legend = FALSE, title = "")
dev.off()
# plot variables against the sPCA components
plotVar(x.pca, comp = c(1, 2), var.names = TRUE,  
        title = 'metaData PCA comp 1 - 2')
# plot samples and variables against the sPCA components
pdf("spls_sPCA_mTgene_abudance_biplot.pdf", 8, 6)
#png("spls_PCA_mTgene_abudance_biplot.png", 
#    width = 6, height = 5, units = "in", res = 400)
biplot(x.pca, cex = 1,  pch.size = 1, ind.names.size = 2.5, cutoff = 0.3,
       xlabs = paste(ncol(X), 1:nrow(X)), #simplify names
       var.arrow.size = 0.01, var.names.size = 0.5, var.arrow.col = 'grey90',
       legend = FALSE, alpha=0.9, var.arrow.length= 0.1,
       group = rownames(X), #multidrug$cell.line$Class,  # colour by sample class
       title = 'mG PCA comp 1 - 2')
dev.off()
##### the end of mT gene abundance #####
################################################################################
############### draw the PCA & sPCA of metabolomes #################
## prepare  metabolome data
rm(list=ls())
xdf=read.csv("all_metabolome.csv", header = T)
xdf=xdf[,c(2,15:68)]

# First column is the factor column
cid <- xdf[[1]]
numeric_data <- xdf[, -1]
# Extract the base date name (before the last underscore) for each numeric column
base_names <- sub("(_[^_]+)$", "", colnames(numeric_data))
# Find unique dates
unique_dates <- unique(base_names)

# Compute row means for each group of replicate columns
aggregated_data <- sapply(unique_dates, function(date) {
  cols <- which(base_names == date)
  rowMeans(numeric_data[, cols, drop = FALSE])
})

# Convert to data frame
aggregated_df <- data.frame(CID = cid, aggregated_data)

# Optional: give cleaner column names
colnames(aggregated_df)[-1] <- unique_dates
rownames(aggregated_df) = aggregated_df[[1]]
aggregated_df = aggregated_df[-1]
# View the result
head(aggregated_df)
df = as.data.frame(t(aggregated_df))
rownames(df)=gsub('X([679])','X0\\1', rownames(df))
rownames(df)=gsub('X','mM', rownames(df))
rownames(df)=gsub('_','\\.', rownames(df))
rownames(df)
# remove NAs from columns
df <- df[, colSums(is.na(df)) < nrow(df)]
df[is.na(df)] <- 0

X = df
X = X[,c(-18631,-18870)]
x.pca = pca(X, center = T, scale = T)
x.pca$variates$X[,1] <- -x.pca$variates$X[,1]
x.pca$variates$X[,2] <- -x.pca$variates$X[,2]
plot(x.pca)

pdf("spls_PCA_mMetabolome_abundance.pdf", 6,4)
# plot final sPCA samples for first two components
plotIndiv(x.pca, comp = c(1, 2), ind.names = TRUE, 
          group = rownames(X),  # use class to colour each sample
          legend = FALSE, title = "")
dev.off()
# plot variables against the sPCA components
plotVar(x.pca, comp = c(1, 2), var.names = TRUE,  
        title = 'metaData PCA comp 1 - 2')
# plot samples and variables against the sPCA components
pdf("spls_sPCA_mMetabolome_abudance_biplot.pdf", 8, 6)
#png("spls_PCA_mTgene_abudance_biplot.png", 
#    width = 6, height = 5, units = "in", res = 400)
biplot(x.pca, cex = 1,  pch.size = 1, ind.names.size = 2, cutoff = 0.4,
       xlabs = paste(ncol(X), 1:nrow(X)), #simplify names
       var.arrow.size = 0.01, var.names.size = 0.5, 
       var.arrow.col = 'grey90', var.arrow.length = 0.1,
       legend = FALSE, alpha=0.9,
       group = rownames(X), #multidrug$cell.line$Class,  # colour by sample class
       title = 'mG PCA comp 1 - 2')
dev.off()
##### the end of metabolome PCA #####
################################################################################
##### metabloome WGCNA #####
rm(list = ls())
library(WGCNA)
library(pheatmap)
library(colorspace)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
ALLOW_WGCNA_THREADS=14
collectGarbage()

# prepare metabolome data
mm0=read.csv("all_metabolome.tsv", sep = "\t", header = T)
mm1=mm0[,c(2, 15:68)]
rownames(mm1)=mm1$ID
mm1=mm1[,-1]
newcols=unique(gsub("_\\d$","", colnames(mm1)))
mm2=as.data.frame(mm1[,1])
rownames(mm2)=rownames(mm1)
for ( i in newcols){
  tdf=mm1[, grepl(i,colnames(mm1))]
  mm2[[i]]=rowMeans(tdf)
}
mm2=mm2[,-1]
colnames(mm2)=gsub("^X","S", colnames(mm2))
colnames(mm2)=gsub("S([679])","S0\\1",colnames(mm2), perl = TRUE)
colnames(mm2)=gsub("_([123])$","-\\1",colnames(mm2),perl = TRUE) 
mm2=mm2[rowSums(mm2) >0,]

otu = as.data.frame(t(mm2))
# start WGCNA

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
abline(h = 500000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 500000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = otu[keepSamples, ] # datExpr = sqrt(otu[keepSamples, ])
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##### 1.c Loading clinical trait data
mdf = read.csv('water_chemistry.csv', header = T, row.names = 1, check.names = F)
dim(mdf) # confirm the dimension of data
rownames(mdf)=gsub("mG19_","S", rownames(mdf))
rownames(mdf)=gsub("\\.","_", rownames(mdf))
rownames(mdf)=gsub("_([123])$","-\\1", rownames(mdf))

all(rownames(mdf)==rownames(datExpr))
datTraits=mdf
collectGarbage()

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
net = blockwiseModules(datExpr, power = 8,
                       #minModuleSize = 20, #TOMType = "unsigned",
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "2019metabolome",
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
     file = "taihu2019_metabolome_networkConstruction-auto.RData");
#load("taihu2019_metabolome_networkConstruction-auto.RData")
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
textMatrix = paste(signif(moduleTraitCor, 2), "|",
                   signif(moduleTraitPvalue, 1), sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(4, 7, 2, 1)+0.1);
# Display the correlation values within a heatmap plot
pdf("Taihu2019_metabolome_wgcna_module-trait_corr_heatmap2.pdf",3,6)
labeledHeatmap(Matrix = moduleTraitCor,
      xLabels = names(datTraits), yLabels = names(MEs),
      ySymbols = names(MEs),      colorLabels = FALSE, 
      xLabelsAngle=90,            colors = diverging_hcl(30, palette='Blue-Red3'), #blueWhiteRed
      textMatrix = textMatrix,    setStdMargins = FALSE,
      cex.text = 0.15,             zlim = c(-1,1),
      cex.lab.x = 0.2,            cex.lab.y = 0.2,
      main = paste("Module-trait relationships"))
pheatmap(moduleTraitCor, border_color = "white",fontsize = 3,
         color = diverging_hcl(30, palette='Blue-Red3'))
dev.off()

#--  analyze the compounds in each significant module
#--  this ko pathway file is needed

# #-- modules associated with temperature
picked_module = c("maroon","darkgreen", "darkturquoise", "mediumorchid", "coral1", "paleturquoise", "lightyellow")
module_dfList <- list()
pdf("wgcan_6module_comp-T_negposCorr-KOpathways_heatmap_4x6.pdf",4,6)
#-- modules associated with temperature
picked_module = c("darkgrey", "maroon","brown")
module_dfList <- list()
pdf("wgcanmM_3module_comp-T_posCorr-KOpathways_heatmap.pdf",8,8)

# #-- modules associated with NH4-N
# picked_module = c("paleturquoise", "white","darkturquoise")
# module_dfList <- list()
# pdf("wgcan_3module_comp-NH4-N_corr_KOpaths_heatmap.pdf",8,8)

# #-- modules associated with NO3_N
# picked_module = c("ivory", "paleturquoise","darkgrey","darkolivegreen","mediumpurple2","lightcyan1")
# module_dfList <- list()
# pdf("wgcan_6module_comp-NO3-N_corr_KOpaths_heatmap.pdf",8,8)

# #-- modules associated with Pi
# picked_module = c("lightcyan", "midnightblue","mediumpurple2")
# module_dfList <- list()
# pdf("wgcan_3module_comp-Pi_corr_KOpaths_heatmap.pdf",8,8)

#-- modules associated with PC
# picked_module = c("violet", "lightcyan1","skyblue4")
# module_dfList <- list()
# pdf("wgcan_3module_comp-PC_corr_KOpaths_heatmap.pdf",8,8)

for (i in 1:length(picked_module)){
  module_comps = names(datExpr)[moduleColors == picked_module[i]]
  module_df=mm0[mm0$ID %in% module_comps, c(2:3, 12,14)]
  rownames(module_df)=module_df$ID
  merge_module_df = merge(module_df, mm2, by = "row.names")
  module_dfList[[i]]=merge_module_df
  pheatmap(log10(merge_module_df[,-1:-5] + 1), cluster_cols=FALSE, fontsize = 3,
         show_rownames = FALSE, main=paste(picked_module[i],"log10value",sep=" "))
}
combined_module_df <- do.call(rbind, module_dfList)
combined_module_df2=combined_module_df[combined_module_df$Metabolite != '-',]
rownames(combined_module_df2)=combined_module_df2$Metabolite
pheatmap(log10(combined_module_df2[,-1:-5] + 1),
         cluster_cols=FALSE, fontsize = 3, show_rownames = TRUE)
#dev.off()

search_terms=combined_module_df2$KEGG.Compound.ID[combined_module_df2$KEGG.Compound.ID != '-']
search_list=unique(unlist(strsplit(search_terms, ";")))

kopath_df=read.csv("Table6_mGmTmM_CnumKOpathway-full.tsv",sep = "\t",header = T)
df_chr <- as.data.frame(lapply(kopath_df, function(x) as.character(x)), stringsAsFactors = FALSE)
df_chr[is.na(df_chr)] <- ""
pattern <- paste0("\\b(", paste(search_list, collapse = "|"), ")\\b")
match_rows <- apply(df_chr, 1, function(row) any(grepl(pattern, row)))
matched_kopath_df <- kopath_df[match_rows, ]

gtgenes=read.csv("Table5-2_mGmT_tpm_bytimeDate.tsv",sep = "\t",header = T)
merge_gtgenes=merge(matched_kopath_df[,c(1,9,10)], gtgenes, by.x = "qseqid", by.y = "Gene")
rownames(merge_gtgenes)=paste(merge_gtgenes$qseqid,merge_gtgenes$ko_symbol, merge_gtgenes$ko_name, sep = "|") 
merge_gtgenes=merge_gtgenes[,-1:-3]
merge_gtgenes=merge_gtgenes[,sort(colnames(merge_gtgenes))]
merge_gtgenes[is.na(merge_gtgenes)]=0
merge_gtgenes=merge_gtgenes[rowSums(merge_gtgenes==0)<15,]
pheatmap(log10(merge_gtgenes[,grepl("^mG",colnames(merge_gtgenes))] +1), cluster_rows = FALSE,
         cluster_cols=FALSE, fontsize = 3, show_rownames = TRUE, main = 'mG_genes_tpm')
pheatmap(log10(merge_gtgenes[,grepl("^mT",colnames(merge_gtgenes))] +1),
         cluster_cols=FALSE, fontsize = 3, show_rownames = TRUE, main = 'mT_genes_tpm')
#dev.off()

#--plotting--KO pathways and enrichment--#

kmdf=matched_kopath_df

library(tidyverse)
library(stringr)
# Step 1: Extract reaction equations
eq_cols <- grep("^rn_equation_", colnames(kmdf), value = TRUE)
# Step 2: Tidy and extract C-numbers from reactions
eq_long <- kmdf %>%
  select(ko_pathway, all_of(eq_cols)) %>%
  pivot_longer(cols = starts_with("rn_equation_"), names_to = "eq_col", values_to = "reaction") %>%
  filter(!is.na(ko_pathway), !is.na(reaction)) %>%
  mutate(C_numbers = str_extract_all(reaction, "C\\d{5}")) %>%
  unnest(C_numbers)
# Step 3: Count total C-numbers and unique reactions per KO-pathway
summary_df <- eq_long %>%
  group_by(ko_pathway) %>%
  summarise(
    total_C = n_distinct(C_numbers),
    total_reactions = n_distinct(reaction)
  )
# Step 4: Count hits from search list
search_list 
hit_df <- eq_long %>%
  filter(C_numbers %in% search_list) %>%
  group_by(ko_pathway) %>%
  summarise(
    hit_C = n_distinct(C_numbers),
    hit_reactions = n_distinct(reaction)
  )
# Step 5: Merge total and hit summaries
k_summary <- summary_df %>%
  left_join(hit_df, by = "ko_pathway") %>%
  replace_na(list(hit_C = 0, hit_reactions = 0))
# Step 6: Hypergeometric test
# Total unique C-numbers and reactions
total_C_all <- length(unique(eq_long$C_numbers))
total_reactions_all <- length(unique(eq_long$reaction))
# Calculate p-values
k_summary <- k_summary %>%
  rowwise() %>%
  mutate(
    pval_C = phyper(hit_C - 1, total_C_all, total_C_all - total_C, total_C, lower.tail = FALSE),
    pval_reactions = phyper(hit_reactions - 1, total_reactions_all, total_reactions_all - total_reactions, total_reactions, lower.tail = FALSE),
    sig_C = ifelse(pval_C < 0.05, "Significant", "Not significant"),
    sig_reaction = ifelse(pval_reactions < 0.05, "Significant", "Not significant")
  ) %>%
  ungroup()
# Step 7: Prepare data for plotting (long format)
# Define a common KO-pathway ordering (based on total_reactions, for example)
ko_order <- k_summary %>%
  arrange(total_reactions) %>%
  pull(ko_pathway)
# Prepare c_total and c_hits
c_total <- k_summary %>%
  select(ko_pathway, value = total_C) %>%
  mutate(ko_pathway = factor(ko_pathway, levels = ko_order))
c_hits <- k_summary %>%
  select(ko_pathway, value = hit_C, significance = sig_C) %>%
  mutate(ko_pathway = factor(ko_pathway, levels = ko_order))
# Prepare r_total and r_hits
r_total <- k_summary %>%
  select(ko_pathway, value = total_reactions) %>%
  mutate(ko_pathway = factor(ko_pathway, levels = ko_order))
r_hits <- k_summary %>%
  select(ko_pathway, value = hit_reactions, significance = sig_reaction) %>%
  mutate(ko_pathway = factor(ko_pathway, levels = ko_order))
# Step 8: Plot C-number bar chart (horizontal, overlay)

#pdf("wgcna_6modul-T_corr_Cnumber-KOpathway.pdf", 8,4)
c_total <- c_total %>%
  mutate(ko_pathway = fct_reorder(ko_pathway, value, .desc = FALSE))
c_hits <- c_hits %>%
  mutate(ko_pathway = factor(ko_pathway, levels = levels(c_total$ko_pathway)))
p1 <- ggplot() +
  geom_col(data = c_total, aes(x = ko_pathway, y = value), fill = "gray80", width = 0.8) +
  geom_col(
    data = c_hits,
    aes(x = ko_pathway, y = value, fill = significance),
    position = "identity",
    #alpha = 0.7,
    width = 0.6
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("Significant" = "coral", "Not significant" = "steelblue"),
    name = "C-number Hits"
  ) +
  labs(
    title = "C-number Hits per KO-Pathway",
    x = "KO-Pathway",
    y = "C-number Count"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
# Step 9: Plot Reaction bar chart (horizontal, overlay)
r_total <- r_total %>%
  mutate(ko_pathway = fct_reorder(ko_pathway, value, .desc = FALSE))
r_hits <- r_hits %>%
  mutate(ko_pathway = factor(ko_pathway, levels = levels(r_total$ko_pathway)))
p2 <- ggplot() +
  geom_col(data = r_total, aes(x = ko_pathway, y = value), fill = "gray80", width = 0.8) +
  geom_col(
    data = r_hits,
    aes(x = ko_pathway, y = value, fill = significance),
    position = "identity",
    #alpha = 0.7,
    width = 0.6
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("Significant" = "coral", "Not significant" = "steelblue"),
    name = "Reaction Hits"
  ) +
  labs(
    title = "Reaction Hits per KO-Pathway",
    x = "KO-Pathway",
    y = "Reaction Count"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
# Step 10: Print plots
print(p1)
print(p2)

dev.off()

# 3.b Gene relationship to trait and important modules: Gene Significance and Module
# Membership
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$T); # ====>>> here to SELECT TRAITs for regression
names(weight) = "T" # ====>>> here to SELECT TRAITs for regression
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
pdf("taihu2019_metabolome-wgcna_modules-T_corr.pdf",6,6)
selec_modnames=c("darkgrey","maroon","darkturquoise", "darkgreen","brown","mediumorchid",
                 "coral1","paleturquoise","mediumpurple2")
for (i in 1:length(selec_modnames)){
  #i=1
  #module = modNames[i] ## ====>>>>> pick the MODULE here
  #column = match(module, modNames);
  module = selec_modnames[i] ## ====>>>>> pick the MODULE here
  column = match(module, selec_modnames);
  moduleGenes = moduleColors==module;
  #sizeGrWindow(7, 7);
  #par(mfrow = c(1,1));
  plott<-verboseScatterplot(geneModuleMembership[moduleGenes, column],
                            geneTraitSignificance[moduleGenes, 1],
                            xlab = paste("Module Membership in", module, "module"),
                            ylab = "metabolome significance for Temperature",
                            main = paste("Module membership vs. T significance\n"),
                            abline = T, abline.color = module,cex.axis = 1.2, 
                            cex = 1.2, pch=16, cex.main = 0.9, cex.lab = 1.2,
                            displayAsZero=0.2, col = module)
}
dev.off()

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

TOM = TOMsimilarityFromExpr(datExpr, power = 7);

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
                                 edgeFile = paste("F-Diel_edges_cutoff37_", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("F-Diel_nodes_cutoff37_", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.37,# 
                                 nodeNames = modProbes,
                                 includeColNames = TRUE,
                                 #altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]
  )
}
# the end of metabolome WGCNA
################################################################################
##### mG gene WGCNA #####
rm(list = ls())
library(WGCNA)
library(pheatmap)
library(colorspace)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
ALLOW_WGCNA_THREADS=14
collectGarbage()

# prepare metagenome data
mg0=read.csv("mG19_merged_genefamily_tpm_sorted.tsv", sep = "\t", header = T)
mg0=mg0[rowSums(is.na(mg0))<7,]
otu = as.data.frame(t(mg0))

# start WGCNA

# hellinger transformation
#sotu=sqrt(otu)
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
abline(h = 12000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 12000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = otu[keepSamples, ] # datExpr = sqrt(otu[keepSamples, ])
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##### 1.c Loading clinical trait data
mdf = read.csv('water_chemistry.csv', header = T, row.names = 1, check.names = F)
dim(mdf) # confirm the dimension of data
mdf=mdf[c(-3,-10),]

rownames(mdf)==rownames(datExpr)
datTraits=mdf
collectGarbage()

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
net = blockwiseModules(datExpr, power = 8,
                       #minModuleSize = 20, #TOMType = "unsigned",
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "2019metabolome",
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
     file = "taihu2019_mG_networkConstruction-auto.RData");

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
textMatrix = paste(signif(moduleTraitCor, 2), "|",
                   signif(moduleTraitPvalue, 1), sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(4, 7, 2, 1)+0.1);
# Display the correlation values within a heatmap plot
pdf("Taihu2019_mG-gene_wgcna_module-trait_corr_heatmap_3x6.pdf",3,6)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits), yLabels = names(MEs),
               ySymbols = names(MEs),      colorLabels = FALSE, 
               xLabelsAngle=90,            colors = diverging_hcl(30, palette='Blue-Red3'), #blueWhiteRed
               textMatrix = textMatrix,    setStdMargins = FALSE,
               cex.text = 0.15,             zlim = c(-1,1),
               cex.lab.x = 0.2,            cex.lab.y = 0.2,
               main = paste("Module-trait relationships"))
pheatmap(moduleTraitCor, border_color = "white",fontsize = 3,
         color = diverging_hcl(30, palette='Blue-Red3'))
dev.off()

#--  analyze the compounds in each significant module
#--  this ko pathway file is needed

# #--  modules of mG genes associated with temperature
picked_module = c("red", "blue","greenyellow", "midnightblue")
module_dfList <- list()
pdf("wgcanMG_4modules-T_corr-KOpathways_heatmap_3x6.pdf",3,6)

# #--  modules of mG genes associated with PC
# picked_module = c("salmon")
# module_dfList <- list()
# pdf("wgcanMG_1modules-PC_corr-KOpathways_heatmap.pdf",8,8)

#--  modules of mG genes associated with NO3-N
# picked_module = c("magenta","purple")
# module_dfList <- list()
# pdf("wgcanMG_1modules-No3-N_corr-KOpathways_heatmap.pdf",8,8)

for (i in 1:length(picked_module)){
  module_genes = names(datExpr)[moduleColors == picked_module[i]]
  module_df=datExpr[,colnames(datExpr) %in% module_genes]
  module_df=as.data.frame(t(module_df))
  module_dfList[[i]]=module_df
  pheatmap(log10(module_df+ 1), cluster_cols=FALSE, fontsize = 5,
           show_rownames = FALSE, main=paste(picked_module[i],"log10value",sep=" "))
}
combined_module_df <- do.call(rbind, module_dfList)
pheatmap(log10(combined_module_df + 1), main = "mGgenes-T_all module log10value",
         cluster_cols=FALSE, fontsize = 5, show_rownames = FALSE)
#dev.off()

search_terms=rownames(combined_module_df)
search_list=unique(search_terms)

# Required libraries
library(dplyr)
library(ggplot2)
library(forcats)

# Example DE gene list (replace with your actual list)
DE_genes <- search_list

# Read your dataframe (if not already in R)
df <- read.csv("Table6_mGmTmM_CnumKOpathway-full.tsv",sep = "\t",header = T)
df=df[,1:11]
# Filter out rows with NA in ko_name
df_clean <- df %>%
  filter(!is.na(ko_name)) %>%               # remove NA pathways
  filter(ko_name != "hypothetical protein") # remove uninformative pathway name
# Define background gene universe and overlap with DE genes
all_genes <- unique(df_clean$qseqid)
DE_genes <- intersect(DE_genes, all_genes)  # restrict DE to genes present in data

# Map each pathway to its gene set
pathway_gene_map <- df_clean %>%
  group_by(ko_name) %>%
  summarize(genes_in_pathway = list(unique(qseqid)), .groups = "drop")

# Perform hypergeometric test
pathway_stats <- pathway_gene_map %>%
  rowwise() %>%
  mutate(
    M = length(all_genes),                        # total background genes
    n = length(DE_genes),                         # number of DE genes
    K = length(genes_in_pathway),                 # number of genes in this pathway
    k = sum(genes_in_pathway %in% DE_genes),      # DE genes in this pathway
    pval = phyper(k - 1, K, M - K, n, lower.tail = FALSE)  # hypergeometric p-value
  ) %>%
  ungroup() %>%
  mutate(
    padj = p.adjust(pval, method = "BH"),         # adjust for multiple testing
    sig_flag = padj < 0.05                        # flag significant pathways
  ) %>%
  filter(k > 0) %>%                               # keep only pathways with DE hits
  arrange(desc(K)) %>%                            # sort by total gene count
  mutate(ko_name = fct_reorder(ko_name, K))       # reorder factor for plotting

# Plot: horizontal bar chart with overlay
ggplot(pathway_stats, aes(x = ko_name)) +
  geom_bar(aes(y = K), stat = "identity", fill = "gray80", width = 0.9) +  # background bar
  geom_bar(
    aes(y = k, fill = sig_flag),
    stat = "identity",
    width = 0.54,  # 60% width
    position = position_nudge(x = 0)
  ) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "coral")) +
  coord_flip() +
  labs(
    x = "Pathway",
    y = "Number of Genes",
    title = "Pathway Enrichment by Hypergeometric Test",
    subtitle = "Gray: total genes; Orange: significant DE genes (padj < 0.05); Blue: non-significant",
    fill = "Significant",
    caption = "Overlay bars scaled to 60% width"
  ) +
  theme_minimal(base_size = 5) +
  theme(
    axis.text.y = element_text(size = 5)  # <- smaller font for y-axis tick labels (pathways)
  )
dev.off()

# 3.b Gene relationship to trait and important modules: Gene Significance and Module
# Membership
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$T); # ====>>> here to SELECT TRAITs for regression
names(weight) = "T" # ====>>> here to SELECT TRAITs for regression
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
pdf("taihu2019_metabolome-wgcna_modules-T_corr.pdf",6,6)
selec_modnames=c("darkgrey","maroon","darkturquoise", "darkgreen","brown","mediumorchid",
                 "coral1","paleturquoise","mediumpurple2")
for (i in 1:length(selec_modnames)){
  #i=1
  #module = modNames[i] ## ====>>>>> pick the MODULE here
  #column = match(module, modNames);
  module = selec_modnames[i] ## ====>>>>> pick the MODULE here
  column = match(module, selec_modnames);
  moduleGenes = moduleColors==module;
  #sizeGrWindow(7, 7);
  #par(mfrow = c(1,1));
  plott<-verboseScatterplot(geneModuleMembership[moduleGenes, column],
                            geneTraitSignificance[moduleGenes, 1],
                            xlab = paste("Module Membership in", module, "module"),
                            ylab = "metabolome significance for Temperature",
                            main = paste("Module membership vs. T significance\n"),
                            abline = T, abline.color = module,cex.axis = 1.2, 
                            cex = 1.2, pch=16, cex.main = 0.9, cex.lab = 1.2,
                            displayAsZero=0.2, col = module)
}
dev.off()

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

TOM = TOMsimilarityFromExpr(datExpr, power = 7);

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
                                 edgeFile = paste("F-Diel_edges_cutoff37_", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("F-Diel_nodes_cutoff37_", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.37,# 
                                 nodeNames = modProbes,
                                 includeColNames = TRUE,
                                 #altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]
  )
}
# the end of mG WGCNA
##################################### WGCNA Analysis ############################################
##### mT gene WGCNA #####
rm(list = ls())
library(WGCNA)
library(pheatmap)
library("colorspace")
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
ALLOW_WGCNA_THREADS=14
collectGarbage()

# prepare metagenome data
mt0=read.csv("Table5_mGmT_geneAb_across-TimeSite.tsv", sep = "\t", header = T)
rownames(mt0)=mt0[,1]
mt0=mt0[,-1:-2]
colnames(mt0)
mt0=mt0[,grepl("_TPM",colnames(mt0))]
mt0=mt0[rowSums(is.na(mt0))<7,]
colnames(mt0)=gsub("X([679])","X0\\1",colnames(mt0), perl = TRUE)
colnames(mt0)=gsub("_TPM","", gsub("X","mT19_", colnames(mt0)))
mt0=mt0[,sort(colnames(mt0))]
otu = as.data.frame(t(mt0))

# start WGCNA

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
abline(h = 1000000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 1000000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = otu[keepSamples, ] # datExpr = sqrt(otu[keepSamples, ])
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##### 1.c Loading clinical trait data
mdf = read.csv('water_chemistry.csv', header = T, row.names = 1, check.names = F)
dim(mdf) # confirm the dimension of data
mdf=mdf[-2,]
rownames(mdf)=gsub("mG","mT", rownames(mdf))

rownames(mdf)==rownames(datExpr)
datTraits=mdf
collectGarbage()

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
net = blockwiseModules(datExpr, power = 8,
                       minModuleSize = 20, TOMType = "signed",
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "2019mT",
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
     file = "taihu2019_mT_networkConstruction-auto.RData");

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
textMatrix = paste(signif(moduleTraitCor, 2), "|",
                   signif(moduleTraitPvalue, 1), sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(4, 7, 2, 1)+0.1);
# Display the correlation values within a heatmap plot
pdf("Taihu2019_mT-gene_wgcna_module-trait_corr_heatmap.pdf",6,6)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits), yLabels = names(MEs),
               ySymbols = names(MEs),      colorLabels = FALSE, 
               xLabelsAngle=90,            colors = diverging_hcl(30, palette='Blue-Red3'), #blueWhiteRed
               textMatrix = textMatrix,    setStdMargins = FALSE,
               cex.text = 0.3,             zlim = c(-1,1),
               cex.lab.x = 0.5,            cex.lab.y = 0.5,
               main = paste("Module-trait relationships"))
pheatmap(moduleTraitCor, border_color = "white",fontsize = 8,
         color = diverging_hcl(30, palette='Blue-Red3'))
dev.off()

# #--  analyze the compounds in each significant module
# # modules associated with PC
# picked_module = c("steelblue","darkgrey","yellowgreen","magenta")
# module_dfList <- list()
# pdf("wgcanMT_4modules-PC_corr-KOpathways_heatmap.pdf",8,8)

# modules associated with PC
picked_module = c("grey60")
module_dfList <- list()
pdf("wgcanMT_4modules-NO3-N_corr-KOpathways_heatmap.pdf",8,8)

for (i in 1:length(picked_module)){
  module_genes = names(datExpr)[moduleColors == picked_module[i]]
  module_df=datExpr[,colnames(datExpr) %in% module_genes]
  module_df=as.data.frame(t(module_df))
  module_dfList[[i]]=module_df
  pheatmap(log10(module_df+ 1), cluster_cols=FALSE, fontsize = 8,
           show_rownames = FALSE, main=paste(picked_module[i],"log10value",sep=" "))
}
combined_module_df <- do.call(rbind, module_dfList)
pheatmap(log10(combined_module_df + 1), main = "mT-T_all module log10value",
         cluster_cols=FALSE, fontsize = 8, show_rownames = FALSE)
#dev.off()

search_terms=rownames(combined_module_df)
search_list=unique(search_terms)

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(forcats)

# Example DE gene list (replace with your actual list)
DE_genes <- search_list

# Read your dataframe (if not already in R)
df <- read.csv("Table6_mGmTmM_CnumKOpathway-full.tsv",sep = "\t",header = T)
df=df[,1:11]
# Step 1: Filter out rows without pathway info or with uninformative label
df_clean <- df %>%
  filter(!is.na(ko_name)) %>%                      # remove NA pathways
  filter(ko_name != "hypothetical protein")        # remove uninformative pathway name

# Step 2: Define the background gene universe
all_genes <- unique(df_clean$qseqid)
DE_genes <- intersect(DE_genes, all_genes)  # ensure DE genes are in data

# Step 3: Group genes by pathway
pathway_gene_map <- df_clean %>%
  group_by(ko_name) %>%
  summarize(genes_in_pathway = list(unique(qseqid)), .groups = "drop")

# Step 4: Perform hypergeometric test
pathway_stats <- pathway_gene_map %>%
  rowwise() %>%
  mutate(
    M = length(all_genes),                        # total genes
    n = length(DE_genes),                         # DE genes
    K = length(genes_in_pathway),                 # genes in pathway
    k = sum(genes_in_pathway %in% DE_genes),      # DE genes in pathway
    pval = phyper(k - 1, K, M - K, n, lower.tail = FALSE)  # enrichment p-value
  ) %>%
  ungroup() %>%
  mutate(
    padj = p.adjust(pval, method = "BH"),
    sig_flag = padj < 0.05
  ) %>%
  filter(k > 0) %>%                      # keep only pathways with at least 1 DE gene
  arrange(desc(K)) %>%                   # sort by total genes per pathway
  mutate(ko_name = fct_reorder(ko_name, K))

# Step 5: Plot
ggplot(pathway_stats, aes(x = ko_name)) +
  geom_bar(aes(y = K), stat = "identity", fill = "gray80", width = 0.9) +  # total genes
  geom_bar(
    aes(y = k, fill = sig_flag),
    stat = "identity",
    width = 0.54,  # 60% width overlay
    position = position_nudge(x = 0)
  ) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "coral")) +
  coord_flip() +
  labs(
    x = "Pathway",
    y = "Number of Genes",
    title = "Pathway Enrichment by Hypergeometric Test",
    subtitle = "Gray: total genes; Orange: significant (padj < 0.05); Blue: not significant",
    fill = "Significant",
    caption = "Overlay bars are scaled to 60% width"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 6)  # smaller y-axis labels
  )

dev.off()

# 3.b Gene relationship to trait and important modules: Gene Significance and Module
# Membership
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$T); # ====>>> here to SELECT TRAITs for regression
names(weight) = "T" # ====>>> here to SELECT TRAITs for regression
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
pdf("taihu2019_metabolome-wgcna_modules-T_corr.pdf",6,6)
selec_modnames=c("darkgrey","maroon","darkturquoise", "darkgreen","brown","mediumorchid",
                 "coral1","paleturquoise","mediumpurple2")
for (i in 1:length(selec_modnames)){
  #i=1
  #module = modNames[i] ## ====>>>>> pick the MODULE here
  #column = match(module, modNames);
  module = selec_modnames[i] ## ====>>>>> pick the MODULE here
  column = match(module, selec_modnames);
  moduleGenes = moduleColors==module;
  #sizeGrWindow(7, 7);
  #par(mfrow = c(1,1));
  plott<-verboseScatterplot(geneModuleMembership[moduleGenes, column],
                            geneTraitSignificance[moduleGenes, 1],
                            xlab = paste("Module Membership in", module, "module"),
                            ylab = "metabolome significance for Temperature",
                            main = paste("Module membership vs. T significance\n"),
                            abline = T, abline.color = module,cex.axis = 1.2, 
                            cex = 1.2, pch=16, cex.main = 0.9, cex.lab = 1.2,
                            displayAsZero=0.2, col = module)
}
dev.off()

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

TOM = TOMsimilarityFromExpr(datExpr, power = 7);

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
                                 edgeFile = paste("F-Diel_edges_cutoff37_", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("F-Diel_nodes_cutoff37_", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.37,# 
                                 nodeNames = modProbes,
                                 includeColNames = TRUE,
                                 #altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]
  )
}
##################################### WGCNA Analysis ############################################



