################################################################################
##################   multiomics integration ####################
rm(list = ls())
library(mixOmics) # import the mixOmics library
set.seed(123) # for reproducibility, remove for normal use
load('mixomics_data/nmt_data_processed.RData')
X1 <- data$rna # select three of the five dataframes to explore
X2 <- data$met_genebody
X3 <- data$acc_genebody
# compile these into a single X object
X <- list(rna = X1, methylation = X2, accessibility = X3) 
lapply(X, dim) # check dimensions

## initial analysis
# select arbitrary values of features to keep
list.keepX = c(25, 25)
list.keepY = c(25, 25)

# generate three pairwise PLS models
pls1 <- spls(X[["rna"]], X[["methylation"]], 
             keepX = list.keepX, keepY = list.keepY)
pls2 <- spls(X[["rna"]], X[["accessibility"]], 
             keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(X[["methylation"]], X[["accessibility"]], 
             keepX = list.keepX, keepY = list.keepY)

# plot features of first PLS
par(mfrow=c(1,3))
plotVar(pls1, cutoff = 0.5, title = "(a) RNA vs Methylation",
        legend = c("RNA", "Methylation"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2,2),
        col = c('darkorchid', 'lightgreen'))

# plot features of second PLS
plotVar(pls2, cutoff = 0.5, title = "(b) RNA vs Accessibility",
        legend = c("RNA", "Accessibility"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2,2),
        col = c('darkorchid', 'lightgreen'))

# plot features of third PLS
plotVar(pls3, cutoff = 0.5, title = "(c) Methylation vs Accessibility",
        legend = c("Methylation", "Accessibility"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2,2),
        col = c('darkorchid', 'lightgreen'))
# calculate correlation of RNA and methylation
cor(pls1$variates$X, pls1$variates$Y)
# calculate correlation of RNA and accessibility
cor(pls2$variates$X, pls2$variates$Y)
# calculate correlation of methylation and accessibility
cor(pls3$variates$X, pls3$variates$Y) 

## Initial Multiblock sPLS model
# for square matrix filled with 0.5s
design = matrix(0.5, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) = 0 # set diagonal to 0s

basic.mbspls.model = block.spls(X, indY = 1, # generate basic model
                                ncomp = 5, 
                                design = design)
print(design)
# set up ncomp and number of features
choice.ncomp <- 3 
choice.keepX <- list(rna = rep(50, 3), # 50 features per each component 
                     methylation = rep(50, 3), 
                     accessibility = rep(50, 3))
## final MB sPLS model
# generate final model using "tuned" parameters
final.mbspls.model = block.spls(X, indY = 1,  
                                ncomp = choice.ncomp, 
                                keepX = choice.keepX,
                                design = design)


## visualization of models: plots
plotIndiv(final.mbspls.model, ind.names = FALSE,
          group = as.factor(cell_metadata$lineage), 
          pch = as.factor(cell_metadata$stage),
          col.per.group = color.mixo(1:8), 
          legend = TRUE, legend.title = 'Lineage', legend.title.pch = 'Stage',
          blocks = 2)
# plot arrows
symbols <- list(rna = 1, methylation = 6, accessibility = 7)
samples <- sample(cell_metadata$lineage, 30, replace = FALSE)
# Original vector
vec <- cell_metadata$lineage  # or your actual vector
selected_indices <- sample(1:826, 10)
samples <- rep(0, 826)
samples[selected_indices] <- vec[selected_indices]
plotArrow(final.mbspls.model, ind.names = FALSE,
          group = as.factor(samples),
          pch = symbols, pch.size = 1)
# correlation circle plot
plotVar(final.mbspls.model, var.names = FALSE,
        legend = TRUE, cutoff = 0.5,
        pch = c(0,1,2))
# circos plot
circosPlot(final.mbspls.model, 
           group = cell_metadata$lineage, 
           cutoff = 0.8,
           Y.name = 'rna')

####### end of multiomics integration. ######
################################################################################
##################   sPLS for mG taxa and metadata ####################
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

## prepare metadata
metadata = read.delim('water_chemistry.csv', header = T, row.names = 1,
                      sep = ',', check.names = F)
dim(metadata) # confirm the dimension of data
rownames(metadata) = gsub('mG19_', 'mG', rownames(metadata))
metadata = metadata[,c(-13:-14)]

rownames(df)==rownames(metadata)

## load library
library(mixOmics) # import the mixOmics library
set.seed(5249) # for reproducibility, remove for normal use

X = df
Y = metadata
dim(X)
dim(Y)


# initial analysis
dev.off()
par(mfrow=c(1,2))
pca.gene <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.gene$variates$X[,1] <- -pca.gene$variates$X[,1]
pca.gene$variates$X[,2] <- -pca.gene$variates$X[,2]

pca.clinic <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.gene)
plot(pca.clinic)
plotIndiv(pca.gene, comp = c(1, 2), 
          group = rownames(X),#  #liver.toxicity$treatment[, 4], 
          ind.names = rownames(X), #liver.toxicity$treatment[, 3], 
          legend = FALSE, title = 'mG genus')
plotIndiv(pca.clinic, comp = c(1, 2), 
          group = rownames(Y), # liver.toxicity$treatment[, 4], 
          ind.names = rownames(Y), # liver.toxicity$treatment[, 3], 
          legend = FALSE, title = 'metadata')

## sPLS tuning does not work for some reason, so it's skipped
# # initial sPLS
# spls.liver <- spls(X = X, Y = Y, ncomp = 5, mode = 'regression')
# ## tuning sPLS
# # set ncomp: repeated CV tuning of component count
# perf.spls.liver <- perf(spls.liver, validation = 'Mfold',
#                         folds = 10, nrepeat = 3) 
# plot(perf.spls.liver, criterion = 'Q2.total')
# # set keepX--the number of variables: 
# # set range of test values for number of variables to use from X dataframe
# list.keepX <- c(seq(2, 18, 2))
# # set range of test values for number of variables to use from Y dataframe
# list.keepY <- c(2:13) 

spls.result <- spls(X, Y, keepX = c(350, 100), keepY = c(3, 1), max.iter=500)  # run the method 
plotIndiv(spls.result) # plot the samples
plotVar(spls.result)   # plot the variables
selectVar(spls.result, comp = 1)$X$name 
# depict weight assigned to each of these variables
plotLoadings(spls.result, method = 'mean', contrib = 'max') 

plotVar(spls.result, cex = c(2,4), var.names = c(TRUE, TRUE))
# network or CIM
dev.off
color.edge <- color.GreenRed(30)  # set the colours of the connecting lines
# X11() # To open a new window for Rstudio
cim(spls.result, comp = 1:2, xlab = "Physiochemical factor", ylab = "Genus",
    row.cex = 0.003, save = "pdf", name.save = "spls_pairwise_mGgenus_metadata_heatmap"
    )

network(spls.result, comp = 1:2,
        cutoff = 0.5, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        #save = 'png', # save as a png to the current working directory
        name.save = 'sPLS Liver Toxicity Case Study Network Plot')
cim(spls.result, comp = 1:2, xlab = "clinic", ylab = "genes", save = "pdf")

############ the end of sPLS for mG taxa and metadata ##################
################################################################################
###########   sPLS for mG genes and metadata ####################
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

## prepare metadata
metadata = read.delim('water_chemistry.csv', header = T, row.names = 1,
                      sep = ',', check.names = F)
dim(metadata) # confirm the dimension of data
rownames(metadata) = gsub('mG19_', 'mG', rownames(metadata))
metadata = metadata[,c(-13:-14)]
metadata = metadata[(rownames(metadata) != 'mG11.3.1'),]
rownames(df)==rownames(metadata)

## load library
library(mixOmics) # import the mixOmics library
set.seed(5249) # for reproducibility, remove for normal use

X = df
Y = metadata
dim(X)
dim(Y)


# initial analysis
dev.off()
par(mfrow=c(1,2))
pca.gene <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.gene$variates$X[,1] <- -pca.gene$variates$X[,1]
#pca.gene$variates$X[,2] <- -pca.gene$variates$X[,2]

pca.clinic <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.gene)
plot(pca.clinic)
plotIndiv(pca.gene, comp = c(1, 2), 
          group = rownames(X),#  #liver.toxicity$treatment[, 4], 
          ind.names = rownames(X), #liver.toxicity$treatment[, 3], 
          legend = FALSE, title = 'mG genus')
plotIndiv(pca.clinic, comp = c(1, 2), 
          group = rownames(Y), # liver.toxicity$treatment[, 4], 
          ind.names = rownames(Y), # liver.toxicity$treatment[, 3], 
          legend = FALSE, title = 'metadata')

## sPLS tuning does not work for some reason, so it's skipped
# # initial sPLS
# spls.liver <- spls(X = X, Y = Y, ncomp = 5, mode = 'regression')
# ## tuning sPLS
# # set ncomp: repeated CV tuning of component count
# perf.spls.liver <- perf(spls.liver, validation = 'Mfold',
#                         folds = 10, nrepeat = 3) 
# plot(perf.spls.liver, criterion = 'Q2.total')
# # set keepX--the number of variables: 
# # set range of test values for number of variables to use from X dataframe
# list.keepX <- c(seq(2, 18, 2))
# # set range of test values for number of variables to use from Y dataframe
# list.keepY <- c(2:13) 

spls.result <- spls(X, Y, keepX = c(5000, 5000), keepY = c(3, 2), max.iter=500)  # run the method 
plotIndiv(spls.result) # plot the samples
plotVar(spls.result)   # plot the variables
selectVar(spls.result, comp = 1)$X$name 
# depict weight assigned to each of these variables
plotLoadings(spls.result, method = 'mean', contrib = 'max') 

plotVar(spls.result, cex = c(2,4), var.names = c(TRUE, TRUE))
# network or CIM
dev.off
color.edge <- color.GreenRed(30)  # set the colours of the connecting lines
# X11() # To open a new window for Rstudio
cim(spls.result, comp = 1:2, xlab = "Physiochemical factor", ylab = "Genus",
    row.cex = 0.003, save = "pdf", name.save = "spls_pairwise_mGgene_metadata_heatmap"
)

network(spls.result, comp = 1:2,
        cutoff = 0.5, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        #save = 'png', # save as a png to the current working directory
        name.save = 'sPLS Liver Toxicity Case Study Network Plot')
cim(spls.result, comp = 1:2, xlab = "clinic", ylab = "genes", save = "pdf")


X = df
############ the end of sPLS for mG genes and metadata ##################
################################################################################
###########   sPLS for metabolome and metadata ####################
## prepare  metabolome data
rm(list=ls())
xdf=read.csv("all_metabolome.csv", header = T)
xdf=xdf[,c(2,15:68)]

# Assume your dataframe is called `df`
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

## prepare metadata
metadata = read.delim('water_chemistry.csv', header = T, row.names = 1,
                      sep = ',', check.names = F)
dim(metadata) # confirm the dimension of data
rownames(metadata) = gsub('mG19_', 'mM', rownames(metadata))
metadata = metadata[,c(-13:-14)]

rownames(df)==rownames(metadata)

## load library
library(mixOmics) # import the mixOmics library
set.seed(5249) # for reproducibility, remove for normal use

X = df
Y = metadata
dim(X)
dim(Y)

X = X[,c(-18631,-18870)]

# initial analysis
dev.off()
par(mfrow=c(1,2))
pca.gene <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.gene$variates$X[,1] <- -pca.gene$variates$X[,1]
pca.gene$variates$X[,2] <- -pca.gene$variates$X[,2]

pca.clinic <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.gene)
plot(pca.clinic)
plotIndiv(pca.gene, comp = c(1, 2), 
          group = rownames(X),#  #liver.toxicity$treatment[, 4], 
          ind.names = rownames(X), #liver.toxicity$treatment[, 3], 
          legend = FALSE, title = 'mT gene, PCA comp 1 - 2')
plotIndiv(pca.clinic, comp = c(1, 2), 
          group = rownames(Y), # liver.toxicity$treatment[, 4], 
          ind.names = rownames(Y), # liver.toxicity$treatment[, 3], 
          legend = FALSE, title = 'metaData, PCA comp 1 - 2')

## sPLS tuning does not work for some reason, so it's skipped
# # initial sPLS
# spls.liver <- spls(X = X, Y = Y, ncomp = 5, mode = 'regression')
# ## tuning sPLS
# # set ncomp: repeated CV tuning of component count
# perf.spls.liver <- perf(spls.liver, validation = 'Mfold',
#                         folds = 10, nrepeat = 3) 
# plot(perf.spls.liver, criterion = 'Q2.total')
# # set keepX--the number of variables: 
# # set range of test values for number of variables to use from X dataframe
# list.keepX <- c(seq(2, 18, 2))
# # set range of test values for number of variables to use from Y dataframe
# list.keepY <- c(2:13) 

spls.result <- spls(X, Y, keepX = c(5000, 5000), keepY = c(4, 4),
                    mode="canonical", max.iter=500)  # run the method 
plotIndiv(spls.result) # plot the samples
plotVar(spls.result)   # plot the variables
selectVar(spls.result, comp = 1)$X$name 
# depict weight assigned to each of these variables
plotLoadings(spls.result, method = 'mean', contrib = 'max') 

plotVar(spls.result, cex = c(3,4), var.names = c(TRUE, TRUE))
# network or CIM
color.edge <- color.GreenRed(50)  # set the colours of the connecting lines
cim(spls.result, comp = 1:2, xlab = "metadata", ylab = "Genus",
    row.cex = 0.03, cutoff = 0.6
  )

# X11() # To open a new window for Rstudio
network(spls.result, comp = 1:2,
        cutoff = 0.5, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        #save = 'png', # save as a png to the current working directory
        name.save = 'sPLS Liver Toxicity Case Study Network Plot')
############ the end of sPLS for mM metabolome and metadata ##############
################################################################################
###########   sPLS for mT and metadata ####################
rm(list=ls())
xdf=read.delim("summary_table_5_across-time.tsv",sep = '\t')
rownames(xdf) =xdf[[1]]
xdf=xdf[, seq(5,73,4)]
colnames(xdf)=gsub("X([679])","X0\\1", colnames(xdf))
colnames(xdf)=gsub("^X","mT", colnames(xdf))
colnames(xdf)=gsub("_TPM","", colnames(xdf))
# remove NAs from columns
xdf <- xdf[, colSums(is.na(xdf)) < nrow(xdf)]
xdf[is.na(xdf)] <- 0
df = as.data.frame(t(xdf))
df = df[order(rownames(df)),]
rownames(df)

## prepare metadata
metadata = read.delim('water_chemistry.csv', header = T, row.names = 1,
                      sep = ',', check.names = F)
dim(metadata) # confirm the dimension of data
rownames(metadata) = gsub('mG19_', 'mT', rownames(metadata))
metadata = metadata[,c(-13:-14)]

rownames(df)==rownames(metadata)

## load library
library(mixOmics) # import the mixOmics library
set.seed(5249) # for reproducibility, remove for normal use

X = df
Y = metadata
dim(X)
dim(Y)

X = X[,-c(17,143,222,412,446,461,487,765,929,1115,1146,1228,1395,2129,2166,2409,2589,2794,2843,3132,3149,3251,3273,3383,3389,3450,3460,3501,3594,3624,3632,3756,3856,3859,4028,4132,4534,4538,4546,4737,4859,4984,5000,5056,5189,5207,5243,5329,5518,5581,5582,5600,5721,5811,5844,5887,5943,5955,6034,6059,6252,6339,6425,6429,6456,6471,6578,6716,6751,6901,6947,7049,7262,7302,7655,7810,7849,7911,8189,8735,8822,9056,9398,9600,9737,9826,9886,9979,9987,10048,10230,10455,10494,10533,10534,10809,11045,11073,11101,11932,11933,12144,12292,12571,12667,13026,13032,13058,13195,13262,13337,13369,13702,13802,13847,13905,14190,14390,14670,14793,15206,15620,15745,15857,16104,16223,16287,16295,16339,16511,16607,16626,16680,16686,16971,17130,17291,17292,17469,17601,17615,17677,17749)]

# initial analysis
dev.off()
par(mfrow=c(1,2))
pca.gene <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.gene$variates$X[,1] <- -pca.gene$variates$X[,1]
pca.gene$variates$X[,2] <- -pca.gene$variates$X[,2]

pca.clinic <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.gene)
plot(pca.clinic)
plotIndiv(pca.gene, comp = c(1, 2), 
          group = rownames(X),#  #liver.toxicity$treatment[, 4], 
          ind.names = rownames(X), #liver.toxicity$treatment[, 3], 
          legend = FALSE, title = 'mT gene, PCA comp 1 - 2')
plotIndiv(pca.clinic, comp = c(1, 2), 
          group = rownames(Y), # liver.toxicity$treatment[, 4], 
          ind.names = rownames(Y), # liver.toxicity$treatment[, 3], 
          legend = FALSE, title = 'metaData, PCA comp 1 - 2')

## sPLS tuning does not work for some reason, so it's skipped
# # initial sPLS
# spls.liver <- spls(X = X, Y = Y, ncomp = 5, mode = 'regression')
# ## tuning sPLS
# # set ncomp: repeated CV tuning of component count
# perf.spls.liver <- perf(spls.liver, validation = 'Mfold',
#                         folds = 10, nrepeat = 3) 
# plot(perf.spls.liver, criterion = 'Q2.total')
# # set keepX--the number of variables: 
# # set range of test values for number of variables to use from X dataframe
# list.keepX <- c(seq(2, 18, 2))
# # set range of test values for number of variables to use from Y dataframe
# list.keepY <- c(2:13) 

spls.result <- spls(X, Y, keepX = c(5000, 5000), keepY = c(4, 4),
                    mode="canonical", max.iter=500)  # run the method 
plotIndiv(spls.result) # plot the samples
plotVar(spls.result)   # plot the variables
selectVar(spls.result, comp = 1)$X$name 
# depict weight assigned to each of these variables
plotLoadings(spls.result, method = 'mean', contrib = 'max') 

plotVar(spls.result, cex = c(3,4), var.names = c(TRUE, TRUE))
# network or CIM
color.edge <- color.GreenRed(50)  # set the colours of the connecting lines
cim(spls.result, comp = 1:2, xlab = "metadata", ylab = "Genus",
    row.cex = 0.02, cutoff = 0.6, save = "pdf",
    name.save = "spls_pairwise_mT_metadata_heatmap"
)

# X11() # To open a new window for Rstudio
network(spls.result, comp = 1:2,
        cutoff = 0.5, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        #save = 'png', # save as a png to the current working directory
        name.save = 'sPLS Liver Toxicity Case Study Network Plot')
############ the end of sPLS for mM metabolome and metadata ##############
################################################################################
