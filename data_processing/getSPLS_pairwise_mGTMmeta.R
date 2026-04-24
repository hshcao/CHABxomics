################################################################################
##################   sPLS for mG taxa and metadata ## ##################
## prepare  mG taxon abundance data
rm(list=ls())
tdf=read.delim("mG_genus_abundance.txt",sep = '\t', row.names = 1)
tdf = tdf[-16:-18,]

## prepare metadata
metadata = read.delim('water_chemistry.csv', header = T, row.names = 1,
                      sep = ',', check.names = F)
dim(metadata) # confirm the dimension of data
rownames(metadata) = gsub('mG19_', 'mG', rownames(metadata))
metadata = metadata[,c(-13:-14)]

rownames(tdf)==rownames(metadata)

## load library
library(mixOmics) # import the mixOmics library
set.seed(5249) # for reproducibility, remove for normal use
library(corrplot)
library(viridis)

X = tdf
Y = metadata
dim(X)
dim(Y)
cor(X$Microcystis, Y$PC)

cor_mat = cor(X)
diag(cor_mat) <- NA
# Find indices of correlations > 0.5
high_cor_idx <- which(abs(cor_mat) > 0.5, arr.ind = TRUE)

# Extract value and variable names
high_cor_values <- data.frame(
  Var1 = rownames(cor_mat)[high_cor_idx[, 1]],
  Var2 = colnames(cor_mat)[high_cor_idx[, 2]],
  Correlation = cor_mat[high_cor_idx]
)

# Optional: Remove duplicate pairs (i.e., Var1 vs Var2 is the same as Var2 vs Var1)
high_cor_values <- high_cor_values[high_cor_values$Var1 < high_cor_values$Var2, ]

highlighted_vars <- c("Limnobacter","Variovorax","Bosea","Pseudoxanthomonas",
                      "Rhizobium","Xanthomonadaceae","Brevundimonas", "Rhizobium", 
                      "Flavobacterium","Comamonadaceae","Sphingopyxis",
                      "Hyphomonas", "Agrobacterium","Pseudomonas","Flavihumibacter",
                      "Microcystis", "Nostoc","Anabaena", "Dolichospermum")
subset_high_cor <- high_cor_values[
  high_cor_values$Var1 %in% highlighted_vars | high_cor_values$Var2 %in% highlighted_vars,
]

library(reshape2)
library(ggplot2)

# Create square correlation matrix from pairwise data
mat <- reshape2::acast(subset_high_cor, Var1 ~ Var2, value.var = "Correlation")

# Convert to long format for ggplot
mat_long <- reshape2::melt(mat, na.rm = TRUE)
dev.off()
# Plot heatmap
ggplot(mat_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0,
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "High Correlations (Subset)", x = "", y = "")


spls.result <- spls(X, Y, keepX = c(350, 130), keepY = c(6, 3), max.iter=500)  # run the method 
plotIndiv(spls.result) # plot the samples
plotVar(spls.result)   # plot the variables
selectVar(spls.result, comp = 1)$X$name 
# depict weight assigned to each of these variables
plotLoadings(spls.result, method = 'mean', contrib = 'max') 

plotVar(spls.result, cex = c(2,4), var.names = c(TRUE, TRUE), )
# network or CIM
color.edge <- color.GreenRed(30)  # set the colours of the connecting lines
# X11() # To open a new window for Rstudio
dev.off()
cim(spls.result, comp = 1:2, xlab = "Metadatar", ylab = "Genus", cutoff = 0.5,
    row.cex = 0.2)
selected_vars <- selectVar(spls.result, comp = 1)$X$name
# Extract from original matrix
top_features <- X[, selected_vars]
labels <- colnames(top_features)
highlighted_vars <- c("Limnobacter","Variovorax","Bosea","Pseudoxanthomonas",
            "Rhizobium","Xanthomonadaceae","Brevundimonas", "Rhizobium", 
            "Flavobacterium","Comamonadaceae","Sphingopyxis",
            "Hyphomonas", "Agrobacterium","Pseudomonas","Flavihumibacter",
            "Microcystis", "Nostoc","Anabaena", "Dolichospermum")
label_colors <- ifelse(labels %in% highlighted_vars, "red", "black")

# Compute pairwise correlation matrix
cor_mat <- cor(top_features)

col=magma(50, alpha = 1, begin = 0, end = 1, direction = 1)
corrplot(cor_mat, method = "color", type = "upper", tl.cex=0.35, 
         tl.col = label_colors, col=col)

############ the end of sPLS for mG taxa and metadata ##################
################################################################################
###########   sPLS for mG genes and metadata ####################
## prepare  mG gene abundance data
rm(list=ls())
xdf=read.delim("Table5_mGmT_geneAb_across-TimeSite.tsv",sep = '\t')
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
    row.cex = 0.003, #save = "pdf", name.save = "spls_pairwise_mGgene_metadata_heatmap"
)

network(spls.result, comp = 1:2,
        cutoff = 0.5, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        #save = 'png', # save as a png to the current working directory
        name.save = 'sPLS Liver Toxicity Case Study Network Plot')
cim(spls.result, comp = 1:2, xlab = "clinic", ylab = "genes")


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

