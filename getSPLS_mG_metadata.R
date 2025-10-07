################################################################################
##################   sPLS for mG and metadata ## ##################
## prepare  mG taxon abundance data
xdf=read.delim("mG19_genus_abundance.txt",sep = '\t')
xdf=xdf[c(-1:-2, (-nrow(xdf)+1),-nrow(xdf)),]
colnames(xdf)=gsub("_([679])\\.","_0\\1\\.", colnames(xdf))
colnames(xdf)=gsub("mG19_","", colnames(xdf))

df=xdf[,-1]
rownames(df)=xdf[,1]
df=df[,order(colnames(df))]
df = as.data.frame(t(df))
colnames(df) <- sapply(strsplit(colnames(df), ";"), function(x) tail(x, 1))

## prepare metadata
metadata = read.csv('water_chemistry.csv', header = T, row.names = 1, check.names = F)
dim(metadata) # confirm the dimension of data
rownames(metadata) = gsub('mG19_', '', rownames(metadata))




################## the end ##################
################################################################################
