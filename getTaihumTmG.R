################################################################################
#################### Question #1 ####################
########## relationships between the abundance measured by mG, mT, Chla, and PC

########## #1.1 Analysis of meta data only ##########
### read in meta data
rm(list=ls())
library(pheatmap)
meta=read.csv("water_chemistry.csv",header = T)
rownames(meta)=gsub("\\.","-",meta[,1])
meta_mg=meta_mt=meta[,-1]                  
meta_mg=meta_mg[!grepl("6-15-3|11-3-2|12-4-1",rownames(meta_mg)),]
rownames(meta_mt)=gsub("mG","mT",rownames(meta_mt))
meta_cm=cor(meta[,-1])
# draw plots abundance over time
pt<-plot(1:tNsamples, pop_abd[,1], type="b", lwd=1.5, las=2, lty=1, col=palette[1], ylim = c(0,20),
         pch=1, cex.axis = 1.1, cex.lab=1.1,ylab="Abundance of four species",xaxt="n", xlab ="",main = sites[i])
axis(1,at = 1:tNsamples,labels = rownames(pop_abd), las=2)
points(1:tNsamples, pop_abd[,2],col=palette[2], pch=18,ylim = c(0,20))
lines(1:tNsamples, pop_abd[,2], lwd=1.5, las=2, lty=1, col=palette[2],
      cex.axis = 1.1,cex.lab=1.1,ylab="",xaxt="n", xlab = "",ylim = c(0,20))
points(1:tNsamples, pop_abd[,3],col=palette[3], pch=2,ylim =  c(0,20))
lines(1:tNsamples, pop_abd[,3], lwd=1.5, las=2, lty=1, col=palette[3],
      cex.axis = 1.1,cex.lab=1.1,ylab="",xaxt="n", xlab = "",ylim =  c(0,20))  
points(1:tNsamples, pop_abd[,4],col=palette[4], pch=16,ylim =  c(0,20))
lines(1:tNsamples, pop_abd[,4], lwd=1.5, las=2, lty=1, col=palette[4],
      cex.axis = 1.1,cex.lab=1.1,ylab="",xaxt="n", xlab = "",ylim =  c(0,20))
points(1:tNsamples, pop_abd[,5],col=palette[5], pch=18,ylim =  c(0,20))
lines(1:tNsamples, pop_abd[,5], lwd=1.5, las=2, lty=1, col=palette[5],
      cex.axis = 1.1,cex.lab=1.1,ylab="",xaxt="n", xlab = "",ylim =  c(0,20))
par(new = T)
pop_abd=as.data.frame(pop_abd)
with(pop_abd, plot(Date, Cylindrospermopsis, pch=1, type = "l", ylim = c(0,0.4),
                   axes=F, xlab=NA, ylab=NA, cex=1.2, col = palette[2]))
axis(side = 4)
mtext(side = 4, line = 2, "Cylindrospermopsis abundance", col = palette[2])
legend("topright", legend = species, col = palette, pch = c(1,0,2,16,18),
       lty = 1, cex = 1.2, merge = TRUE)

#pdf("meta_meta_correlation.pdf",5,5)
pheatmap(as.matrix(meta_cm))
#dev.off()
#pdf("Chl-PC_correlation_meta_data.pdf",8,8)
par(mfrow=c(2,2))
#non-transformed data
plot(Chl ~ PC, data = meta,main="raw Chl~PC")
lines(lowess(meta[,5], meta[,4]), col = "blue")
abline(lm(Chl~PC, data = meta),col="red")
cpc=cor(meta[,4],meta[,5],method = "pearson")
legend("bottomright",legend=paste("R = ",round(cpc,3),sep = ""),bty = "o")
# Chl/PC both LN transformed data
plot(log(Chl) ~ log(PC), data = meta,main="lnChl~lnPC")
lines(lowess(log(meta[,5]), log(meta[,4])), col = "blue")
abline(lm(log(Chl) ~ log(PC), data = meta),col="red")
cpc=cor(log(meta[,5]), log(meta[,4]),method = "pearson")
legend("bottomright",legend=paste("R = ",round(cpc,3),sep = ""),
       bty = "o",col = "red")
#dev.off()

#pdf("Chl_PC-MC_correlation_meta_data.pdf",8,8)
par(mfrow=c(2,2))
#non-transformed data
plot(log(PC) ~ log(MC), data = meta,main="lnPC~lnPC")
lines(lowess(log(meta[,13]), log(meta[,5])), col = "blue")
abline(lm(log(PC) ~ log(MC), data = meta),col="red")
cpc=cor(log(meta[,13]), log(meta[,5]),method = "pearson")
legend("bottomright",legend=paste("R = ",round(cpc,3),sep = ""),
       bty = "o",col = "red")
# Chl/PC both LN transformed data
plot(log(Chl) ~ log(MC), data = meta,main="lnChl~lnMC")
lines(lowess(log(meta[,13]), log(meta[,4])), col = "blue")
abline(lm(log(Chl) ~ log(MC), data = meta),col="red")
cpc=cor(log(meta[,13]), log(meta[,4]),method = "pearson")
legend("bottomright",legend=paste("R = ",round(cpc,3),sep = ""),
       bty = "o",col = "red")
#dev.off()
########## End of #1.1 Analysis of meta data only ##########

########## # 1.2 Correlation in abundance with b/w mG and mT metrics ##########
rm(list=ls())
##### read in files for the number of clean reads in each file of MG and mT
## mG reads file
reads_mg=read.delim("mG19_files_readCount.tsv",sep = "\t",header = F)
colnames(reads_mg)=c("sample","readCount")
reads_mg[,2]=reads_mg[,2]/4
rownames(reads_mg)=reads_mg[,1]
# three samples are removed due to low numbers of reads
reads_mg=reads_mg[!grepl("6-15-3|11-3-2|12-4-1",reads_mg[,1]),]
## mT reads file
reads_mt=read.delim("mT19_files_readCount.tsv",sep = "\t",header = F)
colnames(reads_mt)=c("sample","readCount")
reads_mt[,2]=reads_mt[,2]/4
rownames(reads_mt)=reads_mt[,1]
#################### Microsystis aeruginosa ####################
### read in gene abundance from mG reads
geneAbd_mg=read.delim("Mar_mG19_allSample_geneAbundance_readCount.txt",header = T,sep = "\t")
colnames(geneAbd_mg)=gsub("\\.","-",colnames(geneAbd_mg))
geneAbd_mg=geneAbd_mg[,!grepl("6-15-3|11-3-2|12-4-1",colnames(geneAbd_mg))]
TotReads_persample_mg=colSums(geneAbd_mg[,2:ncol(geneAbd_mg)])
names(TotReads_persample_mg)=colnames(geneAbd_mg)[2:ncol(geneAbd_mg)]
# relative 1: total readCount/total reads in the sample
Rel_TotReads_persample_mg=TotReads_persample_mg/reads_mg[,2]

## read in gene abundance from mT reads
geneAbd_mt=read.delim("Mar_mT19_allSamples_geneAbundance_readCount.txt",header = T,sep = "\t")
colnames(geneAbd_mt)=gsub("\\.","-",colnames(geneAbd_mt))
TotReads_persample_mt=colSums(geneAbd_mt[,2:ncol(geneAbd_mt)])
names(TotReads_persample_mt)=colnames(geneAbd_mt)[2:ncol(geneAbd_mt)]
# relative 1: total readCount/total reads in the sample
Rel_TotReads_persample_mt=TotReads_persample_mt/reads_mt[,2]
Rel_popAbd_mt_mar=Rel_TotReads_persample_mt
names(Rel_popAbd_mt_mar)=gsub("_[679]-","_0\\1-",names(Rel_popAbd_mt_mar))

##### correlation between mG, mT, and meta data
#pdf("MAR_mGmT_correlation_total_relative_readCount.pdf",6,6)
### correlation using total readCount
## total mG vs total mT in replicate 1
toReads_mt_mg=TotReads_persample_mt
names(toReads_mt_mg)=gsub("mT","mG",names(toReads_mt_mg))
toReads_mt_mg_1=toReads_mt_mg[grepl("-1$",names(toReads_mt_mg),perl = T)]
names(toReads_mt_mg_1)=gsub("-1$","",names(toReads_mt_mg_1),perl = T)
toReads_mt_mg_1=toReads_mt_mg_1[match(names(TotReads_persample_mg),names(toReads_mt_mg_1))]
par(mfrow=c(2,2))
plot(toReads_mt_mg_1 ~ TotReads_persample_mg, xlab="Total mG",ylab = "Total mT (Rep 1)")
abline(lm(toReads_mt_mg_1 ~ TotReads_persample_mg),col="red")
cpc=cor(toReads_mt_mg_1, TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")
### total mG vs total mT in replicate 2
toReads_mt_mg_2=toReads_mt_mg[grepl("-2$",names(toReads_mt_mg),perl = T)]
names(toReads_mt_mg_2)=gsub("-2$","",names(toReads_mt_mg_2),perl = T)
toReads_mt_mg_2=toReads_mt_mg_2[match(names(TotReads_persample_mg),names(toReads_mt_mg_2))]
plot(toReads_mt_mg_2 ~ TotReads_persample_mg, xlab="Total mG",ylab = "Total mT (Rep 2)")
abline(lm(toReads_mt_mg_2 ~ TotReads_persample_mg),col="red")
cpc=cor(toReads_mt_mg_2, TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")

##### correlation using relative readCount
### relative mG vs relative mT in replicate 1
relReads_mt_mg=Rel_TotReads_persample_mt
names(relReads_mt_mg)=gsub("mT","mG",names(relReads_mt_mg))
relReads_mt_mg_1=relReads_mt_mg[grepl("-1$",names(relReads_mt_mg),perl = T)]
names(relReads_mt_mg_1)=gsub("-1$","",names(relReads_mt_mg_1),perl = T)
relReads_mt_mg_1=relReads_mt_mg_1[match(names(Rel_TotReads_persample_mg),names(relReads_mt_mg_1))]
plot(relReads_mt_mg_1 ~ Rel_TotReads_persample_mg,xlab="Rel mGreads",ylab = "relative_mT_rep1")
abline(lm(relReads_mt_mg_1 ~ Rel_TotReads_persample_mg),col="red")
cpc=cor(relReads_mt_mg_1, Rel_TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")

### relative mG vs relative mT in replicate 2
relReads_mt_mg_2=relReads_mt_mg[grepl("-2$",names(relReads_mt_mg),perl = T)]
names(relReads_mt_mg_2)=gsub("-2$","",names(relReads_mt_mg_2),perl = T)
relReads_mt_mg_2=relReads_mt_mg_2[match(names(Rel_TotReads_persample_mg),names(relReads_mt_mg_2))]
plot(relReads_mt_mg_2 ~ Rel_TotReads_persample_mg,xlab="Rel mGreads",ylab = "relative_mT_rep2")
abline(lm(relReads_mt_mg_2 ~ Rel_TotReads_persample_mg),col="red")
cpc=cor(relReads_mt_mg_2, Rel_TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")
###
#dev.off()
####################
### combined metadata mG and mT reads
names(Rel_TotReads_persample_mg)==names(toReads_mt_mg_2)
tdf=cbind(TotReads_persample_mg,Rel_TotReads_persample_mg, toReads_mt_mg_1,
          toReads_mt_mg_2,relReads_mt_mg_1,relReads_mt_mg_2)
colnames(tdf)=c("MARmGtotal","MARmGrelative","MARmTtotal_rep1","MARmTtotal_rep2","MARmTrelative_rep1","MARmTrelative_rep2")
tdf <- data.frame(names = row.names(tdf), tdf)
meta[,1]=rownames(meta)
meta[,1]=gsub("_0([679])-","_\\1-",meta[,1])
meta_mTG2=merge(meta_mTG,tdf,by.x = "X",by.y = "names", all.x = T)
meta_mTG2_cm=cor(meta_mTG2[,-1])
#save(meta_mTG2,file = "metadata_MAR-Acir-mGT_combined.txt")
library(pheatmap)
pdf("MAR_Dol-Acir_meta_corr_heatmap.pdf",8,8)
pheatmap(as.matrix(meta_mTG2_cm))
dev.off()

#################### Dolichosperumum sp ACIR310F ##############
### read in gene abundance from mG reads
geneAbd_mg=read.delim("Acir_mG19_allSamples_geneAbundance_readCount.txt",header = T,sep = "\t")
colnames(geneAbd_mg)=gsub("\\.","-",colnames(geneAbd_mg))
colnames(geneAbd_mg)=gsub("Acir-","",colnames(geneAbd_mg))
geneAbd_mg=geneAbd_mg[,!grepl("6-15-3|11-3-2|7-23-2|11-3-1|12-4-1",colnames(geneAbd_mg))]
TotReads_persample_mg=colSums(geneAbd_mg[,2:ncol(geneAbd_mg)])
names(TotReads_persample_mg)=colnames(geneAbd_mg)[2:ncol(geneAbd_mg)]
# relative 1: total readCount/total reads in the sample
total_reads_mg_sub=reads_mg[c(-1,-8),2]
Rel_TotReads_persample_mg=TotReads_persample_mg/total_reads_mg_sub

## read in gene abundance from mT reads
geneAbd_mt=read.delim("Acir_mT19_allSamples_geneAbundance_readCount.txt",header = T,sep = "\t")
colnames(geneAbd_mt)=gsub("\\.","-",colnames(geneAbd_mt))
colnames(geneAbd_mt)=gsub("Acir-","",colnames(geneAbd_mt))
TotReads_persample_mt=colSums(geneAbd_mt[,2:ncol(geneAbd_mt)])
names(TotReads_persample_mt)=colnames(geneAbd_mt)[2:ncol(geneAbd_mt)]
# relative 1: total readCount/total reads in the sample
Rel_TotReads_persample_mt=TotReads_persample_mt/reads_mt[,2]
Rel_popAbd_mt_Acir=Rel_TotReads_persample_mt
names(Rel_popAbd_mt_Acir)=gsub("_[679]-","_0\\1-",names(Rel_popAbd_mt_Acir))

##### correlation between mG, mT, and meta data
#pdf("Acir_mGmT_correlation_total_relative_tpm.pdf",6,6)
### correlation using total readCount
## total mG vs total mT in replicate 1
toReads_mt_mg=TotReads_persample_mt
names(toReads_mt_mg)=gsub("mT","mG",names(toReads_mt_mg))
#names(toReads_mt_mg)=gsub("Acir-mT","mG",names(toReads_mt_mg))
toReads_mt_mg_1=toReads_mt_mg[grepl("-1$",names(toReads_mt_mg),perl = T)]
names(toReads_mt_mg_1)=gsub("-1$","",names(toReads_mt_mg_1),perl = T)
toReads_mt_mg_1=toReads_mt_mg_1[match(names(TotReads_persample_mg),names(toReads_mt_mg_1))]
par(mfrow=c(2,2))
plot(toReads_mt_mg_1 ~ TotReads_persample_mg, xlab="Total mG",ylab = "Total mT (Rep 1)")
abline(lm(toReads_mt_mg_1 ~ TotReads_persample_mg),col="red")
cpc=cor(toReads_mt_mg_1, TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")
### total mG vs total mT in replicate 2
toReads_mt_mg_2=toReads_mt_mg[grepl("-2$",names(toReads_mt_mg),perl = T)]
names(toReads_mt_mg_2)=gsub("-2$","",names(toReads_mt_mg_2),perl = T)
toReads_mt_mg_2=toReads_mt_mg_2[match(names(TotReads_persample_mg),names(toReads_mt_mg_2))]
plot(toReads_mt_mg_2 ~ TotReads_persample_mg, xlab="Total mG",ylab = "Total mT (Rep 2)")
abline(lm(toReads_mt_mg_2 ~ TotReads_persample_mg),col="red")
cpc=cor(toReads_mt_mg_2, TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")

##### correlation using relative readCount
### relative mG vs relative mT in replicate 1
relReads_mt_mg=Rel_TotReads_persample_mt
names(relReads_mt_mg)=gsub("mT","mG",names(relReads_mt_mg))
relReads_mt_mg_1=relReads_mt_mg[grepl("-1$",names(relReads_mt_mg),perl = T)]
names(relReads_mt_mg_1)=gsub("-1$","",names(relReads_mt_mg_1),perl = T)
relReads_mt_mg_1=relReads_mt_mg_1[match(names(Rel_TotReads_persample_mg),names(relReads_mt_mg_1))]
plot(relReads_mt_mg_1 ~ Rel_TotReads_persample_mg,xlab="Rel mGreads",ylab = "relative_mT_rep1")
abline(lm(relReads_mt_mg_1 ~ Rel_TotReads_persample_mg),col="red")
cpc=cor(relReads_mt_mg_1, Rel_TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")

### relative mG vs relative mT in replicate 2
relReads_mt_mg_2=relReads_mt_mg[grepl("-2$",names(relReads_mt_mg),perl = T)]
names(relReads_mt_mg_2)=gsub("-2$","",names(relReads_mt_mg_2),perl = T)
relReads_mt_mg_2=relReads_mt_mg_2[match(names(Rel_TotReads_persample_mg),names(relReads_mt_mg_2))]
plot(relReads_mt_mg_2 ~ Rel_TotReads_persample_mg,xlab="Rel mGreads",ylab = "relative_mT_rep2")
abline(lm(relReads_mt_mg_2 ~ Rel_TotReads_persample_mg),col="red")
cpc=cor(relReads_mt_mg_2, Rel_TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")
###
#dev.off()

### combined metadata mG and mT reads
names(Rel_TotReads_persample_mg)==names(toReads_mt_mg_2)
tdf=cbind(TotReads_persample_mg,Rel_TotReads_persample_mg, toReads_mt_mg_1,
          toReads_mt_mg_2,relReads_mt_mg_1,relReads_mt_mg_2)
colnames(tdf)=c("ACmGtotal","ACmGrelative","ACmTtotal_rep1","ACmTtotal_rep2","ACmTrelative_rep1","ACmTrelative_rep2")
tdf <- data.frame(names = row.names(tdf), tdf)
meta[,1]=rownames(meta)
meta[,1]=gsub("_0([679])-","_\\1-",meta[,1])
meta_mTG=merge(meta,tdf,by.x = "X", by.y = "names", all.y = T)
meta_mTG_cm=cor(meta_mTG[,-1])

library(pheatmap)
pdf("meta_ACmTG_cor.pdf",5,5)
pheatmap(as.matrix(meta_mTG_cm))
dev.off()
##############################
#################### Dolichosperumum sp 1059 ##############
### read in gene abundance from mG reads
geneAbd_mg=read.delim("Doli_mG19_allSamples_geneAbundance_readCount.txt",header = T,sep = "\t")
colnames(geneAbd_mg)=gsub("\\.","-",colnames(geneAbd_mg))
colnames(geneAbd_mg)=gsub("Doli-","",colnames(geneAbd_mg))
geneAbd_mg=geneAbd_mg[,!grepl("6-15-3|11-3-2|12-4-1|7-23-2",colnames(geneAbd_mg))]
TotReads_persample_mg=colSums(geneAbd_mg[,2:ncol(geneAbd_mg)])
names(TotReads_persample_mg)=colnames(geneAbd_mg)[2:ncol(geneAbd_mg)]
# relative 1: total readCount/total reads in the sample
Rel_TotReads_persample_mg=TotReads_persample_mg/reads_mg[-8,2]

## read in gene abundance from mT reads
geneAbd_mt=read.delim("Doli_mT19_allSamples_geneAbundance_readCount.txt",header = T,sep = "\t")
geneAbd_mt=geneAbd_mt[,-2]
colnames(geneAbd_mt)=gsub("\\.","-",colnames(geneAbd_mt))
colnames(geneAbd_mt)=gsub("Doli-","",colnames(geneAbd_mt))
TotReads_persample_mt=colSums(geneAbd_mt[,2:ncol(geneAbd_mt)])
names(TotReads_persample_mt)=colnames(geneAbd_mt)[2:ncol(geneAbd_mt)]
# relative 1: total readCount/total reads in the sample
Rel_TotReads_persample_mt=TotReads_persample_mt/reads_mt[,2]
Rel_popAbd_mt_Doli=Rel_TotReads_persample_mt
names(Rel_popAbd_mt_Doli)=gsub("_[679]-","_0\\1-",names(Rel_popAbd_mt_Doli))
##### correlation between mG, mT, and meta data
#pdf("Doli_mGmT_correlation_total_relative_readCount.pdf",6,6)
### correlation using total readCount
## total mG vs total mT in replicate 1
toReads_mt_mg=TotReads_persample_mt
names(toReads_mt_mg)=gsub("mT","mG",names(toReads_mt_mg))
toReads_mt_mg_1=toReads_mt_mg[grepl("-1$",names(toReads_mt_mg),perl = T)]
names(toReads_mt_mg_1)=gsub("-1$","",names(toReads_mt_mg_1),perl = T)
toReads_mt_mg_1=toReads_mt_mg_1[match(names(TotReads_persample_mg),names(toReads_mt_mg_1))]
par(mfrow=c(2,2))
plot(toReads_mt_mg_1 ~ TotReads_persample_mg, xlab="Total mG",ylab = "Total mT (Rep 1)")
abline(lm(toReads_mt_mg_1 ~ TotReads_persample_mg),col="red")
cpc=cor(toReads_mt_mg_1, TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")
### total mG vs total mT in replicate 2
toReads_mt_mg_2=toReads_mt_mg[grepl("-2$",names(toReads_mt_mg),perl = T)]
names(toReads_mt_mg_2)=gsub("-2$","",names(toReads_mt_mg_2),perl = T)
toReads_mt_mg_2=toReads_mt_mg_2[match(names(TotReads_persample_mg),names(toReads_mt_mg_2))]
plot(toReads_mt_mg_2 ~ TotReads_persample_mg, xlab="Total mG",ylab = "Total mT (Rep 2)")
abline(lm(toReads_mt_mg_2 ~ TotReads_persample_mg),col="red")
cpc=cor(toReads_mt_mg_2, TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")

##### correlation using relative readCount
### relative mG vs relative mT in replicate 1
relReads_mt_mg=Rel_TotReads_persample_mt
names(relReads_mt_mg)=gsub("mT","mG",names(relReads_mt_mg))
relReads_mt_mg_1=relReads_mt_mg[grepl("-1$",names(relReads_mt_mg),perl = T)]
names(relReads_mt_mg_1)=gsub("-1$","",names(relReads_mt_mg_1),perl = T)
relReads_mt_mg_1=relReads_mt_mg_1[match(names(Rel_TotReads_persample_mg),names(relReads_mt_mg_1))]
plot(relReads_mt_mg_1 ~ Rel_TotReads_persample_mg,xlab="Rel mGreads",ylab = "relative_mT_rep1")
abline(lm(relReads_mt_mg_1 ~ Rel_TotReads_persample_mg),col="red")
cpc=cor(relReads_mt_mg_1, Rel_TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")

### relative mG vs relative mT in replicate 2
relReads_mt_mg_2=relReads_mt_mg[grepl("-2$",names(relReads_mt_mg),perl = T)]
names(relReads_mt_mg_2)=gsub("-2$","",names(relReads_mt_mg_2),perl = T)
relReads_mt_mg_2=relReads_mt_mg_2[match(names(Rel_TotReads_persample_mg),names(relReads_mt_mg_2))]
plot(relReads_mt_mg_2 ~ Rel_TotReads_persample_mg,xlab="Rel mGreads",ylab = "relative_mT_rep2")
abline(lm(relReads_mt_mg_2 ~ Rel_TotReads_persample_mg),col="red")
cpc=cor(relReads_mt_mg_2, Rel_TotReads_persample_mg,method = "pearson")
legend("bottomright",legend=paste("r = ",round(cpc,3),sep = ""),bty = "o")
###
#dev.off()

### combined metadata mG and mT reads
names(Rel_TotReads_persample_mg)==names(toReads_mt_mg_2)
tdf=cbind(TotReads_persample_mg,Rel_TotReads_persample_mg, toReads_mt_mg_1,
          toReads_mt_mg_2,relReads_mt_mg_1,relReads_mt_mg_2)
colnames(tdf)=c("DOmGtotal","DOmGrelative","DOmTtotal_rep1","DOmTtotal_rep2","DOmTrelative_rep1","DOmTrelative_rep2")
tdf <- data.frame(names = row.names(tdf), tdf)
meta_mTG3=merge(meta_mTG2,tdf,by.x = "X",by.y = "names", all.x = T)
#save(meta_mTG3,file = "metadata_MAR-Acir-Doli_mGT_combined.txt")
meta_mTG3_cm=cor(meta_mTG3[,-1])

library(pheatmap)
#pdf("MAR_Acir_Doli_meta_mTG_cor.pdf",9,9)
pheatmap(as.matrix(meta_mTG3_cm))
#dev.off()
#################### END// Question #1 // ####################


################################################################################
#################### Question #2 temporal population dynamics ####################
#################################################################################
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
goafiles=file.path("GOAs",list.files("GOAs",pattern = "_GOA"))
##### 2. get query pathway files of all five genomes
#pathwayfiles=file.path("pathfiles",list.files(path = "pathfiles",pattern = "pid_locus_pathway_sort.txt"))
###### 3. get the folder with gene abundance and metadata file

##### 4. loop through each site and the two years for WGCNA # not possible for adjust power
### 4.1 Load water chemistry metadata and combined metadata with mG
meta=read.csv("water_chemistry.csv",header = T)
rownames(meta)=gsub("\\.","-",meta[,1])
rownames(meta)=gsub("mG", "mT", rownames(meta))
meta_mt=meta[,-1]
meta_mt=meta_mt[order(rownames(meta_mt)),]
meta_mt$dates=rownames(meta_mt)
load("metadata_MAR-Acir-Doli_mGT_combined.txt")
write.table(meta_mTG3,"metadata_Mar_Dolilichosperum_Acir310F-1059_mGT_combined.txt",
            quote = F, sep = "\t")
meta_mTG3[,1]=gsub("mG","mT",meta_mTG3[,1])
meta_mTG3[,1]=gsub("_([679])-","_0\\1-",meta_mTG3[,1])

## datTraits file
datTraits=merge(meta_mt,meta_mTG3[,grepl("X|mGrel", colnames(meta_mTG3))], 
                 by.x = "dates", by.y = "X", all.x = TRUE)
rownames(datTraits)=datTraits[,1]
datTraits=datTraits[,-1]
### 4.2 WGCNA analysis for five genome gene expression files at each site/year
# pattern may differ
genefiles=file.path(".",list.files(pattern="mT19_allSamples_geneAbundance_tpm"))
#for (g in 1:5){
g=3
#detach("package:topGO", unload = TRUE)
#detach("package:WGCNA", unload = TRUE)
#library(WGCNA)
# get gene expression profile for each genome
## read in gene abundance from mT reads
geneAbd_mt=read.delim("Mar_mT19_allSamples_geneAbundance_tpm.txt",header = T,sep = "\t", row.names = 1)
colnames(geneAbd_mt)=gsub("\\.","-",colnames(geneAbd_mt))
colnames(geneAbd_mt)=gsub("_([679])-","_0\\1-",colnames(geneAbd_mt))
dates=sort(unique(gsub("-[12]$","", colnames(geneAbd_mt))))

mat=as.data.frame(matrix(0,nrow=length(dates),ncol=nrow(geneAbd_mt)))
rownames(mat)=dates
colnames(mat)=rownames(geneAbd_mt)
for (d in 1:length(dates)) {
  #d=1
  tempi=as.matrix(geneAbd_mt[,grepl(rownames(mat)[d],colnames(geneAbd_mt))])
  tempii=apply(tempi,1,mean)
  mat[d,]=as.numeric(tempii)
}
otu=mat[,colSums(mat>0)>7]
colnames(otu)=gsub("_\\D+\\w+$","",colnames(otu),perl = T)
library(heatmap3)
pdf("Mar_gene_expression_time-series.pdf",8,10)
heatmap3(t(otu), Colv = NA,labRow = "",margins = c(10, 5))
dev.off()
## WGCNA analysis starts
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

## check for outliers
sampleTree = hclust(dist(otu), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,5,2,0))
#pdf("MAR_wgcna_all-sample_cluster.pdf", width = 8, height = 8);
pt<-plot(sampleTree, main = "MAR", sub="", xlab="Sample", 
         cex.lab = 1.8, cex.axis = 1.5, cex.main = 1.2, cex=1.5, ylim=c(0,0.5)
)
#dev.off()

# Plot a line to show the cut
abline(h = 0.02, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 0.02, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = otu[keepSamples, ] # datExpr = sqrt(otu[keepSamples, ])
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Before we continue with network construction and module detection, 
# we visualize how the meta data relate to the sample dendrogram
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf("Mar_wgcna_tpm_meta_combined.pdf",5,5)
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
                       TOMType = "unsigned", minModuleSize = 15,#unsigned
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       #saveTOMFileBase = paste(outbname,"TOM",sep = ""),
                       verbose = 3);
#write.table(table(net$colors),file = "Plastics_taxa-in-modules.txt",sep = "\t")
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
listname= "MAR" #paste(sites[4+f],species[g],sep = "_")
# add module-trait association to list
#list_mod_trait[[listname]]=moduleTraitCor
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
#pdf("Mar_wgcna_modules_metadata_cor_heatmap.pdf",10,10)
par(mar = c(5, 8, 0, 0)+0.1);
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50), #blueWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-metadata relationships"));
#dev.off()

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Temp); # ====>>> here to SELECT TRAITs for regression
names(weight) = "Temp" # ====>>> here to SELECT TRAITs for regression
# names (colors) of the modules to pick the module
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste( "p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# 3.b.2 Intramodular analysis: identifying genes with high GS and MM
modNames
module = modNames[6] ## ====>>>>> pick the MODULE here
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
#pdf(paste("Mar_wgcna",module,"module-gene-membership.pdf",sep="_"),5,5)
plott<-verboseScatterplot(geneModuleMembership[moduleGenes, column],
                          geneTraitSignificance[moduleGenes, 1],
                          xlab = paste("Module Membership in", module, "module"),
                          ylab = paste("Gene significance for",module, "module"),
                          main = paste("Module membership vs. gene significance\n"),
                          abline = T, abline.color = module, 
                          cex = 1.2, pch=16, cex.main = 0.9,
                          cex.lab = 1.2, cex.axis = 1.2,
                          col = module)
#dev.off()

## Export module to cytoscape
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
collectGarbage()
TOM = TOMsimilarityFromExpr(datExpr, power = power);
modNames = substring(names(MEs), 3)
for (m in 1:length(modNames)){
  modules = modNames[m]
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("Mar",modules,"cytos-edge.txt", sep="_"),
                                 nodeFile = paste("Mar",modules,"cytos-nodes.txt", sep="_"),
                                 weighted = TRUE,
                                 threshold = 0.8, # correlation cutoff
                                 nodeNames = modProbes,
                                 includeColNames = TRUE,
                                 #altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]
  )
}

## module functional enrichment and visualization
library(gplots)
library(ggplot2)
library(dplyr)
# Display the correlation values within a heatmap plot
pdf("Mar_wgcna_module_meta_cor2.pdf",10,10)
par(mar = c(4, 6, 1, 1)+0.5)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50), #blueWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
par(mar = c(4, 2, 1, 6)+0.5)
heatmap.2(moduleTraitCor,col= greenred(100), scale="none",trace="none", symkey = T,
          density.info="none", dendrogram=c("both"),  cexRow =1, cexCol = 1.2)
#dev.off()
## enrichment of GO terms for the modules
library(topGO)
#library(UpSetR)
# prepare files for GO enrichment
geneID2GO<-readMappings(file = goafiles[3])
allGenes <- unique(names(geneID2GO))
modNames = substring(names(MEs), 3)
goterms=c("BP","CC","MF")
for (md in 15:length(modNames)){
  #md=14
  modules = modNames[md]
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modGenes = unique(probes[inModule]);
  # GO term enrichment 
  modGenes_factor <- factor(as.integer(allGenes %in% modGenes))
  names(modGenes_factor)=allGenes
  for (gt in 1:length(goterms)){
    #gt=3
    outPDF=paste("Mar_wgcna",modules,goterms[gt],"topGO.pdf",sep = "_")
    outTXT=paste("Mar_wgcna",modules,goterms[gt],"topGO.txt",sep = "_")
    GOdata <- new("topGOdata", ontology = goterms[gt], allGenes = modGenes_factor,
                  nodeSize = 5, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    resultFisher
    topnodes=length(resultFisher@score)
    allRes <- GenTable(GOdata, classicFisher = resultFisher,
                       orderBy = "classicFisher", topNodes = topnodes)
    write.table(allRes, outTXT, quote = F, sep = "\t")
    par(mar = c(0, 0, 0, 0)+0.1)
    pdf(outPDF,8,8)
    print(showSigOfNodes(GOdata, score(resultFisher), 
                   firstSigNodes = 11, useInfo = 'def'))
    print(legend("topleft",legend = paste(outPDF,"ME", modules,goterms[gt],"GOterms",sep = "_"),bty = "n"))
    print(showSigOfNodes(GOdata, score(resultFisher), 
                         firstSigNodes = 11, useInfo = 'all'))
    print(legend("topleft",legend = paste(outPDF,"ME", modules,goterms[gt],"GOterms",sep = "_"),bty = "n"))
    dev.off()
  }
}
print("HAHA")


####################### END of MAR WGCNA and topGO analysis ####################
################################################################################
######## Acir wgcna and module function enrichment in GO terms  ################
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
goafiles=file.path("GOAs",list.files("GOAs",pattern = "_GOA"))
##### 2. get query pathway files of all five genomes
#pathwayfiles=file.path("pathfiles",list.files(path = "pathfiles",pattern = "pid_locus_pathway_sort.txt"))
###### 3. get the folder with gene abundance and metadata file

##### 4. loop through each site and the two years for WGCNA # not possible for adjust power
### 4.1 Load water chemistry metadata and combined metadata with mG
meta=read.csv("water_chemistry.csv",header = T)
rownames(meta)=gsub("\\.","-",meta[,1])
rownames(meta)=gsub("mG", "mT", rownames(meta))
meta_mt=meta[,-1]
meta_mt=meta_mt[order(rownames(meta_mt)),]
meta_mt$dates=rownames(meta_mt)
load("metadata_MAR-Acir-Doli_mGT_combined.txt")
meta_mTG3[,1]=gsub("mG","mT",meta_mTG3[,1])
meta_mTG3[,1]=gsub("_([679])-","_0\\1-",meta_mTG3[,1])

## datTraits file
datTraits=merge(meta_mt,meta_mTG3[,grepl("X|mGrel", colnames(meta_mTG3))], 
                by.x = "dates", by.y = "X", all.x = TRUE)
rownames(datTraits)=datTraits[,1]
datTraits=datTraits[,-1]
### 4.2 WGCNA analysis for five genome gene expression files at each site/year
# pattern may differ
genefiles=file.path(".",list.files(pattern="mT19_allSamples_geneAbundance_tpm"))

geneAbd_mt=read.delim("Acir_mT19_allSamples_geneAbundance_tpm.txt",
                      header = T,sep = "\t", row.names = 1)
colnames(geneAbd_mt)=gsub("\\.","-",colnames(geneAbd_mt))
colnames(geneAbd_mt)=gsub("_([679])-","_0\\1-",colnames(geneAbd_mt))
colnames(geneAbd_mt)=gsub("Acir-","",colnames(geneAbd_mt))
dates=sort(unique(gsub("-[12]$","", colnames(geneAbd_mt))))
mat=as.data.frame(matrix(0,nrow=length(dates),ncol=nrow(geneAbd_mt)))
rownames(mat)=dates
colnames(mat)=rownames(geneAbd_mt)
for (d in 1:length(dates)) {
  #d=1
  tempi=as.matrix(geneAbd_mt[,grepl(rownames(mat)[d],colnames(geneAbd_mt))])
  tempii=apply(tempi,1,mean)
  mat[d,]=as.numeric(tempii)
}
otu=mat[,colSums(mat>0)>3]
colnames(otu)=gsub("_ACIR310F_RS\\d+$","",colnames(otu),perl = T)
colnames(otu)=gsub("_[a-z]\\w+.+$","",colnames(otu),perl = T)
library(heatmap3)

pdf("Dolichospermum_Acir310F_gene_expression_time-series.pdf",8,10)
heatmap3(t(otu), Colv = NA,labRow = "")
dev.off()
##### WGCNA analysis starts
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

## check for outliers
sampleTree = hclust(dist(otu), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,5,2,0))
#pdf("Acir_wgcna_all-sample_cluster.pdf", width = 8, height = 8);
pt<-plot(sampleTree, main = "Acir", sub="", xlab="Sample", 
         cex.lab = 1.8, cex.axis = 1.5, cex.main = 1.2, cex=1.5, ylim=c(0,0.5)
)
#dev.off()

# Plot a line to show the cut
abline(h = 0.02, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 0.02, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = otu[keepSamples, ] # datExpr = sqrt(otu[keepSamples, ])
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Before we continue with network construction and module detection, 
# we visualize how the meta data relate to the sample dendrogram
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf("Acir_wgcna_tpm_meta_combined.pdf",5,5)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",
                    cex.colorLabels = 0.85, cex.dendroLabels = 1, 
                    cex.rowText = 0.8,las=2
)
#dev.off()
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
                       TOMType = "unsigned", minModuleSize = 15,#unsigned
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       #saveTOMFileBase = paste(outbname,"TOM",sep = ""),
                       verbose = 3);
#write.table(table(net$colors),file = "Plastics_taxa-in-modules.txt",sep = "\t")
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
listname= "Acir" #paste(sites[4+f],species[g],sep = "_")
# add module-trait association to list
#list_mod_trait[[listname]]=moduleTraitCor
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
#pdf("Acir_wgcna_modules_metadata_cor_heatmap.pdf",10,10)
par(mar = c(5, 8, 0, 0)+0.1);
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50), #blueWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.0,
               zlim = c(-1,1),
               main = paste("Module-metadata relationships"));
#dev.off()

## Gene relationship to trait and important modules: Gene Significance and Module Membership
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$Temp); # ====>>> here to SELECT TRAITs for regression
names(weight) = "Temp" # ====>>> here to SELECT TRAITs for regression
# names (colors) of the modules to pick the module
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste( "p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
# 3.b.2 Intramodular analysis: identifying genes with high GS and MM
modNames
module = modNames[4] ## ====>>>>> pick the MODULE here
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
#pdf(paste("Mar_wgcna",module,"module-gene-membership.pdf",sep="_"),5,5)
plott<-verboseScatterplot(geneModuleMembership[moduleGenes, column],
                          geneTraitSignificance[moduleGenes, 1],
                          xlab = paste("Module Membership in", module, "module"),
                          ylab = paste("Gene significance for",module, "module"),
                          main = paste("Module membership vs. gene significance\n"),
                          abline = T, abline.color = module, 
                          cex = 1.2, pch=16, cex.main = 0.9,
                          cex.lab = 1.2, cex.axis = 1.2,
                          col = module)
#dev.off()

## Export module to cytoscape
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
collectGarbage()
TOM = TOMsimilarityFromExpr(datExpr, power = power);
modNames = substring(names(MEs), 3)
for (m in 1:length(modNames)){
  modules = modNames[m]
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("Acir",modules,"cytos-edge.txt", sep="_"),
                                 nodeFile = paste("Acir",modules,"cytos-nodes.txt", sep="_"),
                                 weighted = TRUE,
                                 threshold = 0.8, # correlation cutoff
                                 nodeNames = modProbes,
                                 includeColNames = TRUE,
                                 #altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]
  )
}

## module functional enrichment and visualization
library(gplots)
library(ggplot2)
library(dplyr)
# Display the correlation values within a heatmap plot
#pdf("Acir_wgcna_module_meta_cor.pdf",8,8)
par(mar = c(4, 6, 1, 1)+0.5)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50), #blueWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
par(mar = c(4, 2, 1, 6)+0.5)
heatmap.2(moduleTraitCor,col= greenred(100), scale="none",trace="none", symkey = T,
          density.info="none", dendrogram=c("both"),  cexRow =1, cexCol = 1.2)
#dev.off()
## enrichment of GO terms for the modules
library(topGO)
#library(UpSetR)
# prepare files for GO enrichment
geneID2GO<-readMappings(file = goafiles[1])
allGenes <- unique(names(geneID2GO))
modNames = substring(names(MEs), 3)
goterms=c("BP","CC","MF")
for (md in 1:length(modNames)){
  #md=1
  modules = modNames[md]
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modGenes = unique(probes[inModule]);
  # GO term enrichment 
  modGenes_factor <- factor(as.integer(allGenes %in% modGenes))
  names(modGenes_factor)=allGenes
  for (gt in 1:length(goterms)){
    #gt=2
    outPDF=paste("Acir_wgcna",modules,goterms[gt],"topGO.pdf",sep = "_")
    outTXT=paste("Acir_wgcna",modules,goterms[gt],"topGO.txt",sep = "_")
    GOdata <- new("topGOdata", ontology = goterms[gt], allGenes = modGenes_factor,
                  nodeSize = 5, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    resultFisher
    topnodes=length(resultFisher@score)
    allRes <- GenTable(GOdata, classicFisher = resultFisher,
                       orderBy = "classicFisher", topNodes = topnodes)
    write.table(allRes, outTXT, quote = F, sep = "\t")
    par(mar = c(0, 0, 0, 0)+0.1)
    pdf(outPDF,8,8)
    showSigOfNodes(GOdata, score(resultFisher), 
                         firstSigNodes = topnodes, useInfo = 'def')
    print(legend("topleft",legend = paste(outPDF,"ME", modules,goterms[gt],"GOterms",sep = "_"),bty = "n"))
    print(showSigOfNodes(GOdata, score(resultFisher), 
                         firstSigNodes = topnodes, useInfo = 'all'))
    print(legend("topleft",legend = paste(outPDF,"ME", modules,goterms[gt],"GOterms",sep = "_"),bty = "n"))
    dev.off()
  }
}
print("HAHA")

#################################################################################
##### code to used for later #########
  ## Query pathway enrichment
  # read in query pathways
  for (p in 1:length(pathNames)){
    #p=1
    df_allpath_temp=genome_paths[which(genome_paths$shortID != "NoPath"), c(2,4)]
    df_1path_temp=genome_paths[which(genome_paths$shortID==pathNames[p]), c(2,4)]
    hyper_m=nrow(df_1path_temp)
    hyper_n=hyper_N-hyper_m
    hyper_k=length(which((df_allpath_temp$Locus %in% modGenes) == TRUE))
    hyper_x=length(which((df_1path_temp$Locus %in% modGenes) == TRUE))
    p.value <-  phyper(q=hyper_x -1, m=hyper_m, n=hyper_n, k=hyper_k, lower.tail=FALSE)
    df_pathEnrich[p,m]=p.value
  }
}

detach("package:topGO", unload = TRUE)
# plot enriched GO terms for all modules
par(mar = c(2, 2, 1, 1)+0.1)
upset(fromList(modtermList), nsets = length(modtermList),nintersects = NA,order.by = "freq",
      point.size = 3, matrix.color = "grey3", 
      mainbar.y.label = "Number of GO terms in intersections", sets.x.label = "Number of GO terms",
      text.scale = c(1.6, 2, 1.5, 1.5, 1.5, 1.8)
)
# plot enriched GO terms by modules
termTemp_all=unique(unlist(Reduce(c,modxterm_List)))
termTemp_df=data.frame("Terms"=termTemp_all)
for (i in 1:length(modxterm_List)) {
  temp_comp=unique(unlist(modxterm_List[i]))
  termTemp_df[,i+1]=ifelse(termTemp_df[,1] %in% temp_comp, 1,0)
}
rownames(termTemp_df)=termTemp_df[,1]
colnames(termTemp_df)[2:ncol(termTemp_df)]=names(modxterm_List)
termTemp_df=termTemp_df[,-1]
par(mar = c(4, 2, 1, 10)+0.5)
heatmap.2(as.matrix(termTemp_df),col= greenred(100), scale="none",trace="none", key = F,
          cexRow =0.7, cexCol = 1.2, srtCol = 45)
library(gplots)
library(heatmap3)
# write.table(df_pathEnrich, paste(outbname,"modules_queryPath_enrichment.txt",sep = "_"),
#            sep="\t", row.names = TRUE, quote = F, col.names = NA)
# plot enriched query pathways for modules with bubble plot
cmat=df_pathEnrich
ccmat=data.frame("Pathway"=rep(rownames(cmat),ncol(cmat)),
                 "Coefficient"=as.numeric(as.matrix(cmat[,1:ncol(cmat)])),
                 "Module"=rep(colnames(cmat),each=nrow(cmat))
)
ccmat=ccmat[!grepl("Akinete|Glycolysis|TCA",ccmat[,1]),]
ccmat$Correlation=ifelse(ccmat[,2]>0,"Positive","Negative")
ccmat$Color=ifelse(ccmat[,2]>0,"red","blue")
ccmat2=ccmat[with(ccmat, order(Module)),]
par(mar = c(4, 4, 0, 3)+0.1)
ccmat2 %>%
  ggplot(aes(x=Pathway, y=Module, size=-log10(abs(as.numeric(Coefficient))), na.rm = TRUE)) +
  geom_point(alpha=0.9,na.rm = TRUE) +
  scale_size(range = c(0, 8), name="-log10(p value)")+
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(color="black",size = 10, face = "bold"))
dev.off()

## draw plots of 2 year/5 species after pooling the data
# heatmap of one or more species
tempList_df=as.data.frame(matrix(1,1,8))
colnames(tempList_df)=c("TN","NH4+","NO2-","NO3-","TOC","TP","SRP","popAbd")
for (i in 1:length(list_mod_trait)) {
  #i=2
  tempName=names(list_mod_trait)[i]
  tempList=list_mod_trait[grepl(tempName,names(list_mod_trait))]
  tempList_df0=list_mod_trait[[i]]
  rownames(tempList_df0)=gsub("_b1[56]","",
                              paste(names(list_mod_trait)[i], rownames(tempList_df0), sep = "_"),perl=T)
  tempList_df=rbind(tempList_df,tempList_df0)
}
tempList_df=tempList_df[-1,]
tempList_df$color=ifelse(grepl("Pla|Ana",rownames(tempList_df)),"red","black")

tempList_df2=tempList_df[grepl("Ana|Pla",rownames(tempList_df)),]
tempList_df2=tempList_df2[grepl("b16",rownames(tempList_df2)),]
pdf("b16-Ana_Pla_module-trait_heatmap.pdf",8,6)
heatmap.2(as.matrix(tempList_df2[,-ncol(tempList_df2)]),trace = "none",col= greenred(200), scale="none",symkey = T,
          density.info="none", dendrogram=c("both"),cexRow =1, cexCol = 1.2)
dev.off()
# plot Upset and GO terms
siteyear=c("Cyl","Pla","b15","b16")
for (i in 1:4) {
  i=4
  tempName=siteyear[i]
  tempList=list_mod_topGO[grepl(tempName,names(list_mod_topGO))]
  tempList=tempList[grepl("Pla|Cyl",names(tempList))]
  termList=list()
  termvList=list()
  for (j in 1:length(tempList)) {
    listname=names(tempList)[j]
    termList[[listname]]=tempList[[j]][,1]
    termvList[[listname]]=tempList[[j]]
  }
  # plot Upset for GO terms only
  pdf("b16_Cyl_Pla_modules_Upset.pdf",5,5)
  par(mar = c(2, 5, 0, 1)+0.1)
  upset(fromList(termList[]), nsets = length(termList),nintersects = NA,order.by = "freq",
        point.size = 3, matrix.color = "grey3", 
        mainbar.y.label = "Number of GO terms in intersections", sets.x.label = "Number of GO terms",
        text.scale = c(1.6, 1.6, 1.3, 1.5, 1.2, 1.6)
  )
  dev.off()
  # plot heatmap for GO terms with p values
  termTemp_all=data.frame("GO.ID"="xxxx","Term"="yyyy")
  for (i in 1:length(termvList)) {
    termTemp_all=rbind(termTemp_all,termvList[[i]][,1:2])
  }
  termTemp_all=termTemp_all[-1,]
  termTemp_all=termTemp_all[!duplicated(termTemp_all),]
  for (i in 1:length(termvList)) {
    termTemp_all=merge(termTemp_all,termvList[[i]][,c(1,3)], "GO.ID","GO.ID",all.x=T)
    colnames(termTemp_all)[2+i]=names(termvList)[i]
  }
  termTemp_all2=termTemp_all[,-1:-2]
  termTemp_all2 <- as.data.frame(sapply(termTemp_all2, as.numeric))
  rownames(termTemp_all2)=paste(termTemp_all[,1],termTemp_all[,2],sep = "|")
  termTemp_all2[is.na(termTemp_all2)]<-0
  termTemp_all2[termTemp_all2>0]=1
  termTemp_all2=termTemp_all2[,colSums(termTemp_all2)>0]
  pdf("b16_Cyl_Pla_module_GOterms.pdf",10,6)
  par(mar = c(6, 1, 1, 2)+0.5)
  heatmap3(termTemp_all2,col= bluered(100), scale="none", cexRow =0.8, cexCol = 1.2, 
           na.rm = T, 
           legendfun=function() showLegend(legend=c("Presence","Absence"),
                                           col=c("red","blue")))
  dev.off()
}

# get the GO terms for all five species of two years
rm(list = ls())
library(heatmap3)
library(gplots)
load("b15-b16_WGCNA_5sp_mod-trait-topGO.RData")
termsAll=c()
for (i in 1:length(list_mod_topGO)) {
  termsAll=c(termsAll,
             paste(list_mod_topGO[[i]][,1],list_mod_topGO[[i]][,2],sep = "|"))
}
termsAll=unique(termsAll)
termsAll=data.frame("terms"=termsAll)
for (i in 1:length(list_mod_topGO)) {
  temp_terms=list_mod_topGO[[i]]
  temp_terms$cterms=paste(temp_terms[,1],temp_terms[,2],sep = "|")
  termsAll=merge(termsAll, temp_terms[,3:4],by.x = "terms",by.y = "cterms",all.x = T)
  colnames(termsAll)[1+i]=names(list_mod_topGO[i])
}
termsAll_df=apply(termsAll[,-1],2,as.numeric)
rownames(termsAll_df)=termsAll[,1]
termsAll_df=termsAll_df[rowSums(is.na(termsAll_df))<ncol(termsAll_df),]
termsAll_df=termsAll_df[,colSums(is.na(termsAll_df))<nrow(termsAll_df)]
termsAll_df=ifelse(is.na(termsAll_df),0,1)
termsAll_sdf=data.frame("termIDs"=rownames(termsAll_df))
siteyear=apply(expand.grid(c("b15","b16"), c("Ana","Cyl","Mic","Nos","Pla")), 1, paste, collapse="_")
for (i in 1:length(siteyear)) {
  tempdf=as.data.frame(termsAll_df[,grepl(siteyear[i], colnames(termsAll_df))])
  tempdf$presence=ifelse(rowSums(tempdf)>=1,1,0)
  tempV=rownames(tempdf)[tempdf$presence==1]
  termsAll_sdf[,i+1] = ifelse(termsAll_sdf[,1] %in% tempV, 1,0)
  colnames(termsAll_sdf)[1+i]=siteyear[i]
}
termsAll_sdf2=termsAll_sdf[,-1]
rownames(termsAll_sdf2)=termsAll_sdf[,1]
termsAll_sdf2=termsAll_sdf2[,colSums(termsAll_sdf2)>0]
pdf("b15-b16_5sp_wgcna_all-mod-GOterms.pdf",10,6)
par(mar = c(6, 1, 1, 10)+0.5)
heatmap3(termsAll_sdf2[,c(8,9, 1:2, 7,5:6,3:4)],col= bluered(200), scale="none", cexRow =0.3, cexCol = 1.2, 
         na.rm = T, Colv = NA,
         legendfun=function() showLegend(legend=c("Presence","Absence"),
                                         col=c("red","blue")))
heatmap3(termsAll_sdf2[,c(8,1,5,3,9,2,6,4)],col= bluered(200), scale="none", cexRow =0.3, cexCol = 1.2, 
         na.rm = T, Colv = NA, 
         legendfun=function() showLegend(legend=c("Presence","Absence"),
                                         col=c("red","blue")))

dev.off()

#################################################################################################
#################################################################################################
##### sort enrichment of query pathways and GO terms #####
### sort query pathways
rm(list = ls())
pathNames=c("AAT","Akinete","CCM","Glycolysis","MAA","N","NF","Osmosis","OSR",
            "P","PBS","PSI","PSII","PUFA","Sugar","Sulfur","TCA","TEVit","Toxin","Vesicle")
poolpaths=as.data.frame(matrix(2, 20,1))
poolpaths[,1]=pathNames
colnames(poolpaths)[1]="pathName"
QEfiles=list.files(pattern = "Path_enrichment")
site_year=c("b15","b16","e15","e16")
species=c("ana","cyl","mic","nos","pla")
cNames="pathName"
for (gg in 1:5) {
  #gg=1
  Gfiles=QEfiles[grepl(species[gg],QEfiles)]
  for (s in 3:4) {
    #s=2
    sfile=read.delim(Gfiles[s])
    poolpaths=merge(poolpaths,sfile,by.x = "pathName",by.y = "X", all.x = T)
    cNames=c(cNames,paste(species[gg], site_year[s], colnames(sfile)[-1],sep = "_"))
  }
  colnames(poolpaths)=cNames
}
DF=poolpaths
is.num <- sapply(DF, is.numeric)
DF[is.num] <- lapply(DF[is.num], round, 4)
write.table(DF, "e15-16_five_species_modules_queryPath_enrichment.txt",
            sep="\t", row.names = TRUE, quote = F, col.names = NA)
DF2=DF[,-1]
rownames(DF2)=DF[,1]
DF3=DF2[,c()]
## bubble plot
rm(list=ls())
library(ggplot2)
library(dplyr)
cmat=read.delim("b15-16_five_species_modules_queryPath_enrichment.txt",
                sep = "\t",row.names = 1,check.names = FALSE, stringsAsFactors = F)
cmat[cmat=="-"]<-NA
ccmat=data.frame(cmat[,1:2],"Sp-year"=rep(colnames(cmat)[2],20))
colnames(ccmat)=c("Pathway","Coefficient","Species")
for (i in 3:ncol(cmat)) {
  tdf=data.frame(cmat[,c(1,i)],"Species"=rep(colnames(cmat)[i],20))
  colnames(tdf)=c("Pathway","Coefficient","Species")
  ccmat=rbind(ccmat,tdf)
}
ccmat[,2]=as.numeric(ccmat[,2])
ccmat$Sign=ifelse(ccmat[,2]>0,"apositive","negative")
# Most basic bubble plot
ccmat[is.na(ccmat[,4]), 3]<-NA
ccmat=ccmat[!grepl("Akinete",ccmat[,1]),]
ccmat$Pathway=gsub("Glycolysis","zGlycolysis",ccmat$Pathway)
ccmat$Pathway=gsub("TCA","zTCA",ccmat$Pathway)
ccmat[,3]=gsub("-","_",ccmat[,3],perl = T)
ccmat2=ccmat[with(ccmat, order(Species)),]
ccmat2=ccmat2[grepl("b1[56]",ccmat2[,3],perl = TRUE),]
#pdf("SPLS_Query-paths_vs_species-abundance_buoy_2years2.pdf",6.5,3.5)
ccmat2 %>%
  #mutate(Paths = factor(Paths, Paths)) %>%
  ggplot(aes(x=Pathway, y=Species, size=-log10(abs(as.numeric(Coefficient))),color=Sign, na.rm = TRUE)) +
  geom_point(alpha=0.9,na.rm = TRUE) +
  scale_size(range = c(0, 8), name="Coefficient")+
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 12, face = "bold"))
dev.off()

##### sort GO terms with heatmap and Upset
rm(list = ls())
options(stringsAsFactors = FALSE)

# get all GO terms
goFolders=list.files(pattern = "WGCNA_")
goFiles=""
for (i in 1:4) {
  goFiles_temp=file.path(goFolders[i],list.files(path = goFolders[i],pattern = "_go.txt"))
  goFiles=c(goFiles, goFiles_temp)
}
goFiles=goFiles[-1]
goDF=data.frame("GO.ID","Term")
colnames(goDF)=c("GO.ID","Term")
for (i in 1:length(goFiles)) {
  #i=1
  tempDF=read.delim(goFiles[i])
  goDF=rbind(goDF,tempDF[tempDF$classicFisher<0.05,1:2]) # significantly enriched GO terms
}
goDF=goDF[!duplicated(goDF),]
goDF=goDF[-1,]
colnames(goDF)=c("goID","Annotation")

# get all GO terms all sites and years
goterms=goDF
site_year=c("b15","b16","e15","e16")
species=c("ana","cyl","mic","nos","pla")
cName=c("goID","Annotation")
for (i in 1:length(species)) {
  #i=1
  goFiles_temp1=goFiles[grepl(species[i],goFiles)]
  goFiles_temp2=""
  for (j in 1:length(site_year)) {
    #j=1
    goFiles_temp2=c(goFiles_temp2, goFiles_temp1[grepl(site_year[j],goFiles_temp1)])
  }
  goFiles_temp2=goFiles_temp2[-1]
  #read in all files
  for (jj in 1:length(goFiles_temp2)) {
    #jj=1
    goFile_df=read.delim(goFiles_temp2[jj],header = T,sep = "\t", row.names = 1)
    goterms=merge(goterms,goFile_df[,c(1,6)],by.x = "goID",by.y = "GO.ID", all.x = T)
    moduleName1=gsub("WG.+\\/|ME_|_go\\.txt","",goFiles_temp2[jj], perl = T)
    moduleName2=sub("\\w+?_", paste(species[i],"_",sep = ""),moduleName1, perl = T)
    cName=c(cName,moduleName2)
    print(cName)
    colnames(goterms)=cName
  }
}
write.table(goterms,"matrix_all-go-terms_all-samples.txt",sep="\t",
            row.names = TRUE, quote = F, col.names = NA)

# plot these GO terms based site and year
# 2015 and 2016 together
goterms2=goterms[,-1:-2]
rownames(goterms2)=goterms[,1]
#b_terms=c("_b1[56]_","_b15_","_b16_")
b_terms=c("_e1[56]_","_e15_","_e16_")
library(heatmap3)
pdf("efls15-16_5species_module_GOs.pdf", 8,10)
par(mar=c(7,2,0,0))
for (ii in 1:length(b_terms)) {
  #  ii=1
  goterms2_b15_16=goterms2[,grepl(b_terms[ii],colnames(goterms2), perl = T)]
  goterms2_b15_16[goterms2_b15_16>0.05]<-NA
  goterms2_b15_16=goterms2_b15_16[rowSums(is.na(goterms2_b15_16[,1:ncol(goterms2_b15_16)])) < ncol(goterms2_b15_16),]
  Species=c(rep("darkgreen", length(colnames(goterms2_b15_16)[grepl(species[1],colnames(goterms2_b15_16))])),
            rep("brown", length(colnames(goterms2_b15_16)[grepl(species[2],colnames(goterms2_b15_16))])),
            rep("orange", length(colnames(goterms2_b15_16)[grepl(species[3],colnames(goterms2_b15_16))])),
            rep("purple", length(colnames(goterms2_b15_16)[grepl(species[4],colnames(goterms2_b15_16))])),
            rep("red", length(colnames(goterms2_b15_16)[grepl(species[5],colnames(goterms2_b15_16))]))
  )
  heatmap3(as.matrix(goterms2_b15_16), Colv = NA, Rowv = NA, labRow = "", scale = "none",
           ColSideColors=Species, ColAxisColors = 1 )
}
dev.off()

# plot each species between 15 and 16
pdf("efls15-16_1-species_module_GOs.pdf", 5.5,10)
par(mar=c(7,3,1,0))
for (ii in 1:length(species)) {
  goterms2_b15_16=goterms2[,(grepl(species[ii],colnames(goterms2), perl = T) & grepl("_b1",colnames(goterms2), perl = T))]
  goterms2_b15_16[goterms2_b15_16>0.05]<-NA
  goterms2_b15_16=goterms2_b15_16[rowSums(is.na(goterms2_b15_16[,1:ncol(goterms2_b15_16)])) < ncol(goterms2_b15_16),]
  sideLabs=ifelse(grepl("15",colnames(goterms2_b15_16)),"blue","red")
  heatmap3(as.matrix(goterms2_b15_16), Colv = NA, Rowv = NA, labRow = "", scale = "none", 
           margins = c(10, 5), ColSideColors=sideLabs, ColAxisColors = 1 )
}
dev.off()

### sort enrichment of GO terms with UpsetR
rm(list = ls())
options(stringsAsFactors = FALSE)
library(UpSetR)
# example of list input (list of named vectors)
listInput <- list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13),
                  two = c(1, 2, 4, 5, 10),
                  three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))
upset(fromList(listInput), order.by = "freq")
# put all GO terms into a list, one sample as a named vector
listInput=list()
for (jj in 1:length(goFiles)) {
  moduleName1=gsub("WG.+\\/|wgcna_ME_|_go\\.txt","",goFiles[jj], perl = T)
  moduleName2=sub("([a-z])", "ME\\1",moduleName1, perl = T)
  goFile_df=read.delim(goFiles[jj],header = T,sep = "\t", row.names = 1)
  goFile_df=goFile_df[goFile_df$classicFisher<0.05,1:2]
  listInput[[moduleName2]]=paste(goFile_df[,1],goFile_df[,2],sep = "|")
}
pdf("wgcna_modules_goterms_venndiagram_UpsetR.pdf", 12,12)
upset(fromList(listInput), nsets = length(listInput),nintersects = NA,order.by = "freq",
      point.size = 3, matrix.color = "grey3", 
      mainbar.y.label = "GO-term intersection size", sets.x.label = "Number of GO terms",
      text.scale = c(1.6, 2, 1.5, 1.5, 1.5, 1.8)
)
dev.off()

######################################################################################
######################################################################################
###### Plot pop abundance with gene-expr profile clustering tree  #######
rm(list=ls())
library(WGCNA)
library(heatmap3)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
collectGarbage()
#====== 1. get the folder with gene abundance matrix
folders=c("buoy_2015","buoy_2016","efls_2015","efls_2016")
abfiles=file.path("pop_abundance",list.files("pop_abundance",pattern = "newDate.txt"))
sites=c("BUOY","BUOY","EFLS","EFLS","b15","b16","e15","e16")
# loop through each site and the two years
for (f in 1:4){
  #f=1  # f: study site and year
  #====== 1.c Loading water chemistry metadata
  pop_abd=read.delim(abfiles[f],header = T, sep = "\t",row.names = 1,check.names = F)
  pop_abd=t(pop_abd[order(rownames(pop_abd)),])
  #==== 1.e five genome gene expression files at each site/year
  # pattern may differ
  genefiles=file.path(folders[f],list.files(path = folders[f],pattern="gene_abundance_tpm"))
  sizeGrWindow(12,9)
  par(cex = 0.6);
  pdf(paste("Tree_order",sites[4+f],"all-sample_tree_abundance_combined.pdf",sep = "_"), 4,20)
  par(mfrow=c(10,1),mar = c(0.5,5,2,0.5))
  for (g in 1:5){
    #  g=1
    # get gene expression profile for each genome
    tpm=t(read.delim(genefiles[g],header = T, sep = "\t", check.names = F, row.names = 1))
    dates=sort(unique(gsub("[ab]$","", rownames(tpm))))
    mat=as.data.frame(matrix(0, nrow = length(dates), ncol = ncol(tpm)))
    rownames(mat)=dates
    colnames(mat)=colnames(tpm)
    for (d in 1:length(dates)) {
      tempi=as.matrix(tpm[grepl(dates[d],rownames(tpm)),])
      class(tempi)
      if(nrow(tempi)>1){
        tempii=apply(tempi,2,mean)
        mat[d,]=as.numeric(tempii)
      } else {
        mat[d,]=as.numeric(tempi)
      }
    }
    mat.ra=mat/rowSums(mat)
    mat.ra=mat.ra[,colSums(mat.ra>0)>3]
    otu=mat.ra
    outbname=paste(unlist(strsplit(gsub("\\/|_gene.+$","-",genefiles[g]),"-"))[2],sites[f+4],sep = "_")
    #=== 1.f WGCNA analysis
    gsg = goodSamplesGenes(otu, verbose = 3);
    gsg$allOK
    # check for outliers
    sampleTree = hclust(dist(otu), method = "average");
    pt<-plot(sampleTree, main = outbname, sub="", xlab="Sample", 
             cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.2, cex=1.5,las=2
    )
    legend("topright",LETTERS[2*g-1], cex = 1.2, bty = "n", text.font =2, bg = "gray90")
    tdf=data.frame(Date=1:nrow(pop_abd),Abundance=pop_abd[,g][sampleTree$order],
                   Rel_Ab=pop_abd[,g]/rowSums(pop_abd[,1:5])[sampleTree$order])
    xrange <- range(tdf$Date)
    yrange <- range(tdf$Abundance)
    #plot(tdf$Date, tdf$Abundance, type="b", lwd=1.5,las=2,lty=2, col="cyan", cex.axis = 1.2,
    #     cex.lab=1.5, ylab="Abundance", xaxt="n")
    plot(tdf$Date, tdf$Abundance, type="b", lwd=1.5, las=2, lty=1, col="black", cex.axis = 1.2,
         cex=2.5,pch=20,cex.lab=1.5,ylab="Population Abundance", xaxt="n")
    legend("topright",LETTERS[2*g], cex = 1.5, bty = "n", text.font =2, bg = "gray90")
  }
  dev.off()
}
######################################################################################
##################################### SNP analysis #################################
########## SNP analysis 1: GO enrichment
##### plot total stats data #####
rm(list=ls())
options(stringsAsFactors = F)
SNP_folders=file.path("SNP",list.files("SNP",pattern = "snp_sum_"))
species=c("ana","cyl","mar","nos","plk")
#for (i in 1:length(SNP_folders)) {
i=1
fpaths=file.path(SNP_folders[i],list.files(SNP_folders[i]))
SNP_folders=basename(fpaths)
outname=unique(gsub("\\/.+$","",gsub("SNP\\/snp_sum_","",fpaths, perl = T), perl=T))
pdf(paste(outname,"_variants_3sp_stat.pdf",sep = ""),16,20)
par(mfrow=c(5,5),mar=c(7,5,1,0))
for (j in 1:length(species)) {
  #j=5
  temp_fpath=fpaths[grepl(species[j],fpaths)]
  tNsamples=length(file.path(temp_fpath[3],list.files(temp_fpath[3],pattern = "stats.csv")))
  tdf=as.data.frame(matrix(0,9,(3*tNsamples+1)))
  tdf[,1]=c("totalVar","sSNP","nsSNP","Indel","nonCDS","r-sSNP","r-nsSNP","r-Indel","r-nonCDS")
  for (g in 1:length(temp_fpath)) {
    #  g=3
    tSNP_files=file.path(temp_fpath[g],list.files(temp_fpath[g],pattern = "stats.csv"))
    colnames(tdf)[((g-1)*tNsamples+2):(g*tNsamples+1)]=gsub("_stats.csv","",basename(tSNP_files))
    for (k in 1:length(tSNP_files)) {
      #k=1
      x=read.csv(tSNP_files[k],header = T)
      tdf[,((g-1)*tNsamples+1+k)]=unlist(c(x[1,],x[1,2:ncol(x)]/x[1,1]))
    }
  }
  tdf2=as.data.frame(t(tdf[,-1]))
  colnames(tdf2)=tdf[,1]
  #sDate=gsub("(^\\d)","\\1-", gsub("\\w+_0","",rownames(tdf2)[1:tNsamples]),perl = T)
  sDate=gsub("^\\w+?_","",rownames(tdf2)[1:tNsamples], perl = T)
  tdf2$Date=rep(sDate,3)
  for (r in 1:5) {
    #  r=1
    pt<-plot(1:tNsamples, tdf2[1:tNsamples,r], type="b", lwd=1.5, las=2, lty=1, col="blue", 
             pch=1, cex.axis = 1.1, cex.lab=1.1,ylab="Count",xaxt="n", xlab ="",main = species[j])
    legend("topright",colnames(tdf2)[r], cex = 1.2, bty = "n", text.font =2, bg = "gray90")
    axis(1,at = 1:tNsamples,labels = sDate, las=2)
    points(1:tNsamples, tdf2[(tNsamples+1):(2*tNsamples),r],col="red", pch="*",cex=1.6)
    lines(1:tNsamples, tdf2[(tNsamples+1):(2*tNsamples),r], lwd=1.5, las=2, lty=2, col="red",
          cex.axis = 1.1,cex.lab=1.1,ylab="",xaxt="n", xlab = "")
    points(1:tNsamples, tdf2[(2*tNsamples+1):(3*tNsamples),r],col="green", pch=2)
    lines(1:tNsamples, tdf2[(2*tNsamples+1):(3*tNsamples),r], lwd=1.5, las=2, lty=2, col="green",
          cex.axis = 1.1,cex.lab=1.1,ylab="",xaxt="n", xlab = "") 
  }
}
dev.off()
#}
##########
##### SNP plot gene expression levels and GO/query patway enrichement #####
rm(list=ls())
library(topGO)
library(circlize)
library(randomcoloR)
options(stringsAsFactors = FALSE)
##### 1. get the folders with gene exppression files
gene_folders=c("buoy_2015","buoy_2016","efls_2015","efls_2016")
sites=c("BUOY","BUOY","EFLS","EFLS","b15","b16","e15","e16")
gene_files=file.path(gene_folders,list.files(gene_folders, pattern = "gene_abundance"))
##### 2. get SNP folders
SNP_folders=file.path("SNP",list.files("SNP",pattern = "snp_sum_"))
SNP_folders=SNP_folders[!grepl("_stat.pdf",SNP_folders)]
species=c("ana90","cylCS505","marNIES843","Nostoc7107","plkNIEVCYA1268")
palette <- distinctColorPalette(length(species))
##### 3. get GOA files of all five genomes
goafiles=file.path("GOAfiles",list.files("GOAfiles",pattern = "topGO_GOA"))
goafiles=goafiles[!grepl("NIES3756",goafiles)]
snpGOlist=list()
##### 4. Analysis start: lopping through sites and year of SNP profiles
# loop through site with loop through species embedded
#for (i in 1:length(SNP_folders)) {
i=4
temp_list=list()
# read in gene expression file paths
gene_files=file.path(gene_folders[i],list.files(gene_folders[i], pattern = "dance_tpm.txt"))
gene_files=gene_files[!grepl("nies3756",gene_files)]
# read in snp file paths
snp_fpaths=file.path(SNP_folders[i],list.files(SNP_folders[i]))
snp_fpaths=snp_fpaths[grepl("ana90|cylCS505|marNIES843|nos7107|plkNIEVCYA1268",snp_fpaths)]
#pdf(paste(basename(SNP_folders[i]),"_stat.pdf",sep = ""),16,8)
#par(mfrow=c(1,1),mar=c(7,5,1,0))
# loop through species
#par(mfrow=c(4,1))
for (j in 3:4) {
  j=5
  # get geneIDs of expressed genes in all the samples
  geneX=read.delim(gene_files[j], check.names = F)
  colnames(geneX)[1]="locusID"
  geneX=geneX[rowSums(is.na(geneX[,2:ncol(geneX)])) < (ncol(geneX)-1),]
  # get all the locus tag from all the samples
  snp_files=file.path(snp_fpaths[j],list.files(snp_fpaths[j],pattern = "[ABab]\\.csv"))
  Nsamples=length(snp_files)
  DFsnp=data.frame("locus_tag"=geneX[,1])
  for (s in 1:Nsamples) {
    #s=1
    snpX=read.csv(snp_files[s],header = T)
    tdf=snpX[1:length(unique(snpX[,ncol(snpX)])),3:4]
    tdf[,1]=unique(gsub('\\"| |gene_id',"",snpX[,ncol(snpX)], perl = T))
    for (l in 1:nrow(tdf)) {    # l=loucs tag
      #l=1
      ttdf=snpX[grepl(tdf[l,1],snpX[,6]),3:5]
      type=unique(ttdf[,1])
      type[which(type=="S")]="SY"
      variants=paste(unique(ttdf[,2]),collapse = "|")
      tdf[l,2]=paste(paste(type, collapse="|"),variants,sep =  ";")
    }
    DFsnp=merge(DFsnp,tdf,by.x = "locus_tag",by.y = "MUTATION",all.x=T)
    colnames(DFsnp)[s+1]=gsub("^.+?_0*|\\.csv", "",basename(snp_files[s]), perl = T)
  } # all samples are read in
  
  #  ## 4.1 this section plots the figures of the genes with SNP across all the samples
  #  DFsnp2=DFsnp[rowSums(is.na(DFsnp[,2:ncol(DFsnp)])) < (ncol(DFsnp)-1),]
  #  if(j==1){
  #    plot(1:Nsamples, colSums(!is.na(DFsnp2[,2:ncol(DFsnp2)])), type="b",
  #              lwd=1.2, cex.axis = 1.1, cex.lab=1.1, ylim=c(0,500),
  #              col=palette[j], ylab="Count",xaxt="n", xlab ="")
  #     axis(1,at = 1:Nsamples,labels = colnames(DFsnp2)[2:(Nsamples+1)], las=2)
  #  }else{
  #  points(1:Nsamples, colSums(!is.na(DFsnp2[,2:ncol(DFsnp2)])), col=palette[j], pch=j,cex=1.2)
  #     lines(1:Nsamples, colSums(!is.na(DFsnp2[,2:ncol(DFsnp2)])), lwd=1.2, lty=2,
  #           col=palette[j], ylab="",xaxt="n", xlab = "")
  #  }
  #  legend("topleft", legend = species, col = palette, pch = 1:5,
  #         lty = 1, merge = TRUE)
  # }
  # dev.off()
  
  # ## 4.2 this section plot the expresion levels of the genes of SNP types
  # DFsnp2=DFsnp[rowSums(is.na(DFsnp[,2:ncol(DFsnp)])) < (ncol(DFsnp)-1),]
  # pdf(paste(gene_folders[i],species[j],"SNP_gene-tpm.pdf",sep = "_"),6,30)
  # if((ncol(DFsnp2) %% 2)==0){
  #   par(mfrow=c(ncol(DFsnp2)/2,2),mar=c(4.2,4,0,0)+0.1)
  # }else{
  #     par(mfrow=c((ncol(DFsnp2)-1)/2,2),mar=c(4.2,4,0,0)+0.1)
  # }
  # pvalues=c()
  # for (s in 2:ncol(DFsnp2)){
  #   tryCatch({
  #     #s=2
  #     snpXdf=merge(geneX[,c(1,s)],DFsnp2[,c(1,s)],by.x = "locusID",by.y = "locus_tag",all.x = T)
  #     snpXdf=snpXdf[which(snpXdf[,2]>0),]
  #     All=snpXdf[,2]
  #     SNV=snpXdf[(grepl("NS",snpXdf[,3]) | grepl("SY",snpXdf[,3])),2]
  #     NS=snpXdf[grepl("NS",snpXdf[,3]),2]
  #     SY=snpXdf[grepl("SY",snpXdf[,3]),2]
  #     noSNV=snpXdf[is.na(snpXdf[,3]),2]
  #     cdf=data.frame("X"=c(All,SNV,NS,SY,noSNV),
  #                    "Group"=c(rep("All",length(All)),rep("SNV",length(SNV)),
  #                                rep("NS",length(NS)),rep("SY",length(SY)),
  #                                rep("noSNV",length(noSNV)))
  #     )
  #     boxplot(cdf[,1]~cdf[,2],las=2, ylab="Gene level (TPM)",cex.axis=1.2,cex.lab=1.2,
  #               outline=F,boxwex=0.2, cex=1.2)
  #       abline(h=median(snpXdf[,2]), col="blue")
  #       kt=wilcox.test(All, SNV)
  #       legend("topleft",legend = paste(colnames(DFsnp2)[s],"p =",round(kt$p.value,8), sep = " "),
  #              bty = "n",cex = 1.2)
  #   },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # }
  # dev.off()
  # } # species level across all samples
  
  ####### 4.3 GO term enrichment of the genes with SNVs
  #####  4.3.1 get GOA files ready for GO enrichment of genes of each species
  geneID2GO<-readMappings(file = goafiles[j])
  allGenes <- unique(names(geneID2GO))
  DFsnp2=DFsnp[rowSums(is.na(DFsnp[,2:ncol(DFsnp)])) < (ncol(DFsnp)-1),] # remove all na rows
  #outPDF=paste(sites[4+i],species[j], "SNP_GO-enrich.pdf",sep = "_")
  #pdf(outPDF,10,10)
  #for (g in c(2:15,17:31)) { #i2 j1g16
  #for (g in c(2,3,7:8,10,13,18,19, 21:26,28:29)) { #  i2 j5 planktothrix
  #for (g in c(2:7,9:ncol(DFsnp2))) { #i3,j2g8
  #for (g in c(2:7,9:ncol(DFsnp2))) { #i3 j4g8
  #for (g in c(2:7,10,19,21:ncol(DFsnp2))) { #i3 j5 g8g9g20
  for (g in c(3:7,10:19,21:ncol(DFsnp2))) { #i4 j2g8 j4g8 j5g2g8g9g20
    #for (g in c(2:ncol(DFsnp2))) {
    #g=16
    moduleGenes=DFsnp2[!is.na(DFsnp2[,g]),1]
    moduleGenes_factor <- factor(as.integer(allGenes %in% moduleGenes))
    names(moduleGenes_factor)=allGenes
    GOdata <- new("topGOdata", ontology = "BP", allGenes = moduleGenes_factor,
                  nodeSize = 10, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    #resultFisher
    allRes <- GenTable(GOdata, classicFisher = resultFisher,
                       orderBy = "classicFisher", topNodes = 100)
    #allRes
    outlist_name=paste(species[j],colnames(DFsnp2)[g],sep = "_")
    print(paste(g, "is good", sep = " "))
    temp_list[[outlist_name]]=allRes # list of samples with GO enrichment results
    #print(showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all'))
    #print(legend("topleft",paste("SNP-GO",sites[4+i],species[j],colnames(DFsnp2[g]),sep = "_"),bty = "n"))
  }
  #dev.off()
}
snpGOlist[[sites[4+i]]]=temp_list
#}
# save(snpGOlist,file = "snpGO_list.RData")

##### 4.3.2 sort the SNP and WGCNA GO enrichment of genes of each species (a table of all GO temrs were obtained)
#### sort GO enrichment results and do circos plots
## WGCNA GO terms plotting
rm(list = ls())
library(heatmap3)
library(RColorBrewer)
options(stringsAsFactors = FALSE)
## WGCNA GO terms
allGOs=read.delim("matrix_all-go-terms_all-samples.txt",
                  row.names = 1, header = T, sep = "\t")
species=as.vector(outer(c("ana","cyl","mic","nos","pla"),c("b15","b16","e15","e16"), 
                        paste, sep = "_"))
allGOs_mat=allGOs[,1:2]
for (i in 1:length(species)) {
  #i=4
  ttdf1=as.data.frame(allGOs[, grepl(species[i], colnames(allGOs))])
  ttdf1$GOs=allGOs[,1]
  ttdf2=data.frame(goID="GOs",pvalue="pvalue")
  for (j in 1:(ncol(ttdf1)-1)) {
    #j=1
    ttdf3=ttdf1[,c(ncol(ttdf1),j)]
    colnames(ttdf3)=colnames(ttdf2)
    ttdf3=ttdf3[which(ttdf3[,2]<0.05),]
    ttdf3 = ttdf3[order(ttdf3[,1], ttdf3[,2]),] # remove redundant GO terms
    ttdf3 = ttdf3[!duplicated(ttdf3$goID),]
    ttdf2=rbind(ttdf2,ttdf3)
  }
  allGOs_mat=merge(allGOs_mat,ttdf2,by.x = "goID",by.y = "goID",all.x = T)
  colnames(allGOs_mat)[2+i]=species[i]
}
allGOs_mat = allGOs_mat[order(allGOs_mat[,1], allGOs_mat[,2]),] # remove redundant GO terms
allGOs_mat = allGOs_mat[!duplicated(allGOs_mat$goID),]
allGOs_wgcna=data.matrix(allGOs_mat[,-1:-2])
rownames(allGOs_wgcna)=paste(allGOs_mat[,1],allGOs_mat[,2],sep = "|")
allGOs_wgcna[is.na(allGOs_wgcna)]<-0
allGOs_wgcna[allGOs_wgcna>0]<-1
allGOs_wgcna=allGOs_wgcna[rowSums(allGOs_wgcna)>1,]
write.table(allGOs_wgcna,file = "WGCNA_GO-terms_5sp_2sites_sort-p05.txt",quote = F,
            row.names =T, sep = "\t",col.names = NA)
save(allGOs_wgcna,file = "WGCNA_GO-terms_5sp_2sites_sort-p05.Rdata")
# tdf2=data.matrix(allGOs_wgcna[,grepl("_e1[56]", colnames(allGOs_wgcna),perl = T)])
# tdf2=tdf2[rowSums(tdf2)>0,]
# palette2=brewer.pal(n=10, name = "Paired")
# ColSideColors<-cbind(
#   species= rep(brewer.pal(n=5, name = "YlOrRd"),each=2),
#   Year=palette2
# )
#   #pdf("Buoy_15-16_WGCNA_GOs_5sp_allSamples_005_pairwise.pdf",16,16)
#   #pdf("Efls_15-16_WGCNA_GOs_5sp_allSamples_005_pairwise.pdf",16,16)
#   #pdf("Buoy_15-16_WGCNA_GOs_5sp_allSamples_001_pairwise.pdf",16,16)
#   #pdf("Efls_15-16_WGCNA_GOs_5sp_allSamples_001_pairwise.pdf",16,16)
#   par(mfrow=c(1,1))
#   heatmap3(tdf2[,c(1,6,2,7,3,8,4,9,5,10)], Colv = NA, Rowv = T, scale = "none", 
#            cexRow = 0.8, cexCol = 1.2, ColSideColors = ColSideColors, ColAxisColors = 2, 
#            distfun = function(x) dist(x,method = 'euclidean'))
#   dev.off()

### SNP GO terms plotting
rm(list=ls())
library(RColorBrewer)
library(heatmap3)
load("snpGO_list.RData")
load("4sites_SNP_GOterms.RData")
load("WGCNA_GO-terms_5sp_2sites_sort-p05.Rdata")
tdf=goterms[-1,]
colnames(tdf)=c("GOs","Annotation")
for(si in seq_len(length(snpGOlist))){
  #si=1
  templist1=snpGOlist[[si]] # one site with all GO terms of five species with multipe samples
  for (ss in 1:length(templist1)) {
    #ss=1
    templist2=templist1[ss]
    tempdf=templist2[[1]]
    colnames(tempdf)[ncol(tempdf)]=paste(names(snpGOlist)[si],names(templist1[ss]),sep = "_")
    tempdf2=tempdf[as.numeric(tempdf[,ncol(tempdf)]) <=0.01,c(1,ncol(tempdf))]
    tdf=merge(tdf,tempdf2[,c(1,ncol(tempdf2))],by.x="GOs",by.y="GO.ID",all.x=T)
  }
}
tdf=tdf[!duplicated(tdf), ]

tdf1=tdf[,-1:-2]
rownames(tdf1)=paste(tdf[,1],tdf[,2],sep = "|")
#plot: b15 and b16 all species across all samples
tdf2=tdf1[,grepl("^e1[56]_",colnames(tdf1),perl = T)] # select buoy or efls
tdf2[is.na(tdf2)]<-0
tdf2=data.matrix(tdf2)
tdf2[tdf2>0]=1
tdf2=tdf2[rowSums(tdf2) < ncol(tdf2)/4 & rowSums(tdf2) > 0,]

GOcolors=ifelse(rownames(tdf2) %in% rownames(allGOs_wgcna),"red","blue")
palette=brewer.pal(n=5, name = "Set2")
ColSideColors<-cbind(
  Species=c(rep(palette[1],length(colnames(tdf2)[grepl("e15_ana90",colnames(tdf2))])),
            rep(palette[2],length(colnames(tdf2)[grepl("e15_cylCS5",colnames(tdf2))])),
            rep(palette[3],length(colnames(tdf2)[grepl("e15_marNI",colnames(tdf2))])),
            rep(palette[4],length(colnames(tdf2)[grepl("e15_Nostoc",colnames(tdf2))])),
            rep(palette[5],length(colnames(tdf2)[grepl("e15_plkN",colnames(tdf2))])),
            rep(palette[1],length(colnames(tdf2)[grepl("e16_ana90",colnames(tdf2))])),
            rep(palette[2],length(colnames(tdf2)[grepl("e16_cylCS5",colnames(tdf2))])),
            rep(palette[3],length(colnames(tdf2)[grepl("e16_marNI",colnames(tdf2))])),
            rep(palette[4],length(colnames(tdf2)[grepl("e16_Nostoc",colnames(tdf2))])),
            rep(palette[5],length(colnames(tdf2)[grepl("e16_plkN",colnames(tdf2))]))
  ),
  Year=c(
    rep("green", length(colnames(tdf2)[grepl("^e15_",colnames(tdf2), perl = T)])),
    rep("yellow", length(colnames(tdf2)[grepl("^e16_",colnames(tdf2), perl = T)])))
)
#pdf("Buoy_15-16_SNP_GOs_5sp_allSamples.pdf",18,16)
#pdf("Buoy_15-16_SNP_GOs_5sp_allSamples_001-top33pct.pdf",16,6)
#pdf("Buoy_15-16_SNP_GOs_5sp_allSamples_001-bottom25pct.pdf",16,12)
#pdf("Efls_15-16_SNP_GOs_5sp_allSamples.pdf",18,16)
#pdf("Efls_15-16_SNP_GOs_5sp_allSamples_001-top33pct.pdf",16,6)
#pdf("Efls_15-16_SNP_GOs_5sp_allSamples_001-bottom25pct.pdf",16,12)
par(mfrow=c(1,1))
heatmap3(tdf2, Colv = NA, Rowv = T, scale = "none", cexRow = 0.7, cexCol = 0.25,
         ColSideColors = ColSideColors, ColAxisColors = 1,
         RowSideColors = GOcolors, RowAxisColors = 1,
         distfun = function(x) dist(x,method = 'euclidean'))
dev.off()
################### SNP analysis 2: Query pathways enrichment with SNVs #######################
rm(list = ls())
library(heatmap3)
library(RColorBrewer)
options(stringsAsFactors = FALSE)
##### 1: get query pathway files of all five genomes
qpaths=c("Akinete","Vesicle","Toxin","Osmosis","PUFA","AAT","PBS","MAA","NF","CCM",
         "Sugar","Sulfur","MetalR","N","TEVit","P","DrugR","OSR","PSI", "PSII")
pathNames=c(qpaths[order(qpaths)],"Glycolysis","TCA")
pathwayfiles=file.path("pathfiles",list.files("pathfiles",pattern = "pid_locus_pathway_sort"))
sites=c("b15","b16","e15","e16")
##### 2. get SNP folders
SNP_folders=file.path("SNP",list.files("SNP",pattern = "snp_sum_"))
SNP_folders=SNP_folders[!grepl("_stat.pdf",SNP_folders)]
species=c("ana90","cylCS505","marNIES843","nos7107","plkNIEVCYA1268")
qpaths_df=data.frame("queryPath"=pathNames)
##### 3. get query files of all five genomes
for (i in 1:length(SNP_folders)) {
  #i=1
  temp_list=list()
  # read in snp file paths
  snp_fpaths=file.path(SNP_folders[i],list.files(SNP_folders[i]))
  snp_fpaths=snp_fpaths[grepl("ana90|cylCS505|marNIES843|nos7107|plkNIEVCYA1268",snp_fpaths)]
  for (g in 1:length(species)) {
    #g=1
    genome_paths=read.delim(pathwayfiles[g])
    snpfiles=file.path(SNP_folders[i],species[g],list.files(snp_fpaths[g],pattern = "[ABab].csv"))
    df_pathEnrich=as.data.frame(matrix(2,length(pathNames),length(snpfiles)+1))
    samplesNames=gsub("_0","_",gsub("\\.csv","",basename(snpfiles)))
    colnames(df_pathEnrich)=c("paths",paste(sites[i],samplesNames,sep = "_"))
    df_pathEnrich[,1]=pathNames
    hyper_N=length(which(genome_paths$shortID !="NoPath"))
    # hypergeometric test
    for (s in 1:length(snpfiles)){
      #s=1
      snpX=read.csv(snpfiles[s],header = T)
      loci=unique(gsub('\\"| |gene_id',"",snpX[,ncol(snpX)], perl = T))
      temp_pvalues=c()
      for (p in 1:length(pathNames)) {
        #p=1
        df_allpath_temp=genome_paths[which(genome_paths$shortID != "NoPath"), c(2,4)]
        df_1path_temp=genome_paths[which(genome_paths$shortID==pathNames[p]), c(2,4)]
        hyper_m=nrow(df_1path_temp)
        hyper_n=hyper_N-hyper_m
        hyper_k=length(which((df_allpath_temp$Locus %in% loci) == TRUE))
        hyper_x=length(which((df_1path_temp$Locus %in% loci) == TRUE))
        (p.value <-  phyper(q=hyper_x -1, m=hyper_m, n=hyper_n, k=hyper_k, lower.tail=FALSE))
        temp_pvalues=c(temp_pvalues,p.value)
      }
      df_pathEnrich[,s+1]=temp_pvalues
    }
    qpaths_df=merge(qpaths_df, df_pathEnrich,by.x = "queryPath",by.y = "paths",all.x = T)
  }
}
#write.table(qpaths_df,file = "queryPaths_enriched_2sites-2yrs.txt",quote = F,
#            row.names = T, col.names = T, sep = "\t")
#save(qpaths_df,file = "queryPaths_enriched_2sites-2yrs.Rdata")
load("queryPaths_enriched_2sites-2yrs.Rdata")
tdf=qpaths_df[,grepl("e1[56]",colnames(qpaths_df),perl = T)]
rownames(tdf)=qpaths_df[,1]
tdf[tdf>0.05]=0
tdf[tdf<=0.05 & tdf>0]=1
#tdf=tdf[rowSums(tdf)>1,]
palette=brewer.pal(n=5, name = "Set2")
ColSideColors<-cbind(
  Species=c(rep(palette[1],length(colnames(tdf)[grepl("e15_ana90",colnames(tdf))])),
            rep(palette[2],length(colnames(tdf)[grepl("e15_cylCS5",colnames(tdf))])),
            rep(palette[3],length(colnames(tdf)[grepl("e15_marNI",colnames(tdf))])),
            rep(palette[4],length(colnames(tdf)[grepl("e15_nos7107",colnames(tdf))])),
            rep(palette[5],length(colnames(tdf)[grepl("e15_plkNIE",colnames(tdf))])),
            rep(palette[1],length(colnames(tdf)[grepl("e16_ana90",colnames(tdf))])),
            rep(palette[2],length(colnames(tdf)[grepl("e16_cylCS5",colnames(tdf))])),
            rep(palette[3],length(colnames(tdf)[grepl("e16_marNI",colnames(tdf))])),
            rep(palette[4],length(colnames(tdf)[grepl("e16_nos7107",colnames(tdf))])),
            rep(palette[5],length(colnames(tdf)[grepl("e16_plkNIE",colnames(tdf))]))
  ),
  Year=c(
    rep("green", length(colnames(tdf)[grepl("e15_",colnames(tdf), perl = T)])),
    rep("yellow", length(colnames(tdf)[grepl("e16_",colnames(tdf), perl = T)])))
)
#pdf("Buoy_15-16_SNP_querypaths_5sp_allSamples.pdf",18,6)
#pdf("Buoy_15-16_SNP_querypaths_5sp_allSamples_sorted.pdf",16,4)
#pdf("Efls_15-16_SNP_querypaths_5sp_allSamples.pdf",18,6)
#pdf("Efls_15-16_SNP_querypaths_5sp_allSamples_sorted.pdf",16,4)
par(mfrow=c(1,1))
heatmap3(tdf, Colv = NA, Rowv = T, scale = "none", cexRow = 0.9, cexCol = 0.25,
         ColSideColors = ColSideColors, ColAxisColors = 1,
         distfun = function(x) dist(x,method = 'euclidean'),
         legendfun=function() showLegend(legend=c("Enriched","Un-enriched"), col=c("red","blue"))
)
dev.off()

######################################################################################

######################################################################################
################## spls of abundance and pathways of each species   #############
rm(list=ls())
library(heatmap3)
library(spls)
options(stringsAsFactors = F)
folders=c("buoy_2015","buoy_2016","efls_2015","efls_2016")
pathwayfiles=file.path("pathfiles",list.files(path = "pathfiles",pattern = "locus_pathway_sort.txt"))
abfiles=file.path("pop_abundance",list.files("pop_abundance",pattern = "newDate.txt"))
temp=c()
for (i in 1:5){
  tx=read.delim(pathwayfiles[i],header = T, sep = "\t")
  temp=c(unique(tx$shortID),temp)
}
allpaths=unique(temp[!grepl("NoPath",temp, perl = T)])
corr.mat=as.data.frame(matrix("nd",length(allpaths),4*5+1))
corr.mat[,1]=sort(allpaths)
species=c("Ana","Cra","Mae","Nos","Pag")
sets=c("b15","b16","e15","e16")
colnames(corr.mat)=c("P21ID",paste(rep(species,4), rep(sets,each=5),sep = "-"))
colnum=1
# the two loops should be run one by one as K will vary
#for (f in 1:4){
f=3
# abundance of each population at each site/year
pop_abd=read.delim(abfiles[f],header = T, sep = "\t",row.names = 1,check.names = F)
pop_abd=pop_abd[-nrow(pop_abd),]
pop_abd=t(pop_abd[order(rownames(pop_abd)),])
# gene abundance files at each site/year
genefiles=file.path(folders[f],list.files(path = folders[f],pattern="gene_abundance_tpm"))
#for (g in 1:5){
g=1
colnum=(f-1)*5+g+1
tpm=read.delim(genefiles[g],header = T, sep = "\t", check.names = F)
colnames(tpm)[1]="X"
pathways=read.delim(pathwayfiles[g],header = T, sep = "\t")
pathways=pathways[pathways$pathname!= "NoPath",-3]
xmerge=merge(tpm, pathways, by.x = "X", by.y = "Locus")
xmerge=xmerge[,c(ncol(xmerge), 2:(ncol(xmerge)-2))]
xpaths=unique(xmerge$shortID)
if(g==3){
  xpaths=xpaths[!grepl("Akinete",xpaths)]
}
xpaths=xpaths[!grepl("TCA|Glycolysis",xpaths, perl = T)]
mat=matrix(0, nrow = length(xpaths), ncol = ncol(xmerge)-1)
rownames(mat)=xpaths
colnames(mat)=colnames(xmerge)[2:ncol(xmerge)]
# get the mean expression values of the genes in each pathways with a loop
for (j in 1:nrow(mat)){
  temp=apply(xmerge[grepl(xpaths[j], xmerge$shortID),2:ncol(xmerge)],2,mean)
  mat[j,]=temp
}
dates=unique(gsub("[ab]$","", colnames(xmerge)[2:ncol(xmerge)], perl = T))
submat=matrix(0, nrow = length(xpaths), ncol = length(dates))
rownames(submat)=xpaths
colnames(submat)=dates
for (d in 1:length(dates)) {
  tempi=as.matrix(mat[,grep(dates[d],colnames(mat))])
  class(tempi)
  if(ncol(tempi)>1){
    tempii=apply(tempi,1,mean)
    submat[,d]=tempii
  } else {
    submat[,d]=as.numeric(tempi)
  }
}
submat=t(submat[rowSums(submat>0)>2,])
# # standard and angular transformation of the data
# getArc<-function(x){
#   p=(x-min(x))/(max(x)-min(x))
#   p2=asin(sqrt(p))
#   return(p2)
# }
# submat=apply(submat,2, getArc)
# pop_abd=apply(pop_abd,1,getArc)
tcv=cv.spls(submat,pop_abd[,g], eta = seq(0.1,0.9,0.1), K=3:7)
tfit=spls(submat,pop_abd[,g], eta = tcv$eta.opt, K = tcv$K.opt)
print(tfit)
coef.f <- as.data.frame(coef(tfit))
coef.f$pid=rownames(coef.f)
colnames(coef.f)[1]="V1id"
ymerge=merge(corr.mat[,1:2],coef.f,by.x = "P21ID",by.y = "pid",all.x = T)
corr.mat[,colnum]=ymerge[,3]
#}
#}
write.table(corr.mat,"spls_popAbundance-pathway_corr-matrix.txt",
            quote = F,row.names = F, col.names = T, sep = "\t")
outheatmap=paste("plots/","buoy15-16_popAbundance.pdf",sep = "")
pdf("plots/efls15-16_popAbundance.pdf",10,10)
heatmap3(as.matrix(corr.mat[c(-4:-5,-17), 2:11]),Colv = T,Rowv = NA, 
         labRow = corr.mat[c(-4:-5,-17),1])
heatmap3(as.matrix(corr.mat[c(-4:-5,-17), 12:21]),Colv = T,Rowv = NA, 
         labRow = corr.mat[c(-4:-5,-17),1])
heatmap3(as.matrix(corr.mat[c(-4:-5,-17), 2:21]),Colv = T,Rowv = NA, 
         labRow = corr.mat[c(-4:-5,-17),1])
dev.off()
## make bubble plots of the query pathway association
rm(list=ls())
library(ggplot2)
library(dplyr)
cmat=read.delim("HC_spls_popAbundance-pathway_corr-matrix.txt",sep = "\t",check.names = FALSE, stringsAsFactors = F)
cmat[cmat=="-"]<-NA
ccmat=data.frame(cmat[,1:2],"Sp-year"=rep(colnames(cmat)[2],20))
colnames(ccmat)=c("Pathway","Coefficient","Species")
for (i in 3:ncol(cmat)) {
  tdf=data.frame(cmat[,c(1,i)],"Species"=rep(colnames(cmat)[i],20))
  colnames(tdf)=c("Pathway","Coefficient","Species")
  ccmat=rbind(ccmat,tdf)
}
ccmat[,2]=as.numeric(ccmat[,2])
ccmat$Sign=ifelse(ccmat[,2]>0,"apositive","negative")
# Most basic bubble plot
ccmat[is.na(ccmat[,4]), 3]<-NA
ccmat=ccmat[!grepl("Akinete",ccmat[,1]),]
ccmat$Pathway=gsub("Glycolysis","zGlycolysis",ccmat$Pathway)
ccmat$Pathway=gsub("TCA","zTCA",ccmat$Pathway)
ccmat[,3]=gsub("-","_",ccmat[,3],perl = T)
ccmat2=ccmat[with(ccmat, order(Species)),]
ccmat2=ccmat2[grepl("b15",ccmat2[,3],perl = TRUE),]
par(mfrow=c(1,2))
pdf("SPLS_Query-paths_vs_species-abundance_buoy_2015.pdf",6.5,3.5)
ccmat2 %>%
  #mutate(Paths = factor(Paths, Paths)) %>%
  ggplot(aes(x=Pathway, y=Species, size=-log10(abs(as.numeric(Coefficient))),color=Sign, na.rm = TRUE)) +
  geom_point(alpha=0.9,na.rm = TRUE) +
  scale_size(range = c(0, 8), name="Coefficient")+
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 12, face = "bold"))
dev.off()

######################################################################################
######################################################################################
################## spls of pop abundance and meta data   #############
rm(list=ls())
library(heatmap3)
library(spls)
options(stringsAsFactors = F)
folders=c("buoy_2015","buoy_2016","efls_2015","efls_2016")
abfiles=file.path("pop_abundance",list.files("pop_abundance",pattern = "newDate.txt"))
species=c("Ana","Cyl","Mic","Nos","Pla")
sites=c("b15","b16","e15","e16")

metadata=read.table("harsha_metadata.txt",header = T, sep = "\t")
corr.mat=as.data.frame(matrix("1",ncol(metadata)-3,4*5+1))
corr.mat[,1]=sort(colnames(metadata)[4:ncol(metadata)])
colnames(corr.mat)=c("Nutrient",sort(as.vector(outer(c("b15","b16","e15","e16"),c("ana","cyl","mic","nos","pla"), 
                                                     paste, sep = "_"))))
# the two loops should be run one by one as K will vary
#for (i= in 1:4){
i=1
# abundance of each population at each site/year
pop_abd=read.delim(abfiles[i],header = T, sep = "\t",row.names = 1,check.names = F)
pop_abd=pop_abd[-nrow(pop_abd),]
pop_abd=t(pop_abd[order(rownames(pop_abd)),])
# metadata at each site/year
metad=metadata[grepl(sites[i], metadata[,2]),c(-1,-3)]
rownames(metad)=metad[,1]
metad=metad[,-1]
#for (g in 1:5){
g=5
j=(i-1)*5+g+1
tcv=cv.spls(metad,pop_abd[1:nrow(metad),g], eta = seq(0.1,0.9,0.1), K=3:7)
tfit=spls(metad,pop_abd[1:nrow(metad),g], eta = tcv$eta.opt, K = tcv$K.opt)
print(tfit)
coef.f <- as.data.frame(coef(tfit))
coef.f$pid=rownames(coef.f)
colnames(coef.f)[1]="V1id"
ymerge=merge(corr.mat[,1:2],coef.f,by.x = "Nutrient",by.y = "pid",all.x = T)
corr.mat[,j]=ymerge[,3]
#}
#}
write.table(corr.mat,"spls_popAbundance-metaData_corr-matrix.txt",
            quote = F,row.names = F, col.names = T, sep = "\t")
save(corr.mat,file = "spls_popAbundance-metaData_corr-matrix.RData")

## make bubble plots of the pop abundance and metadata spls correlations
rm(list=ls())
library(ggplot2)
load("spls_popAbundance-metaData_corr-matrix.RData")
tdf=data.frame(Nutrient="xx",Coefficient=0,Species="Species")
for (i in 2:ncol(corr.mat)) {
  ttdf=as.data.frame(corr.mat[,c(1,i)])
  ttdf[,3]=rep(colnames(corr.mat)[i],nrow(corr.mat))
  colnames(ttdf)[2:3]=c("Coefficient",Species="Species")
  tdf=rbind(tdf,ttdf)
}
tdf$Sign=ifelse(tdf[,2]>0,"apositive","negative")
ccmat=tdf[-1,]
ccmat[,3]=gsub("([be]1[56])_([a-z]{3})","\\2_\\1",ccmat[,3], perl = T)
ccmat[,3]=gsub("(^[acnp])","\\U\\1",ccmat[,3], perl = T)
ccmat2=ccmat[grepl("_b1[56]",ccmat[,3],perl = TRUE),]
pdf("SPLS_metaData_vs_species-abundance_buoy_2years4.pdf",8,3)
ccmat2 %>%
  #mutate(Paths = factor(Paths, Paths)) %>%
  ggplot(aes(x=Nutrient, y=Species, size=-log10(abs(as.numeric(Coefficient))),color=Sign, na.rm = TRUE)) +
  geom_point(alpha=0.9,na.rm = TRUE) +
  scale_size(range = c(0, 8), name="Coefficient")+
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 12, face = "bold"))
dev.off()

######################################################################################
############ make pop abundane bar plots with query pathways spls regression #########
rm(list=ls())
library(ggplot2)
library(dplyr)
options(stringsAsFactors = F)
folders=c("buoy_2015","buoy_2016","efls_2015","efls_2016")
abfiles=file.path("pop_abundance",list.files("pop_abundance",pattern = "newDate.txt"))
species=c("Ana","Cyl","Mic","Nos","Pla")
sites=c("b15","b16","e15","e16")
# get the abundance of each species/site
abdf=data.frame(Sps="xxxx",Abd=0)
for (i in 1:2) {
  pop_abd=read.delim(abfiles[i],header = T, sep = "\t",row.names = 1,check.names = F)
  pop_abd=pop_abd[-nrow(pop_abd),]
  pop_tdf=data.frame(Sps=paste(sites[i],rownames(pop_abd),sep = "_"),Abd=apply(pop_abd, 1, mean))
  abdf=rbind(abdf,pop_tdf)
}
abdf=abdf[-1,]
abdf=abdf[order(abdf$Abd),]
rownames(abdf)=substr(abdf[,1],start = 1,stop = 7)

abdf15=abdf[grepl("b1[56]",abdf$Sps, perl = T),]
abdf15[,1]=paste(letters[1:5],substr(abdf15[,1],start = 1,stop = 7), sep = "_")
abdf15[,2]
barplot(as.matrix(abdf15[,2]))
p<-ggplot(data=abdf15, aes(x=Sps, y=Abd)) +
  geom_bar(stat="identity", fill="steelblue", width = 0.3)+
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 12, face = "bold"))
p
# Horizontal bar plot
p + coord_flip()

## make bubble plots of the pop abundance and metadata spls correlations
load("spls_popAbundance-metaData_corr-matrix.RData")
tdf=data.frame(Nutrient="xx",Coefficient=0,Species="Species")
for (i in 2:ncol(corr.mat)) {
  ttdf=as.data.frame(corr.mat[,c(1,i)])
  ttdf[,3]=rep(colnames(corr.mat)[i],nrow(corr.mat))
  colnames(ttdf)[2:3]=c("Coefficient",Species="Species")
  tdf=rbind(tdf,ttdf)
}
tdf$Sign=ifelse(tdf[,2]>0,"apositive","negative")
ccmat=tdf[-1,]
ccmat[,3]=gsub("([be]1[56])_([a-z]{3})","\\2_\\1",ccmat[,3], perl = T)
ccmat[,3]=gsub("(^[acnp])","\\U\\1",ccmat[,3], perl = T)
ccmat2=ccmat[grepl("_b1[56]",ccmat[,3],perl = TRUE),]
pdf("SPLS_metaData_vs_species-abundance_buoy_2years4.pdf",8,3)
ccmat2 %>%
  #mutate(Paths = factor(Paths, Paths)) %>%
  ggplot(aes(x=Nutrient, y=Species, size=-log10(abs(as.numeric(Coefficient))),color=Sign, na.rm = TRUE)) +
  geom_point(alpha=0.9,na.rm = TRUE) +
  scale_size(range = c(0, 8), name="Coefficient")+
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 12, face = "bold"))
dev.off()


######################################################################################
################## metaT data analysis  #########################################
## Analyses: heatmap of the query pathways with pop abundance and metadata
rm(list=ls())
library(heatmap3)
library("spls")
folders=c("buoy_2015","buoy_2016","efls_2015","efls_2016")
pathwayfiles=file.path("pathfiles",list.files(path = "pathfiles",pattern = "locus_pathway_sort.txt"))

for (f in 1:4){
  #f=1
  genefiles=file.path(folders[f],list.files(path = folders[f],pattern="gene_abundace_tpm"))
  for (g in 1:5){
    #g=1
    tpm=read.delim(genefiles[g],header = T, sep = "\t", check.names = F)
    colnames(tpm)[1]="X"
    pathways=read.delim(pathwayfiles[g],header = T, sep = "\t")
    pathways=pathways[pathways$pathname!= "NoPath",-3]
    xmerge=merge(tpm, pathways, by.x = "X", by.y = "Locus")
    xmerge=xmerge[,c(ncol(xmerge), 2:(ncol(xmerge)-2))]
    
    xpaths=unique(xmerge$shortID)
    if(g==3){
      xpaths=xpaths[!grepl("Akinete",xpaths)]
    }
    mat=matrix(0, nrow = length(xpaths), ncol = ncol(xmerge)-1)
    rownames(mat)=xpaths
    colnames(mat)=colnames(xmerge)[2:ncol(xmerge)]
    # get the mean expression values of the genes in each pathways with a loop
    for (j in 1:nrow(mat)){
      temp=apply(xmerge[grepl(xpaths[j], xmerge$shortID),2:ncol(xmerge)],2,mean)
      mat[j,]=temp
    }
    dates=unique(gsub("[ab]$","", colnames(xmerge)[2:ncol(xmerge)], perl = T))
    submat=matrix(0, nrow = length(xpaths), ncol = length(dates))
    rownames(submat)=xpaths
    colnames(submat)=dates
    for (d in 1:length(dates)) {
      tempi=as.matrix(mat[,grep(dates[d],colnames(mat))])
      class(tempi)
      if(ncol(tempi)>1){
        tempii=apply(tempi,1,mean)
        submat[,d]=tempii
      } else {
        submat[,d]=as.numeric(tempi)
      }
    }
    submat=submat[rowSums(submat>0)>2,]
    # # standard and angular transformation of the data
    # getArc<-function(x){
    #   p=(x-min(x))/(max(x)-min(x))
    #   p2=asin(sqrt(p))
    #   return(p2)
    # }
    # submat=apply(submat,2, getArc)
    outfile1=gsub("\\/","_",genefiles[g], perl = TRUE)
    # file name with no transformation
    outfile=paste("plots/",gsub("gene_abundace_tpm[.]txt","pathway_NT-Colv_heatmap.pdf",outfile1),sep = "")
    # file name of angular transformation
    #outfile=gsub("gene_abundace_tpm[.]txt","pathway_arcT_heatmap.pdf",outfile1)
    pdf(outfile,18,18)
    heatmap3(submat,cexRow = 1.6,cexCol = 1.6, balanceColor = T,Colv = NA)
    dev.off()
  }
}

######################################################################################
################# reorganize the geneID-query pathway files ##########################
# pathway files for each genome of five genomes
rm(list=ls())
pfiles=list.files(pattern = "locus.pid.pathway")
for (i in 1:5){
  #i=1
  allpaths=read.table(pfiles[i],header = F,sep = "\t")
  colnames(allpaths)=c("ProtID","Locus","pathname")
  allpaths$pathname=gsub(pattern = "Photosynthesis\\|Calvin", "Calvin\\|Calvin", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "Photosynthesis\\|PSB", "PBS\\|PBS", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "PentosePipathname\\|PPP", "PPP\\|PPP", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "CO2CCM", "CCM", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "Vesible", "Vesicle", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "phoQ", "Phosphorus\\|Pi\\phoQ", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "EDpathname", "Glycolysis", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "AAPEP", "AAT", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "IonDrug", "TEVit", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "Nitrogen", "N", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "Nfix", "NF", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "Osmoprotect", "Osmosis", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "PS-I", "PSI", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "PS-II", "PSII", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "Stress", "OSR", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "Toxins", "Toxin", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "UnsatFA", "PUFA", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "PS-II", "PSII", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "PS-II", "PSII", allpaths$pathname)
  allpaths$pathname=gsub(pattern = "Phosphorus", "P", allpaths$pathname)
  allpaths$shortID=ifelse(allpaths$pathname=="NoPath","NoPath",gsub("\\|.+","",allpaths$pathname,perl = TRUE))
  outfile=gsub(".locus.pid.pathway","_pid_locus_pathway_sort.txt",pfiles[i])
  write.table(allpaths,outfile,sep = "\t",quote = F,row.names = F,col.names = T)
}

###########################################################################################
######## Calculat relative metric of pathways among all cyanobacterial genomes ############
rm(list = ls())
par()$oma
op<-par(no.readonly=TRUE)

##### read in dataset for 19 query pathways
x=read.table("all.paths.sta.2", header=T,sep="\t")
xs=x[,c(-2:-5,-9, -16:-24)]
colnames(xs)=c("Strain","Pathway","NP","NC","NPC","LP","FLP","LC","FLC","PC")
bi=read.csv("Blooming_indicence.csv",header = T)
#bi=bi[,c(1,5)]
#colnames(bi)=colnames(x)[c(1,4)]
bi=x[1:113,c(1,4)]
pathlist=unique(xs[,2])
pathlist=pathlist[c(4,19,10,11,1,13,7,8,2,17,18,5,9, 6,12,3,16,15,14)]
pathlist_corr=c("Vesicle","Toxin","Osmosis","PUFA","AAT","PBS","MAA","NF","CCM",
                "Sugar","Sulfur","MetalR","N","TEVit","P","DrugR","OSR","PSI",
                "PSII")
par(oma=c(0,1,0,0))

##### this portion of code standardize and then average all eight metrics for each pathway
##### mean relative metric=(metric-minimum)/(max-min)
mmeans=bi # set up an dataframe for mmeans
mmeans[,3:26]=0
colnames(mmeans)[3:26]=c(pathlist_corr,"Calvin","ETC","Glycolysis","PPP","TCA")
relNP=relNC=relNPC=relPC=mmedian=mmeans2=mmeans1=mmeans

#pdf("query_19_pathways_abs-metric1.pdf",7,7)
#pdf("query_19_pathways_rel-std-metric.pdf",7,7)
for (ii in c(1:19)) {
  #ii=2
  tempath=xs[which(xs[,2]==pathlist[ii]),]
  tempathm=merge(bi,tempath, by.x="Strain",by.y="Strain")
  maxs=apply(tempathm[,4:11],2,max) #comment out for abs abundance
  mins=apply(tempathm[,4:11],2,min) #comment out for abs abundance
  range=maxs-mins                   #comment out for abs abundance
  #maxs[which(maxs==0)]=1 #comment out for abs abundance
  tempathm[,4:11]=t((t(tempathm[,4:11])-mins)/range) #comment out for abs abundance
  tempathm$colors=ifelse(tempathm$BlmFreq>0,"B","NB")
  tempathm=tempathm[,colSums(is.na(tempathm))!=nrow(tempathm)]
  # par(mfrow=c(2,4))
  # for (jj in 4:(ncol(tempathm)-1)) {
  #   boxplot(tempathm[,jj]~tempathm$colors,col=c("blue","grey"), main=pathlist_corr[ii],
  #           ylab=paste("Relative",colnames(tempathm[jj]),sep = " "),
  #           cex.names=1.2,cex.axis=1.3,cex.lab=1.4,outline=FALSE
  #   )
  # }
  #dev.off()
  tempathm=tempathm[,colSums(is.na(tempathm))!=nrow(tempathm)]
  mmeans[,(ii+2)]=apply(tempathm[,4:(ncol(tempathm)-1)],1,mean)
  mmeans1[,(ii+2)]=apply(tempathm[,4:6],1,mean) # starting at 22
  mmeans2[,(ii+2)]=apply(tempathm[,c(4:6,(ncol(tempathm)-1))],1,mean)
  #mmedian[,(ii+2)]=apply(tempathm[,4:11],1,median)
  relNP[,(ii+2)]=tempathm[,4]
  relNC[,(ii+2)]=tempathm[,5]
  relNPC[,(ii+2)]=tempathm[,6]
  #relPC[,(ii+2)]=tempathm[,11]
}

##### read in core pathways

corep=list.files(pattern = "core_path_")
cdf=read.table(corep[1],header = T,sep = "\t")
cdf=cdf[,1:9]
cdf[,5:9]=0
colnames(cdf)[5:9]=c("Calvin","ETC","Glycolysis","PPP","TCA")

for(ii in 1:5){
  temp=read.table(corep[ii],header = T,sep = "\t")
  cdf[,ii+4]=rowSums(temp[,5:ncol(temp)])
}
cdf2=cdf
cdf2[,5:9]=0
for (i in 5:9) {
  temp=(cdf[,i]-min(cdf[,i]))/(max(cdf[,i])-min(cdf[,i]))
  cdf2[,i]=temp
}
mmeans_merge=merge(mmeans[,-22:-27],cdf2[,c(1,5:9)], by.x = "Strain",by.y = "Strain")
par(op)

###########################################################################################

###########################################################################################
##################  reorganize Meta data from Harsha Lake #####################
# sampleID are set to: m-dd-byy
rm(list=ls())
xn=read.delim("Lu_nutrients.csv",sep = ",", header = TRUE)
# add consistent dates to nutrient data
xn$Date=as.character(xn$Date)
for (i in 1:nrow(xn)){
  #i=104  
  if(xn$Date[i] == ""){
    tempv=unlist(strsplit(as.character(xn[i,1]),"\\/"))
    tv1=gsub("0","",unlist(strsplit(as.character(xn[i,1]),"\\/"))[1])
    tv2=as.numeric(unlist(strsplit(as.character(xn[i,1]),"\\/"))[2])
    tv2=ifelse(tv2 <10,paste("0",tv2,sep = ""),tv2)
    tv=paste(tv1,tv2,sep = "-")
    tv=ifelse(grepl("BUOY", xn$Site[i]),paste(tv,"b",sep = "-"),paste(tv,"e",sep = "-"))
    tv=ifelse(grepl("2015", xn$SampleDate[i]),paste(tv,"15",sep = ""),paste(tv,"16",sep = ""))
    xn$Date[i]=tv
  } else {
    tv=ifelse(grepl("BUOY", xn$Site[i]),paste(xn$Date[i],"b",sep = "-"),paste(xn$Date[i],"e",sep = "-"))
    tv=ifelse(grepl("2015", xn$SampleDate[i]),paste(tv,"15",sep = ""),paste(tv,"16",sep = ""))
    xn$Date[i]=tv
  }
}
# reorganize data into 10 columns of dates, sites,and parameters
nutrient=as.data.frame(matrix(0,nrow = length(unique(xn$Date)),ncol = 10), stringsAsFactors=FALSE)
colnames(nutrient)=c("SampleDate","Date","Site","TN","TNH4","TNO2","TNO2-3","TOC","TP","TRP")
nutr=data.frame(xfators=c("TN","TNH4","TNO2","TNO2-3","TOC","TP","TRP"),values=rep("na",7))

for (i in 1:nrow(nutrient)){
  tt=as.data.frame(xn[grepl(unique(xn$Date)[i],xn$Date),c(1,ncol(xn),3, 8:9)],stringsAsFactors=FALSE)
  tt2=apply(tt,2, as.character)
  nutrient[i,1:3]=tt2[1,1:3]
  tt2=as.numeric(t(merge(tt[,4:5],nutr,by.x = "Parameter",by.y = "xfators",all.y=TRUE))[2,])
  nutrient[,4:10]=tt2
}
getRpl<-function(x){
  x[which(x==629)]<-0
  return(x)
}
nutrient[,4:10]<-apply(nutrient[,4:10],2,getRpl)
apply(nutrient[,4:10],2,max)
getMax<-function(x){
  x[which(x==0)]<-55
  return(x)
}
nutrient[,4:10]<-apply(nutrient[,4:10],2,getMax)
write.table(nutrient,"harsha_metadata.txt",quote = F,sep = "\t", row.names = F)
##################################################################################
################## changes sample IDs for species abundance ######################
# sampleID are set to: m-dd-byy
# abundance of all four sites are stored in a list, abdlist, and written to files
rm(list=ls())
xabfiles=list.files(pattern = "Species_abundance_")
abdlist=list()
i=1  # buoy2015
xabf=read.csv(xabfiles[i],header = T)
rownames(xabf)=xabf[,1]
xabf=xabf[,-1]
dates=sort(unique(gsub("[Xab]","",colnames(xabf))))
mat=matrix(0,nrow = nrow(xabf),ncol = length(dates))
colnames(mat)=dates
rownames(mat)=rownames(xabf)
for (k in 1:length(dates)) {
  #k=1
  tempi=as.matrix(xabf[,grep(dates[k],colnames(mat))])
  class(tempi)
  if(ncol(tempi)>1){
    tempii=apply(tempi,1,mean)
    mat[,k]=tempii
  } else {
    mat[,k]=as.numeric(tempi)
  }
}
dates=gsub("2015","-b15", dates)
dates=gsub("^06","6-", dates, perl = T)
dates=gsub("^07","7-", dates, perl = T)
colnames(mat)=dates
abdlist[["buoy2015"]]<-mat
write.table(mat,"species_abundance_buoy2015_newDate.txt",quote = F,sep = "\t",
            row.names = T, col.names = NA)
i=2 # buoy2016
xabf=read.csv(xabfiles[i],header = T)
rownames(xabf)=xabf[,1]
xabf=xabf[,-1]
dates=sort(unique(gsub("[A-Za-z.]","",colnames(xabf))))
mat=matrix(0,nrow = nrow(xabf),ncol = length(dates))
colnames(mat)=dates
rownames(mat)=rownames(xabf)
for (k in 1:length(dates)) {
  #k=1
  tempi=as.matrix(xabf[,grep(dates[k],colnames(mat))])
  class(tempi)
  if(ncol(tempi)>1){
    tempii=apply(tempi,1,mean)
    mat[,k]=tempii
  } else {
    mat[,k]=as.numeric(tempi)
  }
}
dates=gsub("_","-", dates)
dates=paste(dates, "-b16", sep = "")
colnames(mat)=dates
abdlist[["buoy2016"]]<-mat
write.table(mat,"species_abundance_buoy2016_newDate.txt",quote = F,sep = "\t",
            row.names = T, col.names = NA)

i=3 # efls2015
xabf=read.csv(xabfiles[i],header = T)
rownames(xabf)=xabf[,1]
xabf=xabf[,-1]
dates=sort(unique(gsub("[A-Za-z]","",colnames(xabf))))
mat=matrix(0,nrow = nrow(xabf),ncol = length(dates))
colnames(mat)=dates
rownames(mat)=rownames(xabf)
for (k in 1:length(dates)) {
  tempi=as.matrix(xabf[,grep(dates[k],colnames(mat))])
  class(tempi)
  if(ncol(tempi)>1){
    tempii=apply(tempi,1,mean)
    mat[,k]=tempii
  } else {
    mat[,k]=as.numeric(tempi)
  }
}
dates=gsub("[.]","-", dates)
dates=paste(dates, "-e15", sep = "")
colnames(mat)=dates
abdlist[["efls2015"]]<-mat
write.table(mat,"species_abundance_efls2015_newDate.txt",quote = F,sep = "\t",
            row.names = T, col.names = NA)
i=4 # efls2015
xabf=read.csv(xabfiles[i],header = T)
rownames(xabf)=xabf[,1]
xabf=xabf[,-1]
dates=sort(unique(gsub("[A-Za-z]","",colnames(xabf))))
mat=matrix(0,nrow = nrow(xabf),ncol = length(dates))
colnames(mat)=dates
rownames(mat)=rownames(xabf)
for (k in 1:length(dates)) {
  #k=1
  tempi=as.matrix(xabf[,grep(dates[k],colnames(mat))])
  class(tempi)
  if(ncol(tempi)>1){
    tempii=apply(tempi,1,mean)
    mat[,k]=tempii
  } else {
    mat[,k]=as.numeric(tempi)
  }
}
dates=gsub("[.]","-", dates)
dates=paste(dates, "-e16", sep = "")
colnames(mat)=dates
abdlist[["efls2016"]]<-mat
write.table(mat,"species_abundance_efls2016_newDate.txt",quote = F,sep = "\t",
            row.names = T, col.names = NA)
