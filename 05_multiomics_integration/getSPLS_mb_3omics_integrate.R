## ----global_options, include=FALSE------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
library(knitr)
knitr::opts_chunk$set(dpi = 100, 
                      echo= TRUE, 
                      warning=FALSE, 
                      message=FALSE, 
                      fig.show=TRUE, 
                      fig.keep = 'all', 
                      fig.align = "center",
                      out.width = "70%")
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library
set.seed(123) # for reproducibility, remove for normal use
rm(list = ls())
load("spls3omics_integrated.RData")
lapply(X, dim)
# ## ---- fig.show = "hold", out.width = "33%", fig.height = 6, fig.cap = "FIGURE 1: Circle Correlation Plots for pairwise PLS models on the gastrulation data. Only displays the top 25 features for each dimension, subsetting by those with a correlation above 0.5."----
# # select arbitrary values of features to keep
# list.keepX = c(300, 300)
# list.keepY = c(300, 300)
# 
# # generate three pairwise PLS models
# pls1 <- spls(X[["rna"]], X[["methylation"]], keepX = list.keepX, keepY = list.keepY)
# pls2 <- spls(X[["rna"]], X[["accessibility"]], keepX = list.keepX, keepY = list.keepY)
# pls3 <- spls(X[["methylation"]], X[["accessibility"]], keepX = list.keepX, keepY = list.keepY)
# 
# # plot features of first PLS
# plotVar(pls1, cutoff = 0.5, title = "(a) RNA vs Methylation",
#         legend = c("RNA", "Methylation"),
#         var.names = FALSE, style = 'graphics',
#         pch = c(16, 17), cex = c(2,1),
#         col = c('darkorchid', 'lightgreen'))
# 
# # plot features of second PLS
# plotVar(pls2, cutoff = 0.5, title = "(b) RNA vs Accessibility",
#         legend = c("RNA", "Accessibility"),
#         var.names = FALSE, style = 'graphics',
#         pch = c(16, 17), cex = c(1.5,1),
#         col = c('darkorchid', 'lightgreen'))
# 
# # plot features of third PLS
# plotVar(pls3, cutoff = 0.5, title = "(c) Methylation vs Accessibility",
#         legend = c("Methylation", "Accessibility"),
#         var.names = FALSE, style = 'graphics',
#         pch = c(16, 17), cex = c(1.25,1),
#         col = c('darkorchid', 'lightgreen'))
# 
# ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# cor(pls1$variates$X, pls1$variates$Y) # calculate correlation of RNA and methylation
# cor(pls2$variates$X, pls2$variates$Y) # calculate correlation of RNA and accessibility
# cor(pls3$variates$X, pls3$variates$Y) # calculate correlation of methylation and accessibility

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
design = matrix(1, ncol = length(X), nrow = length(X), # for square matrix filled with 0.5s
                dimnames = list(names(X), names(X)))
diag(design) = 0 # set diagonal to 0s
choice.ncomp <- 2 
choice.keepX <- list(mG = rep(2000, choice.ncomp), # 50 features per each component per dataframe
                     mT = rep(2000, choice.ncomp), 
                     mM = rep(2000, choice.ncomp))

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
final.mbspls.model = block.spls(X, indY = 3,  # generate final model using "tuned" parameters
                                ncomp = choice.ncomp, 
                                keepX = choice.keepX,
                                design = design)
res = block.spls(X, indY = 3,  # generate final model using "tuned" parameters
                 ncomp = choice.ncomp, 
                 keepX = choice.keepX,
                 design = design)

# Extract selected variable names for component 1
sel <- selectVar(res, comp = 1)

selected_mG <- sel$mG$name
selected_mM <- sel$mM$name
selected_mT <- sel$mT$name
mG_top <- mG[, selected_mG, drop = FALSE]
mM_top <- mM[, selected_mM, drop = FALSE]
mT_top <- mT[, selected_mT, drop = FALSE]
cor_mG_mM <- cor(mG_top, mM_top)
cor_mG_mT <- cor(mG_top, mT_top)
cor_mM_mT <- cor(mM_top, mT_top)


library(tidyr)
library(dplyr)

# mG vs mM
df_mG_mM <- as.data.frame(as.table(cor_mG_mM)) %>%
  rename(var1 = Var1, var2 = Var2, correlation = Freq) %>%
  mutate(block1 = "mG", block2 = "mM")

# mG vs mT
df_mG_mT <- as.data.frame(as.table(cor_mG_mT)) %>%
  rename(var1 = Var1, var2 = Var2, correlation = Freq) %>%
  mutate(block1 = "mG", block2 = "mT")

# mM vs mT
df_mM_mT <- as.data.frame(as.table(cor_mM_mT)) %>%
  rename(var1 = Var1, var2 = Var2, correlation = Freq) %>%
  mutate(block1 = "mM", block2 = "mT")
# Combine all
all_correlations <- bind_rows(df_mG_mM, df_mG_mT, df_mM_mT)

# Optional: filter by strength
all_correlations_filtered <- all_correlations %>% filter(abs(correlation) > 0.7)

# Export
write.csv(all_correlations_filtered, "pairwise_correlations_top_variables_top2000.csv", row.names = FALSE)

top_vars <- data.frame(
  variable = c(selected_mG, selected_mM, selected_mT),
  block = c(rep("mG", length(selected_mG)),
            rep("mM", length(selected_mM)),
            rep("mT", length(selected_mT)))
)

write.csv(top_vars, "top_variable_ids_top2000.csv", row.names = FALSE)


## ---- out.width = "90%", fig.cap = "FIGURE 2: Sample plot for sPLS2 performed on the gastrulation data. Samples are projected into the space spanned by the components yielded from the RNA dataset."----
plotIndiv(final.mbspls.model, ind.names = FALSE,
          group = as.factor(metadata$time), 
          pch = as.factor(metadata$site),
          col.per.group = color.mixo(1:6), 
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Site',
          blocks = 1)

## ---- echo = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(425)

samples <- sample(1:dim(X$mG)[1], 17)

X1 <- X$mG[samples, ] 
X2 <- X$mT[samples, ] 
X3 <- X$mM[samples, ] 
X.arrow <- list(mG = X1, mT = X2, mM = X3)

final.mbspls.model.arrow = block.spls(X.arrow, indY = 3, 
                                      ncomp = choice.ncomp, 
                                      keepX = choice.keepX,
                                      design = design)

## ---- out.width = "90%", fig.cap = "FIGURE 3: Arrow plot from the sPLS2 performed on the gastrulation data. The star indicates the location of the centroid of a sample across all the three datsets. The tip of each arrow shows the location of that same sample in the space spanned by the components associated to a specific dataset."----
pdf("spls_mb_mGMT_arrowplot6x6.pdf",6,6)
symbols <- list(mG = 1, mT = 6, mM = 10)
Factor = as.factor(metadata$time[samples])
plotArrow(final.mbspls.model.arrow, ind.names = FALSE, arrow.size = 0.15,
          arrow.length = 0.08, arrow.alpha = 0.5,
          group = Factor,box=TRUE, axes.box="all",
          pch = symbols, pch.size = 1)
dev.off()

## ---- out.width = "90%", fig.cap = "FIGURE 4: Correlation circle plot from the sPLS2 performed on the gastrulation data"----------------------------------------------------
pdf("spls_mb_mGMT_3mGMTvaplot4x2.pdf",4,2.25)
plotVar(final.mbspls.model, var.names = FALSE,
        legend = TRUE, cutoff = 0.5, cex = c(0.4, 0.4, 0.4),
        pch = c(0,1,2), )
dev.off()
## ---- echo = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------
X.2 = list(mG = X$mG, 
           mT = X$mT)
list.keepX = list(mG = rep(100, 2), mT = rep(100,2))
list.keepY = c(rep(100, 2))

final.mbspls.model.circos = block.spls(X = X.2, Y = X$mM,
                                       ncomp = 2, keepX = list.keepX,
                                       keepY = list.keepY, design = 'full')

## ---- out.width = "90%", fig.cap = "FIGURE 5: Circos plot from multiblock sPLS performed on the gastrulation data The plot represents the correlations greater than 0.8 between variables of different types, represented on the side quadrants", results = FALSE----
library(circlize)
pdf("mbSPLS_mGMT_circos_bottomHalf.pdf", width = 6, height = 10)
# Set plot region to bottom half (xmin=0, xmax=1, ymin=0, ymax=0.5)
par(fig = c(1, 1, 1, 1),    # bottom half
    mar = c(1, 1, 1, 1),      # small margins
    new = FALSE)

par(fig = c(0, 1, 0, 0.5),           # bottom half of the PDF page
    mar  = c(5.1, 4.1, 4.1, 2.1),    # ← normal/standard R margins
    new  = FALSE)

# Start a fresh plot
selected <- selectVar(final.mbspls.model.circos, comp = 1)

n_genes <- length(selected$X$name)    # X variables (e.g., genes)
n_metabs <- length(selected$Y$name)  # Y variables (e.g., metabolites)

# Safe gap setup
gap_genes <- if (n_genes > 1) rep(2, n_genes - 1) else numeric(0)
gap_metabs <- if (n_metabs > 1) rep(2, n_metabs - 1) else numeric(0)
# Apply to circos
circos.clear()
circos.par(gap.after = c(gap_genes, 10, gap_metabs, 10))

# Use mixOmics function to get input data and layout
circosPlot(final.mbspls.model.circos,
           group = metadata$time,
           legend = FALSE,
           size.labels = 0,     # suppress labels
           cutoff = 0.9,
           comp = 1,
           linkWidth = 0.05,
           size.variables = 0.2,
           Y.name = 'mM'
)
# Then manually add labels again
# You'll need the variable names and positions from your model
labels <- final.mbspls.model.circos$X$names
sectors <- final.mbspls.model.circos$X$var.names

for (i in seq_along(labels)) {
  circos.text(
    sector.index = sectors[i],
    x = 0,
    y = 1.5,
    labels = labels[i],
    facing = "clockwise",
    niceFacing = TRUE,
    cex = 0.6,
    adj = c(0, 0.5)
  )
}
dev.off()

### 
pdf("circos_bottom_half.pdf", width = 6, height = 8)
# Set plot region to bottom half (xmin=0, xmax=1, ymin=0, ymax=0.5)
par(fig = c(0, 1, 0, 0.5),    # bottom half
    mar = c(1, 1, 1, 1),      # small margins
    new = FALSE)
circosPlot(final.mbspls.model.circos, group = metadata$time, legend=T, size.labels=0.5,
           cutoff = 0.9, comp = 1, linkWidth=0.05, size.variables=0.2,
           Y.name = 'mM')
dev.off()

##------------functional inference with mM, MG, and MT --------------------
# mM-enriched KEGG pathways
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

# read in paired mM, mG, and mT data
cdf=read.csv("top_variable_ids_top2000.csv", header = T)

module_dfList <- list()
module_comps = cdf[cdf$block=='mM',1]

module_df=mm0[mm0$ID %in% module_comps, c(2:3, 12,14)]
rownames(module_df)=module_df$ID
merge_module_df = merge(module_df, mm2, by = "row.names")
#dev.off()
#pdf("spls_paired_mMmGmT_mMcomps_keggPathway_enriched.pdf", 8,8)
library(pheatmap)
pheatmap(log10(merge_module_df[,-1:-5] + 1), cluster_cols=FALSE, fontsize = 7,
     show_rownames = FALSE, main="log10value")

merge_module_df2=merge_module_df[merge_module_df$Metabolite != '-',]
rownames(merge_module_df2)=merge_module_df2$Metabolite
#merge_module_df2=merge_module_df2[c(3, 4,6,7,8,19,34, 59,65, 68:70, 74, 76, 96),]
#pdf("mbSPLS_mM_namedMets_heatmap.pdf",5,2.25)
pheatmap(log10(merge_module_df2[,-1:-5] + 1), fontsize_col = 4, fontsize_row = 4,
         cluster_cols=FALSE, fontsize = 4, show_rownames = TRUE)
#dev.off()

search_terms=merge_module_df2$KEGG.Compound.ID[merge_module_df2$KEGG.Compound.ID != '-']
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
         cluster_cols=FALSE, fontsize = 6, show_rownames = TRUE, main = 'mG_genes_tpm')
pheatmap(log10(merge_gtgenes[,grepl("^mT",colnames(merge_gtgenes))] +1),
         cluster_cols=FALSE, fontsize = 6, show_rownames = TRUE, main = 'mT_genes_tpm')
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

write.csv(r_total,'mbSPLS_comp_rxn_kopaths.csv', quote = F, row.names = F)
write.csv(c_total,'mbSPLS_comp_cNo_kopaths.csv', quote = F, row.names = F)

dev.off()

##--------- enrich kegg pathways with mG genes associated with mM ----------
#--  modules of mG and mT genes associated with NO3-N
rm(list = ls())
cdf=read.csv("top_variable_ids_top2000.csv", header = T)

pdf("spls_paired_mGMP_mTgene_keggPathway_enriched.pdf",8,8)

#module_genes = cdf[cdf$block=="mG",1]
module_genes = cdf[cdf$block=="mT",1]
module_genes=gsub("mT","UniRef90",module_genes)
search_list=unique(module_genes)

xdf=read.csv("Table5-2_mGmT_tpm_bytimeDate.tsv",sep = "\t",header = T)
mg_df=xdf[xdf$Gene %in% search_list,1:18]
rownames(mg_df)=mg_df[,1]
mg_df=mg_df[,-1]
mg_df=mg_df[rowSums(is.na(mg_df)) <5,]
log_mg_df=log10(mg_df[]+1)
library(pheatmap)
pheatmap(log_mg_df, cluster_rows = TRUE, cluster_cols=FALSE, fontsize_col = 8, 
         show_rownames = FALSE, main = 'mT genes tpm')
#dev.off()

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

df <- read.csv("Table6_mGmTmM_CnumKOpathway-full.tsv",sep = "\t",header = T)
kopath_df=df[,1:11]
# Step 1: Filter out rows without ko_pathway
filtered_df <- kopath_df %>%
  filter(!is.na(ko_pathway))

# Step 2: Create universe and hits
universe_genes <- unique(filtered_df$qseqid)
input_genes <- intersect(search_list, universe_genes)

# Step 3: Map genes to KO pathways
gene_to_pathway <- filtered_df %>%
  select(qseqid, ko_pathway) %>%
  distinct()

# Step 4: For each pathway, count hits and total genes
pathway_counts <- gene_to_pathway %>%
  group_by(ko_pathway) %>%
  summarise(
    K = n(),  # total genes annotated to this pathway
    k = sum(qseqid %in% input_genes),  # genes in input list
    genes = paste(unique(qseqid[qseqid %in% input_genes]), collapse = ",")
  ) %>%
  ungroup()

# Step 5: Hypergeometric test
N <- length(universe_genes)     # total genes with a ko_pathway
n <- length(input_genes)        # number of genes in search_list

pathway_stats <- pathway_counts %>%
  mutate(
    pval = phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE),
    padj = p.adjust(pval, method = "BH"),
    sig_flag = padj < 0.05
  )

# Step 6: Clean ko_pathway names for better plotting
# Step 6: Clean ko_pathway names for better plotting (no capitalization)
pathway_stats <- pathway_stats %>%
  filter(K > 2) %>%  # optionally filter low-count pathways
  filter(!str_detect(ko_pathway, "Superkiller protein|^Uncharacterized")) %>%
  arrange(desc(K)) %>%
  mutate(
    ko_pathway = factor(ko_pathway, levels = ko_pathway)  # preserve original names and order
  )
pathway_stats=pathway_stats[pathway_stats$k>0,]
#write.csv(pathway_stats, "mbSPLS_mT_gene_kopaths.csv", row.names = F, quote = F)
pathway_stats=pathway_stats[pathway_stats$K>1,]
pathway_stats=pathway_stats[!grepl("^ko03010",pathway_stats$ko_pathway),]
# Step 7: Plot
ggplot(pathway_stats, aes(x = ko_pathway)) +
  geom_bar(aes(y = K), stat = "identity", fill = "gray80", width = 0.9) +  # background bar
  geom_bar(
    aes(y = k, fill = sig_flag),
    stat = "identity",
    width = 0.54,
    position = position_nudge(x = 0)
  ) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "coral")) +
  coord_flip() +
  labs(
    x = "Pathway",
    y = "Number of Genes",
    #title = "Pathway Enrichment by Hypergeometric Test",
    #subtitle = "Gray: total genes; Coral: enriched (FDR < 0.05); Blue: not enriched",
    fill = "Significant",
    caption = "Overlay bars: genes from input list"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(size = 4))

dev.off()
################### 
### draw a venn diagram of ko pathways by mG, mM, and mT
mG=read.csv("mbSPLS_mG_gene_kopaths2.csv", header = T)
mG=gsub(" .*$","", mG$ko_pathway)

mT=read.csv("mbSPLS_mT_gene_kopaths2.csv", header = T)
mT=gsub(" .*$","", mT$ko_pathway)

mM=read.csv("mbSPLS_comp_cNo_kopaths2.csv", header = T)
mM=gsub(" .*$","", mM$ko_pathway)

library(eulerr)

venn_list <- list(
  mG = unique(mG),
  mT = unique(mT),
  mM = unique(mM)
)

# Fit the Euler diagram
fit <- euler(venn_list)

# Extract counts and compute percentages
counts <- round(fit$original.values)
total <- sum(counts)
percentages <- round(100 * counts / total, 1)

# Create custom labels
custom_labels <- paste0(counts, " (", percentages, "%)")

# Plot with custom labels
pdf("mbSPLS_mGMT_euler_counts_percents.pdf", width = 3, height = 2.25)
plot(fit,
     fills = c("coral", "steelblue", "lightgreen"),
     labels = list(font = 1),
     quantities = list(labels = custom_labels, cex = 0.6))
dev.off()



library(ggvenn)
pdf("mbSPLS_mGMT_venn.pdf", 2,4)
venn_list <- list(
  mG = unique(mG),
  mT = unique(mT),
  mM = unique(mM)
)
ggvenn(
  venn_list,
  fill_color = c("coral", "steelblue", "lightgreen"),
  stroke_size = 0.2,
  set_name_size = 3,
  text_size = 1.5, fill_alpha = 0.6
)
dev.off()
#############################################################################
############ below are the codes for processing the input files #########
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('mixomics_Data/nmt_data_processed.RData') # load the gastrulation data
X1 <- data$rna # select three of the five dataframes to explore
X2 <- data$met_genebody
X3 <- data$acc_genebody
X <- list(rna = X1, methylation = X2, accessibility = X3) # compile these into a single X object
lapply(X, dim) # check dimensions

##---- preparing my own data sets
# metadata, metaG, mT, and mM
##--prepare metadata
metadata = read.delim('water_chemistry.csv', header = T, row.names = 1,
                      sep = ',', check.names = F)
dim(metadata) # confirm the dimension of data
rownames(metadata) = gsub('mG19_', 'S', rownames(metadata))
rnames = rownames(metadata)
# as mG misses one sample 11.3.1, so this is removed from all datasets
metadata = metadata[-13,]
metadata <- data.frame(lapply(metadata, as.numeric))
metadata = as.matrix(metadata)
rownames(metadata) = rnames[-13]

embryo = rownames(metadata)
stage = paste("Site", rep(1:3, 6),sep = "")[-13]
lineage = gsub("\\.[1-3]$","",rownames(metadata),perl = T)
lineage = gsub("^S","D",lineage,perl = T)
lineage = lineage
cell_metadata2=list(embryo = embryo, stage = stage, lineage= lineage)

##----- X1=mG: metaGenome preparation
# sample 11.3.1 is missing
xdf=read.delim("summary_table_5_across-time.tsv",sep = '\t')
rownames(xdf) =xdf[[1]]
xdf=xdf[, seq(3,71,4)]
colnames(xdf)=gsub("X([679])","X0\\1", colnames(xdf))
colnames(xdf)=gsub("^X","S", colnames(xdf))
colnames(xdf)=gsub("_abundance","", colnames(xdf))
df = as.data.frame(t(xdf))
df = df[rowSums(is.na(df)) < ncol(df),]
df = df[order(rownames(df)),]
rnames = rownames(df)
df[is.na(df)] <- 0
df <- data.frame(lapply(df, as.numeric))
rownames(df) = rnames
mG = as.matrix(df)

##-----X2=mT: metatranscriptome preparation
# sample 11.3.1 is missing
rm(xdf)
xdf=read.delim("summary_table_5_across-time.tsv",sep = '\t')
rownames(xdf) = gsub("UniRef90","mT",xdf[[1]])
xdf=xdf[, seq(5,73,4)]
colnames(xdf)=gsub("X([679])","X0\\1", colnames(xdf))
colnames(xdf)=gsub("^X","S", colnames(xdf))
colnames(xdf)=gsub("_TPM","", colnames(xdf))
# remove NAs from columns
xdf <- xdf[, colSums(is.na(xdf)) < nrow(xdf)]
xdf[is.na(xdf)] <- 0
df = as.data.frame(t(xdf))
df = df[order(rownames(df)),]
df = df[-13,]
rnames = rownames(df)
df <- data.frame(lapply(df, as.numeric))
rownames(df) = rnames
mT=as.matrix(df)

##----- X3=mM, metabolome preparation
# sample 11.3.1 is missing
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
rownames(df)=gsub('X','S', rownames(df))
rownames(df)=gsub('_','\\.', rownames(df))
df = df[-13,]
df_numeric <- data.frame(lapply(df, as.numeric))
mM = as.matrix(df)

data2 = list(meta = metadata, rna=mG, methylation = mT, accessibility = mM)
X=list(rna = mG, methylation = mT, accessibility = mM)
lapply(X, dim)
## ---- fig.show = "hold", out.width = "33%", fig.height = 6, fig.cap = "FIGURE 1: Circle Correlation Plots for pairwise PLS models on the gastrulation data. Only displays the top 25 features for each dimension, subsetting by those with a correlation above 0.5."----
# select arbitrary values of features to keep
list.keepX = c(300, 300)
list.keepY = c(300, 300)

# generate three pairwise PLS models
pls1 <- spls(X[["rna"]], X[["methylation"]], keepX = list.keepX, keepY = list.keepY)
pls2 <- spls(X[["rna"]], X[["accessibility"]], keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(X[["methylation"]], X[["accessibility"]], keepX = list.keepX, keepY = list.keepY)

# plot features of first PLS
plotVar(pls1, cutoff = 0.5, title = "(a) RNA vs Methylation",
        legend = c("RNA", "Methylation"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2,1),
        col = c('darkorchid', 'lightgreen'))

# plot features of second PLS
plotVar(pls2, cutoff = 0.5, title = "(b) RNA vs Accessibility",
        legend = c("RNA", "Accessibility"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(1.5,1),
        col = c('darkorchid', 'lightgreen'))

# plot features of third PLS
plotVar(pls3, cutoff = 0.5, title = "(c) Methylation vs Accessibility",
        legend = c("Methylation", "Accessibility"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(1.25,1),
        col = c('darkorchid', 'lightgreen'))

## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cor(pls1$variates$X, pls1$variates$Y) # calculate correlation of RNA and methylation
cor(pls2$variates$X, pls2$variates$Y) # calculate correlation of RNA and accessibility
cor(pls3$variates$X, pls3$variates$Y) # calculate correlation of methylation and accessibility


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
design = matrix(0.5, ncol = length(X), nrow = length(X), # for square matrix filled with 0.5s
                dimnames = list(names(X), names(X)))
diag(design) = 0 # set diagonal to 0s

basic.mbspls.model = block.spls(X, indY = 1, # generate basic model
                                ncomp = 5, 
                                design = design)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
choice.ncomp <- 2 
choice.keepX <- list(rna = rep(5000, choice.ncomp), # 50 features per each component per dataframe
                     methylation = rep(5000, choice.ncomp), 
                     accessibility = rep(5000, choice.ncomp))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
final.mbspls.model = block.spls(X, indY = 1,  # generate final model using "tuned" parameters
                                ncomp = choice.ncomp, 
                                keepX = choice.keepX,
                                design = design)


## ---- out.width = "90%", fig.cap = "FIGURE 2: Sample plot for sPLS2 performed on the gastrulation data. Samples are projected into the space spanned by the components yielded from the RNA dataset."----
plotIndiv(final.mbspls.model, ind.names = FALSE,
          group = as.factor(cell_metadata2$lineage), 
          pch = as.factor(cell_metadata2$stage),
          col.per.group = color.mixo(1:6), 
          legend = TRUE, legend.title = 'Lineage', legend.title.pch = 'Stage',
          blocks = 2)

## ---- echo = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(425)

samples <- sample(1:dim(X$rna)[1], 17)

X1 <- X$rna[samples, ] 
X2 <- X$methylation[samples, ] 
X3 <- X$accessibility[samples, ] 
X.arrow <- list(rna = X1, methylation = X2, accessibility = X3)

final.mbspls.model.arrow = block.spls(X.arrow, indY = 1, 
                                ncomp = choice.ncomp, 
                                keepX = choice.keepX,
                                design = design)

## ---- out.width = "90%", fig.cap = "FIGURE 3: Arrow plot from the sPLS2 performed on the gastrulation data. The star indicates the location of the centroid of a sample across all the three datsets. The tip of each arrow shows the location of that same sample in the space spanned by the components associated to a specific dataset."----
symbols <- list(rna = 1, methylation = 6, accessibility = 10)
Factor = as.factor(cell_metadata2$lineage[samples])
plotArrow(final.mbspls.model.arrow, ind.names = FALSE, arrow.size = 0.3,
          arrow.length = 0.15, arrow.alpha = 0.5,
          group = Factor,box=TRUE,
          pch = symbols, pch.size = 2)

## ---- out.width = "90%", fig.cap = "FIGURE 4: Correlation circle plot from the sPLS2 performed on the gastrulation data"----------------------------------------------------
plotVar(final.mbspls.model, var.names = FALSE,
        legend = TRUE, cutoff = 0.5, cex = c(0.5, 0.5, 0.5),
        pch = c(0,1,2), )

## ---- echo = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------
X.2 = list(methylation = X$methylation, 
         accessibility = X$accessibility)
list.keepX = list(methylation = rep(1000, 2), accessibility = rep(1000,2))
list.keepY = c(rep(1000, 2))

final.mbspls.model.circos = block.spls(X = X.2, Y = data2$accessibility,
                                  ncomp = 2, keepX = list.keepX,
                                  keepY = list.keepY, design = 'full')


## ---- out.width = "90%", fig.cap = "FIGURE 5: Circos plot from multiblock sPLS performed on the gastrulation data The plot represents the correlations greater than 0.8 between variables of different types, represented on the side quadrants", results = FALSE----
circosPlot(final.mbspls.model.circos, group = cell_metadata2$lineage,
           cutoff = 0.99, comp = 1, linkWidth=0.05, size.variables=0.45,
           Y.name = 'rna')
#dev.off()
