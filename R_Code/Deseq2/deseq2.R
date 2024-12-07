setwd("path/to/working/directory")


library(DESeq2)
library(tidyverse)

counts_data <- read.csv('Original_data.csv', row.names = 1)
head(counts_data)
which(duplicated(counts_data$gene) == TRUE)

colData <- read.csv('Meta_data_for_deseq2.csv', row.names = 1, header = TRUE)


# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))
# are they in the same order?
all(colnames(counts_data) == rownames(colData))



## Normalization
dds_norm <- vst(dds)
norm.counts <- assay(dds_norm)

# PCA
pca <- prcomp(t(counts_data))
pca.dat <- pca $x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digit =2)
pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2))+
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0("PC1: ", pca.var.percent[1],'%'),
       y = paste0("PC2: ", pca.var.percent[2], '%'))




# construction of DESeqDataSet object

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~Type)
dds
summary(dds)

dds75 <- dds[rowSums(counts(dds) >= 10) >= 85,] ## Remove all genes with counts < 10 in more than 75% of samples
nrow(dds75)
dds75



# set the factor level
dds75$Type <- relevel(dds75$Type, ref='Normal')


# Run DESeq
dds <- DESeq(dds75)
res <- results(dds)
dim(res)
summary(res)

# Write the output
write.csv(res, "deseq_result.csv" )


# Filtering out the DEGs
results1 =read.csv("deseq_result.csv", header=TRUE, row.names = 1)
EnhancedVolcano(res,
                lab = NA,
                x = 'log2FoldChange',
                y = 'pvalue')

library(dplyr)
Upregulated_data <- results1 %>% filter(padj < 0.05 & log2FoldChange > 1)
dim(Upregulated_data)
write.csv(Upregulated_data ,file="upregulated.csv")

downregulated_data <- results1 %>% filter(padj < 0.05 & log2FoldChange < -1)
dim(downregulated_data)
write.csv(downregulated_data ,file="downregulated.csv")



 ### gene id conversation if needed

install.packages('org.Hs.eg.db')


library(biomaRt)
library(annotables)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

# input list of Ensembl ID's
ensembl.ids <- read.delim('list_of_Ensembl_ID.csv', header = T)


# method 1: biomaRt
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)

ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')

attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

x <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
      filters = "ensembl_gene_id",
      values = ensembl.ids$X,
      mart = ensembl.con)

write.csv(x,"up_r_name.csv")

# method 2: annotables
grch38 %>%
  filter(ensgene %in% ensembl.ids$V1)


# method 3: annotation DBs

keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

mapIds(org.Hs.eg.db,
       keys = ensembl.ids$V1,
       keytype = 'ENSEMBL',
       column = 'SYMBOL')

keytypes(EnsDb.Hsapiens.v86)
columns(EnsDb.Hsapiens.v86)

mapIds(EnsDb.Hsapiens.v86,
       keys = ensembl.ids$V1,
       keytype = 'GENEID',
       column = 'SYMBOL')






