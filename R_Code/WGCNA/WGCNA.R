setwd("path/to/working/directory")

library(WGCNA)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(doParallel)
library(DESeq2)
allowWGCNAThreads()


data <- read.csv('Original_data.csv',row.names = 1,header = TRUE)
meta <- read.csv('Meta_data_for_WGCNA.csv', row.names = 1, header = TRUE)        # For the Construction of Module-Trait Relationships
meta_deseq <- read.csv('Meta_data_for_deseq2.csv', row.names = 1, header = TRUE) # For data normalization

# Checking the equality of samples
all(rownames(meta_deseq) %in% colnames(data))
all(rownames(meta_deseq) == colnames(data))

# Good gene extraction
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)
data <- data[gsg$goodGenes == TRUE,]

# Good sample extraction
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

pca <- prcomp(t(data))
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


samples.to.be.excluded <- c('CGGA_P154', 'CGGA_1699', 'CGGA_1635', 'CGGA_1542', 'CGGA_1467', 'CGGA_1494', 'CGGA_1678', 'CGGA_1481', 'CGGA_1418', 'CGGA_1498', 'CGGA_1378', 'CGGA_P116', 'CGGA_1666', 'CGGA_1422', 'CGGA_P25', 'CGGA_1086', 'CGGA_1767', 'CGGA_1486', 'CGGA_P180')
data <- data[,!(colnames(data) %in% samples.to.be.excluded)]
meta <- meta[!(rownames(meta) %in% samples.to.be.excluded),]
meta_deseq <- meta_deseq[!(rownames(meta_deseq) %in% samples.to.be.excluded),]


# Data Normalization
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta_deseq,
                              design = ~1)
dds
dds75 <- dds[rowSums(counts(dds) >= 10) >= 60,] # Remove all genes with counts < 10 in more than 75% of samples
nrow(dds75) 
dds_norm <- vst(dds75)
norm.counts <- assay(dds_norm) %>% 
  t()



## Calculating power

# First method
power_value <- c(c(1:10), seq(from = 12, to = 50, by = 2))
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power_value,
                         networkType = "unsigned",
                         verbose = 5)
sft.data <- sft$fitIndices


a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.89, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
grid.arrange(a1, a2, nrow = 2)





# Second Method
powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft = pickSoftThreshold(t(dds_count), powerVector = powers, networkType = "signed", verbose = 5) 
par(mfrow = c(1,2));
cex1 = 0.8;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


### Calculation Of WGCNA
#power <- sft$powerEstimate

power <- 18

adjacency = adjacency(norm.counts, power = power)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

save(dissTOM, file = "dissTOM.RData")

geneTree = hclust(as.dist(dissTOM), method = "average");
#pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 30);

table(dynamicMods)
length(table(dynamicMods)) 
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file = "4-module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

MEList = moduleEigengenes(norm.counts, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");

sizeGrWindow(7, 6)

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres=0
abline(h=MEDissThres, col = "red")


merge = mergeCloseModules(norm.counts, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs 


plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Unmerged dynamic", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()


## Construction of Module Trait Relationships

nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(mergedMEs, meta, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)


## Ploting Method 1
heatmap.data <- merge(mergedMEs, meta, by = 'row.names')
rownames(heatmap.data) <- heatmap.data[,1]
heatmap.data <- heatmap.data[,-1 ] 

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[5:8],
             y = names(heatmap.data)[1:4],
             col = c("#34C759", "#F7DC6F", "#FFA07A", "#FF3737", "#8B0A1A"))


## Ploting Method 2
textMatrix =  paste(signif(module.trait.corr, 2), "\n(",
                    signif(module.trait.corr.pvals, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.corr)
pdf("module_trait_glioblastoma_2.pdf", width = 10, height = 15)
par(mar = c(6, 9, 2, 1));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module.trait.corr,
               xLabels = colnames(meta),
               yLabels = colnames(mergedMEs),
               ySymbols = colnames(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


# Module Gene Extraction
module.gene.mapping <- as.data.frame(mergedColors)
module.gene.mapping <- cbind(colnames(norm.counts),module.gene.mapping)
gene <-module.gene.mapping %>% 
  filter(`mergedColors` == 'darkorange2') 
#%>% rownames()
write.csv(gene, 'Data/New/short/Grade2/Without_DEGs/Genes/darkorange2.csv', row.names = FALSE)




modNames = substring(names(mergedMEs), 3)
#MET = orderMEs(cbind(MEs, Verru))
geneModuleMembership = as.data.frame(cor(norm.counts, mergedMEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(norm.counts), meta, use = "p");
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(meta), sep="");
names(GSPvalue) = paste("p.GS.", names(meta), sep="");


module = "darkorange2"
# Rename to moduleColors
moduleColors = mergedColors
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
#par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for glioblastoma",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "purple")

GS <- 0.5
MM <- 0.7

abline(h=GS, col = "black")
abline(v=MM, col = "black")
dev.off()

library(dplyr)
#mm_threshold <- 0.85
#gs_threshold <- 0.85
module <-as.data.frame(geneModuleMembership[moduleGenes, column])
significance <- as.data.frame(geneTraitSignificance[moduleGenes, 1])
colnames(significance) <- 'sig'
colnames(module) <- 'membership'

binding <- cbind(gene,module,significance)
filter <-  binding %>% filter(abs(membership) > MM & abs(sig) > GS)

row.names(filter) <- filter$`colnames(norm.counts)`
print(rownames(filter))
rownames(filter)
# Filtered genes based on module membership and gene significance

write.csv(rownames(filter), 'Data/New/short/Grade2/Without_DEGs/Genes/hubdarkorange2.csv', row.names = FALSE)
