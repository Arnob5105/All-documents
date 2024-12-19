#Install from Bioconductor repository
BiocManager::install("maftools")
#Install from GitHub repository
BiocManager::install("PoisonAlien/maftools")
library(maftools)

rm(list = ls())
setwd("F:/WES/maffiles/merge/")
## import data
new <- read.maf(maf= "F:/WES/maffiles/merge/_maftools.maf")

count <- read.csv("F:/WES/maffiles/merge/counts.csv")
a <- unique(count)
## plot VAF

plotVaf(
  new,
  vafCol = NULL,
  #genes = ,
  top = 20,
  orderByMedian = TRUE,
  keepGeneOrder = FALSE,
  flip = FALSE,
  fn = NULL,
  gene_fs = 0.8,
  axis_fs = 0.8,
  height = 5,
  width = 5,
  showN = TRUE,
  color = NULL
)


plotVaf(
  new,
  #vafCol = NULL,
  genes = 'ATRX' )
  #top = 20)
  
  
## summary of the maffile
plotmafSummary(maf = new, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

mafbarplot(maf = new, gene = c("DLX6", "MTCH2", "SSC5D","ANKRD36", "VWF", "FLG", "ANKRD36", "VWF", "FLG", "WDR89", "GOLGA6L2", "GNAQ", "MUC17", "OR10G7", "IRS2","GXYLT1", "TMEM163", "MUC16", "KMT2C", "ZSCAN5B", "KRT18", "OR10G7", "PABPC3", "ZNF208")) 
  ###"ANKRD36", "VWF", "FLG", "WDR89", "GOLGA6L2", "GNAQ", "MUC17", "OR10G7", "IRS2", "GXYLT1", "TMEM163", "MUC16", "KMT2C"))

## somatic mutational landscape
oncoplot(maf = new,gene_mar = 10, top = 20)

## with different color and adding Transition versus transversion
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

oncoplot(maf = new, top = 20, gene_mar = 7, fontSize = 0.9, colors = vc_cols)

### adding pathways

oncoplot(maf = new, pathways = "sigpw",gene_mar = 7,topPathways = 5)

### adding biological process

oncoplot(maf = new, pathways = "smgbp",gene_mar = 14,topPathways = 5)

## plot titv summary ( Transition versus transversion )
new.titv = titv(maf = new, plot = TRUE, useSyn = TRUE)
plotTiTv(res = new.titv)

## rainfall plot (hyper muted genomic region)
rain <- rainfallPlot(maf = new, tsb = NULL , detectChangePoints = TRUE, pointSize = 0.8, ref.build = 'hg38')

## lolipop plot for showing mutation in genes
new <- read.maf(maf= "F:/WES/maffiles/merge/_maftools _new.maf")

a <- lollipopPlot(
  maf = new,
  gene = 'FOXC1',
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
  showDomainLabel = FALSE,
  labPosSize = 1.2,
  printCount = TRUE,
  #labelPos = 103 ,
  labPosAngle = 0,
  repel = FALSE,
  legendTxtSize = 1.4,
  axisTextSize = c(1.3,1.3),
  labelOnlyUniqueDoamins = FALSE,
  titleSize = c(1.4,1.2),
  pointSize = 1.8,
  domainLabelSize = 1.4
  )

write.csv(a, "WDR89.csv")

laml.maf = system.file('extdata', 'protein_domains.csv', package = 'maftools')
print(laml.maf)

laml = read.csv(laml.maf)
prot = .getdomains(geneID = 'MUC3A', refSeqID = refSeqID, proteinID = proteinID)

rd <- readRDS("C:/Users/HP/AppData\Local\R\win-library\4.3\maftools\extdata")
saveRDS(laml, "protein_domains.RDs")



##  Compare mutational load against TCGA cohort 

gbm.mutload = tcgaCompare(maf = new, cohortName = 'GBM_Cohort', logscale = TRUE, capture_size = 50)


## mutually exclusive and inclusive 


# considering 15 out of 20, because top 5 showed no interaction
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c( "gray", "red", "blue", "green", "orange", "purple", "brown", "pink")

library(RColorBrewer)
custom_palette <- RColorBrewer::brewer.pal(n = 12, name = "Greens")

b <- somaticInteractions(maf = new, 
                    #top = 20,
                    genes = c("IRS2", 'GNAQ', "ANKRD36", "FLG", "KMT2C", "OR10G7", "SSC5D", "TMEM163", "VWF", "GXYLT1", "MTCH2", "KMT2C"),
                    countStats = "all",
                    countType = "all",
                    nShiftSymbols = 2,
                    sigSymbolsSize = 3,
                    #geneOrder = c('GNAQ', "IRS2", "ANKRD36", "FLG", "KMT2C", "OR10G7", "SSC5D", "TMEM163", "VWF", "GXYLT1", "MTCH2", "KMT2C", "GOLGA6L2", "MUC16", "DLX6", "WDR89"),
                    fontSize = 0.7,
                    pvalue = c(0.10, 0.20),
                    colPal = "BrBG")



# considering all top 20

somaticInteractions(maf = new, 
                    top = 20,
                    countStats = "all",
                    countType = "all",
                    nShiftSymbols = 2,
                    sigSymbolsSize = 3,
                    #sigSymbolsFontSize = 5,	
                    countsFontSize = 5,
                    fontSize = 0.7,
                    pvalue = c(0.1, 0.15))
                    


# Oncodrive genes
 
gbm.sig = oncodrive(maf = new, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')

plotOncodrive(res = gbm.sig, fdrCutOff = 0.1, useFraction = TRUE, colCode = NULL, bubbleSize = 1, labelSize = 0.8)

write.csv(gbm.sig, "driver.csv")
# pfam domain

gbm.pfam = pfamDomains(maf = new, AACol = 'Protein_Change', top = 5)

## Drug-gene interaction
dgi = drugInteractions(maf = new, fontSize = 0.75)

dnmt3a.dgi = drugInteractions(genes = "MUC16", drugs = TRUE)
dnmt3a.dgi[,.(Gene, interaction_types, drug_name, drug_claim_name)]


### adding pathways

oncoplot(maf = new, pathways = "sigpw",gene_mar = 7,topPathways = 5)

### adding biological process

oncoplot(maf = new, pathways = "smgbp",gene_mar = 14,topPathways = 5)


## pathway
pws = pathways(maf = new, plotType = 'treemap', pathways = "sigpw" )

plotPathways(maf = new, pathlist = pws)

plotSignatures(maf = new, title_size = 1.2, sig_db = "SBS")




laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read.maf(maf = laml.maf)
PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS")

## annovar to maf
var <- annovarToMaf(annovar = 'F:/WES/maffiles/annovar_tumor.txt', refBuild = 'hg38',tsbCol = 'sample_id')
annovar<- read.maf(maf = var)
a <- annovar@data 
plotmafSummary(maf = annovar, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
rainfallPlot(maf = annovar, detectChangePoints = TRUE, pointSize = 0.8)



## mutational signatures

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
laml.tnm = trinucleotideMatrix(maf = new, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")




