setwd("path/to/working/directory")


library(survminer)
library(survival)
library(tidyverse)
library(DESeq2)

# Data Normalization
raw <- read.csv('Original_data.csv', row.names = 1)
meta <- read.csv('Meta_data_for_deseq2.csv', row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = raw,
                              colData = meta,
                              design = ~1)

dds_norm <- vst(dds,blind = TRUE)
dds_count <- assay (dds_norm) 
df <- dds_count %>% as.data.frame() 
write.csv(df, 'vst_normalized_data.csv')


# Survival Analysis
data <- read.csv('LGG/Survival_data_lgg.csv')
surv_object <- Surv(data$time, data$status)
data$HIF1A <- ifelse(data$gene > median(data$gene), "High", "Low")
km_fit <- survfit(surv_object ~ HIF1A, data = data)
p1 <- ggsurvplot(km_fit, data = data, pval = TRUE, risk.table = FALSE, 
           title = "Primary Glioma (LGG)", 
           xlab = "Time (Months)", ylab = "Survival Probability")


data2 <- read.csv('HGG/Meta_survival_hgg.csv')
surv_object <- Surv(data2$time, data2$status)
data2$HIF1A <- ifelse(data2$gene > median(data2$gene), "High", "Low")
km_fit2 <- survfit(surv_object ~ HIF1A, data = data2)
p2 <- ggsurvplot(km_fit2, data = data2, pval = TRUE, risk.table = FALSE, 
           title = "Secondary Glioma (HGG)", 
           xlab = "Time (Months)", ylab = "Survival Probability")

# Making Grid
gridExtra::grid.arrange(p1$plot,p2$plot, ncol=2)

