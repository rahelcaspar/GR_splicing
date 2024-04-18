library("ggplot2")
library("ggrepel")

#new_data (6 samples)
#load deseq2-dataset from nf-core/rnaseq pipeline
load("Studium/Master_Thesis/Ergebnisse_new_data/new_data_Pipeline_Output/deseq2.dds.RData")

#add experimental conditions to deseq2-dataset
raw_counts_dds <- environment(dds@dispersionFunction)[["dds"]]@assays@data@listData[["counts"]]
sample_conditions <- read.csv("Studium/Master_Thesis/Daten_download/new_data_Samples_Uebersicht.csv")

dds2 <- DESeqDataSetFromMatrix(countData = raw_counts_dds, colData = sample_conditions, design = ~condition)
save("dds2", file="Studium/Master_Thesis/Ergebnisse_new_Data/new_data_dds_cond.RData")

#get the variance-stabilized-transformation of the raw count data
vst_data <- vst(dds2, blind = FALSE)

#colored PCA ...
# ... of the top 500 genes
#with point labels
png(filename = "Studium/Master_Thesis/Ergebnisse_new_data/new_data_pca_top500_IDs.png", width = 800, height = 600)
pcaData <- plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition, label=sample_conditions$sample)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("Hepatocytes PCA on vst-data", subtitle = "Top 500 genes") +
  geom_text_repel() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)) + 
  coord_fixed()
dev.off()

# ... of all genes
#with point labels
png(filename = "Studium/Master_Thesis/Ergebnisse_new_data/new_data_pca_IDs.png", width = 800, height = 600)
pcaData <- plotPCA(vst_data, intgroup = "condition", ntop = length(vst_data), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition, label=sample_conditions$sample)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("Hepatocytes PCA on vst-data", subtitle = " ") +
  geom_text_repel() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)) + 
  coord_fixed()
dev.off()
