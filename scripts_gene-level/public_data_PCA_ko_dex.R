library("ggplot2")
library("ggrepel")

#public data 
#load deseq2-dataset from nf-core/rnaseq pipeline
load("Studium/Master_Thesis/Ergebnisse_public_data/public_data_Pipeline_Output/deseq2.dds.RData")

#KO vs WT PCAs  (color – KO/WT, shape - timepoint)
# 16 KO, 17 WT samples

#load conditions table for KO vs WT
sample_conditions <- read.csv("Studium/Master_Thesis/Daten_download/samples_knockout.csv")
raw_counts_dds <- environment(dds@dispersionFunction)[["dds"]]@assays@data@listData[["counts"]]
#filter raw counts
raw_counts_ko <- raw_counts_dds[,colnames(raw_counts_dds) %in% sample_conditions$sample]

dds_ko <- DESeqDataSetFromMatrix(countData = raw_counts_ko, colData = sample_conditions, design = ~ timepoint + knockout)
dds_ko$knockout <- factor(colData(dds_ko)$knockout, levels=c("WT", "KO"))
vst_ko <- vst(dds_ko, blind = FALSE)

#colored PCA KO vs WT over timepoints (with point labels)
# ... of the top 500 genes
png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_pca_top500_ko.png", width = 800, height = 600)
pcaData <- plotPCA(vst_ko, intgroup = c("knockout","timepoint"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=knockout, shape=timepoint, label=sample_conditions$sample)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text_repel() +
  ggtitle("Liver PCA on vst-counts", subtitle = "Top 500 genes") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)) + 
  coord_fixed()
dev.off()

# ... of all genes
png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_pca_ko.png", width = 800, height = 600)
pcaData <- plotPCA(vst_ko, intgroup = c("knockout","timepoint"), ntop = length(vst_ko), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=knockout, shape=timepoint, label=sample_conditions$sample)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text_repel() +
  ggtitle("Liver PCA on vst-counts", subtitle = " ") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)) + 
  coord_fixed()
dev.off()


#dex vs control PCAs (color – dex/control, shape - timepoint)
# 5 dex, 4 control samples

#load conditions table for dex vs control
sample_conditions <- read.csv("Studium/Master_Thesis/Daten_download/samples_dex.csv")
raw_counts_dds <- environment(dds@dispersionFunction)[["dds"]]@assays@data@listData[["counts"]]
#filter raw counts
raw_counts_dex <- raw_counts_dds[,colnames(raw_counts_dds) %in% sample_conditions$sample]

dds_dex <- DESeqDataSetFromMatrix(countData = raw_counts_dex, colData = sample_conditions, design = ~ timepoint + dex)
vst_dex <- vst(dds_dex, blind = FALSE)

#colored PCA Dex vs control over timepoints (with point labels)
# ... of the top 500 genes
png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_pca_top500_dex.png", width = 800, height = 600)
pcaData <- plotPCA(vst_dex, intgroup = c("dex","timepoint"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=dex, shape=timepoint, label=sample_conditions$sample)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text_repel() +
  ggtitle("Liver PCA on vst-counts", subtitle = "Top 500 genes") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)) + 
  coord_fixed()
dev.off()

# ... of all genes
png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_pca_dex.png", width = 800, height = 600)
pcaData <- plotPCA(vst_dex, intgroup = c("dex","timepoint"), ntop = length(vst_dex), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=dex, shape=timepoint, label=sample_conditions$sample)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text_repel() +
  ggtitle("Liver PCA on vst-counts", subtitle = " ") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)) + 
  coord_fixed()
dev.off()

