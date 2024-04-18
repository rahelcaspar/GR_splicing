#public data: dex vs control DESeq2 (5 vs 4 samples)
library("DESeq2")
library("ggplot2")
library("ggrepel")

#load deseq2-dataset from nf-core/rnaseq pipeline
load("Studium/Master_Thesis/Ergebnisse_public_data/public_data_Pipeline_Output/deseq2.dds.RData")

sample_conditions <- read.csv("Studium/Master_Thesis/Daten_download/samples_dex.csv")
raw_counts_dds <- environment(dds@dispersionFunction)[["dds"]]@assays@data@listData[["counts"]]
raw_counts_dex <- raw_counts_dds[,colnames(raw_counts_dds) %in% sample_conditions$sample]


#ZT00 (3 dex vs 2 control samples)
raw_counts_dex1 <- raw_counts_dex[,(sample_conditions$timepoint == "ZT00")]
sample_conditions1 <- sample_conditions[(sample_conditions$timepoint == "ZT00"),]

#create deseq2-dataset
dds_dex1 <- DESeqDataSetFromMatrix(countData = raw_counts_dex1, colData = sample_conditions1, design = ~ dex)

#run DESeq2
dds_dex1 <- DESeq(dds_dex1)
res_dex1 <- results(dds_dex1)
#sort results by adjusted p-value
res_dex1 <- res_dex1[order(res_dex1$padj),]

#save results from DESeq2-run
save(res_dex1, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_res_dex_ZT00_deseq2.RData")
write.csv(res_dex1, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_res_dex_ZT00_deseq2.csv", quote = FALSE)

#volcano plot  
png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_dex_ZT00_volcano-plot.png")
with(res_dex1, plot(log2FoldChange, -log10(padj), pch = 20, main = "Liver Volcano Plot Dex vs Control (ZT00)"))
# Add colored points: blue if padj<0.05, red if |log2FC|>1 and padj<0.05
with(subset(res_dex1, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res_dex1, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
dev.off()


#ZT12 (2 dex vs 2 control samples)
raw_counts_dex2 <- raw_counts_dex[,(sample_conditions$timepoint == "ZT12")]
sample_conditions2 <- sample_conditions[(sample_conditions$timepoint == "ZT12"),]

#create deseq2-dataset
dds_dex2 <- DESeqDataSetFromMatrix(countData = raw_counts_dex2, colData = sample_conditions2, design = ~ dex)

#run DESeq2
dds_dex2 <- DESeq(dds_dex2)
res_dex2 <- results(dds_dex2)
#sort results by adjusted p-value
res_dex2 <- res_dex2[order(res_dex2$padj),]

#save results from DESeq2-run
save(res_dex2, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_res_dex_ZT12_deseq2.RData")
write.csv(res_dex2, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_res_dex_ZT12_deseq2.csv", quote = FALSE)

#volcano plot  
png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_dex_ZT12_volcano-plot.png")
with(res_dex2, plot(log2FoldChange, -log10(padj), pch = 20, main = "Liver Volcano Plot Dex vs Control (ZT12)"))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05
with(subset(res_dex2, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res_dex2, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
dev.off()


#both timepoints (5 dex vs 4 control samples)
#create deseq2-dataset
dds_dex <- DESeqDataSetFromMatrix(countData = raw_counts_dex, colData = sample_conditions, design = ~ timepoint + dex)

#run DESeq2
dds_dex <- DESeq(dds_dex)
res_dex <- results(dds_dex)
#sort results by adjusted p-value
res_dex <- res_dex[order(res_dex$padj),]

#save results from DESeq2-run
save(res_dex, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_res_dex_deseq2.RData")
write.csv(res_dex, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_res_dex_deseq2.csv", quote = FALSE)

#volcano plot  
png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_dex_volcano-plot.png")
with(res_dex, plot(log2FoldChange, -log10(padj), pch = 20, main = "Liver Volcano Plot Dex vs Control"))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05
with(subset(res_dex, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res_dex, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
dev.off()
