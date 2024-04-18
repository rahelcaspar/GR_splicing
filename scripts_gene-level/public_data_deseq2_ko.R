#public data: KO vs WT DESeq2 (16 vs 17 samples)
library("DESeq2")
library("ggplot2")
library("ggrepel")

#load deseq2-dataset from nf-core/rnaseq pipeline
load("Studium/Master_Thesis/Ergebnisse_public_data/public_data_Pipeline_Output/deseq2.dds.RData")


#create deseq2-dataset
sample_conditions <- read.csv("Studium/Master_Thesis/Daten_download/samples_knockout.csv")
raw_counts_dds <- environment(dds@dispersionFunction)[["dds"]]@assays@data@listData[["counts"]]
raw_counts_ko <- raw_counts_dds[,colnames(raw_counts_dds) %in% sample_conditions$sample]

dds_ko <- DESeqDataSetFromMatrix(countData = raw_counts_ko, colData = sample_conditions, design = ~ timepoint + knockout)
dds_ko$knockout <- factor(colData(dds_ko)$knockout, levels=c("WT", "KO"))

#run DESeq2
dds_ko <- DESeq(dds_ko)
res_ko <- results(dds_ko)
#sort results by adjusted p-value
res_ko <- res_ko[order(res_ko$padj),]

#save results from DESeq2-run
save(res_ko, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_res_ko_deseq2.RData")
write.csv(res_ko, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_res_ko_deseq2.csv", quote = FALSE)


#volcano plot  
png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_ko_volcano-plot.png")
with(res_ko, plot(log2FoldChange, -log10(padj), pch = 20, main = "Liver Volcano Plot KO_WT"))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05
with(subset(res_ko, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res_ko, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
#with(subset(res_ko, -log10(padj) > 30 & abs(log2FoldChange)>1), text(log2FoldChange, -log10(padj), labels = rownames(res_ko)))
dev.off()
