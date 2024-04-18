#new_data (6 samples)
library("DESeq2")
#load deseq2-dataset from nf-core/rnaseq pipeline
load("Studium/Master_Thesis/Ergebnisse_new_data/new_data_Pipeline_Output/deseq2.dds.RData")


#create deseq2-dataset
raw_counts_dds <- environment(dds@dispersionFunction)[["dds"]]@assays@data@listData[["counts"]]
sample_conditions <- read.csv("Studium/Master_Thesis/Daten_download/new_data_Samples_Uebersicht.csv")

dds2 <- DESeqDataSetFromMatrix(countData = raw_counts_dds, colData = sample_conditions, design = ~condition)

#run DESeq2
dds2 <- DESeq(dds2)
res <- results(dds2)
#sort results by adjusted p-value
res <- res[order(res$padj),]


#save results from DESeq2-run
save(res, file = "Studium/Master_Thesis/Ergebnisse_new_data/new_data_res_deseq2.RData")
write.csv(res, file = "Studium/Master_Thesis/Ergebnisse_new_data/new_data_res_deseq2.csv", quote = FALSE)

#volcano plot
png(filename = "Studium/Master_Thesis/Ergebnisse_new_data/new_data_volcano-plot.png")
with(res, plot(log2FoldChange, -log10(padj), pch = 20, main = "Hepatocytes Volcano Plot"))
# Add colored points: blue if padj<0.05, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))
dev.off()

