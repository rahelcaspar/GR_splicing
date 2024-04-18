library("gprofiler2")

#Dex vs control ZT00
#get deseq2 results
deseq2_res_dex1 <- read.csv("Studium/Master_Thesis/Ergebnisse_public_data/deseq2/public_data_res_dex_ZT00_deseq2.csv")
#sort by log2FC
deseq2_res_dex1 <- deseq2_res_dex1[order(deseq2_res_dex1$log2FoldChange, decreasing = TRUE),]

#Dex vs control Gene Set Enrichment Analysis (ZT00)
#only for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
geneIDs <- deseq2_res_dex1[!is.na(deseq2_res_dex1$pvalue),]

geneIDs_enriched <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange > 1)),1]
geneIDs_enriched <- geneIDs_enriched[!is.na(geneIDs_enriched)]

geneIDs_depleted <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange < -1)),1]
geneIDs_depleted <- geneIDs_depleted[!is.na(geneIDs_depleted)]

#save list of enriched and list of depleted genes
write.csv(geneIDs_enriched, file = "Studium/Master_Thesis/Ergebnisse_public_data/deseq2_enriched_genes_dex00.csv", quote = FALSE, row.names = FALSE)
write.csv(geneIDs_depleted, file = "Studium/Master_Thesis/Ergebnisse_public_data/deseq2_depleted_genes_dex00.csv", quote = FALSE, row.names = FALSE)

gostresult <- gost(query = list("genes_enriched" = geneIDs_enriched, "genes_depleted" = geneIDs_depleted), organism = "mmusculus", ordered_query = TRUE)

#save results
result_enrichment_dex1 <- apply(gostresult$result, 2, as.character)
write.csv(result_enrichment_dex1, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_analysis_results_dex_ZT00.csv", quote = FALSE, row.names = FALSE)

#visualization
interactive_plot <- gostplot(gostresult)
htmlwidgets::saveWidget(interactive_plot, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_dex_ZT00.html")

png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_dex_ZT00.png", height = 710, width = 1603)
gostplot(gostresult, interactive = FALSE)
dev.off()


#Dex vs control Overrepresentation Analysis (ZT00)
#only for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
genes <- deseq2_res_dex1[!is.na(deseq2_res_dex1$pvalue),]

geneIDs_enriched <- genes[((genes$padj < 0.05) & (genes$log2FoldChange > 1)),1]
geneIDs_enriched <- geneIDs_enriched[!is.na(geneIDs_enriched)]

geneIDs_depleted <- genes[((genes$padj < 0.05) & (genes$log2FoldChange < -1)),1]
geneIDs_depleted <- geneIDs_depleted[!is.na(geneIDs_depleted)]

multi_gostresult <- gost(query = list("genes_enriched" = geneIDs_enriched, "genes_depleted" = geneIDs_depleted), organism = "mmusculus", ordered_query = FALSE)

#save results
result_overrepresentation_dex1 <- apply(multi_gostresult$result, 2, as.character)
write.csv(result_overrepresentation_dex1, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_results_dex_ZT00.csv", quote = FALSE, row.names = FALSE)

#visualization
interactive_plot2 <- gostplot(multi_gostresult)
htmlwidgets::saveWidget(interactive_plot2, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_dex_ZT00.html")

png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_dex_ZT00.png", height = 710, width = 1603)
gostplot(multi_gostresult, interactive = FALSE)
dev.off()


#Dex vs control ZT12
#get deseq2 results
deseq2_res_dex2 <- read.csv("Studium/Master_Thesis/Ergebnisse_public_data/deseq2/public_data_res_dex_ZT12_deseq2.csv")
#sort by log2FC
deseq2_res_dex2 <- deseq2_res_dex2[order(deseq2_res_dex2$log2FoldChange, decreasing = TRUE),]

#Dex vs control Gene Set Enrichment Analysis (ZT12)
#only for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
geneIDs <- deseq2_res_dex2[!is.na(deseq2_res_dex2$pvalue),]

geneIDs_enriched <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange > 1)),1]
geneIDs_enriched <- geneIDs_enriched[!is.na(geneIDs_enriched)]

geneIDs_depleted <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange < -1)),1]
geneIDs_depleted <- geneIDs_depleted[!is.na(geneIDs_depleted)]

#save list of enriched and list of depleted genes
write.csv(geneIDs_enriched, file = "Studium/Master_Thesis/Ergebnisse_public_data/deseq2_enriched_genes_dex12.csv", quote = FALSE, row.names = FALSE)
write.csv(geneIDs_depleted, file = "Studium/Master_Thesis/Ergebnisse_public_data/deseq2_depleted_genes_dex12.csv", quote = FALSE, row.names = FALSE)

gostresult <- gost(query = list("genes_enriched" = geneIDs_enriched, "genes_depleted" = geneIDs_depleted), organism = "mmusculus", ordered_query = TRUE)

#save results
result_enrichment_dex2 <- apply(gostresult$result, 2, as.character)
write.csv(result_enrichment_dex2, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_analysis_results_dex_ZT12.csv", quote = FALSE, row.names = FALSE)

#visualization
interactive_plot <- gostplot(gostresult)
htmlwidgets::saveWidget(interactive_plot, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_dex_ZT12.html")

png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_dex_ZT12.png", height = 710, width = 1603)
gostplot(gostresult, interactive = FALSE)
dev.off()


#Dex vs control Overrepresentation Analysis (ZT12)
#only for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
genes <- deseq2_res_dex2[!is.na(deseq2_res_dex2$pvalue),]

geneIDs_enriched <- genes[((genes$padj < 0.05) & (genes$log2FoldChange > 1)),1]
geneIDs_enriched <- geneIDs_enriched[!is.na(geneIDs_enriched)]

geneIDs_depleted <- genes[((genes$padj < 0.05) & (genes$log2FoldChange < -1)),1]
geneIDs_depleted <- geneIDs_depleted[!is.na(geneIDs_depleted)]

multi_gostresult <- gost(query = list("genes_enriched" = geneIDs_enriched, "genes_depleted" = geneIDs_depleted), organism = "mmusculus", ordered_query = FALSE)

#save results
result_overrepresentation_dex2 <- apply(multi_gostresult$result, 2, as.character)
write.csv(result_overrepresentation_dex2, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_results_dex_ZT12.csv", quote = FALSE, row.names = FALSE)

#visualization
interactive_plot2 <- gostplot(multi_gostresult)
htmlwidgets::saveWidget(interactive_plot2, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_dex_ZT12.html")

png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_dex_ZT12.png", height = 710, width = 1603)
gostplot(multi_gostresult, interactive = FALSE)
dev.off()


#Dex vs control (both timepoints)
#get deseq2 results
deseq2_res_dex <- read.csv("Studium/Master_Thesis/Ergebnisse_public_data/deseq2/public_data_res_dex_deseq2.csv")
#sort by log2FC
deseq2_res_dex <- deseq2_res_dex[order(deseq2_res_dex$log2FoldChange, decreasing = TRUE),]

#Dex vs control Gene Set Enrichment Analysis (both timepoints)
#only for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
geneIDs <- deseq2_res_dex[!is.na(deseq2_res_dex$pvalue),]

geneIDs_enriched <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange > 1)),1]
geneIDs_enriched <- geneIDs_enriched[!is.na(geneIDs_enriched)]

geneIDs_depleted <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange < -1)),1]
geneIDs_depleted <- geneIDs_depleted[!is.na(geneIDs_depleted)]

#save list of enriched and list of depleted genes
write.csv(geneIDs_enriched, file = "Studium/Master_Thesis/Ergebnisse_public_data/deseq2_enriched_genes_dex.csv", quote = FALSE, row.names = FALSE)
write.csv(geneIDs_depleted, file = "Studium/Master_Thesis/Ergebnisse_public_data/deseq2_depleted_genes_dex.csv", quote = FALSE, row.names = FALSE)

gostresult <- gost(query = list("genes_enriched" = geneIDs_enriched, "genes_depleted" = geneIDs_depleted), organism = "mmusculus", ordered_query = TRUE)

#save results
result_enrichment_dex <- apply(gostresult$result, 2, as.character)
write.csv(result_enrichment_dex, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_analysis_results_dex.csv", quote = FALSE, row.names = FALSE)

#visualization
interactive_plot <- gostplot(gostresult)
htmlwidgets::saveWidget(interactive_plot, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_dex.html")

png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_dex.png", height = 710, width = 1603)
gostplot(gostresult, interactive = FALSE)
dev.off()


#Dex vs control Overrepresentation Analysis (both timepoints)
#only for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
genes <- deseq2_res_dex[!is.na(deseq2_res_dex$pvalue),]

geneIDs_enriched <- genes[((genes$padj < 0.05) & (genes$log2FoldChange > 1)),1]
geneIDs_enriched <- geneIDs_enriched[!is.na(geneIDs_enriched)]

geneIDs_depleted <- genes[((genes$padj < 0.05) & (genes$log2FoldChange < -1)),1]
geneIDs_depleted <- geneIDs_depleted[!is.na(geneIDs_depleted)]

multi_gostresult <- gost(query = list("genes_enriched" = geneIDs_enriched, "genes_depleted" = geneIDs_depleted), organism = "mmusculus", ordered_query = FALSE)

#save results
result_overrepresentation_dex <- apply(multi_gostresult$result, 2, as.character)
write.csv(result_overrepresentation_dex, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_results_dex.csv", quote = FALSE, row.names = FALSE)

#visualization
interactive_plot2 <- gostplot(multi_gostresult)
htmlwidgets::saveWidget(interactive_plot2, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_dex.html")

png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_dex.png", height = 710, width = 1603)
gostplot(multi_gostresult, interactive = FALSE)
dev.off()

