library("gprofiler2")

#get deseq2 results
deseq2_res_ko <- read.csv("Studium/Master_Thesis/Ergebnisse_public_data/deseq2/public_data_res_ko_deseq2.csv")
#sort by log2FC
deseq2_res_ko <- deseq2_res_ko[order(deseq2_res_ko$log2FoldChange, decreasing = TRUE),]


#KO vs WT Gene Set Enrichment Analysis
#only for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
geneIDs <- deseq2_res_ko[!is.na(deseq2_res_ko$pvalue),]

geneIDs_enriched <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange > 1)),1]
geneIDs_enriched <- geneIDs_enriched[!is.na(geneIDs_enriched)]

geneIDs_depleted <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange < -1)),1]
geneIDs_depleted <- geneIDs_depleted[!is.na(geneIDs_depleted)]

#save list of enriched and list of depleted genes
write.csv(geneIDs_enriched, file = "Studium/Master_Thesis/Ergebnisse_public_data/deseq2_enriched_genes_ko.csv", quote = FALSE, row.names = FALSE)
write.csv(geneIDs_depleted, file = "Studium/Master_Thesis/Ergebnisse_public_data/deseq2_depleted_genes_ko.csv", quote = FALSE, row.names = FALSE)

gostresult <- gost(query = list("genes_enriched" = geneIDs_enriched, "genes_depleted" = geneIDs_depleted), organism = "mmusculus", ordered_query = TRUE)

#save results
result_enrichment_ko <- apply(gostresult$result, 2, as.character)
write.csv(result_enrichment_ko, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_analysis_results_ko.csv", quote = FALSE, row.names = FALSE)

#visualization
interactive_plot <- gostplot(gostresult)
htmlwidgets::saveWidget(interactive_plot, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_ko.html")

png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_gs_enrichment_ko.png", height = 710, width = 1603)
gostplot(gostresult, interactive = FALSE)
dev.off()



#KO vs WT Overrepresentation Analysis
#only for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
genes <- deseq2_res_ko[!is.na(deseq2_res_ko$pvalue),]

geneIDs_enriched <- genes[((genes$padj < 0.05) & (genes$log2FoldChange > 1)),1]
geneIDs_enriched <- geneIDs_enriched[!is.na(geneIDs_enriched)]

geneIDs_depleted <- genes[((genes$padj < 0.05) & (genes$log2FoldChange < -1)),1]
geneIDs_depleted <- geneIDs_depleted[!is.na(geneIDs_depleted)]

multi_gostresult <- gost(query = list("genes_enriched" = geneIDs_enriched, "genes_depleted" = geneIDs_depleted), organism = "mmusculus", ordered_query = FALSE)

#save results
result_overrepresentation_ko <- apply(multi_gostresult$result, 2, as.character)
write.csv(result_overrepresentation_ko, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_results_ko.csv", quote = FALSE, row.names = FALSE)

#visualization
interactive_plot2 <- gostplot(multi_gostresult)
htmlwidgets::saveWidget(interactive_plot2, file = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_ko.html")

png(filename = "Studium/Master_Thesis/Ergebnisse_public_data/public_data_overrepresentation_ko.png", height = 710, width = 1603)
gostplot(multi_gostresult, interactive = FALSE)
dev.off()

