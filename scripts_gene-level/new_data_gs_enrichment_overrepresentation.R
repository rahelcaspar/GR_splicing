library("gprofiler2")

#get deseq2 results 
deseq2_res <- read.csv("Studium/Master_Thesis/Ergebnisse_new_data/deseq2/new_data_res_deseq2.csv")
#sort by log2FC
deseq2_res <- deseq2_res[order(deseq2_res$log2FoldChange, decreasing = TRUE),]


#Gene Set Enrichment Analysis with gprofiler2
#only for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)

#geneIDs <- deseq2_res[[1]]
geneIDs <- deseq2_res[!is.na(deseq2_res$pvalue),]

geneIDs_enriched <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange > 1)),1]
geneIDs_enriched <- geneIDs_enriched[!is.na(geneIDs_enriched)]

geneIDs_depleted <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange < -1)),1]
geneIDs_depleted <- geneIDs_depleted[!is.na(geneIDs_depleted)]

#save list of enriched and list of depleted genes
write.csv(geneIDs_enriched, file = "Studium/Master_Thesis/Ergebnisse_new_data/deseq2_enriched_genes.csv", quote = FALSE, row.names = FALSE)
write.csv(geneIDs_depleted, file = "Studium/Master_Thesis/Ergebnisse_new_data/deseq2_depleted_genes.csv", quote = FALSE, row.names = FALSE)

gostresult <- gost(query = list("genes_enriched" = geneIDs_enriched, "genes_depleted" = geneIDs_depleted), organism = "mmusculus", ordered_query = TRUE)

#save results
result_enrichment <- apply(gostresult$result, 2, as.character)
write.csv(result_enrichment, file = "Studium/Master_Thesis/Ergebnisse_new_data/gs_enrichment_analysis_results.csv", quote = FALSE, row.names = FALSE)

#visualization
interactive_plot <- gostplot(gostresult)
htmlwidgets::saveWidget(interactive_plot, file = "Studium/Master_Thesis/Ergebnisse_new_data/gs_enrichment.html")

png(filename = "Studium/Master_Thesis/Ergebnisse_new_data/gs_enrichment.png", height = 710, width = 1603)
gostplot(gostresult, interactive = FALSE)
dev.off()


#overrepresentation analysis with gprofiler2
#only for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
genes <- deseq2_res[!is.na(deseq2_res$pvalue),]

geneIDs_enriched <- genes[((genes$padj < 0.05) & (genes$log2FoldChange > 1)),1]
geneIDs_enriched <- geneIDs_enriched[!is.na(geneIDs_enriched)]

geneIDs_depleted <- genes[((genes$padj < 0.05) & (genes$log2FoldChange < -1)),1]
geneIDs_depleted <- geneIDs_depleted[!is.na(geneIDs_depleted)]

multi_gostresult <- gost(query = list("genes_enriched" = geneIDs_enriched, "genes_depleted" = geneIDs_depleted), organism = "mmusculus", ordered_query = FALSE)

#save results
result_overrepresentation <- apply(multi_gostresult$result, 2, as.character)
write.csv(result_overrepresentation, file = "Studium/Master_Thesis/Ergebnisse_new_data/overrepresentation_results.csv", quote = FALSE, row.names = FALSE)

#visualization
interactive_plot2 <- gostplot(multi_gostresult)
htmlwidgets::saveWidget(interactive_plot2, file = "Studium/Master_Thesis/Ergebnisse_new_data/overrepresentation.html")

png(filename = "Studium/Master_Thesis/Ergebnisse_new_data/overrepresentation.png", height = 710, width = 1603)
gostplot(multi_gostresult, interactive = FALSE)
dev.off()
