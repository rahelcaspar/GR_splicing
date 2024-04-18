library(VennDiagram)
#public dex & new data - (dex vs control) - venn diagram

#public dex vs control ZT00
#get deseq2 results
deseq2_res_dex1 <- read.csv("Studium/Master_Thesis/Ergebnisse_public_data/deseq2/public_data_res_dex_ZT00_deseq2.csv")

#filter for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
geneIDs00 <- deseq2_res_dex1[!is.na(deseq2_res_dex1$pvalue),]

geneIDs_enriched00 <- geneIDs00[((geneIDs00$padj < 0.05) & (geneIDs00$log2FoldChange > 1)),1]
geneIDs_enriched00 <- geneIDs_enriched00[!is.na(geneIDs_enriched00)]
write.csv(as.character(geneIDs_enriched00), file = "Studium/Master_Thesis/Ergebnisse_public_data/public_dex_venn_diagram/public_dex_geneIDs_enriched00.csv", quote = FALSE, row.names = FALSE)

geneIDs_depleted00 <- geneIDs00[((geneIDs00$padj < 0.05) & (geneIDs00$log2FoldChange < -1)),1]
geneIDs_depleted00 <- geneIDs_depleted00[!is.na(geneIDs_depleted00)]
write.csv(as.character(geneIDs_depleted00), file = "Studium/Master_Thesis/Ergebnisse_public_data/public_dex_venn_diagram/public_dex_geneIDs_depleted00.csv", quote = FALSE, row.names = FALSE)


#public dex vs control ZT12
#get deseq2 results
deseq2_res_dex2 <- read.csv("Studium/Master_Thesis/Ergebnisse_public_data/deseq2/public_data_res_dex_ZT12_deseq2.csv")

#filter for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
geneIDs12 <- deseq2_res_dex2[!is.na(deseq2_res_dex2$pvalue),]

geneIDs_enriched12 <- geneIDs12[((geneIDs12$padj < 0.05) & (geneIDs12$log2FoldChange > 1)),1]
geneIDs_enriched12 <- geneIDs_enriched12[!is.na(geneIDs_enriched12)]
write.csv(as.character(geneIDs_enriched12), file = "Studium/Master_Thesis/Ergebnisse_public_data/public_dex_venn_diagram/public_dex_geneIDs_enriched12.csv", quote = FALSE, row.names = FALSE)

geneIDs_depleted12 <- geneIDs12[((geneIDs12$padj < 0.05) & (geneIDs12$log2FoldChange < -1)),1]
geneIDs_depleted12 <- geneIDs_depleted12[!is.na(geneIDs_depleted12)]
write.csv(as.character(geneIDs_depleted12), file = "Studium/Master_Thesis/Ergebnisse_public_data/public_dex_venn_diagram/public_dex_geneIDs_depleted12.csv", quote = FALSE, row.names = FALSE)


#public dex vs control (both timepoints)
#get deseq2 results
deseq2_res_dex <- read.csv("Studium/Master_Thesis/Ergebnisse_public_data/deseq2/public_data_res_dex_deseq2.csv")

#filter for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
geneIDs <- deseq2_res_dex[!is.na(deseq2_res_dex$pvalue),]

geneIDs_enriched <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange > 1)),1]
geneIDs_enriched <- geneIDs_enriched[!is.na(geneIDs_enriched)]
write.csv(as.character(geneIDs_enriched), file = "Studium/Master_Thesis/Ergebnisse_public_data/public_dex_venn_diagram/public_dex_geneIDs_enriched.csv", quote = FALSE, row.names = FALSE)

geneIDs_depleted <- geneIDs[((geneIDs$padj < 0.05) & (geneIDs$log2FoldChange < -1)),1]
geneIDs_depleted <- geneIDs_depleted[!is.na(geneIDs_depleted)]
write.csv(as.character(geneIDs_depleted), file = "Studium/Master_Thesis/Ergebnisse_public_data/public_dex_venn_diagram/public_dex_geneIDs_depleted.csv", quote = FALSE, row.names = FALSE)


#new data (hepatocytes): dex vs control 
#get deseq2 results
deseq2_res_dex_h <- read.csv("Studium/Master_Thesis/Ergebnisse_new_data/deseq2/new_data_res_deseq2.csv")

#filter for genes with padj < 0.05 and log2FoldChange > 1 (enriched) or log2FoldChange < 1 (depleted)
geneIDs_h <- deseq2_res_dex_h[!is.na(deseq2_res_dex_h$pvalue),]

geneIDs_enriched_h <- geneIDs_h[((geneIDs_h$padj < 0.05) & (geneIDs_h$log2FoldChange > 1)),1]
geneIDs_enriched_h <- geneIDs_enriched_h[!is.na(geneIDs_enriched_h)]
write.csv(as.character(geneIDs_enriched_h), file = "Studium/Master_Thesis/Ergebnisse_new_data/venn_diagram/new_data_geneIDs_enriched.csv", quote = FALSE, row.names = FALSE)

geneIDs_depleted_h <- geneIDs_h[((geneIDs_h$padj < 0.05) & (geneIDs_h$log2FoldChange < -1)),1]
geneIDs_depleted_h <- geneIDs_depleted_h[!is.na(geneIDs_depleted_h)]
write.csv(as.character(geneIDs_depleted_h), file = "Studium/Master_Thesis/Ergebnisse_new_data/venn_diagram/new_data_geneIDs_depleted.csv", quote = FALSE, row.names = FALSE)


#venn diagram for upregulated genes public dex (ZT00, ZT12, both timepoints) & new data (hepatocytes)
venn.diagram(list(geneIDs_enriched00, geneIDs_enriched12, geneIDs_enriched, geneIDs_enriched_h), imagetype = "png", 
             category.names = c("ZT00","ZT12","ZT00 & ZT12", "Hepatocytes"), disable.logging = TRUE, 
             filename = "Studium/Master_Thesis/Ergebnisse_new_data/venn_diagram/public+new_data_dex_venn_up.png", 
             main = "Dex vs Control - Upregulated Genes", main.fontface = "bold", main.fontfamily = "sans", 
             main.pos = c(0.5, 0.95), 
             fill = c("aquamarine", "sandybrown", "lightskyblue1", "plum2"), 
             #circles
             lwd = 0.5,  
             #numbers
             cex = .8, fontface = "bold", fontfamily = "sans", 
             #set names
             cat.cex = 0.8, cat.fontface = "bold", cat.default.pos = "outer", cat.fontfamily = "sans", 
             cat.pos = c(-30,30,-25,23), cat.dist = c(0.035,0.035,0.03,0.03))


#venn diagram for downregulated genes public dex (ZT00, ZT12, both timepoints) & new data (hepatocytes)
venn.diagram(list(geneIDs_depleted00, geneIDs_depleted12, geneIDs_depleted, geneIDs_depleted_h), imagetype = "png", 
             category.names = c("ZT00","ZT12","ZT00 & ZT12", "Hepatocytes"), disable.logging = TRUE, 
             filename = "Studium/Master_Thesis/Ergebnisse_new_data/venn_diagram/public+new_data_dex_venn_down.png", 
             main = "Dex vs Control - Downregulated Genes", main.fontface = "bold", main.fontfamily = "sans", 
             main.pos = c(0.5, 0.95), 
             fill = c("aquamarine", "sandybrown", "lightskyblue1", "plum2"), 
             #circles
             lwd = 0.5, 
             #numbers
             cex = .8, fontface = "bold", fontfamily = "sans", 
             #set names
             cat.cex = 0.8, cat.fontface = "bold", cat.default.pos = "outer", cat.fontfamily = "sans", 
             cat.pos = c(-30,30,-25,23), cat.dist = c(0.035,0.035,0.03,0.03))
