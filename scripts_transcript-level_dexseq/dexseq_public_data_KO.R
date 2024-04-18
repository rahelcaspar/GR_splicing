#DEXSeq (public data - KO vs WT)
library("GenomicFeatures")
library("Rsamtools")
library("GenomicAlignments")
library("DEXSeq")


#convert gtf to txbd
txdb = makeTxDbFromGFF("/nfs/data/references/ensembl110_GRCm39/Mus_musculus.GRCm39.110.gtf")

#preparing the annotation
flattenedAnnotation = exonicParts(txdb, linked.to.single.gene.only = TRUE)
names(flattenedAnnotation) = sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)

#read public_data samples-table (knockout vs wildtype)
sample_conditions <- read.csv("/nfs/home/students/rahelcaspar/MA/Samples_Tabellen/samples_knockout.csv")

#read bam-files (one for each sample)
bamFiles <- c()
for(i in sample_conditions$sample){
  bamFiles <- c(bamFiles, BamFile(paste0("/nfs/proj/gr_splicing/star_rsem/",i,".markdup.sorted.bam")))
  
}
bamFiles = BamFileList(bamFiles)

se = summarizeOverlaps(flattenedAnnotation, bamFiles, singleEnd=FALSE, fragments=TRUE, ignore.strand=FALSE)

#create DEXSeqDataSet
colData(se)$condition = factor(sample_conditions$knockout)
colData(se)$condition <- factor(colData(se)$condition, levels=c("WT", "KO"))
colData(se)$libType = factor(rep("paired-end", times=length(sample_conditions$sample)))
colData(se)$timepoint = factor(sample_conditions$timepoint)
colData(se)@rownames = substring(colData(se)@rownames, first = 1, last = 10)
dxd = DEXSeqDataSetFromSE(se, design= ~ sample + exon + timepoint:exon + condition:exon)

#run DEXSeq 
dxresult = DEXSeq(dxd)

#save results
save("dxd", file = "/nfs/proj/gr_splicing/Ergebnisse_public_data/dexseq/public_data_dexseq-dataset_ko.RData")
save("dxresult", file = "/nfs/proj/gr_splicing/Ergebnisse_public_data/dexseq/public_data_dexseq_result_ko.RData")


#visualization
#create html report
DEXSeqHTML(dxresult, FDR = 0.05, color=c("#FF000080", "#0000FF80"), path = "/nfs/proj/gr_splicing/Ergebnisse_public_data/dexseq/DEXSeqReport_ko")


#MA plot
png("/nfs/proj/gr_splicing/Ergebnisse_public_data/dexseq/plotMA_ko.png")
plotMA(dxresult, alpha = 0.05)
dev.off()
