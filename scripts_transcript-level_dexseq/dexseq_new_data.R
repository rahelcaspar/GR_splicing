#DEXSeq (new data)
library("GenomicFeatures")
library("Rsamtools")
library("GenomicAlignments")
library("DEXSeq")


#convert gtf to txbd
txdb = makeTxDbFromGFF("/nfs/data/references/ensembl110_GRCm39/Mus_musculus.GRCm39.110.gtf")

#preparing the annotation
flattenedAnnotation = exonicParts(txdb, linked.to.single.gene.only = TRUE)
names(flattenedAnnotation) = sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)

#read new_data_Samples_Uebersicht
sample_conditions <- read.csv("/nfs/home/students/rahelcaspar/MA/Samples_Tabellen/new_data_Samples_Uebersicht.csv")

#read bam-files (one for each sample)
bamFiles <- c()
for(i in substring(sample_conditions$sample, first = 2)){
  bamFiles <- c(bamFiles, BamFile(paste0("/nfs/proj/gr_splicing/new_data/star_rsem/",i,".markdup.sorted.bam")))
  
}
bamFiles = BamFileList(bamFiles)

se = summarizeOverlaps(flattenedAnnotation, bamFiles, singleEnd=FALSE, fragments=TRUE, ignore.strand=FALSE)

#create DEXSeqDataSet
colData(se)$condition = factor(sample_conditions$condition)
colData(se)$libType = factor(rep("paired-end", times=length(sample_conditions$sample)))
colData(se)@rownames = substring(colData(se)@rownames, first = 1, last = 9)
dxd = DEXSeqDataSetFromSE(se, design= ~ sample + exon + condition:exon)
sample_names <- c()
for(i in 1:length(dxd@colData@listData$sample)){
  sample_names <- c(sample_names, substring(dxd@colData@listData$sample[[i]], first = 1, last = 9))
}
sample_names <- factor(sample_names)
dxd@colData@listData$sample <- sample_names

#run DEXSeq
dxresult = DEXSeq(dxd)

#save results
save("dxd", file = "/nfs/proj/gr_splicing/Ergebnisse_new_data/dexseq/new_data_dexseq-dataset2.RData")
save("dxresult", file = "/nfs/proj/gr_splicing/Ergebnisse_new_data/dexseq/new_data_dexseq_result2.RData")


#visualization
#create html report
DEXSeqHTML(dxresult, FDR = 0.05, color=c("#FF000080", "#0000FF80"), path = "/nfs/proj/gr_splicing/Ergebnisse_new_data/dexseq/DEXSeqReport")


#MA plot
png("/nfs/proj/gr_splicing/Ergebnisse_new_data/dexseq/plotMA.png")
plotMA(dxresult, alpha = 0.05)
dev.off()

