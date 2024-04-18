library(ggplot2)
#visualize whippet results - public data - dex

#input: /nfs/proj/gr_splicing/whippet-files/output_quant/[new_data/public_ko_wt/public_dex]/[sample-ID].psi.gz
#columns: Type, Psi  (tab-separated)
#input is per sample: (1) save separately, (2) combine samples into conditions (ko/wt, dex/control)
AS_types <- c("AF", "TS", "RI", "AA", "CE", "AD", "TE", "AL")

## read public_dex samples-table: 
dataset_samples <- read.csv("MA/Samples_Tabellen/samples_dex.csv")

#get the AS information from the psi.gz files
sample_events <- list()
cond_stats <- list(dex = c(), control = c())
cond_stats[["dex"]] <- list(CE = c(), AA = c(), AD = c(), RI = c(), TS = c(), TE = c(), AF = c(), AL = c())
cond_stats[["control"]] <- list(CE = c(), AA = c(), AD = c(), RI = c(), TS = c(), TE = c(), AF = c(), AL = c())
for(sampleX in dataset_samples[,"sample"]){
  print(sampleX)
  file_sampleX <- read.csv(gzfile(paste0("/nfs/proj/gr_splicing/whippet-files/output_quant/public_dex/",sampleX,".psi.gz")), sep = "\t")
  file_sampleX <- file_sampleX[,c("Type","Psi")]
  
  #remove NAs, PSI == 0 & PSI == 1 entries
  file_sampleX <- file_sampleX[complete.cases(file_sampleX),]
  file_sampleX[,"Psi"] <- replace(file_sampleX[,"Psi"], (file_sampleX[,"Psi"]==0), NA)
  file_sampleX <- file_sampleX[complete.cases(file_sampleX),]
  file_sampleX[,"Psi"] <- replace(file_sampleX[,"Psi"], (file_sampleX[,"Psi"]==1), NA)
  file_sampleX <- file_sampleX[complete.cases(file_sampleX),]
  
  sample_events[[sampleX]] = file_sampleX
  
  #collect AS-type PSI values per condition
  condX <- dataset_samples[dataset_samples$sample==sampleX,"dex"]
  for(ASx in AS_types){
    cond_stats[[condX]][[ASx]] <- c(cond_stats[[condX]][[ASx]], file_sampleX[file_sampleX$Type==ASx,"Psi"])
  }
}

#per condition: How often which node type (AS type)? (grouped Bar Plot)
AS_Type <- rep(AS_types, each = 2)
condition <- rep(c("dex", "control"), 8)
value <- c()
for(i in 1:16){
  value <- c(value, length(cond_stats[[condition[i]]][[AS_Type[i]]]))
}
data <- data.frame(AS_Type,condition,value)
data$AS_Type <- replace(data$AS_Type, data$AS_Type=="CE", "ES")
AS_types <- replace(AS_types, AS_types=="CE", "ES")
data$AS_Type <- factor(data$AS_Type, levels = AS_types)

#grouped barplot 
png("/nfs/proj/gr_splicing/Ergebnisse_public_data/whippet-plots/public_dex_AS_types_barplot.png", width = 700)
ggplot(data, aes(fill=condition, y=value, x=AS_Type)) + 
  geom_bar(position="dodge", stat="identity") + 
  labs(title = "Frequency of AS Types in Liver", x = "Alternative Splicing Type", y = "Frecuency") + 
  labs(fill = "Condition") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#per condition: per node type: which distribution of PSI values? (grouped Violin Plot)
psi_AS_Type <- c()
psi_condition <- c()
psi <- c()
for(i in 1:16){
  psi <- c(psi, cond_stats[[condition[i]]][[AS_Type[i]]])
  psi_AS_Type <- c(psi_AS_Type, rep(AS_Type[i], length(cond_stats[[condition[i]]][[AS_Type[i]]])))
  psi_condition <- c(psi_condition, rep(condition[i], length(cond_stats[[condition[i]]][[AS_Type[i]]])))
}
psi_data <- data.frame(psi_AS_Type,psi_condition,psi)
psi_data$psi_AS_Type <- replace(psi_data$psi_AS_Type, psi_data$psi_AS_Type=="CE", "ES")
psi_data$psi_AS_Type <- factor(psi_data$psi_AS_Type, levels = AS_types)

#grouped violin plot
png("/nfs/proj/gr_splicing/Ergebnisse_public_data/whippet-plots/public_dex_AS_types_violinplot.png", width = 1500)
ggplot(psi_data, aes(fill=psi_condition, y=psi, x=psi_AS_Type)) + 
  geom_violin(position="dodge") +
  labs(title = "Distribution of PSI Values in Liver", x = "Alternative Splicing Type", y = "PSI") + 
  labs(fill = "Condition") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
