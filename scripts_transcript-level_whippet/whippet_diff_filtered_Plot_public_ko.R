library(ggplot2)

#load output data from whippet delta (public ko)
file_diff <- read.table(gzfile("/nfs/proj/gr_splicing/whippet-files/output_delta/whippet-delta_public_ko_wt.diff.gz"), header = TRUE)

#filter whippet-delta output (public ko) for |DeltaPsi| > 0.1 and Probability > 0.9
file_diff <- file_diff[(abs(file_diff$DeltaPsi) > 0.1),]
file_diff <- file_diff[(file_diff$Probability > 0.9),]

#save in file
write.table(file_diff, file = "/nfs/proj/gr_splicing/whippet-files/output_delta/whippet-delta-filtered_public_ko_wt.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

#check distribution of AS types in filtered data (grouped barplot)
AS_types <- c("AF", "TS", "RI", "AA", "CE", "AD", "TE", "AL")
#columns: Type, Psi_A, Psi_B, DeltaPsi  (tab-separated)
# A = KO, B = WT


#remove NAs, PSI == 0 & PSI == 1 entries from Psi_A
file_diff <- file_diff[complete.cases(file_diff),]
file_diff[,"Psi_A"] <- replace(file_diff[,"Psi_A"], (file_diff[,"Psi_A"]==0), NA)
file_diff <- file_diff[complete.cases(file_diff),]
file_diff[,"Psi_A"] <- replace(file_diff[,"Psi_A"], (file_diff[,"Psi_A"]==1), NA)
file_diff <- file_diff[complete.cases(file_diff),]

#remove NAs, PSI == 0 & PSI == 1 entries from Psi_B
file_diff <- file_diff[complete.cases(file_diff),]
file_diff[,"Psi_B"] <- replace(file_diff[,"Psi_B"], (file_diff[,"Psi_B"]==0), NA)
file_diff <- file_diff[complete.cases(file_diff),]
file_diff[,"Psi_B"] <- replace(file_diff[,"Psi_B"], (file_diff[,"Psi_B"]==1), NA)
file_diff <- file_diff[complete.cases(file_diff),]


# A = KO, B = WT
#per condition: How often which node type (AS type)? (grouped Bar Plot)
AS_Type <- rep(AS_types, each = 2)
condition <- rep(c("KO", "WT"), 8)
value <- c()
for(i in 1:16){
  value <- c(value, length(which(file_diff$Type==AS_Type[i])))
}
data <- data.frame(AS_Type,condition,value)
data$AS_Type <- replace(data$AS_Type, data$AS_Type=="CE", "ES")
AS_types <- replace(AS_types, AS_types=="CE", "ES")
data$AS_Type <- factor(data$AS_Type, levels = AS_types)

#grouped barplot 
png("/nfs/proj/gr_splicing/Ergebnisse_public_data/whippet-plots/public_ko_diff_AS_types_barplot.png", width = 700)
ggplot(data, aes(fill=condition, y=value, x=AS_Type)) + 
  geom_bar(position="dodge", stat="identity") + 
  labs(title = "Frequency of AS Types in Liver", x = "Alternative Splicing Type", y = "Frecuency") + 
  labs(subtitle = "from AS Events with significant Change between Conditions (Mean over all Samples)") +
  labs(fill = "Condition") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.subtitle = element_text(hjust = 0.5))
dev.off()

