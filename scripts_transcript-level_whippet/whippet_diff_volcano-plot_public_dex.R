#load output data from whippet delta (public dex)  |  A = dex, B = control
file_diff <- read.table(gzfile("/nfs/proj/gr_splicing/whippet-files/output_delta/whippet-delta_public_dex.diff.gz"), header = TRUE)

#volcano plot with x = DeltaPsi, y = Probability (from filtered whippet-delta output)
png(filename = "/nfs/proj/gr_splicing/Ergebnisse_public_data/whippet-plots/public_dex_diff_volcano-plot.png", width = 1000, height = 800)
par(mar=c(5, 5, 4, 8), xpd=TRUE)

plot(file_diff$DeltaPsi, file_diff$Probability, pch = 20, main = "Liver (Dex vs Control) AS Events Volcano Plot", xlab = "DeltaPsi", ylab = "Probability", cex.main=1.5, cex.lab=1.3)
#add values to the x-axis
axis(1, at = c(-0.8, -0.7, -0.6, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8))
#different colors for AS types
AS_types <- c("TE", "TS", "CE", "AD", "AA", "RI", "AF", "AL")
AS_colors <- c("purple", "gold1", "skyblue", "blue", "deeppink", "green", "red", "brown")
for(i in 1:8){
  file_diff_x <- file_diff[(file_diff$Type == AS_types[i]),]
  points(file_diff_x$DeltaPsi, file_diff_x$Probability, pch=20, col=AS_colors[i])
}
AS_types <- replace(AS_types, AS_types=="CE", "ES")
legend(x = "topright", inset = c((-0.1), 0), fill = AS_colors, legend = AS_types, title = "AS Type", title.cex = 1.2)
#insert lines (x,y) for significant values (x = 0.1 & x = - 0.1; y = 0.9)
abline(v = 0.1, col="grey30", xpd=FALSE)
abline(v = -0.1, col="grey30", xpd=FALSE)
abline(h = 0.9, col="grey30", xpd=FALSE)

dev.off()
