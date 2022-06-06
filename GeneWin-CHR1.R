library(GenWin)
library(dplyr)
data <- read.table("1Kkb.windowed.weir.fst", header=T)
#subset one region/chromosome
Z <- subset(data,CHROM == "CM036346.1")
fstsubset<-Z[complete.cases(Z),]
test3 <- subset(fstsubset, WEIGHTED_FST >= 0)
Y <- test3$WEIGHTED_FST
map <- test3$BIN_START
#give variables
pdf("spline-Wstat-CM036346.1.pdf")
ZW <- splineAnalyze(Y, map, smoothness = 100,
plotRaw = T, plotWindows = T, method = 4)
dev.off()
#write the table to plot all together in a manhattan-like plot
R <- as.data.frame(ZW$windowData)
write.table(R,"CM036346.1-Fst-GeneWin.txt", quote=F,row=F)
