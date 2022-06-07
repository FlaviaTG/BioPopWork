#to plot GeneWin results all regions together
library(data.table)
library(qqman)
library(filesstrings)
library(ggplot2)
library(gghighlight)

#
files <- list.files(pattern = "*GeneWin.txt")
DT <- rbindlist(sapply(files, fread, simplify = FALSE),
                use.names = TRUE, idcol = "FileName", fill=TRUE)
#give the desire name for each chromosomes from a list of the chromosomes names in: files-MATCH-CHR.txt
names <- read.table("files-MATCH-CHR-sinZW.txt")
names(names) <- c('file','chr')
test <- merge(names, DT, by.x="file", by.y="FileName")
test2 <- test[order(test$chr),]
#
fstsubset<-test2[complete.cases(test2),]
SNP<-c(1:(nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)
#
#with qqman
quantile(mydf$MeanY, c(0.975, 0.995), na.rm = T)
# identify the 95% percentile
my_threshold <- quantile(mydf$MeanY, c(0.975), na.rm = T)
#make character vector of the SNPs to higligt that are on the 95% percentile
snpsOfInterest <- subset(mydf, MeanY > 0.4068114)
snpsOfInterest2 <- snpsOfInterest$SNP

write.table(snpsOfInterest, "Green-jay-Fst-genewin.txt",row.names=F, quote=F)
#save plot. Manhattan plot
#setEPS()
#postscript("Green-jay-Fst-genewin-quantiles-new.eps")
manhattan(mydf,chr="chr",bp="WindowStart",p="MeanY",logp=F, ylab = "Fst", highlight = snpsOfInterest2,chrlabs = c(1:41),suggestiveline = F, genomewideline = F, cex = 1, cex.axis = 0.8)
#dev.off()
#
#
mydf$size <- (mydf$WindowStop-mydf$WindowStart)
mydf$KB <- (mydf$size/1000)
mydf$chr <- as.factor(mydf$chr)

#plot the size of the spline window per chromosomes
#jpeg("spline-size-window-r.jpg", width = 1300, height = 400)
p <- ggplot(data=mydf, aes(x=chr, y=KB, fill=chr)) 
plot <- p + geom_boxplot(show.legend=F) +labs(x="Chromosomes", y="spline window size kb")+ theme_bw() +
  theme(axis.text.x=element_text(size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"), axis.line = element_line(colour = "black"),
	axis.text.y = element_text(size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
	axis.title.y = element_text(size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"), 
	axis.title=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
#dev.off()
