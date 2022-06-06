library(SNPRelate)
library(gdsfmt)
library("devtools")
library(ggplot2)
#### Make the .gds format. Depending on the size of the vcf file it can take some tome. Run an R script to generate the gds file. See example of R script `gds-Format.R` and to run the job see example `Job-gdsFormat.sh`.
snpgdsVCF2GDS("newH-genolike1-greenjay-UNPLACED-10Kkb.vcf.recode.vcf", "genolike1-greenjay-UNPLACED-PRUNE.gds", method="biallelic.only")
genofile <- snpgdsOpen("genolike1-greenjay-UNPLACED-PRUNE.gds", readonly = FALSE)
#annotate population, sex and species category. The order of samples is the same as in your .vcf file header
samp.annot1 <- data.frame(pop.group = c("north","north","south","north","north","south","north","north","south","south","south","north","north","north","north","south","south","north","south","north"),
                         sub.species = c("north","north","south","north","north","south","north","north","south","south","south","north","north","north","north","south","south","north","south","north"))

add.gdsn(genofile, "sample.annot1", samp.annot1)
pop_code <- read.gdsn(index.gdsn(genofile, path="sample.annot1/pop.group"))

#### Run the principal Component Analysis (PCA) on genotypes
pca <- snpgdsPCA(genofile, autosome.only=FALSE)
#### Get variance proportion (%)
See the general values of variance
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#### Get sample id, pop info and sex or other info
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot1/pop.group"))
subspecies_code <- read.gdsn(index.gdsn(genofile, "sample.annot1/sub.species"))
#### Check your sample ID and population code
head(cbind(sample.id, pop_code))
#### Make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(pop_code)[match(pca$sample.id, sample.id)],
    sp = factor(subspecies_code)[match(pca$sample.id, sample.id)],
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
#### See the generated data frame. This eigenvect values can be use for regressions or other analysis
head(tab)
#### Indicate the colors to give to the groups sex, sp
colors <- c("#E69F00","#999999")
colors <- colors[as.numeric(tab$sp)]
#### Define shapes for populations 
shapes = c(0,1,2,3,4,5,6,7,8,9,10) 
shapes <- shapes[as.numeric(tab$pop)]
#### Define margins
par(xpd = T, mar = par()$mar + c(0,0,0,10))
#### Run the plot and specify legend
plot(tab$EV2, tab$EV1, col=colors,pch=shapes, xlab="eigenvector 2 (%)", ylab="eigenvector 1 (%)",cex=1.5,yaxt = "n",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
legend(0.2, 0.2, legend=levels(tab$pop),col=c("#E69F00","#999999"), cex=1.2, pch=c(0,1,2,3,4,5,6,7,8,9,10),
       lwd = 1, lty = 1)
axis(side = 1, labels = FALSE,cex=1.5,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
#### Draw the y-axis.
axis(side = 2,
#### Rotate the labels.
     las = 2,cex=1.5,
#### Adjust the label position.
     mgp = c(3, 0.75, 0))
par(mar=c(7, 5, 5, 6) + 0.1)
#### all egenvectors plots
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], labels=lbls, col=colors, pch=shapes)
#### to save the plots you could do
ggsave(plot=plot,"PCA-greenjays.eps",device=cairo_ps)
pdf("PCA-greenjays.pdf")
plot
dev.off()
setEPS()
postscript("PCA-greenjays2.eps")
plot
dev.off()

#### Relatedness 

ibd.robust <- snpgdsIBDKING(genofile, sample.id=sample.id, num.thread=2, autosome.only=FALSE)
data <- snpgdsIBDSelection(ibd.robust)
#### identify the thresholds
for (i in 1:nrow(data)) {
        if (data$kinship[i]>= 0.354) data$color[i] <- "gray"
        if (data$kinship[i]>= 0.177 & data$kinship[i]< 0.354) data$color[i] <- "red"
        if (data$kinship[i]>= 0.0884 & data$kinship[i]< 0.177) data$color[i] <- "blue"
        if (data$kinship[i]>= 0.0442 & data$kinship[i]< 0.0884) data$color[i] <- "green"
        if (data$kinship[i]>= 0 & data$kinship[i]< 0.0442) data$color[i] <- "cyan"
        if (data$kinship[i]< 0) data$color[i] <- "black"
}
#### save the plot in a pdf
pdf("kingship_greenjays.pdf")
#### Full plot
plot(data$IBS0, data$kinship, xlab="Proportion of Zero IBS",
     ylab="Estimated Kinship Coefficient (KING-robust)",
     pch = 21, col=1, bg = data$color, cex = 2)
abline(h=0.354, lty = 2)
abline(h=0.177, lty = 2)
abline(h=0.0884, lty = 2)
abline(h=0.0442, lty = 2)
abline(h=0, lty = 2)
legend("topright", inset = 0.02, legend = c("1st degree", "2nd degree", "3rd degree", "unrelated"),
       pch = 21, col = 1, pt.bg = c("red","blue","green", "black"), pt.cex = 2, bg = "white", y.intersp =2)
dev.off()
#### Dendrogram
set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2, autosome.only=FALSE))
rv <- snpgdsCutTree(ibs.hc)
pdf("dendogram-greenjays2.pdf")
plot(rv$dendrogram, leaflab="none", main="HapMap Phase II")
table(rv$samp.group)
rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))
#Create 10 groups.
plot(rv2$dendrogram, leaflab="none", main="HapMap Phase II")
#legend("topright", legend=levels(subspecies_code), col=1:nlevels(subspecies_code), pch=19, ncol=3)
legend("bottomright", legend=levels(pop_code), col=1:nlevels(pop_code), pch=19, ncol=3)
dev.off()
#############

