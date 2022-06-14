library(ggplot2)
library("cowplot")
library(reshape2)
library(gghighlight)
#
archivos <- list.files(pattern = "*txt")
for (file in archivos){
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=F, sep="\t")
    filename <- gsub(pattern = "\\.txt$", "",file)
    temp_dataset$ind <- filename
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  if (!exists("dataset")){
    dataset <- read.table(file, header=F, sep="\t")
    filename <- gsub(pattern = "\\.txt$", "",file)
    dataset$ind <- filename
  }
}
#the plot with all
library(data.table)
library(qqman)
library(filesstrings)
#
p<- ggplot(dataset, aes(x =V1, y = V2, group=ind,colour=ind))
#png("PSMC-all-individuals-minimus_bicknell.png")
p+ geom_step() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + scale_x_log10(oob = scales::squish_infinite) +theme(legend.position="none") + ylab("Effective population size(x10^4)") + xlab("Years")+ gghighlight(ind=='sample.CTCACG.minimus'|ind=='sample.CAGCTT.minimus'|ind=='sample.TCCTAG.minimus'|ind=='sample.AGATCC.minimus'|ind=='sample.AGCTTT.minimus'|ind=='sample.TAATGT.minimus'|ind=='sample.CAGTGT.minimus'|ind=='sample.GACTCA.minimus'|ind=='sample.AGACCA.minimus'|ind=='sample.TGAGCC.minimus'|ind=='sample.ACGGTC.minimus'|ind=='sample.TTCGAA.minimus'|ind=='sample.TATCAG.minimus'|ind=='sample.CCATGT.minimus'|ind=='sample.CATTTT.minimus', keep_scales = TRUE, unhighlighted_params = list(size = 2, colour = alpha("black", 0.4)))
###
ggsave("B-PSMC-bicknell-hl.png",plot=p, width=30, height=15, dpi=300)
ggsave("B-PSMC-minimus-hl.png",plot=p, width=30, height=15, dpi=300)
#higligting reference genomes
p+ geom_step() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + scale_x_log10(oob = scales::squish_infinite) +theme(legend.position="none") + ylab("Effective population size(x10^4)") + xlab("Years")+ gghighlight(ind=='sample.CGATGT.bicknell'|ind=='sample.CATTTT.minimus', keep_scales = TRUE)
