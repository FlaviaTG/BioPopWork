setwd("~/Desktop/Desktop/Exercises/output")
library(ggplot2)  
library(reshape2)  
library(plyr)
#### plot reads
ReadCount <- read.table("MappedReads.txt", header=F)
sample <- read.table("Green-jay-20-sampleID.txt", header=F)
# a bit of polishing to remove total (not needed for plotting)  
ReadCountSmall <- data.frame(BAM = ReadCount$V1, Total = ReadCount$V2, Mapped = ReadCount$V3, Unmapped = ReadCount$V4, Good = ReadCount$V5)
# ggplot needs data in a specific layout  
MeltedReadCount = melt(ReadCountSmall, id=c('BAM'))  
names(MeltedReadCount) <- c('BAM', 'Mapping', 'Reads')
MeltedReadCount$SampleID <- sample$V2
# sort the data frame and add fraction to the data frame  
to_graph <- cbind(arrange(MeltedReadCount, SampleID))

# Now all we have to do is plot the data  
gp <- ggplot(data=to_graph, aes(x=SampleID, y=Reads, fill=Mapping)) +  
  geom_bar(stat="identity",position=position_dodge()) +
  scale_y_continuous(labels = comma) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
#
ggsave("read-mapped-green-jay.pdf", plot=gp)
###############
# plot coverage
cover <- read.table("coverage.txt", header=F)
sample <- read.table("Green-jay-20-sampleID.txt", header=F)
cover$SampleID <- sample$V2
names(cover) <- c('BAM', 'aver', 'equal', 'X', 'total', 'SampleID')
cover$realCOV <- cover$X/1.2
#
gp2 <- ggplot(data=cover, aes(x=reorder(SampleID,-realCOV), y=realCOV)) +  
  geom_bar(stat="identity",position=position_dodge()) +
  scale_y_continuous(labels = comma) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
ggsave("genome-coverage-green-jay.pdf", plot=gp2)

