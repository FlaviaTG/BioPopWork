library(SNPRelate)
library(gdsfmt)
library("devtools")
library(ggplot2)
#### Make the .gds format. Depending on the size of the vcf file it can take some tome. Run an R script to generate the gds file. See example of R script `gds-Format.R` and to run the job see example `Job-gdsFormat.sh`.
snpgdsVCF2GDS("newH-genolike1-greenjay-UNPLACED-10Kkb.vcf.recode.vcf", "genolike1-greenjay-UNPLACED-PRUNE.gds", method="biallelic.only")
