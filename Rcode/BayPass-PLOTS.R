library(corrplot)
library(ape)
library(geigen)
library(mvtnorm)
#
#source the baypass R functions (check PATH)
source("/scratch/08752/ftermig/programs/baypass_2.3/utils/baypass_utils.R")
#upload estimate of omega
omegaCV=as.matrix(read.table("Scale-var-Bio1_mat_omega.out"))
pop.names=c("1","2")
dimnames(omegaCV)=list(pop.names,pop.names)
#Compute and visualize the correlation matrix
pdf('correlation-matrix-BIOCLIM1-green-jay.pdf')
cor.mat=cov2cor(omegaCV)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
main=expression("Correlation map based on"~hat(Omega)))
#Visualize the correlation matrix as hierarchical clustering tree
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
plot(bta14.tree,type="p",main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
dev.off()
######Aux model plots##
covaux.snp.res=read.table("anacovaux_summary_betai.out",h=T)
covaux.snp.xtx=read.table("anacovaux_summary_pi_xtx.out",h=T)$M_XtX
graphics.off()
layout(matrix(1:3,3,1))
plot(covaux.snp.res$BF.dB.,xlab="SNP",ylab="BFmc (in dB)")
plot(covaux.snp.res$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"))
plot(covaux.snp.xtx,xlab="SNP",ylab="XtX corrected for SMS")
####The resulting estimates of the posterior mean of the each auxiliary variable δ i under both models (AUX model with no SNP spatial dependency
#and AUX model Bayes Factor, the underlying regression coefficients (posterior mean) and the corrected XtX might be plotted as follows:
covauxisb1.snp.res=read.table("anacovauxisb1_summary_betai.out",h=T)
graphics.off()
layout(matrix(1:2,2,1))
plot(covaux.snp.res$M_Delta,xlab="SNP",ylab=expression(delta[i]),main="AUX model")
plot(covauxisb1.snp.res$M_Delta,xlab="SNP",ylab=expression(delta[i]),
main="AUX model with isb=1")
####################################################
####visualizando resultados del analisis final BF###CALIBRADO
####################################################
covis2.snp.res=read.table("podcovisFINAL_summary_betai_reg.out",h=T)
pdf('BF-POD-BIOCLIM.pdf')
layout(matrix(1:3,3,1))
plot(covis2.snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
plot(covis2.snp.res$eBPis,xlab="SNP",ylab="eBPis")
plot(covis2.snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()
