plot_bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,name_highlighted=F,add_text=F)
{
  if (is.character(res))
    res=read.table(res)
  
  colfstat=5
  colq=colfstat-2
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0.05)
    ok_outliers=FALSE;
  
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
  }
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
  
  return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}
#
library("boa")
##########3
pdf("Bayescan-original-plot.pdf")
mydata2=read.table("newH-genolike2-greenjay_ZW-test.NOmissing.THIN.G_fst.txt",colClasses="numeric", header=TRUE)
plot_bayescan(mydata2,FDR=0.05)
results <-plot_bayescan("newH-genolike2-greenjay_ZW-test.NOmissing.THIN.G_fst.txt",0, FDR=0.05)
results$outliers
results$nb_outliers
###Fst distribution######
mydata=read.table("newH-genolike2-greenjay_ZW-test.NOmissing.THIN.G.sel",colClasses="numeric", header=TRUE)
plot(density(mydata[["Fst1"]]),xlab="Fst1",main=paste("Fst1","posterior distribution"))
plot(density(mydata2[["alpha"]]),xlab="alpha",main=paste("alpha","posterior distribution"))
plot(density(mydata[["logL"]]),xlab="logL",main=paste("logL" ,"posterior distribution"))
boa.hpd(mydata[["Fst1"]],0.05)
boa.hpd(mydata[["Fst2"]],0.05)
###########
plot(density(mydata2[["alpha"]]),xlab="alpha",main=paste("alpha","posterior distribution"))
plot(density(mydata2[["prob"]]),xlab="prob",main=paste("prob","posterior distribution"))
plot(density(mydata2[["log10.PO."]]),xlab="log10.PO.",main=paste("log10.PO.","posterior distribution"))
###
dev.off()
