vocano_plot<-function(x,out="qqplot.pdf"){
  library(calibrate)
  diffdata = read.table(x,header = T,stringsAsFactors = FALSE)
  diffdata=diffdata[!is.na(diffdata$pvalue),]
  diffdata$logpvalue<- -log10(diffdata$pvalue)
  diffdata$logpvalue<- ifelse(diffdata$logpvalue>130, 130, diffdata$logpvalue)
  gene=row.names(diffdata)
  res=cbind(gene,diffdata)
  out_fn = paste(out,".pdf",sep="")
  pdf(out_fn)
  with(res, plot(log2FoldChange, logpvalue, pch=20, main="Volcano plot",col="gray",xlim=c(-10,10),ylim=c(0,130)))
  with(subset(res, padj<.05 & log2FoldChange> 1), points(log2FoldChange, logpvalue, pch=20, col="red4"))
  with(subset(res, padj<.05 & log2FoldChange< -1), points(log2FoldChange, logpvalue, pch=20, col="blue"))
  dev.off()
}
