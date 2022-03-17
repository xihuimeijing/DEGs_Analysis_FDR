#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
usage=function(){
  cat('Description:DEGs identification using mutiple methods. 
Author: Yumei Li, 03/14/2021.
Revision: Add Wilcoxon test using edgeR calculated CPM. Yumei Li, 04/08/2021
Revision: Add Limma-voom analysis. Yumei Li, 05/13/2021
Revision: Add NOISeq analysis. Yumei Li, 06/04/2021
Revision: Add dearseq analysis. Yumei Li, 08/27/2021
Revision: Perform Wilcoxon test using matrixTests package. Yumei Li, 08/30/2021\n',file=stderr())
  cat('Usage:DEGs.R -g=gene.readCount.matrix.tsv -c=Conditions.txt -o=result.tsv\n',file=stderr())
  cat('\t-g\t\tFILE\tGene read count matrix, first column is gene ID and first row is sample name.\n',file=stderr())
  cat('\t-c\t\tFILE\tConditions file with one line separated by tab for the samples.\n',file=stderr())
  cat('\t-f\t\tFLOAT\tFDR cutoff for DEGs [default: 0.05].\n',file=stderr())
  cat('\t-s\t\tSTRING\tSoftware for DEGs calling [d for DESeq2, e for edgeR, w for Wilcoxon test, l for Limma-voom, n for NOISeq, r for dearseq].\n',file=stderr())
  cat('\t-m\t\tSTRING\tTest method for dearseq [asymptotic(default) or permutation]\n',file=stderr())
  cat('\t-e\t\tLOGICAL\tIf use edgeR to normalize and filter raw read count for dearseq analysis [default:F]\n',file=stderr())
  cat('\t-o\t\tSTRING\tOutput file name.\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}

fdrcutoff=0.05
dearseqTest="asymptotic"
edgeRNorm=F

if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-g=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      expFile=arg.split[2]
    }
    if(grepl('^-c=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      conFile=arg.split[2]
    }
    if(grepl('^-f=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      fdrcutoff=as.numeric(arg.split[2])
    }
    if(grepl('^-s=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      software=arg.split[2]
    }
    if(grepl('^-m=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      dearseqTest=arg.split[2]
    }
    if(grepl('^-e=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      edgeRNorm=as.logical(arg.split[2])
    }
    if(grepl('^-o=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      out=arg.split[2]
    }
    if(grepl('^-h', arg)){
      usage()
    }
  }
}else{
  usage()
}

readCount<-read.table(file=expFile, header = T, row.names = 1, stringsAsFactors = F,check.names = F)
conditions<-read.table(file=conFile, header = F)
conditions<-factor(t(conditions))

### DESeq2
if(software == "d"){
  suppressMessages(library("DESeq2"))
  ddsCount = DESeqDataSetFromMatrix( readCount, DataFrame(conditions),~conditions)
  dds<-DESeq(ddsCount)
  res<-results(dds)
  res<-na.omit(res)
  res<-res[res$padj<fdrcutoff,]
  write.table(res, file=out,sep="\t", quote=F,row.names = T,col.names = T)
### edgeR
}else if(software == "e"){
  suppressMessages(library("edgeR"))
  y <- DGEList(counts=readCount,group=conditions)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~conditions)
  y <- estimateDisp(y,design) # Considering robust=T when you think your data has potential outlier issue.
  #perform quasi-likelihood F-tests:
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  res<-topTags(qlf, n=nrow(readCount), adjust.method = "BH", sort.by = "PValue", p.value = fdrcutoff)
  write.table(res, file=out,sep="\t", quote=F,row.names = T,col.names = T)
### Wilcoxon rank-sum test
}else if(software == "w"){
  suppressMessages(library("edgeR"))
  suppressMessages(library(matrixTests))
  y <- DGEList(counts=readCount,group=conditions)
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  ##Method 1
  # res=matrix(nrow=nrow(count_norm),ncol=2)
  # for(i in 1:nrow(count_norm)){
  #   data<-cbind(t(count_norm[i,]),conditions)
  #   res[i,1]=colnames(data)[1]
  #   colnames(data)[1]="gene"
  #   res[i,2]=wilcox.test(gene~conditions, data)$p.value
  # }
  # fdr=p.adjust(res[,2],method = "fdr")
  # outputRst=as.data.frame(cbind(res,fdr),stringsAsFactors=F)
  # outputRst[,2]=as.numeric(outputRst[,2])
  # outputRst[,3]=as.numeric(outputRst[,3])
  # outputRst=na.omit(outputRst)
  # write.table(outputRst[outputRst[,3]<fdrcutoff,], file=out,sep="\t", quote=F,row.names = F,col.names = c("Gene","p-value","FDR"))
  ##Method 2
  # pvalues <- sapply(1:nrow(count_norm),function(i){
  #   data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  #   p=wilcox.test(gene~conditions, data)$p.value
  #   return(p)
  # })
  # fdr=p.adjust(pvalues,method = "fdr")
  ##Method 3
  conditionsLevel<-levels(conditions)
  dataMem1=count_norm[,c(which(conditions==conditionsLevel[1]))]
  dataMem2=count_norm[,c(which(conditions==conditionsLevel[2]))]
  pvalue<-row_wilcoxon_twosample(dataMem1,dataMem2)$pvalue
  fdr=p.adjust(pvalue,method = "BH")
  outputRst=na.omit(as.data.frame(cbind(row.names(count_norm)[which(fdr<fdrcutoff)],pvalue[which(fdr<fdrcutoff)],fdr[which(fdr<fdrcutoff)]),stringsAsFactors=F))
  write.table(outputRst, file=out,sep="\t", quote=F,row.names = F,col.names = c("Gene","p-value","FDR"))
### Limma-voom
}else if(software == "l"){
  suppressMessages(library("edgeR"))
  y <- DGEList(counts=readCount,group=conditions)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~conditions)
  v <- voom(y, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  res=topTable(fit, coef=2, sort.by="P", number=Inf, p.value = fdrcutoff)
  write.table(res, file=out,sep="\t", quote=F,row.names = T,col.names = T)
### NOISeq
}else if(software == "n"){
  suppressMessages(library(NOISeq))
  data<-NOISeq::readData(data=readCount, factors=as.data.frame(conditions))
  res=noiseqbio(data, k=0.5, norm="tmm", factor="conditions",random.seed = 12345, filter = 1, cv.cutoff = 100, cpm = 1)
  outputRst=degenes(res, q=1-fdrcutoff, M=NULL)
  write.table(outputRst, file=out,sep="\t", quote=F,row.names = T,col.names = T)
### dearseq
}else if(software == "r"){
  library(dearseq)
  if(edgeRNorm==TRUE){
    suppressMessages(library("edgeR"))
    y <- DGEList(counts=readCount,group=conditions)
    keep <- filterByExpr(y)
    y <- y[keep,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y,method="TMM")
    count_norm=cpm(y, log=TRUE)
    conditions<-matrix(as.numeric(conditions),ncol=1)
    res=dearseq::dear_seq(exprmat=count_norm, variables2test=conditions, which_test=dearseqTest, parallel_comp=F, preprocessed=T)
    res<-res$pvals
  }else{
    conditions<-matrix(as.numeric(conditions),ncol=1)
    res=dearseq::dear_seq(exprmat=as.matrix(readCount), variables2test=conditions, which_test=dearseqTest, parallel_comp=F, preprocessed=F)
    res<-res$pvals
  }
  write.table(res[res$adjPval<fdrcutoff,], file=out,sep="\t", quote=F,row.names = T,col.names = T)
}else{
  cat('Please provide the correct software parameter: d, e, l, n, r, t or w.\n',file=stderr())
  q(save="no")
}
