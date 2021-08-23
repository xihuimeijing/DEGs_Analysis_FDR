#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
usage=function(){
  cat('Description:DEGs identification using DESeq2, edgeR, Clipper. 
Author: Yumei Li, 03/14/2021.
Revision: Add Wilcoxon test using edgeR calculated CPM. Yumei Li, 04/08/2021
Revision: Add Limma-voom analysis. Yumei Li, 05/13/2021
Revision: Add NOISeq analysis. Yumei Li, 06/04/2021\n',file=stderr())
  cat('Usage:DEGs.R -g=gene.readCount.matrix.tsv -c=Conditions.txt -o=result.tsv\n',file=stderr())
  cat('\t-g\t\tFILE\tGene read count matrix, first column is gene ID and first row is sample name.\n',file=stderr())
  cat('\t-c\t\tFILE\tConditions file with one line separated by tab for the samples.\n',file=stderr())
  cat('\t-f\t\tFLOAT\tFDR cutoff for DEGs [default: 0.05].\n',file=stderr())
  cat('\t-s\t\tSTRING\tSoftware for DEGs calling [d for DESeq2, e for edgeR, c for Clipper, w for Wilcoxon test, l for Limma-voom, n for NOISeq].\n',file=stderr())
  cat('\t-t\t\tSTRING\tContrast score for Clipper [max(default) or diff]\n',file=stderr())
  cat('\t-o\t\tSTRING\tOutput file name.\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}

fdrcutoff=0.05
contrastScore="max"
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
    if(grepl('^-t=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      contrastScore=arg.split[2]
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

if(software == "d"){
  suppressMessages(library("DESeq2"))
  ddsCount = DESeqDataSetFromMatrix( readCount, DataFrame(conditions),~conditions)
  dds<-DESeq(ddsCount)
  res<-results(dds)
  res<-na.omit(res)
  res<-res[res$padj<fdrcutoff,]
  write.table(res, file=out,sep="\t", quote=F,row.names = T,col.names = T)
}else if(software == "e"){
  suppressMessages(library("edgeR"))
  y <- DGEList(counts=readCount,group=conditions)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~conditions)
  y <- estimateDisp(y,design)
  #perform quasi-likelihood F-tests:
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  res<-topTags(qlf, n=nrow(readCount), adjust.method = "BH", sort.by = "PValue", p.value = fdrcutoff)
  write.table(res, file=out,sep="\t", quote=F,row.names = T,col.names = T)
}else if(software == "c"){
  suppressMessages(library("edgeR"))
  suppressMessages(library(Clipper))
  y <- DGEList(counts=readCount,group=conditions)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  count_norm=cpm(y)
  conLevels<-levels(conditions)
  re_clipper <- Clipper(score.exp=log2(count_norm[,which(conditions==conLevels[1])]+1), score.back=log2(count_norm[,which(conditions==conLevels[2])]+1), analysis = "differential", contrast.score=contrastScore, FDR = fdrcutoff)
  res=cbind(row.names(count_norm)[re_clipper$discoveries[[1]]],re_clipper$contrast.score.value[re_clipper$discoveries[[1]]])
  write.table(res, file=out,sep="\t", quote=F,row.names = F,col.names = c("Gene","ContrastScore"))
}else if(software == "w"){
  suppressMessages(library("edgeR"))
  y <- DGEList(counts=readCount,group=conditions)
  keep <- filterByExpr(y)
  y <- y[keep,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y,method="TMM")
  count_norm=cpm(y)
  count_norm<-as.data.frame(count_norm)
  res=matrix(nrow=nrow(count_norm),ncol=2)
  for(i in 1:nrow(count_norm)){
    data<-cbind(t(count_norm[i,]),conditions)
    res[i,1]=colnames(data)[1]
    colnames(data)[1]="gene"
    res[i,2]=wilcox.test(gene~conditions, data)$p.value
  }
  fdr=p.adjust(res[,2],method = "fdr")
  outputRst=as.data.frame(cbind(res,fdr),stringsAsFactors=F)
  outputRst[,2]=as.numeric(outputRst[,2])
  outputRst[,3]=as.numeric(outputRst[,3])
  outputRst=na.omit(outputRst)
  write.table(outputRst[outputRst[,3]<fdrcutoff,], file=out,sep="\t", quote=F,row.names = F,col.names = c("Gene","p-value","FDR"))
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
}else if(software == "n"){
  suppressMessages(library(NOISeq))
  data<-NOISeq::readData(data=readCount, factors=as.data.frame(conditions))
  res=noiseqbio(data, k=0.5, norm="tmm", factor="conditions",random.seed = 12345, filter = 1, cv.cutoff = 100, cpm = 1)
  outputRst=degenes(res, q=1-fdrcutoff, M=NULL)
  write.table(outputRst, file=out,sep="\t", quote=F,row.names = T,col.names = T)
}else{
  cat('Please provide the correct software parameter: d, e, c or w.\n',file=stderr())
  q(save="no")
}


