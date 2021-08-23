#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
usage=function(){
  cat('Description:This script will perform functional enrichment analysis for the input gene list.\n',file=stderr())
  cat('Usage:geneEnrichAnalysis.R peak.bed annotationOut.tsv annoPie.pdf GO.dotplot.pdf\n',file=stderr())
  cat('\t-i\t\tFILE\tInput gene list (can read from STDIN)\n',file=stderr())
  cat('\t-k\t\tSTRING\tKeytype of input gene ["SYMBOL"(default),ENSEMBL","ENTREZID"]\n',file=stderr())
  cat('\t-s\t\tFILE\tThe genome assembly for the input species[mm10,hg19,hg38]\n',file=stderr())
  cat('\t-d\t\tLOGIC\tDatabase to perform enrichment analysis [GO (default), KEGG]\n',file=stderr())
  cat('\t-n\t\tINT\tThe number of ploted enriched pathways/GO terms.[default: 20]\n',file=stderr())
  cat('\t-q\t\tINT\tThe q-value cutoff.[default: 0.05]\n',file=stderr())
  cat('\t-p\t\tSTRING\tPlot type ["b" for barplot (default), "d" for dotplot, "a" for both]\n',file=stderr())
  cat('\t-o\t\tFILE\tOutput prefix for all the output files\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}
if(length(args)==0 || args[1]=="-h"){
  usage()
}
database="GO"
plotN=20
keytype="SYMBOL"
plottype="b"
qvalue=0.05
if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      geneFile=arg.split[2]
    }
    if(grepl('^-k=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      keytype=arg.split[2]
    }
    if(grepl('^-s=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      species=arg.split[2]
    }
    if(grepl('^-d=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      database=arg.split[2]
    }
    if(grepl('^-n=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      plotN=as.numeric(arg.split[2])
    }
    if(grepl('^-q=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      qvalue=as.numeric(arg.split[2])
    }
    if(grepl('^-p=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      plottype=arg.split[2]
    }
    if(grepl('^-o=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      prefix=arg.split[2]
    }
    if(grepl('^-h',arg)){usage()}
  }
}else{
  usage()
}
suppressMessages(library("ggplot2"))
suppressMessages(library("clusterProfiler"))
#1, Read gene list
if(exists("geneFile")){
  geneList<-read.table(file=geneFile, header = F, stringsAsFactors = F)
}else{
  geneList<-read.table(file('stdin'), header = F, stringsAsFactors = F)
}
#2. Get OrgDb
if(species == "mm10"){
  speciesKegg="mmu"
  suppressMessages(library("org.Mm.eg.db"))
  orgdb<-org.Mm.eg.db
}else if(species == "hg19"){
  speciesKegg="hsa"
  suppressMessages(library("org.Hs.eg.db"))
  orgdb<-org.Hs.eg.db
}else if(species == "hg38"){
  speciesKegg="hsa"
  suppressMessages(library("org.Hs.eg.db"))
  orgdb<-org.Hs.eg.db
}else{
  cat('Cannot process the input species',file=stderr())
  q(save="no")
}

#3. Run enrichment analysis
if (database == "GO"){
  if(keytype == "SYMBOL"){
    erst<-enrichGO(gene = geneList$V1,
                   keyType       = keytype,
                   OrgDb         = orgdb,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = qvalue,)
  }else{
    erst<-enrichGO(gene = geneList$V1,
                   OrgDb         = orgdb,
                   keyType       = keytype,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = qvalue,
                   readable      = TRUE)
  }
}else if (database == "KEGG"){
  if(keytype == "SYMBOL"){
    geneID<-bitr(geneList$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
  }else if(keytype == "ENSEMBL"){
    geneID<-bitr(geneList$V1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
  }else if(keytype == "ENTREZID"){
    geneID=geneList$V1
  }
  erst<-enrichKEGG(gene         = geneID,
                   organism     = speciesKegg,
                   keyType      = "ncbi-geneid",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = qvalue)
}else{
  cat('Please set the right database parameter\n',file=stderr())
  q(save="no")
}

#4. Custom plot
EGO<-as.data.frame(erst)
write.table(EGO,file=paste(prefix,".enrichment.tsv",sep=""),quote=F,sep="\t",row.names = F)
EGO$p.adjust<-(-log10(EGO$p.adjust))
ratio<-function(x){as.numeric(strsplit(x,'/')[[1]])[1]/as.numeric(strsplit(x,'/')[[1]])[2]*100}
EGO$GeneRatio<-mapply(FUN=ratio,EGO$GeneRatio,USE.NAMES=F)
subset<-EGO[c(1:plotN),]
subset$Description<-factor(subset$Description,levels=subset$Description)
pdf(file=paste(prefix,".enrichment.pdf",sep=""), width = 6, height = 3)
if(plottype == "b"){
  p<-ggplot(subset,aes(x=p.adjust,y=Description,fill=p.adjust))+geom_bar(stat='identity')+scale_fill_gradient(high = "#132B43", low = "#56B1F7")+
    labs(x="-log10(p.adjust)")+theme(legend.position = "none",axis.text = element_text(color="black"),text=element_text(size=12))
  print(p)
}else if(plottype == "d"){
  p<-ggplot(subset,aes(x=p.adjust,y=Description,color=GeneRatio,size=Count))+geom_point(stat='identity')+
    labs(x="-log10(p.adjust)")+scale_colour_gradient(high = "#132B43", low = "#56B1F7")+theme(axis.text = element_text(color="black"),text=element_text(size=12))
  print(p)
}else if(plottype == "a"){
  p1<-ggplot(subset,aes(x=p.adjust,y=Description,fill=p.adjust))+geom_bar(stat='identity')+
    labs(x="-log10(p.adjust)")+scale_fill_gradient(high = "#132B43", low = "#56B1F7")+theme(legend.position = "none",axis.text = element_text(color="black"),text=element_text(size=15))
  p2<-ggplot(subset,aes(x=p.adjust,y=Description,color=GeneRatio,size=Count))+geom_point(stat='identity')+
    labs(x="-log10(p.adjust)")+scale_colour_gradient(high = "#132B43", low = "#56B1F7")+theme(axis.text = element_text(color="black"),text=element_text(size=12))
  print(p1)
  print(p2)
}
dev.off()
