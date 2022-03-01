########## Load libraries ########## 
library(ggplot2)
library(plyr)
library(dplyr)
library(grid)
library(reshape2)
library(gridExtra)
library(scales)
library(cowplot)

########## Prepare data ##########
# Use setwd() to change the working directory to the directory containing the data for plot
samples=c(list.dirs("./TCGA",recursive = F),list.dirs("./GTEx",recursive = F))
GTEx_gof=readRDS(file="gtex_goodness_of_fit.rds")
TCGA_gof=readRDS(file="tcga_goodness_of_fit.rds")
gof=c(TCGA_gof, GTEx_gof)
samples<-samples[c(1:10,12,11)]

########## Define functions ##########
log2_pseudo1<-trans_new(name="log2_pseudo1",transform = function(x){log2(x+1)},inverse = function(x){2^x-1})
log10_pseudo<-trans_new(name="log10_pseudo",transform = function(x){log10(x+1e-6)},inverse = function(x){10^x-1e-6})

########## Set theme for the plots ##########
theme_update(axis.text = element_text(color="black"),text=element_text(size=12), plot.title=element_text(hjust=0.5), 
             panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA), 
             panel.grid.major = element_blank(), panel.grid.minor = element_blank())

######### Fig. S6-S17 ##########
### Set scale factor for figure A and C
TCGA=604.83
GTEx=562.00
for(i in 1:12){
  ### Read DEGs based on real data
  deseq2=read.delim(file=list.files(samples[i],pattern = "DESeq2.rst.tsv$|DESeq2.0.05.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
  deseq2=deseq2[deseq2$padj<0.01,]
  edger=read.delim(file=list.files(samples[i],pattern = "edgeR.rst.tsv$|edgeR.0.05.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
  edger=edger[edger$FDR<0.01,]
  wilcoxon=read.delim(file=list.files(samples[i],pattern = "WilcoxonTest.rst.tsv$|WilcoxonTest.0.05.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
  wilcoxon=wilcoxon[wilcoxon$FDR<0.01,]
  limma=read.delim(file=list.files(samples[i],pattern = "limma.rst.tsv$|limma.0.05.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
  limma=limma[limma$adj.P.Val<0.01,]
  noiseq=read.delim(file=list.files(samples[i],pattern = paste0("NOISeq.0.01.rst.tsv$"),full.names = TRUE), header = TRUE, stringsAsFactors = F)
  dearseq=read.delim(file=list.files(samples[i],pattern = "dearseqPreprocessed.rst.tsv$|dearseqPreprocessed.0.05.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
  dearseq=dearseq[dearseq$adjPval<0.01,]
  ### A & C: Plot for DEGs number
  shuffle_d<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.DESeq2.1000Perm.number.txt"),header = F)
  shuffle_e<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.edgeR.1000Perm.number.txt"),header = F)
  shuffle_l<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.limma.1000Perm.number.txt"),header = F)
  shuffle_n<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.NOISeq.1000Perm.number.txt"),header = F)
  shuffle_r<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.dearseqPreprocessed.1000Perm.number.txt"),header = F)
  shuffle_w<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.WilcoxonTest.1000Perm.number.txt"),header = F)
  data<-cbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_r,shuffle_w)
  colnames(data)<-c("DESeq2","edgeR","limma-voom","NOISeq","dearseq","Wilcoxon")
  data<-melt(data)
  sumData1<-ddply(data, ~variable, summarise, mean=mean(value),sd=sd(value))
  sumData1$errorU=sumData1$mean+sumData1$sd
  sumData1$errorD=sumData1$mean-sumData1$sd
  sumData1[sumData1$errorD<0,5]=0
  #Barplot for overlapped DEGs
  shuffle_d<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.DESeq2.1000Perm.overlapDEGs.txt"),header = F)
  shuffle_e<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.edgeR.1000Perm.overlapDEGs.txt"),header = F)
  shuffle_l<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.limma.1000Perm.overlapDEGs.txt"),header = F)
  shuffle_n<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.NOISeq.1000Perm.overlapDEGs.txt"),header = F)
  shuffle_r<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.dearseqPreprocessed.1000Perm.overlapDEGs.txt"),header = F)
  shuffle_w<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.WilcoxonTest.1000Perm.overlapDEGs.txt"),header = F)
  data<-cbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_r,shuffle_w)
  colnames(data)<-c("DESeq2","edgeR","limma-voom","NOISeq","dearseq","Wilcoxon")
  data<-melt(data)
  sumData2<-ddply(data, ~variable, summarise, mean=mean(value),sd=sd(value))
  sumData2$errorU=sumData2$mean+sumData2$sd
  sumData2$errorD=sumData2$mean-sumData2$sd
  sumData2[sumData2$errorD<0,5]=0
  #ymin=min(sumData1$mean-sumData1$sd,sumData2$mean-sumData2$sd)
  #if(ymin>0){ymin=0}
  pointsReal<-data.frame(px=c("DESeq2","edgeR","limma-voom","NOISeq","dearseq","Wilcoxon"),py=c(nrow(deseq2),nrow(edger),nrow(limma),nrow(noiseq),nrow(dearseq),nrow(wilcoxon)))
  totalNum=get(strsplit(samples[i],"/")[[1]][2])
  pA<-ggplot(sumData1,aes(x=variable,y=mean))+geom_bar(stat="identity", color="black")+geom_errorbar(aes(ymin=errorD, ymax=errorU),width=0.2,position=position_dodge(0.9))+
    geom_point(data=pointsReal,aes(x=px,y=py, fill="red"), color="red", shape=23, size=1.5)+labs(y="# of identified DEGs\nfrom permuted data",x="")+
    scale_y_continuous(label=comma_format(accuracy = 1),trans=log2_pseudo1, breaks=trans_breaks(function(x){log2(x+1)},function(x){2^x-1},n=6),limits = c(0,max(pointsReal$py)), 
                       sec.axis = sec_axis(trans=~./totalNum, breaks=c(0,0.1,1,10,20,40,60),name="% of identified DEGs\nfrom permuted data"))+
    scale_fill_manual(values="red", labels="# of identified DEGs from the original data", name="")+
    theme(legend.position = "bottom", legend.box.spacing = unit(-0.2,"inch"),legend.key = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
  pC<-ggplot(sumData2,aes(x=variable,y=mean))+geom_bar(stat="identity", color="black")+geom_errorbar(aes(ymin=errorD, ymax=errorU),width=0.2,position=position_dodge(0.9))+
    geom_point(data=pointsReal,aes(x=px,y=py,fill="red"), color="red", shape=23, size=1.5)+
    scale_y_continuous(label=comma_format(accuracy = 1),trans=log2_pseudo1, breaks=trans_breaks(function(x){log2(x+1)},function(x){2^x-1},n=6),limits = c(0,max(pointsReal$py)),
                       sec.axis = sec_axis(trans=~./get(strsplit(samples[i],"/")[[1]][2]), breaks=c(0,0.1,1,10,20,40,60), name="% of identified DEGs from\nboth original and permuted data"))+
    labs(y="# of DEGs identified from\nboth original and permuted data",x="")+
    scale_fill_manual(values="red", labels="# of identified DEGs from the original data", name="")+
    theme(legend.position = "bottom", legend.box.spacing = unit(-0.2,"inch"), legend.key = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
  
  ### B: Distribution for repeated times of DEGs 
  shuffle_d<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.DESeq2.1000Perm.uniqGene.txt"),header = F)
  shuffle_d$V2="DESeq2"
  shuffle_e<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.edgeR.1000Perm.uniqGene.txt"),header = F)
  shuffle_e$V2="edgeR"
  shuffle_n<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.NOISeq.1000Perm.uniqGene.txt"),header = F)
  shuffle_n$V2="NOISeq"
  shuffle_r<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.dearseqPreprocessed.uniqGene.txt"),header = F)
  shuffle_r$V2="dearseq"
  file_l=paste0(samples[i],"/FDR-0.01/shuffledLabel.limma.1000Perm.uniqGene.txt")
  file_w=paste0(samples[i],"/FDR-0.01/shuffledLabel.WilcoxonTest.1000Perm.uniqGene.txt")
  if(file.info(file_l)$size==0){
    shuffle_l<-data.frame(V1=NA,V2=NA)
  }else{
    shuffle_l<-read.delim(file=file_l,header = F)
  }
  shuffle_l$V2="limma-voom"
  if(file.info(file_w)$size==0){
    shuffle_w<-data.frame(V1=NA,V2=NA)
  }else{
    shuffle_w<-read.delim(file=file_w,header = F) 
  }
  shuffle_w$V2="Wilcoxon"
  data<-rbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_r,shuffle_w)
  data$V2=factor(data$V2, levels = c("DESeq2","edgeR","limma-voom","NOISeq","dearseq","Wilcoxon"))
  #data$V1=(data$V1/shuffleNum)*100
  xlabels=paste(seq(1,1000,by=200),paste("(",seq(1,1000,by=200)/10,"%)",sep=""),sep = "\n")
  pB<-ggplot(data,aes(x=V1))+geom_histogram(fill="white",colour="black",bins=30)+facet_wrap(~V2)+
    scale_y_continuous(label=comma_format(accuracy = 1), trans=log2_pseudo1, breaks=trans_breaks(function(x){log2(x+1)},function(x){2^x-1},n=6))+
    scale_x_continuous(breaks = seq(1,1000,by=200), labels = xlabels)+
    labs(y="# of identified DEGs",x="# of permuted datasets where a gene is wrongly identified as a DEG")+
    theme(strip.background = element_blank(), strip.text = element_text(size=12), axis.text.x = element_text(size=8))
  
  ### D: Shuffled DEGs repeat time across real DEGs
  shuffle_d<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.DESeq2.1000Perm.uniqGene.txt"),header = F)
  shuffle_d$V1=shuffle_d$V1/10
  deseq2[,7]<-rownames(deseq2)
  binFactor<-as.factor(ceiling(c(1:nrow(deseq2))/100))
  data<-left_join(deseq2, shuffle_d,by=c("V7"="V2"))
  data[is.na(data$V1),8] = 0
  sortData<-as.data.frame(cbind(c(1:nrow(data)),data[order(abs(data$log2FoldChange),decreasing = TRUE),c(2,8)]))
  subdata<-as.data.frame(cbind(sortData[,-1],binFactor))
  subdata$log2FoldChange<-abs(subdata$log2FoldChange)
  subdata$binFactor<-as.factor(subdata$binFactor)
  finaldata<-ddply(subdata, ~binFactor, summarise, meanP=mean(V1),meanF=mean(log2FoldChange))
  finaldata<-cbind(finaldata,"DESeq2")
  colnames(finaldata)[4]="method"
  shuffle_e<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.edgeR.1000Perm.uniqGene.txt"),header = F)
  edger[,6]<-rownames(edger)
  shuffle_e$V1=shuffle_e$V1/10
  binFactor<-as.factor(ceiling(c(1:nrow(edger))/100))
  data<-left_join(edger, shuffle_e,by=c("V6"="V2"))
  data[is.na(data$V1),7] = 0
  sortData<-as.data.frame(cbind(c(1:nrow(data)),data[order(abs(data$logFC),decreasing = TRUE),c(1,7)]))
  subdata<-as.data.frame(cbind(sortData[,-1],binFactor))
  subdata$binFactor<-as.factor(subdata$binFactor)
  subdata$logFC<-abs(subdata$logFC)
  subdata<-ddply(subdata, ~binFactor, summarise, meanP=mean(V1),meanF=mean(logFC))
  subdata<-cbind(subdata,"edgeR")
  colnames(subdata)[4]="method"
  finaldata<-rbind(finaldata,subdata)
  finaldata$binFactor=as.numeric(as.character(finaldata$binFactor))
  rho1=round(cor(finaldata[finaldata$method=="DESeq2",1], finaldata[finaldata$method=="DESeq2",2],method = "spearman"),digits = 2)
  rho2=round(cor(finaldata[finaldata$method=="edgeR",1], finaldata[finaldata$method=="edgeR",2],method = "spearman"),digits = 2)
  anno_text <- data.frame(
    rho = c(paste0("Spearman's rho = ",rho1),paste0("Spearman's rho = ",rho2)),
    method = c("DESeq2","edgeR")
  )
  xbreaks=seq(1,max(finaldata$binFactor), by=50)
  tmp<-finaldata[finaldata$method=="DESeq2",]
  deseqFC=paste(round(na.omit(tmp[xbreaks,])[,3],digits = 2), collapse = " ")
  tmp<-finaldata[finaldata$method=="edgeR",]
  edgerFC=paste(round(na.omit(tmp[xbreaks,])[,3], digits = 2), collapse = " ")
  pD<-ggplot(finaldata,aes(x=binFactor,y=meanP))+geom_point()+geom_smooth()+ facet_wrap(~method, scales = "free_x") +
    geom_text(data = anno_text,mapping= aes(x=0,y= max(finaldata$meanP),label=rho), hjust=-0.1)+ 
    scale_y_continuous(trans=log10_pseudo, breaks=trans_breaks(function(x){log10(x+1e-6)},function(x){10^x-1e-6},n=6), labels=function(x){signif(x, digits = 1)}, limits = c(0,max(finaldata$meanP)))+
    scale_x_continuous(breaks = xbreaks)+
    labs(y="% of permuted datasets where\na gene is wrongly identified as a DEG",caption = paste0(deseqFC,"\n",edgerFC),
         x="Rank of abs. log2(fold-change)\n(from high to low) of DEGs from the original data")+theme(strip.background = element_blank(), strip.text = element_text(size=12))
  
  ### E: GO enrichment barplot
  file_d=paste0(samples[i],"/FDR-0.01/shuffledLabel.DESeq2.DEGs.gt10percent.GO.1000Perm.enrichment.tsv")
  if(file.info(file_d)$size==1){
    data<-data.frame(ID=NA,Description=NA,GeneRatio=NA,BgRatio=NA,pvalue=NA,p.adjust=NA,qvalue=NA,geneID=NA,Count=NA)
    data=na.omit(data)
  }else{
    data<-read.delim(file=file_d,header = TRUE)
  }
  data<-data[order(data$p.adjust,decreasing = F),]
  data$p.adjust<-(-log10(data$p.adjust))
  subdata<-data[c(1:5),]
  subdata$Description<-factor(subdata$Description,levels=subdata$Description)
  pE1<-ggplot(subdata,aes(x=p.adjust,y=Description,fill=p.adjust))+geom_bar(stat='identity',  width=0.8, fill="#2C5E87")+
    labs(y="",x="-log10(p.adjust)",title="DESeq2")+theme(plot.title = element_text(size=12))
  file_e=paste0(samples[i],"/FDR-0.01/shuffledLabel.edgeR.DEGs.gt10percent.GO.1000Perm.enrichment.tsv")
  if(file.info(file_e)$size==1){
    data<-data.frame(ID=NA,Description=NA,GeneRatio=NA,BgRatio=NA,pvalue=NA,p.adjust=NA,qvalue=NA,geneID=NA,Count=NA)
    data=na.omit(data)
  }else{
    data<-read.delim(file=file_e,header = TRUE)
  }
  data<-data[order(data$p.adjust,decreasing = F),]
  data$p.adjust<-(-log10(data$p.adjust))
  subdata<-data[c(1:5),]
  subdata$Description<-factor(subdata$Description,levels=subdata$Description)
  pE2<-ggplot(subdata,aes(x=p.adjust,y=Description,fill=p.adjust))+geom_bar(stat='identity',  width=0.8, fill="#2C5E87")+
    labs(y="",x="-log10(p.adjust)",title="edgeR")+theme(plot.title = element_text(size=12))
  
  ### F: Boxplots for Goodness-of-fit test
  #### DESeq2
  pvalues1 <- data.frame(count = gof[[i]][[5]])
  pvalues2 <- data.frame(count = gof[[i]][[6]])
  shuffle <- c(pvalues1[,1], pvalues1[,2])
  nonshuffle <- c(pvalues2[,1], pvalues2[,2])
  data <- data.frame(pvalues=c(shuffle, nonshuffle), cond=c(rep("Shuffled DEGs", length(shuffle)),
                                                            rep("Other genes", length(nonshuffle))))
  ymax=max(boxplot.stats(-log10(data[data$cond=="Shuffled DEGs",1]))$stats[5],boxplot.stats(-log10(data[data$cond=="Other genes",1]))$stats[5])
  pF1 <- ggplot(data, aes(x=cond, y=-log10(pvalues))) + geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(0,ymax)) +
    scale_x_discrete(labels=c("Other\ngenes","Genes wrongly\nidentified as DEGs\nfrom any permuted\ndatasets"))+labs(x="",y="Poorness of fit",title="DESeq2")+ 
    annotate(geom="text",x=1.5,y= ymax,label=paste0("p=",signif(wilcox.test(shuffle, nonshuffle)$p.value,digits = 2)))+
    theme(legend.position = "none", axis.text.x = element_text(size = 9),plot.title = element_text(size=12))
  #### edgeR
  pvalues1 <- data.frame(count = gof[[i]][[2]])
  pvalues2 <- data.frame(count = gof[[i]][[3]])
  shuffle <- c(pvalues1[,1], pvalues1[,2])
  nonshuffle <- c(pvalues2[,1], pvalues2[,2])
  data <- data.frame(pvalues=c(shuffle, nonshuffle), cond=c(rep("Shuffled DEGs", length(shuffle)),
                                                            rep("Other genes", length(nonshuffle))))
  ymax=max(boxplot.stats(-log10(data[data$cond=="Shuffled DEGs",1]))$stats[5],boxplot.stats(-log10(data[data$cond=="Other genes",1]))$stats[5])
  pF2 <- ggplot(data, aes(x=cond, y=-log10(pvalues))) + geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(0,ymax)) +
    scale_x_discrete(labels=c("Other\ngenes","Genes wrongly\nidentified as DEGs\nfrom any permuted\ndatasets"))+labs(x="",y="Poorness of fit",title="edgeR")+ 
    annotate(geom="text",x=1.5,y= ymax,label=paste0("p=",signif(wilcox.test(shuffle, nonshuffle)$p.value,digits = 2)))+
    theme(legend.position = "none", axis.text.x = element_text(size = 9), plot.title = element_text(size=12))
  
  ### Plot all in one figure
  prefix=paste0(strsplit(samples[i],'/')[[1]][2],"_",strsplit(samples[i],'/')[[1]][3],"_0.01")
  pdf(file=paste0(prefix,"_suppFig.pdf"),height = 11,width = 8)
  p<-ggdraw()+draw_label(prefix, x=0.5,y=0.98, fontface = "bold")+
    draw_plot(pA,0.02,0.66,0.38,0.3)+draw_plot(pB,0.4,0.66,0.6,0.3)+
    draw_plot(pC,0.02,0.36,0.38,0.3)+draw_plot(pD,0.4,0.36,0.6,0.3)+
    draw_plot(pE1,0,0.18,0.7,0.18)+draw_plot(pE2,0,0,0.7,0.18)+
    draw_plot(pF1,0.7,0.18,0.3,0.18)+draw_plot(pF2,0.7,0,0.3,0.18)+
    draw_plot_label(c("A","B","C","D","E","F"),x=c(0,0.4,0,0.4,0,0.7),y=c(0.96,0.96,0.66,0.66,0.36,0.36),size=15)
  print(p)
  dev.off()
}