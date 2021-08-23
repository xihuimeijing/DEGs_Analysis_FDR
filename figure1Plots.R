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
# Use setwd() to change the working directory to the directory containing count matrices
samples=c(list.dirs("./TCGA",recursive = F),list.dirs("./GTEx",recursive = F))
GTEx_gof=readRDS(file="gtex_goodness_of_fit.rds")
TCGA_gof=readRDS(file="tcga_goodness_of_fit.rds")
gof=c(TCGA_gof, GTEx_gof)
samples<-samples[c(1:10,12,11)]
immuno_gof=readRDS(file="ImmuneTherapy_goodness_of_fit.rds")

########## Define functions ##########
log2_pseudo1<-trans_new(name="log2_pseudo1",transform = function(x){log2(x+1)},inverse = function(x){2^x-1})
log10_pseudo<-trans_new(name="log10_pseudo",transform = function(x){log10(x+1e-6)},inverse = function(x){10^x-1e-6})

########## Set theme for the plots ##########
theme_update(axis.text = element_text(color="black"),text=element_text(size=12), plot.title=element_text(hjust=0.5), 
             panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA), 
             panel.grid.major = element_blank(), panel.grid.minor = element_blank())

########## Figure 1 ##########
immunoSample="./ImmunotherapyData/GSE91061"
### Read DEGs based on real data
deseq2=read.delim(file=list.files(immunoSample,pattern = "DESeq2.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
edger=read.delim(file=list.files(immunoSample,pattern = "edgeR.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
wilcoxon=read.delim(file=list.files(immunoSample,pattern = "WilcoxonTest.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
noiseq=read.delim(file=list.files(immunoSample,pattern = "NOISeq.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
limma=read.delim(file=list.files(immunoSample,pattern = "limma.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
### A & C: Plot for DEGs number
shuffle_d<-read.delim(file=paste0(immunoSample,"/shuffledLabel.DESeq2.1st-3rd.number.txt"),header = F)
shuffle_e<-read.delim(file=paste0(immunoSample,"/shuffledLabel.edgeR.1st-3rd.number.txt"),header = F)
shuffle_l<-read.delim(file=paste0(immunoSample,"/shuffledLabel.limma.1st-3rd.number.txt"),header = F)
shuffle_n<-read.delim(file=paste0(immunoSample,"/shuffledLabel.NOISeq.1st-3rd.number.txt"),header = F)
shuffle_w<-read.delim(file=paste0(immunoSample,"/shuffledLabel.WilcoxonTest.1st-3rd.number.txt"),header = F)
data<-cbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_w)
colnames(data)<-c("DESeq2","edgeR","limma-voom","NOISeq","Wilcoxon")  
data<-melt(data)
sumData1<-ddply(data, ~variable, summarise, mean=mean(value),sd=sd(value))  
sumData1$errorU=sumData1$mean+sumData1$sd
sumData1$errorD=sumData1$mean-sumData1$sd
sumData1[sumData1$errorD<0,5]=0
###Barplot for overlapped DEGs
shuffle_d<-read.delim(file=paste0(immunoSample,"/shuffledLabel.DESeq2.1st-3rd.overlapDEGs.txt"),header = F)
shuffle_e<-read.delim(file=paste0(immunoSample,"/shuffledLabel.edgeR.1st-3rd.overlapDEGs.txt"),header = F)
shuffle_l<-read.delim(file=paste0(immunoSample,"/shuffledLabel.limma.1st-3rd.overlapDEGs.txt"),header = F)
shuffle_n<-read.delim(file=paste0(immunoSample,"/shuffledLabel.NOISeq.1st-3rd.overlapDEGs.txt"),header = F)
shuffle_w<-read.delim(file=paste0(immunoSample,"/shuffledLabel.WilcoxonTest.1st-3rd.overlapDEGs.txt"),header = F)
data<-cbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_w)
colnames(data)<-c("DESeq2","edgeR","limma-voom","NOISeq","Wilcoxon")
data<-melt(data)
sumData2<-ddply(data, ~variable, summarise, mean=mean(value),sd=sd(value))
sumData2$errorU=sumData2$mean+sumData2$sd
sumData2$errorD=sumData2$mean-sumData2$sd
sumData2[sumData2$errorD<0,5]=0
ymax=max(sumData1$errorU,sumData2$errorU)
pointsReal<-data.frame(px=c("DESeq2","edgeR","limma-voom","NOISeq","Wilcoxon"),py=c(nrow(deseq2),nrow(edger),nrow(limma),nrow(noiseq),nrow(wilcoxon)))
pA<-ggplot(sumData1,aes(x=variable,y=mean))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=errorD,ymax=errorU),width=0.2,position=position_dodge(0.9))+
  geom_point(data=pointsReal,aes(x=px,y=py, fill="red"), color="red", shape=23, size=1.5)+ylim(0,ymax)+labs(y="# of identified DEGs",x="")+
  scale_fill_manual(values="red", labels="# of identified DEGs from the original data", name="")+
  theme(legend.position = "bottom", legend.box.spacing = unit(-0.2,"inch"),legend.key = element_blank(), axis.text.x=element_text(angle=45, hjust = 1))
pC<-ggplot(sumData2,aes(x=variable,y=mean))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=errorD,ymax=errorU),width=0.2,position=position_dodge(0.9))+
  geom_point(data=pointsReal,aes(x=px,y=py,fill="red"), color="red", shape=23, size=1.5)+ylim(0,ymax)+
  labs(y="# of DEGs identified from\nboth original and permuted data",x="")+
  scale_fill_manual(values="red", labels="# of identified DEGs from the original data", name="")+
  theme(legend.position = "bottom", legend.box.spacing = unit(-0.2,"inch"), legend.key = element_blank(),axis.text.x=element_text(angle=45, hjust = 1))
### B: Distribution for percentage
shuffleNum=nrow(read.delim(file=paste0(immunoSample,"/shuffledLabel.realRatio.1st-3rd.selected.txt"),header = F))
shuffle_d<-read.delim(file=paste0(immunoSample,"/shuffledLabel.DESeq2.1st-3rd.uniqGene.txt"),header = F)
shuffle_d$V2="DESeq2"
shuffle_e<-read.delim(file=paste0(immunoSample,"/shuffledLabel.edgeR.1st-3rd.uniqGene.txt"),header = F)
shuffle_e$V2="edgeR"
shuffle_n<-read.delim(file=paste0(immunoSample,"/shuffledLabel.NOISeq.1st-3rd.uniqGene.txt"),header = F)
shuffle_n$V2="NOISeq"
shuffle_l<-read.delim(file=paste0(immunoSample,"/shuffledLabel.limma.1st-3rd.uniqGene.txt"),header = F)
shuffle_l$V2="limma-voom"
shuffle_w<-read.delim(file=paste0(immunoSample,"/shuffledLabel.WilcoxonTest.1st-3rd.uniqGene.txt"),header = F) 
shuffle_w$V2="Wilcoxon"
data<-rbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_w)
data$V2=factor(data$V2, levels = c("DESeq2","edgeR","limma-voom","NOISeq","Wilcoxon"))
data$V1=(data$V1/shuffleNum)*100
pB<-ggplot(data,aes(x=V1))+geom_histogram(fill="white",colour="black",bins=30)+facet_wrap(~V2)+xlim(0,100)+
  scale_y_continuous(trans=log2_pseudo1, breaks=trans_breaks(function(x){log2(x+1)},function(x){2^x-1},n=6))+
  labs(y="# of identified DEGs",x="% of permuted datasets where a gene is wrongly identified as a DEG")+
  theme(strip.background = element_blank())
### D: Shuffled DEGs identified percentage across real DEGs
shuffle_d<-read.delim(file=paste0(immunoSample,"/shuffledLabel.DESeq2.1st-3rd.uniqGene.txt"),header = F)
shuffle_d$V1=shuffle_d$V1/shuffleNum
deseq2[,7]<-rownames(deseq2)
data<-left_join(deseq2, shuffle_d,by=c("V7"="V2"))
data[is.na(data$V1),8] = 0
sortData<-as.data.frame(cbind(c(1:nrow(data)),data[order(abs(data$log2FoldChange),decreasing = TRUE),8]))
finaldata<-cbind(sortData,"DESeq2")
colnames(finaldata)[3]="method"
shuffle_e<-read.delim(file=paste0(immunoSample,"/shuffledLabel.edgeR.1st-3rd.uniqGene.txt"),header = F)
shuffle_e$V1=shuffle_e$V1/shuffleNum
edger[,6]<-rownames(edger)
data<-left_join(edger, shuffle_e,by=c("V6"="V2"))
data[is.na(data$V1),7] = 0
sortData<-as.data.frame(cbind(c(1:nrow(data)),data[order(abs(data$logFC),decreasing = TRUE),7]))
sortData<-cbind(sortData,"edgeR")
colnames(sortData)[3]="method"
finaldata<-rbind(finaldata,sortData)
rho1=round(cor(finaldata[finaldata$method=="DESeq2",1], finaldata[finaldata$method=="DESeq2",2],method = "spearman"),digits = 2)
rho2=round(cor(finaldata[finaldata$method=="edgeR",1], finaldata[finaldata$method=="edgeR",2],method = "spearman"),digits = 2)
anno_text <- data.frame(
  rho = c(paste0("Spearman's rho = ",rho1),paste0("Spearman's rho = ",rho2)),
  method = c("DESeq2","edgeR")
)
pD<-ggplot(finaldata,aes(x=V1,y=V2))+geom_point(alpha=0.6)+geom_smooth()+ facet_wrap(~method, scales = "free_x") + ylim(0,1.2)+
    geom_text(data = anno_text,mapping= aes(x=0,y= 1.18,label=rho), hjust=-0.1)+
    labs(y="% of permuted datasets\nwhere a gene is wrongly identified as a DEG",x="Rank of abs. log2(fold-change)\n(from high to low) of DEGs from the original data")+theme(strip.background = element_blank())
### E: GO enrichment barplot
data<-read.delim(file=paste0(immunoSample,"/shuffledLabel.DESeq2.1st-3rd.DEGs.gt10percent.GO.enrichment.tsv"),header = TRUE)
data<-data[order(data$p.adjust,decreasing = F),]
data$p.adjust<-(-log10(data$p.adjust))
subdata<-data[c(1:5),]
subdata$Description<-factor(subdata$Description,levels=subdata$Description)
pE1<-ggplot(subdata,aes(x=p.adjust,y=Description))+geom_bar(stat='identity', width=0.8, fill="#2C5E87")+labs(y="",x="-log10(p.adjust)",title="DESeq2")
data<-read.delim(file=paste0(immunoSample,"/shuffledLabel.edgeR.1st-3rd.DEGs.gt10percent.GO.enrichment.tsv"),header = TRUE)
data<-data[order(data$p.adjust,decreasing = F),]
data$p.adjust<-(-log10(data$p.adjust))
subdata<-data[c(1:5),]
subdata$Description<-factor(subdata$Description,levels=subdata$Description)
pE2<-ggplot(subdata,aes(x=p.adjust,y=Description))+geom_bar(stat='identity', width=0.8, fill="#2C5E87")+labs(y="",x="-log10(p.adjust)",title="edgeR")
### F: Boxplots for Goodness-of-fit test
  #### DESeq2
pvalues1 <- data.frame(count = immuno_gof[[5]])
pvalues2 <- data.frame(count = immuno_gof[[6]])
shuffle <- c(pvalues1[,1], pvalues1[,2])
nonshuffle <- c(pvalues2[,1], pvalues2[,2])
data <- data.frame(pvalues=c(shuffle, nonshuffle), cond=c(rep("Shuffled DEGs", length(shuffle)),
                                                          rep("Other genes", length(nonshuffle))))
log10ps=-log10(data[,1])
pF1 <- ggplot(data, aes(x=cond, y=-log10(pvalues), fill=cond)) + geom_boxplot(outlier.shape = NA) + ylim(0,150)+
    scale_x_discrete(labels=c("Other\ngenes","Genes wrongly\nidentified as DEGs\nfrom any permuted\ndatasets"))+labs(x="",y="Poorness of fit",title="DESeq2")+ 
    annotate(geom="text",x=1.5,y= 140,label="p<2.2e-16")+theme(legend.position = "none")
    #### edgeR
pvalues1 <- data.frame(count = immuno_gof[[2]])
pvalues2 <- data.frame(count = immuno_gof[[3]])
shuffle <- c(pvalues1[,1], pvalues1[,2])
nonshuffle <- c(pvalues2[,1], pvalues2[,2])
data <- data.frame(pvalues=c(shuffle, nonshuffle), cond=c(rep("Shuffled DEGs", length(shuffle)),
                                                          rep("Other genes", length(nonshuffle))))
log10ps=-log10(data[,1])
pF2 <- ggplot(data, aes(x=cond, y=-log10(pvalues), fill=cond)) + geom_boxplot(outlier.shape = NA) + ylim(0,25) +
    scale_x_discrete(labels=c("Other\ngenes","Genes wrongly\nidentified as DEGs\nfrom any permuted\ndatasets"))+labs(x="",y="Poorness of fit",title="edgeR")+ 
    annotate(geom="text",x=1.5,y= 24,label="p<2.2e-16")+theme(legend.position = "none")
### Plot all in one figure
pdf(file="Figure 1.pdf",height = 11,width = 8)
p<-ggdraw()+
  draw_plot(pA,0.02,0.66,0.38,0.3)+draw_plot(pB,0.4,0.66,0.6,0.3)+
  draw_plot(pC,0.02,0.36,0.38,0.3)+draw_plot(pD,0.4,0.36,0.6,0.3)+
  draw_plot(pE1,0,0.18,0.7,0.18)+draw_plot(pE2,0,0,0.7,0.18)+
  draw_plot(pF1,0.7,0.18,0.3,0.18)+draw_plot(pF2,0.7,0,0.3,0.18)+
  draw_plot_label(c("A","B","C","D","E","F"),x=c(0,0.4,0,0.4,0,0.7),y=c(0.96,0.96,0.66,0.66,0.36,0.36),size=15)
print(p)
dev.off()

########## Supplementary Figure 1 ##########
for(i in 1:12){
  deseq2=read.delim(file=list.files(samples[i],pattern = "DESeq2.rst.tsv|DESeq2.0.05.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
  deseq2=deseq2[deseq2$padj<0.01,]
  edger=read.delim(file=list.files(samples[i],pattern = "edgeR.rst.tsv|edgeR.0.05.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
  edger=edger[edger$FDR<0.01,]
  data=list(DESeq2=row.names(deseq2),edgeR=row.names(edger))
  print(upset(fromList(data), sets=c("DESeq2","edgeR"), order.by = "degree", keep.order = TRUE, 
              sets.x.label = "# of DEGs", mainbar.y.label = "Gene Intersections",text.scale = c(1.5,1.5,1.5,1.5,1.5,1.5)))
  print(grid.text(samples[i],x=0.5,y=0.98, gp=gpar(fontsize=15)))
}
deseq2=read.delim(file="./ImmunotherapyData/GSE91061/GSE91061_BMS038109Sample.real.DESeq2.rst.tsv", header = TRUE, stringsAsFactors = F)
edger=read.delim(file="./ImmunotherapyData/GSE91061/GSE91061_BMS038109Sample.real.edgeR.rst.tsv", header = TRUE, stringsAsFactors = F)
data=list(DESeq2=row.names(deseq2),edgeR=row.names(edger))
print(upset(fromList(data), sets=c("DESeq2","edgeR"), order.by = "degree", keep.order = TRUE, 
            sets.x.label = "# of DEGs", mainbar.y.label = "Gene Intersections",text.scale = c(1.5,1.5,1.5,1.5,1.5,1.5)))
print(grid.text("Immunotherapy",x=0.5,y=0.98, gp=gpar(fontsize=15)))
dev.off()

########## Supplementary Figure 2 ##########
deseq2=read.delim(file="./ImmunotherapyData/GSE91061/GSE91061_BMS038109Sample.real.DESeq2.rst.tsv", header = TRUE, stringsAsFactors = F)
p1<-ggplot(deseq2,aes(x=abs(log2FoldChange)))+geom_histogram(bins=20, fill="gray", color="black")+labs(y="# of identified DEGs",x="Abs. log2(fold-change)",title="DESeq2")
edger=read.delim(file="./ImmunotherapyData/GSE91061/GSE91061_BMS038109Sample.real.edgeR.rst.tsv", header = TRUE, stringsAsFactors = F)
p2<-ggplot(edger,aes(x=abs(logFC)))+geom_histogram(bins=20, fill="gray", color="black")+labs(y="# of identified DEGs",x="Abs. log2(fold-change)",title="edgeR")
pdf(file="Supplementary Figure 2.pdf", height = 4, width = 8)
grid.arrange(p1,p2,nrow=1)
dev.off()

########## Supplementary Figure 4-15 ##########
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
  noiseq=read.delim(file=list.files(samples[i],pattern = "NOISeq.0.01.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
    
  ### A & C: Plot for DEGs number
  shuffle_d<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.DESeq2.number.txt"),header = F)
  shuffle_e<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.edgeR.number.txt"),header = F)
  shuffle_l<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.limma.number.txt"),header = F)
  shuffle_n<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.NOISeq.number.txt"),header = F)
  shuffle_w<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.WilcoxonTest.number.txt"),header = F)
  data<-cbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_w)
  colnames(data)<-c("DESeq2","edgeR","limma-voom","NOISeq","Wilcoxon")
  data<-melt(data)
  sumData1<-ddply(data, ~variable, summarise, mean=mean(value),sd=sd(value))
  sumData1$errorU=sumData1$mean+sumData1$sd
  sumData1$errorD=sumData1$mean-sumData1$sd
  sumData1[sumData1$errorD<0,5]=0
  #Barplot for overlapped DEGs
  shuffle_d<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.DESeq2.overlapDEGs.txt"),header = F)
  shuffle_e<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.edgeR.overlapDEGs.txt"),header = F)
  shuffle_l<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.limma.overlapDEGs.txt"),header = F)
  shuffle_n<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.NOISeq.overlapDEGs.txt"),header = F)
  shuffle_w<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.WilcoxonTest.overlapDEGs.txt"),header = F)
  data<-cbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_w)
  colnames(data)<-c("DESeq2","edgeR","limma-voom","NOISeq","Wilcoxon")
  data<-melt(data)
  sumData2<-ddply(data, ~variable, summarise, mean=mean(value),sd=sd(value))
  sumData2$errorU=sumData2$mean+sumData2$sd
  sumData2$errorD=sumData2$mean-sumData2$sd
  sumData2[sumData2$errorD<0,5]=0
  pointsReal<-data.frame(px=c("DESeq2","edgeR","limma-voom","NOISeq","Wilcoxon"),py=c(nrow(deseq2),nrow(edger),nrow(limma),nrow(noiseq),nrow(wilcoxon)))
  pA<-ggplot(sumData1,aes(x=variable,y=mean))+geom_bar(stat="identity", color="black")+geom_errorbar(aes(ymin=errorD, ymax=errorU),width=0.2,position=position_dodge(0.9))+
    geom_point(data=pointsReal,aes(x=px,y=py, fill="red"), color="red", shape=23, size=1.5)+labs(y="# of identified DEGs",x="")+
    scale_y_continuous(trans=log2_pseudo1, breaks=trans_breaks(function(x){log2(x+1)},function(x){2^x-1},n=6),limits = c(0,max(pointsReal$py)))+
    scale_fill_manual(values="red", labels="# of identified DEGs from the original data", name="")+
    theme(legend.position = "bottom", legend.box.spacing = unit(-0.2,"inch"),legend.key = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
  pC<-ggplot(sumData2,aes(x=variable,y=mean))+geom_bar(stat="identity", color="black")+geom_errorbar(aes(ymin=errorD, ymax=errorU),width=0.2,position=position_dodge(0.9))+
    geom_point(data=pointsReal,aes(x=px,y=py,fill="red"), color="red", shape=23, size=1.5)+
    scale_y_continuous(trans=log2_pseudo1, breaks=trans_breaks(function(x){log2(x+1)},function(x){2^x-1},n=6),limits = c(0,max(pointsReal$py)))+
    labs(y="# of DEGs identified from\nboth original and permuted data",x="")+
    scale_fill_manual(values="red", labels="# of identified DEGs from the original data", name="")+
    theme(legend.position = "bottom", legend.box.spacing = unit(-0.2,"inch"), legend.key = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
    
  ### B: Distribution for repeated times of DEGs 
  shuffleNum=nrow(read.delim(file=paste0(samples[i],"/shuffledLabel.realRatio.1st-3rd.selected.txt"),header = F))
  shuffle_d<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.DESeq2.uniqGene.txt"),header = F)
  shuffle_d$V2="DESeq2"
  shuffle_e<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.edgeR.uniqGene.txt"),header = F)
  shuffle_e$V2="edgeR"
  shuffle_n<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.NOISeq.uniqGene.txt"),header = F)
  shuffle_n$V2="NOISeq"
  file_l=paste0(samples[i],"/FDR-0.01/shuffledLabel.limma.uniqGene.txt")
  file_w=paste0(samples[i],"/FDR-0.01/shuffledLabel.WilcoxonTest.uniqGene.txt")
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
  data<-rbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_w)
  data$V2=factor(data$V2, levels = c("DESeq2","edgeR","limma-voom","NOISeq","Wilcoxon"))
  data$V1=(data$V1/shuffleNum)*100
  pB<-ggplot(data,aes(x=V1))+geom_histogram(fill="white",colour="black",bins=30)+facet_wrap(~V2)+xlim(0,100)+
    scale_y_continuous(trans=log2_pseudo1, breaks=trans_breaks(function(x){log2(x+1)},function(x){2^x-1},n=6))+
    labs(y="# of identified DEGs",x="% of permuted datasets where a gene is wrongly identified as a DEG")+
    theme(strip.background = element_blank())
    
  ### D: Shuffled DEGs repeat time across real DEGs
  shuffle_d<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.DESeq2.uniqGene.txt"),header = F)
  shuffle_d$V1=shuffle_d$V1/shuffleNum
  deseq2[,7]<-rownames(deseq2)
  binFactor<-as.factor(ceiling(c(1:nrow(deseq2))/100))
  data<-left_join(deseq2, shuffle_d,by=c("V7"="V2"))
  data[is.na(data$V1),8] = 0
  sortData<-as.data.frame(cbind(c(1:nrow(data)),data[order(abs(data$log2FoldChange),decreasing = TRUE),8]))
  subdata<-as.data.frame(cbind(sortData[,2],binFactor))
  subdata$binFactor<-as.factor(subdata$binFactor)
  finaldata<-ddply(subdata, ~binFactor, summarise, mean=mean(V1),sd=sd(V1))
  finaldata<-cbind(finaldata,"DESeq2")
  colnames(finaldata)[4]="method"
  shuffle_e<-read.delim(file=paste0(samples[i],"/FDR-0.01/shuffledLabel.edgeR.uniqGene.txt"),header = F)
  edger[,6]<-rownames(edger)
  shuffle_e$V1=shuffle_e$V1/shuffleNum
  binFactor<-as.factor(ceiling(c(1:nrow(edger))/100))
  data<-left_join(edger, shuffle_e,by=c("V6"="V2"))
  data[is.na(data$V1),7] = 0
  sortData<-as.data.frame(cbind(c(1:nrow(data)),data[order(abs(data$logFC),decreasing = TRUE),7]))
  subdata<-as.data.frame(cbind(sortData[,2],binFactor))
  subdata$binFactor<-as.factor(subdata$binFactor)
  subdata<-ddply(subdata, ~binFactor, summarise, mean=mean(V1),sd=sd(V1))
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
  pD<-ggplot(finaldata,aes(x=binFactor,y=mean))+geom_point()+geom_smooth()+ facet_wrap(~method, scales = "free_x") +
    geom_text(data = anno_text,mapping= aes(x=0,y= max(finaldata$mean),label=rho), hjust=-0.1)+
    scale_y_continuous(trans=log10_pseudo, breaks=trans_breaks(function(x){log10(x+1e-6)},function(x){10^x-1e-6},n=6),labels = function(x){signif(x*100, digits = 1)})+
    labs(y="% of permuted datasets\nwhere a gene is wrongly identified as a DEG",x="Rank of abs. log2(fold-change)\n(from high to low) of DEGs from the original data")+theme(strip.background = element_blank())
    
  ### E: GO enrichment barplot
  file_d=paste0(samples[i],"/FDR-0.01/shuffledLabel.DESeq2.DEGs.gt10percent.GO.enrichment.tsv")
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
    labs(y="",x="-log10(p.adjust)",title="DESeq2")
  file_e=paste0(samples[i],"/FDR-0.01/shuffledLabel.edgeR.DEGs.gt10percent.GO.enrichment.tsv")
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
    labs(y="",x="-log10(p.adjust)",title="edgeR")
    
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
    theme(legend.position = "none", axis.text.x = element_text(size = 9))
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
    theme(legend.position = "none", axis.text.x = element_text(size = 9))
    
  ### Plot all in one figure
  prefix=paste0("Supplementary Figure ",i+3,strsplit(samples[i],'/')[[1]][2],"_",strsplit(samples[i],'/')[[1]][3],"_0.01")
  pdf(file=paste0(prefix,".pdf"),height = 11,width = 8)
  p<-ggdraw()+
    draw_plot(pA,0.02,0.66,0.38,0.3)+draw_plot(pB,0.4,0.66,0.6,0.3)+
    draw_plot(pC,0.02,0.36,0.38,0.3)+draw_plot(pD,0.4,0.36,0.6,0.3)+
    draw_plot(pE1,0,0.18,0.7,0.18)+draw_plot(pE2,0,0,0.7,0.18)+
    draw_plot(pF1,0.7,0.18,0.3,0.18)+draw_plot(pF2,0.7,0,0.3,0.18)+
    draw_plot_label(c("A","B","C","D","E","F"),x=c(0,0.4,0,0.4,0,0.7),y=c(0.96,0.96,0.66,0.66,0.36,0.36),size=15)
  print(p)
  dev.off()
}



