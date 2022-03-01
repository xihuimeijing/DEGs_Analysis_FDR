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
# Use setwd() to change the working directory to the directory containing the data fro plot
immuno_gof=readRDS(file="ImmuneTherapy_goodness_of_fit.rds")

########## Define functions ##########
log2_pseudo1<-trans_new(name="log2_pseudo1",transform = function(x){log2(x+1)},inverse = function(x){2^x-1})

########## Set theme for the plots ##########
theme_update(axis.text = element_text(color="black"),text=element_text(size=12), plot.title=element_text(hjust=0.5), 
             panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA), 
             panel.grid.major = element_blank(), panel.grid.minor = element_blank())

########## Figure 1 ##########
setwd("./ImmunotherapyData")
### Read DEGs based on real data
deseq2=read.delim(file="GSE91061_BMS038109Sample.real.DESeq2.rst.tsv", header = TRUE, stringsAsFactors = F)
edger=read.delim(file="GSE91061_BMS038109Sample.real.edgeR.rst.tsv", header = TRUE, stringsAsFactors = F)
wilcoxon=read.delim(file="GSE91061_BMS038109Sample.real.WilcoxonTest.rst.tsv", header = TRUE, stringsAsFactors = F)
noiseq=read.delim(file="GSE91061_BMS038109Sample.real.NOISeq.rst.tsv", header = TRUE, stringsAsFactors = F)
dearseq=read.delim(file="GSE91061_BMS038109Sample.real.dearseqPreprocessed.rst.tsv", header = TRUE, stringsAsFactors = F)
file_l="GSE91061_BMS038109Sample.real.limma.rst.tsv"
if(file.info(file_l)$size<=1){
  limma<-data.frame()
}else{
  limma=read.delim(file=file_l, header = TRUE, stringsAsFactors = F)
}

### A & C: Plot for DEGs number
shuffle_d<-read.delim(file="shuffledLabel.DESeq2.number.txt",header = F)
shuffle_e<-read.delim(file="shuffledLabel.edgeR.number.txt",header = F)
shuffle_l<-read.delim(file="shuffledLabel.limma.number.txt",header = F)
shuffle_n<-read.delim(file="shuffledLabel.NOISeq.number.txt",header = F)
shuffle_r<-read.delim(file="shuffledLabel.dearseqPreprocessed.number.txt",header = F)
shuffle_w<-read.delim(file="shuffledLabel.WilcoxonTest.number.txt",header = F)
data<-cbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_r,shuffle_w)
colnames(data)<-c("DESeq2","edgeR","limma-voom","NOISeq","dearseq","Wilcoxon")
data<-melt(data)
sumData1<-ddply(data, ~variable, summarise, mean=mean(value),sd=sd(value))
sumData1$errorU=sumData1$mean+sumData1$sd
sumData1$errorD=sumData1$mean-sumData1$sd
sumData1[sumData1$errorD<0,5]=0
#Barplot for overlapped DEGs
shuffle_d<-read.delim(file="shuffledLabel.DESeq2.overlapDEGs.txt",header = F)
shuffle_e<-read.delim(file="shuffledLabel.edgeR.overlapDEGs.txt",header = F)
shuffle_l<-read.delim(file="shuffledLabel.limma.overlapDEGs.txt",header = F)
shuffle_n<-read.delim(file="shuffledLabel.NOISeq.overlapDEGs.txt",header = F)
shuffle_r<-read.delim(file="shuffledLabel.dearseqPreprocessed.overlapDEGs.txt",header = F)
shuffle_w<-read.delim(file="shuffledLabel.WilcoxonTest.overlapDEGs.txt",header = F)
data<-cbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_r,shuffle_w)
colnames(data)<-c("DESeq2","edgeR","limma-voom","NOISeq","dearseq","Wilcoxon")
data<-melt(data)
sumData2<-ddply(data, ~variable, summarise, mean=mean(value),sd=sd(value))
sumData2$errorU=sumData2$mean+sumData2$sd
sumData2$errorD=sumData2$mean-sumData2$sd
sumData2[sumData2$errorD<0,5]=0
ymax=max(sumData1$errorU,sumData2$errorU)
pointsReal<-data.frame(px=c("DESeq2","edgeR","limma-voom","NOISeq","dearseq","Wilcoxon"),py=c(nrow(deseq2),nrow(edger),nrow(limma),nrow(noiseq),nrow(dearseq),nrow(wilcoxon)))
pA<-ggplot(sumData1,aes(x=variable,y=mean))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=errorD,ymax=errorU),width=0.2,position=position_dodge(0.9))+
  geom_point(data=pointsReal,aes(x=px,y=py, fill="red"), color="red", shape=23, size=1.5)+
  scale_y_continuous(limits=c(0,ymax),sec.axis = sec_axis(trans=~./220.67, name="% of identified DEGs\nfrom permuted data"))+
  labs(y="# of identified DEGs\nfrom permuted data",x="")+scale_fill_manual(values="red", labels="# of identified DEGs from the original data", name="")+
  theme(legend.position = "bottom", legend.box.spacing = unit(-0.2,"inch"),legend.key = element_blank(), axis.text.x=element_text(angle=45, hjust = 1))
pC<-ggplot(sumData2,aes(x=variable,y=mean))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=errorD,ymax=errorU),width=0.2,position=position_dodge(0.9))+
  geom_point(data=pointsReal,aes(x=px,y=py,fill="red"), color="red", shape=23, size=1.5)+
  scale_y_continuous(limits=c(0,ymax),sec.axis = sec_axis(trans=~./220.67, name="% of DEGs identified from\nboth original and permuted data"))+
  labs(y="# of DEGs identified from\nboth original and permuted data",x="")+scale_fill_manual(values="red", labels="# of identified DEGs from the original data", name="")+
  theme(legend.position = "bottom", legend.box.spacing = unit(-0.2,"inch"), legend.key = element_blank(),axis.text.x=element_text(angle=45, hjust = 1))

### B: Distribution for repeated times of DEGs 
shuffle_d<-read.table(file="shuffledLabel.DESeq2.uniqGene.txt",header = F) 
shuffle_d$V2="DESeq2"
shuffle_e<-read.table(file="shuffledLabel.edgeR.uniqGene.txt",header = F)
shuffle_e$V2="edgeR"
shuffle_l<-read.table(file="shuffledLabel.limma.uniqGene.txt",header = F)
shuffle_l$V2="limma-voom"
shuffle_n<-read.table(file="shuffledLabel.NOISeq.uniqGene.txt",header = F)
shuffle_n$V2="NOISeq"
shuffle_r<-read.table(file="shuffledLabel.dearseqPreprocessed.uniqGene.txt",header = F)
shuffle_r$V2="dearseq"
shuffle_w<-read.table(file="shuffledLabel.WilcoxonTest.uniqGene.txt",header = F)
shuffle_w$V2="Wilcoxon"
data<-rbind(shuffle_d,shuffle_e,shuffle_l,shuffle_n,shuffle_r,shuffle_w)
data$V2=factor(data$V2, levels = c("DESeq2","edgeR","limma-voom","NOISeq","dearseq","Wilcoxon"))
xlabels=paste(seq(1,1000,by=200),paste("(",seq(1,1000,by=200)/10,"%)",sep=""),sep = "\n")
pB<-ggplot(data,aes(x=V1))+geom_histogram(fill="white",colour="black",breaks=seq(1,1000,length.out = 31))+facet_wrap(~V2)+
  scale_y_continuous(label=comma_format(accuracy = 1), trans=log2_pseudo1, breaks=trans_breaks(function(x){log2(x+1)},function(x){2^x-1},n=6))+
  scale_x_continuous(breaks = seq(1,1000,by=200), labels = xlabels)+
  labs(y="# of identified DEGs",x="# of permuted datasets where a gene is wrongly identified as a DEG")+
  theme(strip.background = element_blank(),strip.text = element_text(size=12),axis.text.x=element_text(size=8))

### D: Shuffled DEGs repeat time across real DEGs
shuffle_d<-read.table(file="shuffledLabel.DESeq2.uniqGene.txt",header = F)
shuffle_d$V1=shuffle_d$V1/10
deseq2[,7]<-rownames(deseq2)
data<-left_join(deseq2, shuffle_d,by=c("V7"="V2"))
data[is.na(data$V1),8] = 0
sortData<-as.data.frame(cbind(c(1:nrow(data)),data[order(abs(data$log2FoldChange),decreasing = TRUE),c(2,8)]))
finaldata<-cbind(sortData,"DESeq2")
colnames(finaldata)=c("rank","log2FC","percent","method")
shuffle_e<-read.table(file="shuffledLabel.edgeR.uniqGene.txt",header = F)
shuffle_e$V1=shuffle_e$V1/10
edger[,6]<-rownames(edger)
data<-left_join(edger, shuffle_e,by=c("V6"="V2"))
data[is.na(data$V1),7] = 0
sortData<-as.data.frame(cbind(c(1:nrow(data)),data[order(abs(data$logFC),decreasing = TRUE),c(1,7)]))
sortData<-cbind(sortData,"edgeR")
colnames(sortData)=c("rank","log2FC","percent","method")
finaldata<-rbind(finaldata,sortData)
rho1=round(cor(finaldata[finaldata$method=="DESeq2",1], finaldata[finaldata$method=="DESeq2",3],method = "spearman"),digits = 2)
rho2=round(cor(finaldata[finaldata$method=="edgeR",1], finaldata[finaldata$method=="edgeR",3],method = "spearman"),digits = 2)
anno_text <- data.frame(
  rho = c(paste0("Spearman's rho = ",rho1),paste0("Spearman's rho = ",rho2)),
  method = c("DESeq2","edgeR")
)
xbreaks=seq(1,max(finaldata$rank), by=50)
tmp<-finaldata[finaldata$method=="DESeq2",]
deseqFC=paste(round(abs(na.omit(tmp[xbreaks,])[,2]),digits = 2), collapse = " ")
tmp<-finaldata[finaldata$method=="edgeR",]
edgerFC=paste(round(abs(na.omit(tmp[xbreaks,])[,2]), digits = 2), collapse = " ")
pD<-ggplot(finaldata,aes(x=rank,y=percent))+geom_point(alpha=0.6)+geom_smooth()+ facet_wrap(~method, scales = "free_x") + ylim(0,120)+
  geom_text(data = anno_text,mapping= aes(x=0,y= 118,label=rho), hjust=-0.1)+ scale_x_continuous(breaks = xbreaks)+
  labs(y="% of permuted datasets where\na gene is wrongly identified as a DEG", caption = paste0(deseqFC,"\n",edgerFC), 
       x="Rank of abs. log2(fold-change)\n(from high to low) of DEGs from the original data")+theme(strip.background = element_blank(), strip.text = element_text(size=12))

### E: GO enrichment barplot
data<-read.delim(file="shuffledLabel.DESeq2.DEGs.gt10percent.GO.enrichment.tsv",header = TRUE)
data<-data[order(data$p.adjust,decreasing = F),]
data$p.adjust<-(-log10(data$p.adjust))
subdata<-data[c(1:5),]
subdata$Description<-factor(subdata$Description,levels=subdata$Description)
pE1<-ggplot(subdata,aes(x=p.adjust,y=Description))+geom_bar(stat='identity', width=0.8, fill="#2C5E87")+labs(y="",x="-log10(p.adjust)",title="DESeq2")+
  theme(plot.title = element_text(size=12))
data<-read.delim(file="shuffledLabel.edgeR.DEGs.gt10percent.GO.enrichment.tsv",header = TRUE)
data<-data[order(data$p.adjust,decreasing = F),]
data$p.adjust<-(-log10(data$p.adjust))
subdata<-data[c(1:5),]
subdata$Description<-factor(subdata$Description,levels=subdata$Description)
pE2<-ggplot(subdata,aes(x=p.adjust,y=Description))+geom_bar(stat='identity', width=0.8, fill="#2C5E87")+labs(y="",x="-log10(p.adjust)",title="edgeR")+
  theme(plot.title = element_text(size=12))

### F: Boxplots for Goodness-of-fit test
#### DESeq2
pvalues1 <- data.frame(count = immuno_gof[[5]])
pvalues2 <- data.frame(count = immuno_gof[[6]])
shuffle <- c(pvalues1[,1], pvalues1[,2])
nonshuffle <- c(pvalues2[,1], pvalues2[,2])
data <- data.frame(pvalues=c(shuffle, nonshuffle), cond=c(rep("Shuffled DEGs", length(shuffle)),
                                                          rep("Other genes", length(nonshuffle))))
log10ps=-log10(data[,1])
pF1 <- ggplot(data, aes(x=cond, y=-log10(pvalues))) + geom_boxplot(outlier.shape = NA) + ylim(0,150)+
  scale_x_discrete(labels=c("Other\ngenes","Genes wrongly\nidentified as DEGs\nfrom any permuted\ndatasets"))+labs(x="",y="Poorness of fit",title="DESeq2")+ 
  annotate(geom="text",x=1.5,y= 140,label="p<2.2e-16")+theme(legend.position = "none", axis.text.x=element_text(size=9), plot.title = element_text(size=12))
#### edgeR
pvalues1 <- data.frame(count = immuno_gof[[2]])
pvalues2 <- data.frame(count = immuno_gof[[3]])
shuffle <- c(pvalues1[,1], pvalues1[,2])
nonshuffle <- c(pvalues2[,1], pvalues2[,2])
data <- data.frame(pvalues=c(shuffle, nonshuffle), cond=c(rep("Shuffled DEGs", length(shuffle)),
                                                          rep("Other genes", length(nonshuffle))))
log10ps=-log10(data[,1])
pF2 <- ggplot(data, aes(x=cond, y=-log10(pvalues))) + geom_boxplot(outlier.shape = NA) + ylim(0,25) +
  scale_x_discrete(labels=c("Other\ngenes","Genes wrongly\nidentified as DEGs\nfrom any permuted\ndatasets"))+labs(x="",y="Poorness of fit",title="edgeR")+ 
  annotate(geom="text",x=1.5,y= 24,label="p<2.2e-16")+theme(legend.position = "none", axis.text.x=element_text(size=9), plot.title = element_text(size=12))

### Plot all in one figure
pdf(file="ImmunotherapyData_Fig1.pdf",height = 11,width = 8)
p<-ggdraw()+draw_label("Immunotherapy", x=0.5, y=0.98, fontface = "bold")+
  draw_plot(pA,0.02,0.66,0.38,0.3)+draw_plot(pB,0.4,0.66,0.6,0.3)+
  draw_plot(pC,0.02,0.36,0.38,0.3)+draw_plot(pD,0.4,0.36,0.6,0.3)+
  draw_plot(pE1,0,0.18,0.7,0.18)+draw_plot(pE2,0,0,0.7,0.18)+
  draw_plot(pF1,0.7,0.18,0.3,0.18)+draw_plot(pF2,0.7,0,0.3,0.18)+
  draw_plot_label(c("A","B","C","D","E","F"),x=c(0,0.4,0,0.4,0,0.7),y=c(0.96,0.96,0.66,0.66,0.36,0.36),size=15)
print(p)
dev.off()