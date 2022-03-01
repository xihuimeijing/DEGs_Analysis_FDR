library(UpSetR)
library(grid)

########## Prepare data ##########
# Use setwd() to change the working directory to the directory containing the data for plot
samples=c(list.dirs("./TCGA",recursive = F),list.dirs("./GTEx",recursive = F))

pdf(file="allSamples_DESeq2-edgeR_realDEGs_upsetPlot.pdf", width = 5, height = 3.5)
for(i in 1:12){
  deseq2=read.delim(file=list.files(samples[i],pattern = "DESeq2.rst.tsv|DESeq2.0.05.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
  deseq2=deseq2[deseq2$padj<0.01,]
  edger=read.delim(file=list.files(samples[i],pattern = "edgeR.rst.tsv|edgeR.0.05.rst.tsv$",full.names = TRUE), header = TRUE, stringsAsFactors = F)
  edger=edger[edger$FDR<0.01,]
  data=list(DESeq2=row.names(deseq2),edgeR=row.names(edger))
  print(upset(fromList(data), sets=c("DESeq2","edgeR"), order.by = "degree", keep.order = TRUE, 
              sets.x.label = "Number of DEGs", mainbar.y.label = "Gene Intersections",text.scale = c(1.5,1.5,1.5,1.5,1.5,1.5)))
  print(grid.text(samples[i],x=0.5,y=0.98, gp=gpar(fontsize=15)))
}
deseq2=read.delim(file="./ImmunotherapyData/GSE91061_BMS038109Sample.real.DESeq2.rst.tsv", header = TRUE, stringsAsFactors = F)
edger=read.delim(file="./ImmunotherapyData/GSE91061_BMS038109Sample.real.edgeR.rst.tsv", header = TRUE, stringsAsFactors = F)
data=list(DESeq2=row.names(deseq2),edgeR=row.names(edger))
print(upset(fromList(data), sets=c("DESeq2","edgeR"), order.by = "degree", keep.order = TRUE, 
            sets.x.label = "Number of DEGs", mainbar.y.label = "Gene Intersections",text.scale = c(1.5,1.5,1.5,1.5,1.5,1.5)))
print(grid.text("Immunotherapy",x=0.5,y=0.98, gp=gpar(fontsize=15)))
dev.off()