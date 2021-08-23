#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
usage=function(){
  cat('Description: Summarize the real condition ratio in the shuffled lables. 
Author: Yumei Li, 04/19/2021.\n',file=stderr())
  cat('Usage:shuffledLable_ratio.R -r=real.condition.tsv -s=shuffledLable.tsv -o=result.tsv\n',file=stderr())
  cat('\t-r\t\tFILE\tConditions in one line for real data.\n',file=stderr())
  cat('\t-s\t\tFILE\tShuffled conditions, one line for each sample.\n',file=stderr())
  cat('\t-o\t\tSTRING\tOutput file name.\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}

if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-r=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      realCon=arg.split[2]
    }
    if(grepl('^-s=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      shuffledCon=arg.split[2]
    }
    if(grepl('^-o=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      outFile=arg.split[2]
    }
    if(grepl('^-h', arg)){
      usage()
    }
  }
}else{
  usage()
}

realData<-read.table(file=realCon,header = F)
realData<-factor(realData)
conLevels<-levels(realData)
condition1<-which(realData==conLevels[1])
condition2<-which(realData==conLevels[2])

shuffledData<-read.table(file=shuffledCon, header = F)
outMatrix=matrix(nrow = nrow(shuffledData),ncol=2)
for(i in 1:nrow(shuffledData)){
  subdata<-shuffledData[i,condition1]
  outMatrix[i,1]=round(length(subdata[1,which(subdata[1,]==conLevels[1])])/length(condition1),digits=4)
  subdata<-shuffledData[i,condition2]
  outMatrix[i,2]=round(length(subdata[1,which(subdata[1,]==conLevels[2])])/length(condition2),digits=4)
}
write.table(outMatrix, file=outFile,sep="\t", quote=F,row.names = F,col.names = F)
