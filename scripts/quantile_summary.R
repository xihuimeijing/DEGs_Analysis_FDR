#!/usr/bin/env Rscript


args <- commandArgs(TRUE)

usage=function(){
  cat('Usage:quantile_summary.R -i=<file> -c=[1]',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\t\tFILE\tInput file without header line.If not provided, read from stdin\n',file=stderr())
  cat('\t-p\tINT\tThe number of quantiles you want to calculate.[defalut use summary function]\n',file=stderr())
  cat('\t-c\tINT\tThe number of column that you want to summary[1]\n',file=stderr())
  q(save="no")
}

col=1
if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -i')
      }else{
        inFile=arg.split[2]
      }
    }
    if(grepl('^-c=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      col=as.numeric(arg.split[2])
    }
    if(grepl('^-p=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      prob=as.numeric(arg.split[2])
    }
	if(grepl('^-h',arg)){usage()}
  }
}else{
  usage()	
}

if(exists("inFile")){
	data<-read.delim(file=inFile,header=F,comment.char='#')
}else{
	data<-read.delim(file('stdin'), header = F)
}
if(exists("prob")){
  quantile(data[,col],probs=seq(0,1,length.out=prob+1))
}else{
  summary(data[,col])
}
