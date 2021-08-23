#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
usage=function(){
  cat('Description:This script select and output columns based on column names.\n',file=stderr())
  cat('Usage:columnSelect.R \n',file=stderr())
  cat('\t-f\t\tFILE\tInput file to be selected from\n',file=stderr())
  cat('\t-i\t\tFILE\tThe file contained all the selected column names(one name per column) the same with -f file.\n',file=stderr())
  cat('\t-s\t\tSTRING\tThe seperator for the -f file["tab" (default) or "whitespace"]\n',file=stderr())
  cat('\t-n\t\tINT\tA number represents the column number from the first column to be retained.[1]\n',file=stderr())
  cat('\t-o\t\tFILE\tOutput file name\n',file=stderr())
  cat('\t-h\t\tPrint this help information.\n',file=stderr())
  q(save="no")
}
if(length(args)==0 || args[1]=="-h"){
  usage()
}
seperator="tab"
outN=1
if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-f=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      inFile=arg.split[2]
    }
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      columnOut=arg.split[2]
    }
    if(grepl('^-s=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      seperator=arg.split[2]
    }
    if(grepl('^-n=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      outN=as.numeric(arg.split[2])
    }
    if(grepl('^-o=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      outFile=arg.split[2]
    }
    if(grepl('^-h',arg)){usage()}
  }
}else{
  usage()
}
if(seperator == "tab"){
  data<-read.delim(file=inFile,header=T,check.names=F)
}else if(seperator == "whitespace"){
  data<-read.table(file=inFile,header=T,sep=" ",check.names=F)
}else{
  cat('Please set the right seperator parameter\n',file=stderr())
  q(save="no")
}
columnNames<-colnames(data)
columnsOut<-as.vector(as.matrix(read.delim(file=columnOut,header=F)))
subdata<-data[,is.element(columnNames,columnsOut)]
write.table(cbind(data[,c(1:outN)],subdata),file=outFile,quote=F,sep="\t",row.names =F,col.names = c(columnNames[1:outN], colnames(subdata)))
