#!/bin/sh

## Required softwares and scripts
module load R/4.0.2
module load gcc/8.4.0
scriptPath=/dfs5/weil21-lab/yumeil1/scripts

## Change to the working directory
mkdir TCGA GTEx Immunotherapy

#1. at Immunotherapy/ directory
##1.1 Prepare read count matrix
	wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061%5FBMS038109Sample%2Ehg19KnownGene%2Eraw%2Ecsv%2Egz
	sed 's/,/\t/g' GSE91061_BMS038109Sample.hg19KnownGene.raw.csv >GSE91061_BMS038109Sample.hg19KnownGene.raw.tsv
	echo -e "library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	geneid=keys(TxDb.Hsapiens.UCSC.hg19.knownGene, "GENEID")
	idSymbol<-select(org.Hs.eg.db, geneid, c("SYMBOL"))
	write.table(idSymbol, file="geneID.symbol.map.tsv", sep="\t", quote=F, row.names=F, col.names=T)" | Rscript - 
	head -n1 GSE91061_BMS038109Sample.hg19KnownGene.raw.tsv >GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv
	perl ${scriptPath}/join.pl -i2 GSE91061_BMS038109Sample.hg19KnownGene.raw.tsv -f2 1 -i1 geneID.symbol.map.tsv -f1 1 |cut -f2,4- >>GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv
	awk '$1!="NA"' GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv >tmp;mv tmp GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv
	head -n1 GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv |cut -f2-|awk -v OFS="" -v ORS="" '{for(i=1;i<=NF;i++){split($i,a,"_");print a[2]"\t"}print "\n"}'|sed 's/\t$//' >GSE91061_BMS038109Sample.conditions.tsv

##1.2 Run DEG analysis for real datasets
	Rscript ${scriptPath}/DEGs.R -g=GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv -c=GSE91061_BMS038109Sample.conditions.tsv -s=d -o=GSE91061_BMS038109Sample.real.DESeq2.rst.tsv
	Rscript ${scriptPath}/DEGs.R -g=GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv -c=GSE91061_BMS038109Sample.conditions.tsv -s=e -o=GSE91061_BMS038109Sample.real.edgeR.rst.tsv
	Rscript ${scriptPath}/DEGs.R -g=GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv -c=GSE91061_BMS038109Sample.conditions.tsv -s=l -o=GSE91061_BMS038109Sample.real.limma.rst.tsv
	Rscript ${scriptPath}/DEGs.R -g=GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv -c=GSE91061_BMS038109Sample.conditions.tsv -s=n -o=GSE91061_BMS038109Sample.real.NOISeq.rst.tsv
	Rscript ${scriptPath}/DEGs.R -g=GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv -c=GSE91061_BMS038109Sample.conditions.tsv -s=w -o=GSE91061_BMS038109Sample.real.WilcoxonTest.rst.tsv

##1.3 Prepare permutated datasets and run DEG analysis
	echo -e "outMatrix=matrix(nrow=1000,ncol=109)
	set.seed(100)
	index=1
	for(i in 1:1000){
		x=rep(NA,109)
		x[sample(1:109,58)]="On"
		x[is.na(x)]="Pre"
		outMatrix[index,]=x
		index=index+1
	}
	write.table(outMatrix,file="shuffledLabel.basedReal.tsv",sep="\t", quote=F,row.names = F,col.names = F)" | Rscript -
	###Split to run parallelly
	split -l 100 -d shuffledLabel.basedReal.tsv shuffledLabel.basedReal.
	ls shuffledLabel.basedReal.0*|while read file;do
		part=$(echo $file|cut -f3 -d '.');
		mkdir $part
		mv $file $part
		###Run DEG analysis for each part
		cd $part
		index=1
		cat shuffledLabel.basedReal.0*|while read line;do
			echo $line >tmp.conditions.tsv
			Rscript $scriptPath/DEGs.R -g=../GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv -c=tmp.conditions.tsv -s=d -o=shuffledLabel.${index}.DESeq2.rst.tsv
			Rscript $scriptPath/DEGs.R -g=../GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv -c=tmp.conditions.tsv -s=e -o=shuffledLabel.${index}.edgeR.rst.tsv
			Rscript $scriptPath/DEGs.R -g=../GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv -c=tmp.conditions.tsv -s=w -o=shuffledLabel.${index}.WilcoxonTest.rst.tsv
			Rscript $scriptPath/DEGs.R -g=../GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv -c=tmp.conditions.tsv -s=l -o=shuffledLabel.${index}.limma.rst.tsv
			Rscript $scriptPath/DEGs.R -g=../GSE91061_BMS038109Sample.hg19KnownGene.symbol.tsv -c=tmp.conditions.tsv -s=n -o=shuffledLabel.${index}.NOISeq.rst.tsv
			index=$((index+1));
		done;
		cd ../
	done;
	Rscript ${scriptPath}/shuffledLabel_ratio.R -r=GSE91061_BMS038109Sample.conditions.tsv -s=shuffledLabel.basedReal.tsv -o=shuffledLabel.realRatio.tsv
	ls -d */|while read dir;do
		index=1;
		for N in {0..99}; do
			echo -e "${dir}shuffledLabel.${index}" >>shuffledLabel.rstFile.txt;
			index=$((index+1));
		done;
	done;
	paste shuffledLabel.rstFile.txt shuffledLabel.realRatio.tsv >shuffledLabel.realRatio.wtFileName.tsv
	firstQ=$(cut -f2 shuffledLabel.realRatio.wtFileName.tsv |Rscript ${scriptPath}/quantile_summary.R -c=1|tail -n1|awk '{print $2}')
	thirdQ=$(cut -f2 shuffledLabel.realRatio.wtFileName.tsv |Rscript ${scriptPath}/quantile_summary.R -c=1|tail -n1|awk '{print $5}')
	awk -v first=$firstQ -v third=$thirdQ '$2>=first && $2<=third{print $1}' shuffledLabel.realRatio.wtFileName.tsv >shuffledLabel.realRatio.1st-3rd.selected.txt
	###Summary results for permutated datasets
	##DESeq2
	awk '{print $0".DESeq2.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
		wc -l $file|awk '{print $1-1}' >>shuffledLabel.DESeq2.1st-3rd.number.txt
		cut -f1 $file |sed -n '2,$p' >>shuffledLabel.DESeq2.Gene.txt
		cat $file GSE91061_BMS038109Sample.real.DESeq2.rst.tsv |grep -v "baseMean"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>shuffledLabel.DESeq2.1st-3rd.overlapDEGs.txt
	done
	sort shuffledLabel.DESeq2.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >shuffledLabel.DESeq2.1st-3rd.uniqGene.txt
	##edgeR
  awk '{print $0".edgeR.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
      wc -l $file|awk '{print $1-1}' >>shuffledLabel.edgeR.1st-3rd.number.txt
      cut -f1 $file |sed -n '2,$p' >>shuffledLabel.edgeR.Gene.txt
      cat $file GSE91061_BMS038109Sample.real.edgeR.rst.tsv |grep -v "logFC"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>shuffledLabel.edgeR.1st-3rd.overlapDEGs.txt
  done
  sort shuffledLabel.edgeR.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >shuffledLabel.edgeR.1st-3rd.uniqGene.txt
  ##limma
  awk '{print $0".limma.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
    	wc -l $file|awk '{print $1-1}' >>shuffledLabel.limma.1st-3rd.number.txt
     grep -vE "^$|logFC" $file |cut -f1 >>shuffledLabel.limma.Gene.txt
     cat $file GSE91061_BMS038109Sample.real.limma.rst.tsv |grep -vE "^$|logFC"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>shuffledLabel.limma.1st-3rd.overlapDEGs.txt
  done
  sort shuffledLabel.limma.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >shuffledLabel.limma.1st-3rd.uniqGene.txt
	##NOISeq
	awk '{print $0".NOISeq.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
      wc -l $file|awk '{print $1-1}' >>shuffledLabel.NOISeq.1st-3rd.number.txt
      grep -vE "^$|log2FC" $file |cut -f1 >>shuffledLabel.NOISeq.Gene.txt
      cat $file GSE91061_BMS038109Sample.real.NOISeq.rst.tsv |grep -vE "^$|log2FC"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>shuffledLabel.NOISeq.1st-3rd.overlapDEGs.txt
  done
  sort shuffledLabel.NOISeq.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >shuffledLabel.NOISeq.1st-3rd.uniqGene.txt
  ##Wilcoxon test
  awk '{print $0"."WilcoxonTest.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
      wc -l $file|awk '{print $1-1}' >>shuffledLabel.WilcoxonTest.1st-3rd.number.txt
      cut -f1 $file |sed -n '2,$p' >>shuffledLabel.WilcoxonTest.Gene.txt
      cat $file GSE91061_BMS038109Sample.real.WilcoxonTest.rst.tsv |grep -v "Gene"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>shuffledLabel.WilcoxonTest.1st-3rd.overlapDEGs.txt
  done
  sort shuffledLabel.WilcoxonTest.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >shuffledLabel.WilcoxonTest.1st-3rd.uniqGene.txt

	###GO enrichment analysis
	awk '$1>66{print $2}' shuffledLabel.DESeq2.1st-3rd.uniqGene.txt|Rscript /dfs5/weil21-lab/yumeil1/scripts/geneEnrichAnalysis.R -s=hg19 -d=KEGG -n=10 -o=shuffledLabel.DESeq2.1st-3rd.DEGs.gt10percent.KEGG
	awk '$1>66{print $2}' shuffledLabel.DESeq2.1st-3rd.uniqGene.txt|Rscript /dfs5/weil21-lab/yumeil1/scripts/geneEnrichAnalysis.R -s=hg19 -d=GO -n=10 -o=shuffledLabel.DESeq2.1st-3rd.DEGs.gt10percent.GO
	awk '$1>66{print $2}' shuffledLabel.edgeR.1st-3rd.uniqGene.txt|Rscript /dfs5/weil21-lab/yumeil1/scripts/geneEnrichAnalysis.R -s=hg19 -d=KEGG -n=10 -o=shuffledLabel.edgeR.1st-3rd.DEGs.gt10percent.KEGG
  awk '$1>66{print $2}' shuffledLabel.edgeR.1st-3rd.uniqGene.txt|Rscript /dfs5/weil21-lab/yumeil1/scripts/geneEnrichAnalysis.R -s=hg19 -d=GO -n=10 -o=shuffledLabel.edgeR.1st-3rd.DEGs.gt10percent.GO

#2. at GTEx/ directory
GTExPath=/dfs5/weil21-lab/yumeil1/data/GTEx
##2.1 Prepare read count matrix and permutated conditions
	mkdir adipose brainAmygdala-spinalcord cellsEVBlym-MSgland heart prostate-brain wholeBlood-muscle
	cd adipose
	paste $GTExPath/subTissues/Adipose-Subcutaneous.wtGeno.gene_readCount.tsv $GTExPath/subTissues/Adipose-Visceral_Omentum.wtGeno.gene_readCount.tsv |cut -f1-583,586-| \
		awk -v OFS="" -v ORS="" '{print $1"_"$2;for(i=3;i<=NF;i++){print "\t"$i}print "\n"}' >AdiposeSubcutaneous.AdiposeVisceralOmentum.gene_readCount.tsv
	Rscript ../conditions.R 581 469 sub vis 100
	cd ../brainAmygdala-spinalcord
	paste $GTExPath/subTissues/Brain-Amygdala.wtGeno.gene_readCount.tsv $GTExPath/subTissues/Brain-Spinalcord_cervicalc-1.wtGeno.gene_readCount.tsv |cut -f1-131,134- | \
		awk -v OFS="" -v ORS="" '{print $1"_"$2;for(i=3;i<=NF;i++){print "\t"$i}print "\n"}' >BrainAmygdala.BrainSpinalcordCervicalc.gene_readCount.tsv
	Rscript ../conditions.R 129 126 amy spi 100
	cd ../cellsEVBlym-MSgland
	paste $GTExPath/subTissues/Cells-EBV-transformedlymphocytes.wtGeno.gene_readCount.tsv $GTExPath/subTissues/MinorSalivaryGland.wtGeno.gene_readCount.tsv|cut -f1-149,152-| \
		awk -v OFS="" -v ORS="" '{print $1"_"$2;for(i=3;i<=NF;i++){print "\t"$i}print "\n"}' >CellsEVBtransformedLymphocytes.MinorSalivaryGland.gene_readCount.tsv
	Rscript ../conditions.R 147 144 cell gland 100
	cd ../heart
	paste $GTExPath/subTissues/Heart-AtrialAppendage.wtGeno.gene_readCount.tsv $GTExPath/subTissues/Heart-LeftVentricle.wtGeno.gene_readCount.tsv | cut -f1-374,377-| \
		awk -v OFS="" -v ORS="" '{print $1"_"$2;for(i=3;i<=NF;i++){print "\t"$i}print "\n"}' >HeartAtrialAppendage.HeartLeftVentricle.gene_readCount.tsv
	Rscript ../conditions.R 372 386 atr ven 100
	cd ../wholeBlood-muscle
	paste $GTExPath/subTissues/WholeBlood.wtGeno.gene_readCount.tsv $GTExPath/subTissues/Muscle-Skeletal.wtGeno.gene_readCount.tsv|cut -f1-672,675-| \
		awk -v OFS="" -v ORS="" '{print $1"_"$2;for(i=3;i<=NF;i++){print "\t"$i}print "\n"}' >WholeBlood.MuscleSkeletal.gene_readCount.tsv
	Rscript ../conditions.R 670 706 blood muscle 100

##2.2 Run DEG analysis for real datasets
	ls -d */|while read dir;do
		cd $dir;
		prefix=$(ls *gene_readCount.tsv|sed 's/.gene_readCount.tsv//');
		Rscript $scriptPath/DEGs.R -g=${prefix}.gene_readCount.tsv -c=real.conditions.tsv -f=0.01 -s=d -o=${prefix}.DESeq2.0.01.rst.tsv
		Rscript $scriptPath/DEGs.R -g=${prefix}.gene_readCount.tsv -c=real.conditions.tsv -f=0.01 -s=e -o=${prefix}.edgeR.0.01.rst.tsv
		Rscript $scriptPath/DEGs.R -g=${prefix}.gene_readCount.tsv -c=real.conditions.tsv -f=0.01 -s=w -o=${prefix}.WilcoxonTest.0.01.rst.tsv
		Rscript $scriptPath/DEGs.R -g=${prefix}.gene_readCount.tsv -c=real.conditions.tsv -f=0.01 -s=l -o=${prefix}.limma.0.01.rst.tsv
		Rscript $scriptPath/DEGs.R -g=${prefix}.gene_readCount.tsv -c=real.conditions.tsv -f=0.01 -s=n -o=${prefix}.NOISeq.0.01.rst.tsv
		cd../
	done;
	
##2.3 Prepare permutated datasets and run DEG analysis
	ls -d */|while read dir;do
		cd $dir;
		split -l 100 -d shuffledLabel.basedReal.tsv shuffledLabel.basedReal.
		ls shuffledLabel.basedReal.0*|while read file;do
			part=$(echo $file|cut -f3 -d '.');
			mkdir $part
			mv $file $part
			###Run DEG analysis for permutated datasets
			cd $part
			index=1
			cat shuffledLabel.basedReal.0*|while read line;do
				echo $line >tmp.conditions.tsv
				readCount=$(ls ../*.gene_readCount.tsv);
				Rscript $scriptPath/DEGs.R -g=$readCount -c=tmp.conditions.tsv -f=0.01 -s=d -o=shuffledLabel.${index}.DESeq2.0.01.rst.tsv
				Rscript $scriptPath/DEGs.R -g=$readCount -c=tmp.conditions.tsv -f=0.01 -s=e -o=shuffledLabel.${index}.edgeR.0.01.rst.tsv
				Rscript $scriptPath/DEGs.R -g=$readCount -c=tmp.conditions.tsv -f=0.01 -s=w -o=shuffledLabel.${index}.WilcoxonTest.0.01.rst.tsv
				Rscript $scriptPath/DEGs.R -g=$readCount -c=tmp.conditions.tsv -f=0.01 -s=l -o=shuffledLabel.${index}.limma.0.01.rst.tsv
				Rscript $scriptPath/DEGs.R -g=$readCount -c=tmp.conditions.tsv -f=0.01 -s=n -o=shuffledLabel.${index}.NOISeq.0.01.rst.tsv
				index=$((index+1));
			done;
			cd ../
		done;
		cd ../
	done;
	
	ls -d */|while read tissue;do
		cd $tissue;
		ls -d */|while read dir;do
			index=1;
			for N in {0..99}; do
				echo -e "${dir}shuffledLabel.${index}" >>shuffledLabel.rstFile.txt;
				index=$((index+1));
			done;
		done;
		paste shuffledLabel.rstFile.txt shuffledLabel.realRatio.tsv >shuffledLabel.realRatio.wtFileName.tsv
		firstQ1=$(cut -f2 shuffledLabel.realRatio.wtFileName.tsv |Rscript /dfs5/weil21-lab/yumeil1/scripts/quantile_summary.R -c=1|tail -n1|awk '{print $2}')
		thirdQ1=$(cut -f2 shuffledLabel.realRatio.wtFileName.tsv |Rscript /dfs5/weil21-lab/yumeil1/scripts/quantile_summary.R -c=1|tail -n1|awk '{print $5}')
		firstQ2=$(cut -f3 shuffledLabel.realRatio.wtFileName.tsv |Rscript /dfs5/weil21-lab/yumeil1/scripts/quantile_summary.R -c=1|tail -n1|awk '{print $2}')
    thirdQ2=$(cut -f3 shuffledLabel.realRatio.wtFileName.tsv |Rscript /dfs5/weil21-lab/yumeil1/scripts/quantile_summary.R -c=1|tail -n1|awk '{print $5}')
		awk -v first1=$firstQ1 -v third1=$thirdQ1 -v first2=$firstQ2 -v third2=$thirdQ2 '$2>=first1 && $2<=third1 && $3>=first2 && $3<=third2{print $1}' shuffledLabel.realRatio.wtFileName.tsv >shuffledLabel.realRatio.1st-3rd.selected.txt
		###Summary results
		prefix=$(ls *.DESeq2.0.01.rst.tsv|sed 's/.DESeq2.0.01.rst.tsv//');
		fdrs=0.01
		mkdir FDR-$fdrs
		##DESeq2
    awk '{print $0".DESeq2.0.01.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
      awk -v fdr=$fdrs '$7<fdr' $file |wc -l |awk '{print $1-1}' >>FDR-$fdrs/shuffledLabel.DESeq2.number.txt
      awk -v fdr=$fdrs '$7<fdr{print $1}' $file |grep "ENSG" >>FDR-$fdrs/shuffledLabel.DESeq2.Gene.txt
      awk -v fdr=$fdrs '$7<fdr' $file |grep "ENSG" >FDR-$fdrs/tmp.DEGs
      cat FDR-$fdrs/tmp.DEGs ${prefix}.DESeq2.${fdrs}.rst.tsv |grep "^ENSG"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>FDR-$fdrs/shuffledLabel.DESeq2.overlapDEGs.txt
    done
    sort FDR-$fdrs/shuffledLabel.DESeq2.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >FDR-$fdrs/shuffledLabel.DESeq2.uniqGene.txt
		##edgeR
    awk '{print $0".edgeR.0.01.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
      awk -v fdr=$fdrs '$6<fdr' $file |wc -l |awk '{print $1-1}' >>FDR-$fdrs/shuffledLabel.edgeR.number.txt
      awk -v fdr=$fdrs '$6<fdr{print $1}' $file |grep "ENSG" >>FDR-$fdrs/shuffledLabel.edgeR.Gene.txt
      awk -v fdr=$fdrs '$6<fdr' $file |grep "ENSG" >FDR-$fdrs/tmp.DEGs
      cat FDR-$fdrs/tmp.DEGs ${prefix}.edgeR.${fdrs}.rst.tsv |grep "^ENSG"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>FDR-$fdrs/shuffledLabel.edgeR.overlapDEGs.txt
    done
    sort FDR-$fdrs/shuffledLabel.edgeR.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >FDR-$fdrs/shuffledLabel.edgeR.uniqGene.txt
    ##NOISeq
		awk '{print $0".NOISeq.0.01.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
      wc -l $file|awk '{print $1-1}' >>FDR-$fdrs/shuffledLabel.NOISeq.number.txt
      grep "ENSG" $file|cut -f1 >>FDR-$fdrs/shuffledLabel.NOISeq.Gene.txt
      cat $file ${prefix}.NOISeq.${fdrs}.rst.tsv |grep "^ENSG"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>FDR-${fdrs}/shuffledLabel.NOISeq.overlapDEGs.txt
    done
    sort FDR-${fdrs}/shuffledLabel.NOISeq.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >FDR-${fdrs}/shuffledLabel.NOISeq.uniqGene.txt
		##Wilcoxon
    awk '{print $0".WilcoxonTest.0.01.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
      awk -v fdr=$fdrs '$3<fdr' $file |wc -l >>FDR-$fdrs/shuffledLabel.WilcoxonTest.number.txt
      awk -v fdr=$fdrs '$3<fdr{print $1}' $file |grep "ENSG" >>FDR-$fdrs/shuffledLabel.WilcoxonTest.Gene.txt
      awk -v fdr=$fdrs '$3<fdr' $file |grep "ENSG" >FDR-$fdrs/tmp.DEGs
      cat FDR-$fdrs/tmp.DEGs ${prefix}.WilcoxonTest.$fdrs.rst.tsv |grep "^ENSG"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>FDR-$fdrs/shuffledLabel.WilcoxonTest.overlapDEGs.txt
		done
		sort FDR-$fdrs/shuffledLabel.WilcoxonTest.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >FDR-$fdrs/shuffledLabel.WilcoxonTest.uniqGene.txt
		##Limma
		awk '{print $0".limma.0.01.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
			grep "^ENSG" $file|awk -v fdr=$fdrs '$6<fdr'|wc -l >>FDR-$fdrs/shuffledLabel.limma.number.txt
			awk -v fdr=$fdrs '$6<fdr{print $1}' $file |grep "ENSG" >>FDR-$fdrs/shuffledLabel.limma.Gene.txt
      awk -v fdr=$fdrs '$6<fdr' $file |grep "ENSG" >FDR-$fdrs/tmp.DEGs
      cat FDR-$fdrs/tmp.DEGs ${prefix}.limma.$fdrs.rst.tsv |grep "^ENSG"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>FDR-$fdrs/shuffledLabel.limma.overlapDEGs.txt
		done
    sort FDR-$fdrs/shuffledLabel.limma.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >FDR-$fdrs/shuffledLabel.limma.uniqGene.txt
		###GO enrichment analysis
    awk -v OFS="\t" '{print $1,int($2*0.1)}' shuffledLabel.realRatio.1st-3rd.number.tsv|while read sample minNum;do
      awk -v number=$minNum '$1>number{print $2}' FDR-$fdrs/shuffledLabel.edgeR.uniqGene.txt|cut -f2 -d '_'|Rscript $scriptPath/geneEnrichAnalysis.R -s=hg19 -d=GO -n=10 -o=FDR-$fdrs/shuffledLabel.edgeR.DEGs.gt10percent.GO
      awk -v number=$minNum '$1>number{print $2}' FDR-$fdrs/shuffledLabel.DESeq2.uniqGene.txt|cut -f2 -d '_'|Rscript $scriptPath/geneEnrichAnalysis.R -s=hg19 -d=GO -n=10 -o=FDR-$fdrs/shuffledLabel.DESeq2.DEGs.gt10percent.GO
    done
		cd ../
	done;
						
#3. at TCGA/ directory
##3.1 Prepare read count matrix and permutated conditions
	echo -e "BRCA\nKIRC\nLIHC\nLUAD\nPRAD\nTHCA" >sampleGt50.txt
	cat sampleGt50.txt |while read dir;do
		cd $dir;
		wget https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-${dir}.htseq_counts.tsv.gz
		zcat *htseq_counts.tsv.gz |head -n1|sed 's/\t/\n/g'|sed -n '2,$p'|awk '$0~"11A$|11B$|11C$"' >TCGA-${dir}.11.samples.txt;
		cut -f1-3 -d '-' TCGA-${dir}.11.samples.txt >tmp;
		zcat *htseq_counts.tsv.gz |head -n1|cut -f2-|awk '{for(i=1;i<=NF;i++){print $i}}'|grep -f tmp >TCGA-${dir}.normal-tumor.pair.txt
		cd ../
	done;
	###Remove technical replicates for some samples
	grep -vE "TCGA-44-6144|TCGA-55-6969|01B$" LUAD/TCGA-LUAD.normal-tumor.pair.txt >tmp;mv tmp LUAD/TCGA-LUAD.normal-tumor.pair.txt
	grep -v "01B$" PRAD/TCGA-PRAD.normal-tumor.pair.txt >tmp;mv tmp PRAD/TCGA-PRAD.normal-tumor.pair.txt
	cat sampleGt50.txt |while read dir;do
		cd $dir;
		zcat *htseq_counts.tsv.gz >tmp
		Rscript $scriptPath/columnSelect.R -f=tmp -i=TCGA-${dir}.normal-tumor.pair.txt -n=1 -o=TCGA-${dir}.normal-tumor.pair.Count.tsv
		Rscript ../rawCount.R TCGA-${dir}.normal-tumor.pair.Count.tsv TCGA-${dir}.normal-tumor.pair.rawCount.tsv
		###Remove Meta tag from the read count file
		grep -v "^__" TCGA-${dir}.normal-tumor.pair.rawCount.tsv >tmp;mv tmp TCGA-${dir}.normal-tumor.pair.rawCount.tsv
		head -n1 TCGA-${dir}.normal-tumor.pair.rawCount.tsv|cut -f2-|awk -v ORS="" -v OFS="" '{
			for(i=1;i<NF;i++){
				split($i,a,"-");
				if(a[4]~"01"){print "tumor\t"}else{print "normal\t"}
			}
			split($NF,b,"-");
			if(b[4]~"01"){print "tumor\n"}else{print "normal\n"}
		}' >TCGA-${dir}.conditions.tsv
		sampleNum=$(wc -l TCGA-${dir}.normal-tumor.pair.txt|cut -f1 -d ' ');
		sampleNum=$(($sampleNum/2))
		Rscript ../conditions.R $sampleNum $sampleNum normal tumor 200
		cd ../;
	done;

##3.2 Run DEG analysis for real datasets
	cat sampleGt50.txt |while read dir;do
		cd $dir;
		Rscript $scriptPath/DEGs.R -g=TCGA-${dir}.normal-tumor.pair.rawCount.tsv -c=TCGA-${dir}.conditions.tsv -f=0.01 -s=d -o=TCGA-${dir}.DESeq2.0.01.rst.tsv
		Rscript $scriptPath/DEGs.R -g=TCGA-${dir}.normal-tumor.pair.rawCount.tsv -c=TCGA-${dir}.conditions.tsv -f=0.01 -s=e -o=TCGA-${dir}.edgeR.0.01.rst.tsv
		Rscript $scriptPath/DEGs.R -g=TCGA-${dir}.normal-tumor.pair.rawCount.tsv -c=TCGA-${dir}.conditions.tsv -f=0.01 -s=l -o=TCGA-${dir}.limma.0.01.rst.tsv
		Rscript $scriptPath/DEGs.R -g=TCGA-${dir}.normal-tumor.pair.rawCount.tsv -c=TCGA-${dir}.conditions.tsv -f=0.01 -s=l -o=TCGA-${dir}.NOISeq.0.01.rst.tsv
		Rscript $scriptPath/DEGs.R -g=TCGA-${dir}.normal-tumor.pair.rawCount.tsv -c=TCGA-${dir}.conditions.tsv -f=0.01 -s=w -o=TCGA-${dir}.WilcoxonTest.0.01.rst.tsv
		cd ../
	done;
	
##3.3 Prepare permutated datasets and run DEG analysis
	cat sampleGt50.txt|while read dir;do
		cd $dir;
		split -l 100 -d shuffledLabel.basedReal.tsv shuffledLabel.basedReal.
		ls shuffledLabel.basedReal.0*|while read file;do
			part=$(echo $file|cut -f3 -d '.');
			mkdir $part
			mv $file $part
			cd $part
			index=1
			cat shuffledLabel.basedReal.0*|while read line;do
				echo $line >tmp.conditions.tsv
				readCount=$(ls ../*.rawCount.tsv);
				Rscript $scriptPath/DEGs.R -g=$readCount -c=tmp.conditions.tsv -f=0.01 -s=d -o=shuffledLabel.${index}.DESeq2.0.01.rst.tsv
				Rscript $scriptPath/DEGs.R -g=$readCount -c=tmp.conditions.tsv -f=0.01 -s=e -o=shuffledLabel.${index}.edgeR.0.01.rst.tsv
				Rscript $scriptPath/DEGs.R -g=$readCount -c=tmp.conditions.tsv -f=0.01 -s=l -o=shuffledLabel.${index}.limma.0.01.rst.tsv
				Rscript $scriptPath/DEGs.R -g=$readCount -c=tmp.conditions.tsv -f=0.01 -s=w -o=shuffledLabel.${index}.WilcoxonTest.0.01.rst.tsv
				Rscript $scriptPath/DEGs.R -g=$readCount -c=tmp.conditions.tsv -f=0.01 -s=n -o=shuffledLabel.${index}.NOISeq.0.01.rst.tsv
				index=$((index+1));
			done;
			cd ../
		done;
		cd ../
	done;
	###Summary results of permutated datasets
	cat sampleGt50.txt|while read tissue;do
		cd $tissue;
		ls -d */|while read dir;do
			index=1;
			for N in {0..99}; do
				echo -e "${dir}shuffledLabel.${index}" >>shuffledLabel.rstFile.txt;
				index=$((index+1));
			done;
		done;
		paste shuffledLabel.rstFile.txt shuffledLabel.realRatio.tsv >shuffledLabel.realRatio.wtFileName.tsv
		firstQ=$(cut -f2 shuffledLabel.realRatio.wtFileName.tsv |Rscript /dfs5/weil21-lab/yumeil1/scripts/quantile_summary.R -c=1|tail -n1|awk '{print $2}')
		thirdQ=$(cut -f2 shuffledLabel.realRatio.wtFileName.tsv |Rscript /dfs5/weil21-lab/yumeil1/scripts/quantile_summary.R -c=1|tail -n1|awk '{print $5}')
		awk -v first=$firstQ -v third=$thirdQ '$2>=first && $2<=third{print $1}' shuffledLabel.realRatio.wtFileName.tsv >shuffledLabel.realRatio.1st-3rd.selected.txt
		fdrs=0.01
		mkdir FDR-$fdrs;
		##DESeq2
		awk '{print $0".DESeq2.0.01.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
			awk -v fdr=$fdrs '$7<fdr' $file |wc -l |awk '{print $1-1}' >>FDR-$fdrs/shuffledLabel.DESeq2.number.txt
			awk -v fdr=$fdrs '$7<fdr{print $1}' $file |grep "ENSG" >>FDR-$fdrs/shuffledLabel.DESeq2.Gene.txt
			awk -v fdr=$fdrs '$7<fdr' $file |grep "ENSG" >FDR-$fdrs/tmp.DEGs
			cat FDR-$fdrs/tmp.DEGs TCGA-${tissue}.DESeq2.$fdrs.rst.tsv |grep "^ENSG"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>FDR-$fdrs/shuffledLabel.DESeq2.overlapDEGs.txt
		done
		sort FDR-$fdrs/shuffledLabel.DESeq2.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >FDR-$fdrs/shuffledLabel.DESeq2.uniqGene.txt
		##edgeR
		awk '{print $0".edgeR.0.01.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
			awk -v fdr=$fdrs '$6<fdr' $file |wc -l |awk '{print $1-1}' >>FDR-$fdrs/shuffledLabel.edgeR.number.txt
      awk -v fdr=$fdrs '$6<fdr{print $1}' $file |grep "ENSG" >>FDR-$fdrs/shuffledLabel.edgeR.Gene.txt
      awk -v fdr=$fdrs '$6<fdr' $file |grep "ENSG" >FDR-$fdrs/tmp.DEGs
			cat FDR-$fdrs/tmp.DEGs TCGA-${tissue}.edgeR.$fdrs.rst.tsv |grep "^ENSG"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>FDR-$fdrs/shuffledLabel.edgeR.overlapDEGs.txt
		done
		sort FDR-$fdrs/shuffledLabel.edgeR.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >FDR-$fdrs/shuffledLabel.edgeR.uniqGene.txt
		##NOIseq
		awk -v fdr=$fdrs '{print $0".NOISeq.0.01.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
      wc -l $file|awk '{print $1-1}' >>FDR-$fdrs/shuffledLabel.NOISeq.number.txt
      grep "ENSG" $file|cut -f1 >>FDR-$fdrs/shuffledLabel.NOISeq.Gene.txt
      cat $file TCGA-${tissue}.NOISeq.$fdrs.rst.tsv |grep "^ENSG"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>FDR-$fdrs/shuffledLabel.NOISeq.overlapDEGs.txt
    done
    sort FDR-$fdrs/shuffledLabel.NOISeq.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >FDR-$fdrs/shuffledLabel.NOISeq.uniqGene.txt
		##Wilcoxon
		awk '{print $0".WilcoxonTest.0.01.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
      awk -v fdr=$fdrs '$3<fdr' $file |wc -l >>FDR-$fdrs/shuffledLabel.WilcoxonTest.number.txt
      awk -v fdr=$fdrs '$3<fdr{print $1}' $file |grep "ENSG" >>FDR-$fdrs/shuffledLabel.WilcoxonTest.Gene.txt
      awk -v fdr=$fdrs '$3<fdr' $file |grep "ENSG" >FDR-$fdrs/tmp.DEGs
      cat FDR-$fdrs/tmp.DEGs TCGA-${tissue}.WilcoxonTest.$fdrs.rst.tsv |grep "^ENSG"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>FDR-$fdrs/shuffledLabel.WilcoxonTest.overlapDEGs.txt
    done
		sort FDR-$fdrs/shuffledLabel.WilcoxonTest.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >FDR-$fdrs/shuffledLabel.WilcoxonTest.uniqGene.txt
		##limma
		awk '{print $0".limma.rst.tsv"}' shuffledLabel.realRatio.1st-3rd.selected.txt |while read file;do
      awk -v fdr=$fdrs '$6<fdr' $file |wc -l |awk '{print $1-1}' >>FDR-$fdrs/shuffledLabel.limma.number.txt
      awk -v fdr=$fdrs '$6<fdr{print $1}' $file |grep "ENSG" >>FDR-$fdrs/shuffledLabel.limma.Gene.txt
      awk -v fdr=$fdrs '$6<fdr' $file |grep "ENSG" >FDR-$fdrs/tmp.DEGs
      cat FDR-$fdrs/tmp.DEGs TCGA-${tissue}.limma.$fdrs.rst.tsv |grep "^ENSG"|cut -f1|sort|uniq -c|awk '$1>1'|wc -l >>FDR-$fdrs/shuffledLabel.limma.overlapDEGs.txt
    done
    sort FDR-$fdrs/shuffledLabel.limma.Gene.txt|uniq -c |sort -k1,1nr|awk '{print $1"\t"$2}' >FDR-$fdrs/shuffledLabel.limma.uniqGene.txt
		cd ../
	done;
	##GO enrichemnt analysis
	awk -v OFS="\t" '{print $1,int($2*0.1)}' shuffledLabel.realRatio.1st-3rd.number.tsv|while read sample minNum;do
		awk -v number=$minNum '$1>number{print $2}' $sample/FDR-0.01/shuffledLabel.edgeR.uniqGene.txt|cut -f1 -d '.'|Rscript $scriptPath/geneEnrichAnalysis.R -k=ENSEMBL -s=hg19 -d=GO -n=10 -o=$sample/FDR-0.01/shuffledLabel.edgeR.DEGs.gt10percent.GO
		awk -v number=$minNum '$1>number{print $2}' $sample/FDR-0.01/shuffledLabel.DESeq2.uniqGene.txt|cut -f1 -d '.'|Rscript $scriptPath/geneEnrichAnalysis.R -k=ENSEMBL -s=hg19 -d=GO -n=10 -o=$sample/FDR-0.01/shuffledLabel.DESeq2.DEGs.gt10percent.GO
	done

