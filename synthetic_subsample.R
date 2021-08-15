library(DESeq2)
library(edgeR)
library(parallel)
library(NOISeq)
myedger = function(count1, count2, q = 0.05){
  n1 = ncol(count1)
  n2 = ncol(count2)
  cond_idx = rep(2, n1 + n2)
  cond_idx[1:n1] = 1
  dat = cbind(count1, count2)
  if(!is.factor(cond_idx)){
    cond_idx = factor(cond_idx)
  }
  y <- DGEList(counts=dat,group=cond_idx)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  count_norm = cpm(y)
  design <- model.matrix(~cond_idx)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  qlf_i = topTags(qlf, n = nrow(dat), p.value = 1)@.Data[[1]]
  pvalues = qlf_i[match(rownames(count1), rownames(qlf_i)),"PValue"]
  qlf <- topTags(qlf, n = nrow(dat), p.value = q)
  discovery = rownames(qlf)
  output <- list(discovery, pvalues)
  return(output)
}

mylimma = function(count1, count2, q = 0.05){
  n1 = ncol(count1)
  n2 = ncol(count2)
  cond_idx = rep(2, n1 + n2)
  cond_idx[1:n1] = 1
  dat = cbind(count1, count2)
  if(!is.factor(cond_idx)){
    cond_idx = factor(cond_idx)
  }
  y <- DGEList(counts=dat,group=cond_idx)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  count_norm = cpm(y)
  design <- model.matrix(~cond_idx)
  v <- voom(y, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  res=topTable(fit, coef=2, sort.by="P", number=Inf, p.value = 1)
  pvalues = res[match(rownames(count1), rownames(res)),"P.Value"]
  discovery = rownames(res)[which(p.adjust(res$P.Value, method = "BH")< q)]
  output <- list(discovery, pvalues)
  return(output)
}

mydeseq2 = function(count1, count2, q=0.05){
  
  n1 = ncol(count1)
  n2 = ncol(count2)
  cond_idx = rep(2, n1 + n2)
  cond_idx[1:n1] = 1
  dat = cbind(count1, count2)
  rownames(dat) <- 1:nrow(dat)
  if(!is.factor(cond_idx)){
    cond_idx = factor(cond_idx)
  }
  dds <- DESeqDataSetFromMatrix(dat, DataFrame(cond_idx), ~ cond_idx)
  dds <- DESeq(dds)
  res <- results(dds, alpha = q)
  pvalues <- res$padj
  discovery = rownames(res)[which(res$padj <= q)]
  return(list(discovery, pvalues))
}
my.wilcoxon <- function(score_exp, score_back){
  p <- sapply(1:nrow(score_exp), function(j){
    return(wilcox.test(as.numeric(score_exp[j,]), as.numeric(score_back[j,]), alternative = "two.sided")$p.value)
  })
}
my.NOISeq <- function(readCount, conditions){
  data<-NOISeq::readData(data=readCount, factors=as.data.frame(conditions))
  res=noiseqbio(data, k=0.5, norm="tmm", factor="conditions",random.seed = 12345, filter = 1, cv.cutoff = 100, cpm = 1)
  outputRst=degenes(res, M=NULL)
  q = rep(NA, nrow(readCount))
  q[as.numeric(rownames(outputRst))] = 1 - outputRst$prob
  return(q)
}
my.ttest <- function(score_exp, score_back){
  p <- sapply(1:nrow(score_exp), function(j){
    return(t.test(as.numeric(score_exp[j,]), as.numeric(score_back[j,]), alternative = "two.sided")$p.value)
  })
}

find_discovery_wpval = function(pval, q){
  pval.adj = p.adjust(pval, method = 'BH')
  discovery_ls = lapply(q, function(q_i){
    discovery = which(pval.adj <= q_i)
  })
  return(discovery_ls)
}
compute_fdppow = function(discovery_ls, trueidx){
  
  fdp = sum(!discovery_ls %in% trueidx )/max(length(discovery_ls),1)
  pow = sum(discovery_ls %in% trueidx)/length(trueidx)
  return(c(fdp, pow))
}
compute_fdppow2 = function(pvalues, trueidx, q){
  index.p <- which(!is.na(pvalues))
  is.trueDE <- 1:length(pvalues) %in% trueidx
  n.dis <- cumsum(is.trueDE[order(pvalues)])
  FDR.i = 1 - n.dis/1:length(pvalues)
  index.i <- max(which(FDR.i<=q))
  fdp <- FDR.i[index.i]
  pow <- n.dis[index.i]/length(trueidx)
  return(c(fdp, pow))
}


celltype.pairs <- c("AdiposeSubcutaneous.AdiposeVisceralOmentum", "BrainAmygdala.BrainSpinalcordCervicalc",
                    "CellsEVBtransformedLymphocytes.MinorSalivaryGland", "HeartAtrialAppendage.HeartLeftVentricle", 
                    "WholeBlood.MuscleSkeletal", "Brain.Prostate")
tcga.pairs <- c('TCGA-BRCA', 'TCGA-KIRC', 'TCGA-LIHC', 'TCGA-LUAD', 'TCGA-PRAD', 'TCGA-THCA')

#n.sample <- c(5, 10, 15, 20, 25, 30, 35, 40)

n.sample <- c(2:10, 12, 15, 20, 25, 30, 35, 40)
for(j in 1:6){
  count.mat <- read.table(paste0(tcga.pairs[j], '.normal-tumor.pair.rawCount.tsv.1'), sep = '\t', header = TRUE)
  conditions <- read.table(paste0(tcga.pairs[j],'.conditions.tsv'), sep = '\t')
  
  conditions <- substring(conditions[1,], 1, 1)
  condition.names <- unique(conditions)
  count1 <- count.mat[, which(conditions == condition.names[1])+1]
  count2 <- count.mat[, which(conditions == condition.names[2])+1]
  
  dataset <- cbind(count1, count2)
  conditions2 <- c(rep("1", ncol(count1)),rep("2", ncol(count2)))
  
  # y <- DGEList(dataset)
  # y=calcNormFactors(y)
  # cpmtmm=cpm(y)
  # score_exp <- log(cpmtmm[,1:(ncol(count1))]+1)
  # score_back <- log(cpmtmm[,-(1:(ncol(count1)))]+1)
  # re_edger <- myedger(count1, count2, q = 0.000001)
  # dis_edger <- re_edger[[1]]
  # re_limma <- mylimma(count1, count2, q = 0.000001)
  # dis_limma <- re_limma[[1]]
  # re_deseq2 <- mydeseq2(count1, count2, q= 0.000001)
  # dis_deseq2 <- re_deseq2[[1]]
  # re_wilxocon <- my.wilcoxon(score_exp, score_back)
  # re_wilxocon <- find_discovery_wpval(re_wilxocon, q = 0.000001)[[1]]
  # re_noiseq <- my.NOISeq(dataset, conditions2)
  # dis_noiseq <- which(re_noiseq <= 0.000001)
  # dis0 <- list(dis_edger, dis_limma, dis_deseq2, re_wilxocon, dis_noiseq)
  # saveRDS(dis0, paste0(tcga.pairs[j], '_dis.rds'))
  dis0 <- readRDS(paste0(tcga.pairs[j], '_dis.rds'))
  trueDE = intersect(intersect(intersect(intersect(dis0[[1]], dis0[[2]]), dis0[[3]]), dis0[[4]]), dis0[[5]])
  
  ls_subsample <- lapply(n.sample, function(n){
    print(n)
    ls_onepair <- mclapply(1:25, function(it){
      set.seed(it)
      dataset0 <- dataset
      index0 <- c(sample(1:ncol(count1), n), sample((ncol(count1)+1):ncol(dataset0), n))
      dataset0 <- dataset[,index0]
      trueDE0 <- sample(trueDE, length(trueDE)/2)
      fdppow <- list()
      fdppow.i <- list()
      dataset0[-trueDE0,] <- t(apply(dataset0[-trueDE0, ], 1, sample))
      
      #dataset0 <- t(apply(dataset0, 1, sample))
      
      #dataset0 <- dataset0[, sample(1:ncol(dataset), ncol(dataset))]
      count.perm1 <- dataset0[,1:n]
      count.perm2 <- dataset0[,(n+1):(2*n)]
      conditions2 <- c(rep("1", n),rep("2", n))
      thres <- c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001)
      
      re_edger <- myedger(count.perm1, count.perm2)
      dis_edger <- re_edger[[1]]
      p_edger <- re_edger[[2]]
      fdppow[[1]] <- sapply(1:6, function(j){
        dis0 <- which(p.adjust(p_edger, method = "BH")<=thres[j])
        compute_fdppow(dis0, trueDE0)
      })
      fdppow.i[[1]] <- sapply(1:6, function(j){
        compute_fdppow2(p_edger, trueDE0, thres[j])
      })
      
      re_limma <- mylimma(count.perm1, count.perm2)
      dis_limma <- re_limma[[1]]
      p_limma <- re_limma[[2]]
      fdppow[[2]] <- sapply(1:6, function(j){
        dis0 <- which(p.adjust(p_limma, method = "BH")<=thres[j])
        compute_fdppow(dis0, trueDE0)
      })
      fdppow.i[[2]] <- sapply(1:6, function(j){
        compute_fdppow2(p_limma, trueDE0, thres[j])
      })
      
      re_deseq2 <- suppressMessages(try(mydeseq2(count.perm1, count.perm2)))
      if (class(re_deseq2) == "try-error"){
        fdppow[[3]] <- NA
        fdppow.i[[3]] <- NA
      }else{
        fdppow[[3]] <-sapply(1:6, function(j){
          dis0 <- which(re_deseq2[[2]] <=thres[j])
          compute_fdppow(dis0, trueDE0)
        })
        fdppow.i[[3]] <- sapply(1:6, function(j){
          compute_fdppow2(re_deseq2[[2]], trueDE0, thres[j])
        })
      }
      re_wilxocon <- my.wilcoxon(count.perm1, count.perm2)
      fdppow[[4]] <- sapply(1:6, function(j){
        dis0 <- which(p.adjust(re_wilxocon, method = "BH")<=thres[j])
        compute_fdppow(dis0, trueDE0)
      })
      fdppow.i[[4]] <- sapply(1:6, function(j){
        compute_fdppow2(re_wilxocon, trueDE0, thres[j])
      })
      re_noiseq <- suppressMessages(my.NOISeq(dataset0, conditions2))
      fdppow[[5]] <- sapply(1:6, function(j){
        dis0 <- which(re_noiseq<=thres[j])
        compute_fdppow(dis0, trueDE0)
      })
      fdppow.i[[5]] <- sapply(1:6, function(j){
        compute_fdppow2(re_noiseq, trueDE0, thres[j])
      })
      return(list(fdppow, fdppow.i))
    }, mc.cores = 25)
    return(ls_onepair)
  })
  saveRDS(ls_subsample, paste0(tcga.pairs[j], "subsample.rds"))
}



  
n.sample <- c(2:10, 12, 15, 20, 40, 60, 80, 100)
for(j in 1:6){
  if (j <=5){
    count.mat <- read.table(paste0(celltype.pairs[j], '.gene_readCount.tsv'), sep = '\t', header = TRUE)
    
    conditions <- read.table(paste0(celltype.pairs[j],'.conditions.tsv'), sep = '\t')
    
    conditions <- substring(conditions[1,], 1, 1)
    condition.names <- unique(conditions)
    count1 <- count.mat[, which(conditions == condition.names[1])+1]
    count2 <- count.mat[, which(conditions == condition.names[2])+1]
    
  }else{
    Prostate <- read.table('Prostate.wtGeno.gene_readCount.tsv', sep = '\t', header = TRUE)
    Brain <- read.table('Brain-Cortex.wtGeno.gene_readCount.tsv', sep = '\t', header = TRUE)
    count1 <- Prostate[,-c(1:2)]
    count2 <- Brain[,-c(1:2)]
    conditions <- rep(2, ncol(count1) + ncol(count2))
    conditions[1:ncol(count1)] = 1
    conditions = factor(conditions)
  }
  dataset <- cbind(count1, count2)
  dis0 <- readRDS(paste0(celltype.pairs[j], '_dis.rds'))
  trueDE = intersect(intersect(intersect(intersect(dis0[[1]], dis0[[2]]), dis0[[3]]), dis0[[4]]), dis0[[5]])
  ls_subsample <- lapply(n.sample, function(n){
    
    ls_onepair <- mclapply(1:25, function(it){
      set.seed(it)
      dataset0 <- dataset
      index0 <- c(sample(1:ncol(count1), n), sample((ncol(count1)+1):ncol(dataset0), n))
      dataset0 <- dataset[,index0]
      trueDE0 <- sample(trueDE, length(trueDE)/2)
      fdppow <- list()
      fdppow.i <- list()
      dataset0[-trueDE0,] <- t(apply(dataset0[-trueDE0, ], 1, sample))
      
      #dataset0 <- t(apply(dataset0, 1, sample))
      
      #dataset0 <- dataset0[, sample(1:ncol(dataset), ncol(dataset))]
      count.perm1 <- dataset0[,1:n]
      count.perm2 <- dataset0[,(n+1):(2*n)]
      thres <- c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001)
      conditions2 <- c(rep("1", n),rep("2", n))
      
      re_edger <- myedger(count.perm1, count.perm2)
      dis_edger <- re_edger[[1]]
      p_edger <- re_edger[[2]]
      fdppow[[1]] <- sapply(1:6, function(j){
        dis0 <- which(p.adjust(p_edger, method = "BH")<=thres[j])
        compute_fdppow(dis0, trueDE0)
      })
      fdppow.i[[1]] <- sapply(1:6, function(j){
        compute_fdppow2(p_edger, trueDE0, thres[j])
      })
      
      re_limma <- mylimma(count.perm1, count.perm2)
      dis_limma <- re_limma[[1]]
      p_limma <- re_limma[[2]]
      fdppow[[2]] <- sapply(1:6, function(j){
        dis0 <- which(p.adjust(p_limma, method = "BH")<=thres[j])
        compute_fdppow(dis0, trueDE0)
      })
      fdppow.i[[2]] <- sapply(1:6, function(j){
        compute_fdppow2(p_limma, trueDE0, thres[j])
      })
      
      re_deseq2 <- try(mydeseq2(count.perm1, count.perm2))
      if (class(re_deseq2) == "try-error"){
        fdppow[[2]] <- NA
        fdppow.i[[2]] <- NA
      }else{
        fdppow[[3]] <-sapply(1:6, function(j){
          dis0 <- which(re_deseq2[[2]] <=thres[j])
          compute_fdppow(dis0, trueDE0)
        })
        fdppow.i[[3]] <- sapply(1:6, function(j){
          compute_fdppow2(re_deseq2[[2]], trueDE0, thres[j])
        })
      }
      re_wilxocon <- my.wilcoxon(count.perm1, count.perm2)
      fdppow[[4]] <- sapply(1:6, function(j){
        dis0 <- which(p.adjust(re_wilxocon, method = "BH")<=thres[j])
        compute_fdppow(dis0, trueDE0)
      })
      fdppow.i[[4]] <- sapply(1:6, function(j){
        compute_fdppow2(re_wilxocon, trueDE0, thres[j])
      })
      re_noiseq <- my.NOISeq(dataset0, conditions2)
      fdppow[[5]] <- sapply(1:6, function(j){
        dis0 <- which(re_noiseq<=thres[j])
        compute_fdppow(dis0, trueDE0)
      })
      fdppow.i[[5]] <- sapply(1:6, function(j){
        compute_fdppow2(re_noiseq, trueDE0, thres[j])
      })
      return(list(fdppow, fdppow.i))
    }, mc.cores = 25)
    return(ls_onepair)
  })
  saveRDS(ls_subsample, paste0(celltype.pairs[j], "subsample.rds"))
}

