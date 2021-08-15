library(parallel)
library(locfdr)
library(scales)
library(truncdist)
library(ggplot2)
library(wesanderson)
library(egg)
library(grid)
library(Rgb)
library(scales)
library(stringr)
penal_size_width = unit(1.32, "in")
penal_size_height = unit(1.32, "in")

method_tot = c("edgeR", "Limma", "DESeq2", "Wilcoxon", "NOISeq")
col_wilcoxon = rgb(0, 104, 178, max = 255) # blue
col_edger = rgb(200, 140, 210, max = 255) # dark brown
col_deseq2 = rgb(226, 186, 73, max = 255) # orange
col_limma = rgb(200, 215, 144, max = 255) # ginger 
col_noiseq = rgb(130, 209, 169, max = 255) # ginger 
col_palette = c(col_edger, col_limma, col_deseq2, col_wilcoxon, col_noiseq)
names(col_palette) = method_tot
lty_palette = c("dashed", "dashed", "dashed", "solid", "solid")
names(lty_palette) = method_tot
  
mydarkgrey = "#252525"
FDR  = 0.05

axissize = 9
textsize = 9
titletextsize = 15
barwidth = 1
ylim_fdr = 0.5
label_height = 0.01
ylim_pow = 1
lwd = 0.5
pt_size = 1
pt_shape = 1:5

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

transformp = function(x){
  log10(x + 1e-6)
}

transformpow = function(x){
  log10(x + 1)
}

transformn = function(x){
  log2(x)
}

plabels = c(0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.2)
pticks = transformp(plabels)
fdrlabels = c("0", ".001", ".01", ".1", "1", "5", "20")


nlabels = c(2, 4, 8, 20, 40, 100)
nticks = transformn(nlabels)

powlabels = c(0, 0.05, 0.1, 0.2, 0.4, 1)
powticks = transformpow(plabels)


format_text = function(x){
  re = sapply(x, function(x_i){
    if(x_i < 10){
      return(as.character(round(x_i, 1)))
    }else{
      return(as.character(round(x_i)))
    }
  })
  return(re)
}
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

dat2plot <- function(dat, ngene){
  p = ggplot(dat, 
             aes(x = transformp(q), y = transformp(fdp), group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text( colour = 1, size = axissize),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = 'none',
          # axis.title.x = element_blank(),
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA)) +
    scale_x_continuous(limits = c(-6.1, 0),
                       breaks = pticks,
                       expand = c(0,0),
                       labels = fdrlabels) +
    scale_y_continuous(limits = c(-6.1, 0),
                       breaks = pticks,
                       expand = c(0,0),
                       labels = fdrlabels) +
    # coord_cartesian( ylim=c(0, 1.1)) +
    geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(x = transformp(q), y = transformp(fdp), col = methods, linetype = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    scale_linetype_manual(values=lty_palette) +
    geom_point(aes(x = transformp(q), y = transformp(fdp), shape = methods, col = methods), stroke = 1,size = pt_size) +
    scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)") + ylab('Actual FDR (%)')
  p1 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  
  p = ggplot(dat, 
             aes(x = transformp(q), y = power, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text( colour = 1, size = axissize),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          legend.position = 'none',
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA))  +
    scale_x_continuous(limits = c(-6.1, 0),
                       breaks = pticks,
                       expand = c(0,0),
                       labels = fdrlabels) +
    scale_y_continuous(limits = c(0, 1), 
                       breaks = c(0.25, 0.5, 0.75,1), 
                       expand = c(0,0), 
                       labels = c(0.25, 0.5, 0.75,1)*100,
                       sec.axis = sec_axis(~ . * ngene, name = "True Positive DEGs")) + 
    # coord_cartesian( ylim=c(0, 1.1)) +
    geom_line(aes(x = transformp(q), y = power, col = methods, linetype = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    scale_linetype_manual(values=lty_palette) +
    geom_point(aes(x = transformp(q), y = power, shape = methods, col = methods), stroke = 1,size = pt_size) +
    scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)") + ylab('Power (%)')
  p2 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  
  p = ggplot(dat, 
             aes(x = transformp(q), y = power2, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text( colour = 1, size = axissize),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = 'none',
          # axis.title.x = element_blank(),
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA))   +
    scale_x_continuous(limits = c(-6.1, 0),
                       breaks = pticks,
                       expand = c(0,0),
                       labels = fdrlabels) +
    scale_y_continuous(limits = c(0, 1), 
                       breaks = c(0.25, 0.5, 0.75,1), 
                       expand = c(0,0), 
                       labels = c(0.25, 0.5, 0.75,1)*100,
                       sec.axis = sec_axis(~ . * ngene, name = "True Positive DEGs")) + 
    # coord_cartesian( ylim=c(0, 1.1)) +
    geom_line(aes(x = transformp(q), y = power2, col = methods, linetype = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    scale_linetype_manual(values=lty_palette)+
    geom_point(aes(x = transformp(q), y = power2, shape = methods, col = methods), stroke = 1,size = pt_size) +
    scale_shape_manual(values = pt_shape) + xlab("Actual FDR (%)") + ylab('Power (%)')
  p3 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  
  p = p + theme(legend.position = 'right',
                legend.text = element_text(size = textsize),
                legend.key.size = unit(0.4,'in'))
  legend = get_legend(p)
  p = p + theme(legend.position = 'none')
  
  p_ls = list(p1, p2, p3)
  return(list(p_ls, legend))
}

dat2plot_n <- function(dat, thre, ngene){
  dat = dat[dat$q == thre,]

  p = ggplot(dat, 
             aes(x = transformn(n), y = transformp(fdp), group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = 'none',
          # axis.title.x = element_blank(),
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA)) +
    scale_x_continuous(breaks = nticks,
                       expand = c(0,0),
                       labels = nlabels)  +
    scale_y_continuous(limits = c(-6.1, 0),
                       breaks = pticks,
                       expand = c(0,0),
                       labels = fdrlabels) +
    # coord_cartesian( ylim=c(0, 1.1)) +
    geom_hline(yintercept = transformp(thre), lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(transformn(n), y = transformp(fdp), col = methods, linetype = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    scale_linetype_manual(values=lty_palette) + 
    geom_point(aes(transformn(n), y = transformp(fdp), shape = methods, col = methods), stroke = 1,size = pt_size) +
    scale_shape_manual(values = pt_shape) + xlab("Sample Size") + ylab('Actual FDR (%)') + ggtitle(paste0('Target FDR = ', 100*thre, '%'))
  p1 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  
  p = ggplot(dat, 
             aes(x = transformn(n), y = power, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          legend.position = 'none',
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA))  +
    scale_x_continuous(breaks = nticks,
                       expand = c(0,0),
                       labels = nlabels)+
    scale_y_continuous(limits = c(0, 1), 
                       breaks = c(0.25, 0.5, 0.75,1), 
                       expand = c(0,0), 
                       labels = c(0.25, 0.5, 0.75,1)*100,
                       sec.axis = sec_axis(~ . * ngene, name = "True Positive DEGs")) + 
    # coord_cartesian( ylim=c(0, 1.1)) +
    geom_line(aes(x = transformn(n), y = power, col = methods, linetype = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    scale_linetype_manual(values=lty_palette) +
    geom_point(aes(x = transformn(n), y = power, shape = methods, col = methods), stroke = 1,size = pt_size) +
    scale_shape_manual(values = pt_shape) + xlab("Sample Size") + ylab(paste0('Power (%)')) + ggtitle(paste0('Target FDR = ', 100*thre, '%'))
  p2 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  
  p = ggplot(dat, 
             aes(x = transformn(n), y = power2, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = 'none',
          # axis.title.x = element_blank(),
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA))   +
    scale_x_continuous(breaks = nticks,
                       expand = c(0,0),
                       labels = nlabels) +
    scale_y_continuous(limits = c(0, 1), 
                       breaks = c(0.25, 0.5, 0.75,1), 
                       expand = c(0,0), 
                       labels = c(0.25, 0.5, 0.75,1)*100,
                       sec.axis = sec_axis(~ . * ngene, name = "True Positive DEGs")) + 
    # coord_cartesian( ylim=c(0, 1.1)) +
    geom_line(aes(x = transformn(n), y = power2, col = methods, linetype = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    scale_linetype_manual(values=lty_palette) +
    geom_point(aes(x = transformn(n), y = power2, shape = methods, col = methods), stroke = 1,size = pt_size) +
    scale_shape_manual(values = pt_shape) + xlab("Sample Size") + ylab(paste0('Power (%)')) + ggtitle(paste0('Actual FDR = ', 100*thre, '%'))
  p3 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  
  p_ls = list(p1, p2, p3)
  return(p_ls)
}

dat2plot_new <- function(dat, thre, ngene){
  dat = dat[dat$q == thre,]
  
  p = ggplot(dat, 
             aes(x = transformn(n), y = transformp(fdp), group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text( colour = 1, size = axissize),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = 'none',
          plot.title=element_text(hjust=0.5),
          # axis.title.x = element_blank(),
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA)) +
    scale_x_continuous(breaks = nticks,
                       expand = c(0,0),
                       labels = nlabels)  +
    scale_y_continuous(limits = c(-6.1, 0),
                       breaks = pticks,
                       expand = c(0,0),
                       labels = fdrlabels) +
    # coord_cartesian( ylim=c(0, 1.1)) +
    geom_hline(yintercept = transformp(thre), lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(transformn(n), y = transformp(fdp), col = methods, linetype = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    scale_linetype_manual(values=lty_palette) + xlab("Sample Size") + ylab('Actual FDR (%)') + ggtitle(paste0('Target FDR = ', 100*thre, '%'))
  p1 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  
  p = ggplot(dat, 
             aes(x = transformn(n), y = power, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text( colour = 1, size = axissize),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          plot.title=element_text(hjust=0.5),
          # axis.title.x = element_blank(),
          legend.position = 'none',
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA))  +
    scale_x_continuous(breaks = nticks,
                       expand = c(0,0),
                       labels = nlabels)+
    scale_y_continuous(limits = c(0, 1), 
                       breaks = c(0.25, 0.5, 0.75,1), 
                       expand = c(0,0), 
                       labels = c(0.25, 0.5, 0.75,1)*100,
                       sec.axis = sec_axis(~ . * ngene, name = "True Positive DEGs")) + 
    # coord_cartesian( ylim=c(0, 1.1)) +
    geom_line(aes(x = transformn(n), y = power,col = methods, linetype = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    scale_linetype_manual(values=lty_palette) + xlab("Sample Size") + ylab(paste0('Power (%)')) + ggtitle(paste0('Target FDR = ', 100*thre, '%'))
  p2 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  
  p = ggplot(dat, 
             aes(x = transformn(n), y = power2, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.title=element_text(hjust=0.5),
          axis.text = element_text( colour = 1, size = axissize),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = 'none',
          # axis.title.x = element_blank(),
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA))   +
    scale_x_continuous(breaks = nticks,
                       expand = c(0,0),
                       labels = nlabels) +
    scale_y_continuous(limits = c(0, 1), 
                       breaks = c(0.25, 0.5, 0.75,1), 
                       expand = c(0,0), 
                       labels = c(0.25, 0.5, 0.75,1)*100,
                       sec.axis = sec_axis(~ . * ngene, name = "True Positive DEGs")) + 
    # coord_cartesian( ylim=c(0, 1.1)) +
    geom_line(aes(x = transformn(n), y = power2, col = methods, linetype = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    scale_linetype_manual(values=lty_palette) + xlab("Sample Size") + ylab(paste0('Power (%)')) + ggtitle(paste0('Actual FDR = ', 100*thre, '%'))
  p3 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  
  p_ls = list(p1, p2, p3)
  return(p_ls)
}

celltype.pairs <- c("AdiposeSubcutaneous.AdiposeVisceralOmentum", "BrainAmygdala.BrainSpinalcordCervicalc",
                    "CellsEVBtransformedLymphocytes.MinorSalivaryGland", "HeartAtrialAppendage.HeartLeftVentricle", 
                    "WholeBlood.MuscleSkeletal", "Brain.Prostate")
tcga.pairs <- c('BRCA', 'KIRC', 'LIHC', 'LUAD', 'PRAD', 'THCA')

# ls_0<- readRDS("fdppow_plot_0518.rds")
# ls_sub_gtex <- readRDS("fdppow_subsample_gtex.rds")
# ls_sub_tcga <- readRDS("fdppow_subsample_tcga.rds")
setwd("~/Study/research/pvalue/Yumei/070121_nDEG/")
ls_0 <- readRDS("../fdppow_plot_noiseq.rds")
ls_sub_gtex <- readRDS("fdppow_subsample_gtex.rds")
ls_sub_tcga <- readRDS("fdppow_subsample_tcga.rds")

num.trueDE <- readRDS("../num_trueDE_noiseq.rds")

for (i in 1:6){
  dat = ls_0[[1]][[i]]
  ngene = num.trueDE$GTEx[[i]]
  p_ls1 = dat2plot(dat, ngene  = ngene)
  dat = ls_sub_gtex[[i]]
  dat[is.na(dat)] = 0
 
  p_ls2 = dat2plot_new(dat, thre = 0.1, ngene  = ngene)
  p_ls4 = dat2plot_new(dat, thre = 0.01, ngene  = ngene)
  pdf(file =  paste0(celltype.pairs[i], '_fdppow.pdf'), height = unit(7.2, "in"), width = unit(8, "in"))

  p = arrangeGrob(
    ncol = 1,nrow =5, widths = c(unit(7, "in")), heights = c(0.2,2.4,0.2,2.4,2.4),
    textGrob(label = 'A', hjust = 20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call( 'cbind', p_ls1[[1]]), 
    textGrob(label = 'B', hjust =20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call('cbind', p_ls2),
    do.call('cbind', p_ls4)
  )
  grid.arrange(p,p_ls1[[2]],  ncol = 2,nrow = 1, widths = c(unit(7, "in"), unit(1, "in")))
  dev.off()
}

for (i in 1:6){
  dat = ls_0[[1]][[i]]
  ngene = num.trueDE$TCGA[[i]]
  p_ls1 = dat2plot(dat, ngene  = ngene)
  dat = ls_sub_tcga[[i]]
  dat[is.na(dat)] = 0
  
  p_ls2 = dat2plot_new(dat, thre = 0.1, ngene  = ngene)
  p_ls3 = dat2plot_new(dat, thre = 0.05, ngene  = ngene)
  p_ls4 = dat2plot_new(dat, thre = 0.01, ngene  = ngene)
  p_ls5 = dat2plot_new(dat, thre = 0.001, ngene  = ngene)
  pdf(file =  paste0(tcga.pairs[i], '_fdppow.pdf'), height = unit(13, "in"), width = unit(8, "in"))
  p = arrangeGrob(
    ncol = 1,nrow =7, widths = c(unit(7, "in")), heights = c(0.2,2.4,0.2,2.4,2.4,2.4,2.4),
    textGrob(label = 'A', hjust = 20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call( 'cbind', p_ls1[[1]]), 
    textGrob(label = 'B', hjust =20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call('cbind', p_ls2),
    do.call('cbind', p_ls3),
    do.call('cbind', p_ls4),
    do.call('cbind', p_ls5)
  )
  grid.arrange(p,p_ls1[[2]],  ncol = 2,nrow = 1, widths = c(unit(7, "in"), unit(1, "in")))
  dev.off()
}



ls_sub_gtex <- readRDS("fdppow_subsample_gtex_0001.rds")
ls_sub_tcga <- readRDS("fdppow_subsample_tcga_0001.rds")

num.trueDE <- readRDS("num_trueDE_thre0001.rds")

for (i in 1:6){
  dat = ls_0[[1]][[i]]
  ngene = num.trueDE$GTEx[[i]]
  p_ls1 = dat2plot(dat, ngene  = ngene)
  dat = ls_sub_gtex[[i]]
  dat[is.na(dat)] = 0
  
  p_ls2 = dat2plot_new(dat, thre = 0.1, ngene  = ngene)
  p_ls3 = dat2plot_new(dat, thre = 0.05, ngene  = ngene)
  p_ls4 = dat2plot_new(dat, thre = 0.01, ngene  = ngene)
  p_ls5 = dat2plot_new(dat, thre = 0.001, ngene  = ngene)
  pdf(file =  paste0('thre0001/',celltype.pairs[i], '_fdppow.pdf'), height = unit(13, "in"), width = unit(8, "in"))
  
  p = arrangeGrob(
    ncol = 1,nrow =7, widths = c(unit(7, "in")), heights = c(0.2,2.4,0.2,2.4,2.4,2.4,2.4),
    textGrob(label = 'A', hjust = 20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call( 'cbind', p_ls1[[1]]), 
    textGrob(label = 'B', hjust =20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call('cbind', p_ls2),
    do.call('cbind', p_ls3),
    do.call('cbind', p_ls4),
    do.call('cbind', p_ls5)
  )
  grid.arrange(p,p_ls1[[2]],  ncol = 2,nrow = 1, widths = c(unit(7, "in"), unit(1, "in")))
  dev.off()
}

for (i in 1:6){
  dat = ls_0[[1]][[i]]
  ngene = num.trueDE$TCGA[[i]]
  p_ls1 = dat2plot(dat, ngene  = ngene)
  dat = ls_sub_tcga[[i]]
  dat[is.na(dat)] = 0
  
  p_ls2 = dat2plot_new(dat, thre = 0.1, ngene  = ngene)
  p_ls3 = dat2plot_new(dat, thre = 0.05, ngene  = ngene)
  p_ls4 = dat2plot_new(dat, thre = 0.01, ngene  = ngene)
  p_ls5 = dat2plot_new(dat, thre = 0.001, ngene  = ngene)
  pdf(file =  paste0('thre0001/',tcga.pairs[i], '_fdppow.pdf'), height = unit(13, "in"), width = unit(8, "in"))
  p = arrangeGrob(
    ncol = 1,nrow =7, widths = c(unit(7, "in")), heights = c(0.2,2.4,0.2,2.4,2.4,2.4,2.4),
    textGrob(label = 'A', hjust = 20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call( 'cbind', p_ls1[[1]]), 
    textGrob(label = 'B', hjust =20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call('cbind', p_ls2),
    do.call('cbind', p_ls3),
    do.call('cbind', p_ls4),
    do.call('cbind', p_ls5)
  )
  grid.arrange(p,p_ls1[[2]],  ncol = 2,nrow = 1, widths = c(unit(7, "in"), unit(1, "in")))
  dev.off()
}

ls_sub_gtex <- readRDS("fdppow_subsample_gtex_00000001.rds")
ls_sub_tcga <- readRDS("fdppow_subsample_tcga_00000001.rds")

num.trueDE <- readRDS("num_trueDE_thre00000001.rds")

for (i in 1:6){
  dat = ls_0[[1]][[i]]
  ngene = num.trueDE$GTEx[[i]]
  p_ls1 = dat2plot(dat, ngene  = ngene)
  dat = ls_sub_gtex[[i]]
  dat[is.na(dat)] = 0
  
  p_ls2 = dat2plot_new(dat, thre = 0.1, ngene  = ngene)
  p_ls3 = dat2plot_new(dat, thre = 0.05, ngene  = ngene)
  p_ls4 = dat2plot_new(dat, thre = 0.01, ngene  = ngene)
  p_ls5 = dat2plot_new(dat, thre = 0.001, ngene  = ngene)
  pdf(file =  paste0('thre00000001/', celltype.pairs[i], '_fdppow.pdf'), height = unit(13, "in"), width = unit(8, "in"))
  
  p = arrangeGrob(
    ncol = 1,nrow =7, widths = c(unit(7, "in")), heights = c(0.2,2.4,0.2,2.4,2.4,2.4,2.4),
    textGrob(label = 'A', hjust = 20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call( 'cbind', p_ls1[[1]]), 
    textGrob(label = 'B', hjust =20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call('cbind', p_ls2),
    do.call('cbind', p_ls3),
    do.call('cbind', p_ls4),
    do.call('cbind', p_ls5)
  )
  grid.arrange(p,p_ls1[[2]],  ncol = 2,nrow = 1, widths = c(unit(7, "in"), unit(1, "in")))
  dev.off()
}

for (i in 1:6){
  dat = ls_0[[1]][[i]]
  ngene = num.trueDE$TCGA[[i]]
  p_ls1 = dat2plot(dat, ngene  = ngene)
  dat = ls_sub_tcga[[i]]
  dat[is.na(dat)] = 0
  
  p_ls2 = dat2plot_new(dat, thre = 0.1, ngene  = ngene)
  p_ls3 = dat2plot_new(dat, thre = 0.05, ngene  = ngene)
  p_ls4 = dat2plot_new(dat, thre = 0.01, ngene  = ngene)
  p_ls5 = dat2plot_new(dat, thre = 0.001, ngene  = ngene)
  pdf(file =  paste0('thre00000001/', tcga.pairs[i], '_fdppow.pdf'), height = unit(13, "in"), width = unit(8, "in"))
  p = arrangeGrob(
    ncol = 1,nrow =7, widths = c(unit(7, "in")), heights = c(0.2,2.4,0.2,2.4,2.4,2.4,2.4),
    textGrob(label = 'A', hjust = 20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call( 'cbind', p_ls1[[1]]), 
    textGrob(label = 'B', hjust =20, gp = gpar(fontsize= titletextsize, fontface =2)),
    do.call('cbind', p_ls2),
    do.call('cbind', p_ls3),
    do.call('cbind', p_ls4),
    do.call('cbind', p_ls5)
  )
  grid.arrange(p,p_ls1[[2]],  ncol = 2,nrow = 1, widths = c(unit(7, "in"), unit(1, "in")))
  dev.off()
}
### number of trueDE genes ##########
num.trueDE <- list()
num.trueDE[[1]] <- lapply(1:6, function(j){
  dis0 <- readRDS(paste0(celltype.pairs[j], '_dis.rds'))
  trueDE = intersect(intersect(intersect(intersect(dis0[[1]], dis0[[2]]), dis0[[3]]), dis0[[4]]),dis0[[5]])
  return(length(trueDE)%/%2)
})

num.trueDE[[2]] <- lapply(1:6, function(j){
  dis0 <- readRDS(paste0(tcga.pairs[j], '_dis.rds'))
  trueDE = intersect(intersect(intersect(intersect(dis0[[1]], dis0[[2]]), dis0[[3]]), dis0[[4]]),dis0[[5]])
  return(length(trueDE)%/%2)
})
names(num.trueDE) <- c("GTEx", "TCGA")
saveRDS(num.trueDE, "num_trueDE.rds")


num.trueDE <- list()
num.trueDE[[1]] <- lapply(1:6, function(j){
  dis0 <- readRDS(paste0('thre0001/', celltype.pairs[j], '_dis.rds'))
  trueDE = intersect(intersect(intersect(intersect(dis0[[1]], dis0[[2]]), dis0[[3]]), dis0[[4]]),dis0[[5]])
  return(length(trueDE)%/%2)
})

num.trueDE[[2]] <- lapply(1:6, function(j){
  dis0 <- readRDS(paste0('thre0001/thre0001TCGA-', tcga.pairs[j], '_dis.rds'))
  trueDE = intersect(intersect(intersect(intersect(dis0[[1]], dis0[[2]]), dis0[[3]]), dis0[[4]]),dis0[[5]])
  return(length(trueDE)%/%2)
})
names(num.trueDE) <- c("GTEx", "TCGA")
saveRDS(num.trueDE, "num_trueDE_thre0001.rds")

num.trueDE <- list()
num.trueDE[[1]] <- lapply(1:3, function(j){
  dis0 <- readRDS(paste0('thre00000001/', celltype.pairs[j], '_dis.rds'))
  trueDE = intersect(intersect(intersect(intersect(dis0[[1]], dis0[[2]]), dis0[[3]]), dis0[[4]]),dis0[[5]])
  return(length(trueDE)%/%2)
})

num.trueDE[[2]] <- lapply(1:6, function(j){
  dis0 <- readRDS(paste0('thre00000001/TCGA-', tcga.pairs[j], '_dis.rds'))
  trueDE = intersect(intersect(intersect(intersect(dis0[[1]], dis0[[2]]), dis0[[3]]), dis0[[4]]),dis0[[5]])
  return(length(trueDE)%/%2)
})
names(num.trueDE) <- c("GTEx", "TCGA")
saveRDS(num.trueDE, "num_trueDE_thre00000001.rds")
