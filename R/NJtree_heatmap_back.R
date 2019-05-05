#### Usage: Rscript Draw_Mut_sort_heatmap.R <Mutation_binaray_matrix.input> [Outdir]
library(reshape2)
library(ape)
library(ggplot2)
require(extrafont)
font_import(pattern="[A/a]rial", prompt=FALSE)

##read mutation matrix, trans to 0-1 binary matrix（0 represents mutation absent）
##add a new column(0) representing normal sample
snp_ccf_sort <- function(M){
  index_row_col <- snp_binary_sort(M, is.ccf = TRUE)

  names <- c("snp", "sample", "CCF")
  colnames(M) <- names
  snp_samples <- dcast(M, snp~sample, value.var = "CCF")[,-1]
  snp_samples[snp_samples == 0 ] <- 2
  snp_samples[is.na(snp_samples)] <- 0
  snp_samples$normal <- 0

  snp_samples<- apply(snp_samples,2,as.numeric)
  snp_sort <- snp_samples[index_row_col[[1]],index_row_col[[2]]]
  snp_melt <- melt(snp_sort)
  colnames(snp_melt) <- names
  return(snp_melt)
}


##sort the matrix; return row order and coloum order
snp_binary_sort <- function(M, is.ccf = FALSE) {
  names <- c("snp", "sample", "Mutation")
  colnames(M) <- names
  snp_samples <- dcast(M[,-3], snp~sample)[,-1]
  snp_samples[!is.na(snp_samples)] <- 1
  snp_samples[is.na(snp_samples)] <- 0
  snp_samples$normal <- 0
  snp_binary<- apply(snp_samples,2,as.numeric)

  snp_binary <- t(snp_binary)
  sampleOrder <- sort(rowSums(snp_binary), decreasing = TRUE, index.return = TRUE)$ix
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i)
      }
    }
    return(score)
  }

  scores <- apply(snp_binary[sampleOrder, ], 2, scoreCol)
  geneOrder <- sort(scores, decreasing = TRUE, index.return = TRUE)$ix
  index_row_col <- list(geneOrder, sampleOrder)

  snp_binary_sort <- t(snp_binary[sampleOrder,geneOrder])
  snp_binary_melt <- melt(snp_binary_sort)
  colnames(snp_binary_melt) <- names

  if(is.ccf){
    return(index_row_col)
  }
  else{
    return(snp_binary_melt)
  }
}



heatmap_input <- function(maf, use.indel = FALSE, use.ccf = FALSE, ccf = data.frame(), ccf.mutation.id, ccf.mutation.sep){
  if(use.ccf){
    mut_dat <- mut_ccf_sort(snp, ccf = ccf, use.indel = use.indel, ccf.mutation.id, ccf.mutation.sep)
  }
  else{
    mut_dat <- mut_binary_sort(maf, use.indel = use.indel)
  }
  return(mut_dat)
}


mut.heatmap <- function(snp, is.ccf = FALSE){
  type  <- "Mutation"
  if(is.ccf){type <- "CCF"}
  plot_dat <- heatmap_input(snp, is.ccf = is.ccf)
  p_basic <- ggplot(plot_dat, aes(sample, snp)) +
    labs(x = "", y = "") + theme_bw() +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text.y = element_blank()) +  ## 删去y轴刻度标签
    theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 10, hjust = 1, face = "bold")) +
    theme(axis.ticks = element_blank()) +   ## 删去所有刻度线
    theme(legend.title = element_text(size = 12, face = "bold", color = "black")) +   ##标签名字体
    theme(legend.text = element_text(size = 10, face = "bold", color = "black")) +    ##标签字体
    labs(fill = type) + scale_y_continuous(expand = c(0,0))

  if(is.ccf){
    p <- p_basic + geom_tile(aes(fill = plot_dat$CCF)) + scale_fill_gradient(low = "grey", high = "red", na.value="black", limit=c(0, 1))
    ggsave(paste(patientID,"_mut_CCF.pdf",sep = ""), p, width = 4.5, height = 6.5)
    }
  if(!is.ccf){
    p <- p_basic + geom_tile(aes(fill = as.factor(plot_dat$Mutation))) + scale_fill_manual(values = c("grey","red"))
    ggsave(paste(patientID,"_mut.pdf",sep = ""), p, width = 4.5, height = 6.5)
    }
}

mut.heatmap(snp,is.ccf = T)
setwd("C:/Users/HP-1/Desktop")
patientID <- args[1]
snp<-read.table(paste(patientID,".CCF.txt", sep = ""), header = FALSE, quote = "\t")



