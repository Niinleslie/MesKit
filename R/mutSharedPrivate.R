#' Use R code to find the intersect mutations and their types in several samples of one patient
#'
#' @param maf Maf object return from read.Maf()
#' @param show.num a logic parameter to determine whether to show the number of each mutations in the stack plot
#' @param savePlot if save plot of result
#'
#' @return mutSharedPrivate
#' @examples
#' maf.File <- system.file("extdata/multi_lesion/maf", "311252.maf", package = "Meskit")
#' sampleInfo.File <- system.file("extdata/multi_lesion", "sample_info.txt", package = "Meskit")
#' pyCloneCluster <- system.file("extdata/multi_lesion/ccf", "311252.cluster.tsv", package = "Meskit")
#' pyCloneLoci <- system.file("extdata/multi_lesion/ccf", "311252.loci.tsv", package = "Meskit")
#' maf <- readMaf(patientID = "311252", mafFile = maf.File, sampleInfo = sampleInfo.File, refBuild = "hg19")
#' mutSharedPrivate(maf)

mutSharedPrivate <- function(maf, show.num = FALSE, savePlot = FALSE){
  maf.dat <- dplyr::select(maf@data,c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,  Variant_Type,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,Ref_allele_depth,Alt_allele_depth,VAF,CDS_Change,Protein_Change,Tumor_Sample_Barcode))
  patientID <- maf@patientID
  analysisMutSharedPrivate(maf.dat, patientID = patientID, show.num = show.num, savePlot = savePlot)
}
analysisMutSharedPrivate <- function(df = NULL, patientID = patientID, show.num = FALSE, savePlot){
  primitiveLength <- length(df$Hugo_Symbol)
  #separate by quote
  df$Hugo_Symbol <- as.character(df$Hugo_Symbol)
  newGene <- c()
  for (n in 1:length(df$Hugo_Symbol)){
    newGene <- c(newGene, na.omit(strsplit(df$Hugo_Symbol[n], split = ",")[[1]][-1]))
  }
  n=1
  while (n>0 & n<primitiveLength + 1){
    geneNumber <- as.numeric(length(strsplit(as.character(df$Hugo_Symbol[n]),split = ",")[[1]]))
    if (geneNumber > 1){
      df[n,1] <- strsplit(df$Hugo_Symbol[n],split = ",")[[1]][1]
      while (geneNumber>1){
        l=length(df$Hugo_Symbol)
        df$Hugo_Symbol <- as.character(df$Hugo_Symbol)
        newRow <- df[n, ]
        df <- rbind(df,newRow) 
        df[l+1, 2:9] <- df[n, 2:9]
        df[l+1, 1] <- as.character(strsplit(df$Hugo_Symbol[n],split = ",")[[1]][geneNumber])
        geneNumber = geneNumber - 1
      }
    }
    n=n+1
  }
  newGeneNumber <- length(df$Hugo_Symbol) - primitiveLength
  for (n in 1 : newGeneNumber){
    df[primitiveLength + n, 1] <- newGene[n]
  }
  #discrete the different data
  for (n in 1:length(df$Hugo_Symbol)){
    if (length(strsplit(as.character(df[n, 5]),split = "\t")[[1]])>2){
      df[n, ] <- NA
    }
  }
  df <- df[, c(1:5,7,9,15)]
  df <- na.omit(df)
  #organize dataframe
  df$Tumor_Sample_Barcode <- factor(df$Tumor_Sample_Barcode,levels = unique(df$Tumor_Sample_Barcode))
  df$Variant_Classification <- factor(df$Variant_Classification,levels = unique(df$Variant_Classification))
  df <- df[order(df$Tumor_Sample_Barcode), ]
  #get the intersection genes
  combination <- function(vector){
    n <- length(vector)
    num <- 0
    cycle <- 1
    for (i in 1:n){
      num = num + choose(n,i)
    }
    result <- list()
    for (j in 1:n){
      oneresult <- list(combn(vector, j))
      result = c(result,oneresult)
      cycle = cycle + 1
      if (cycle == num)
        break
    }
    return(result)
  }
  df$Reference_Allele <- as.character(df$Reference_Allele)
  df$Tumor_Seq_Allele2 <- as.character(df$Tumor_Seq_Allele2)
  sample.combination <- combination(c(levels(df$Tumor_Sample_Barcode)))
  neededData <- select(tidyr::unite(df, "neededData", Hugo_Symbol, Chromosome, Start_Position, End_Position ,Reference_Allele, Tumor_Seq_Allele2, sep = "_"), neededData)
  df <- cbind(df,neededData)
  #get intersection
  allIntersect <- list()
  getIntersect <- function(sample,information){
    combination_number <- length(sample)/length(sample[, 1])
    sample_number <- length(sample[, 1])
    n <- 1
    while (n > 0 & n < combination_number + 1){
      intersect.n <- list()
      a=1
      while (a < sample_number + 1){
        intersect.n <- c(intersect.n, 
                         list(information$neededData[which(information$Tumor_Sample_Barcode == sample[, n][a])]))
        a <- a+1
      }
      result.n <- Reduce(intersect,intersect.n)
      allIntersect <<- c(allIntersect,list(length(result.n)))
      allIntersect <<- c(allIntersect,list(result.n))
      n <- n+1
    }
  }
  process <- lapply(sample.combination, getIntersect, df)
  #get the dataframe only includes genes
  o <- 1
  all.intersect.only <- list()
  while(o > 0 & o < length(allIntersect) + 1){
    if (o%%2 == 0){
      all.intersect.only <- c(all.intersect.only, list(allIntersect[o][[1]]))
    }
    o <- o+1
  }
  #create the final dataframe
  all.final.frame <- data.frame(0,0,0)
  colnames(all.final.frame) <- c("Sample","Type","Number")
  all.final.frame <- all.final.frame[-1, ]
  #get the first dataframe with one sample
  sample.number <- length(levels(df$Tumor_Sample_Barcode))
  sample.represent <- c(1 : length(levels(df$Tumor_Sample_Barcode)))
  final.frame0 <- data.frame(0,0,0)
  colnames(final.frame0) <- c("Sample","Type","Number")
  final.frame0 <- final.frame0[-1, ]
  n=1
  while(n > 0 & n<sample.number + 1){
    sub.n <- setdiff(sample.represent,sample.represent[n])
    seperate.n <- setdiff(unlist(all.intersect.only[n]),unlist(all.intersect.only[sub.n]))
    final.seperate.n <- c()
    a = 1
    while (a>0 & a<length(seperate.n)+1){
      final.seperate.n <- c(final.seperate.n,as.character(unique(df$Variant_Classification[which (df$neededData==seperate.n[a])])))
      a = a + 1
    }
    try.seperate.n <- data.frame(0,0,0)
    colnames(try.seperate.n) <- c("Sample","Type","Number")
    try.seperate.n[1:length(names(table(final.seperate.n))), 2] <- names(table(final.seperate.n))
    try.seperate.n[1:length(names(table(final.seperate.n))), 3] <- as.character(table(final.seperate.n))
    try.seperate.n$Sample <- sample.combination[[1]][n]
    final.frame0 <- rbind(final.frame0, try.seperate.n)
    n = n + 1
  }
  final.frame0$Type <- factor(final.frame0$Type,levels = unique(final.frame0$Type))
  final.frame0$Sample <- factor(final.frame0$Sample,levels = unique(final.frame0$Sample))
  final.frame0$Number <- as.integer(final.frame0$Number)
  #create orders
  sample.order0 <- data.frame(0,0)
  colnames(sample.order0) <- c("Sample","Number")
  sample.order0 <- sample.order0[-1,]
  sample.order0[1:sample.number,1] <- levels(df$Tumor_Sample_Barcode)[1:sample.number]
  for (n in 1:sample.number){
    sample.order0[n,2] <- sum(as.integer(final.frame0$Number[which(final.frame0$Sample==levels(final.frame0$Sample)[n])]))
  }
  sample.order0 <- sample.order0[order(sample.order0$Number,decreasing =TRUE),]
  sample.order0$Sample <- factor(sample.order0$Sample,levels=sample.order0$Sample)
  order.sample0 <- as.character(sample.order0$Sample)
  #organize the final dataframe with one sample
  final.frame0$Type <- factor(final.frame0$Type,levels = unique(final.frame0$Type))
  final.frame0$Sample <- factor(final.frame0$Sample,levels = order.sample0)
  final.frame0$Number <- as.integer(final.frame0$Number)
  all.final.frame <- rbind(all.final.frame,final.frame0)
  #get other dataframes
  CNM <- function(a, b){
    return(factorial(a)/(factorial(b)*factorial(a-b)))
  }
  sumCNM <- function(a, b){
    sum.final = 0
    while(b>0){
      sum.final <- sum.final+CNM(a,b)
      b=b-1
    }
    return(sum.final)
  }
  q <- 2
  while(q > 1 & q < sample.number + 1){
    p <- sumCNM(sample.number,(q - 1)) + 1
    final.frame1 <- data.frame(0,0,0)
    colnames(final.frame1) <- c("Sample","Type","Number")
    final.frame1 <- final.frame1[-1, ]
    while(p>sumCNM(sample.number, q - 1) & p<sumCNM(sample.number,q) + 1){
      n=1
      types.n <- c()
      while (n > 0 & n<(length(all.intersect.only[p][[1]]) + 1)){
        types.n <- c(types.n,as.character(unique(df$Variant_Classification[which (df$neededData==all.intersect.only[p][[1]][n])])))
        n=n+1
      }
      if(length(types.n) != 0){
        try.n <- data.frame(0,0,0)
        colnames(try.n) <- c("Sample","Type","Number")
        try.n[1:length(names(table(types.n))), 2] <- names(table(types.n))
        try.n[1:length(names(table(types.n))), 3] <- as.character(table(types.n))
        try.n$Sample[1:length(names(table(types.n)))] <- paste(sample.combination[[q]][,(p-sumCNM(sample.number,q-1))], collapse = "∩")
        try.n$Type <- factor(try.n$Type,levels = levels(df$Variant_Classification))
        final.frame1 <- rbind(final.frame1,try.n)
      }
      p=p+1
    }
    r <- sumCNM(sample.number,q - 1) + 1
    sample.order.q <- data.frame(0,0)
    colnames(sample.order.q) <- c("Sample","Number")
    sample.order.q <- sample.order.q[-1, ]
    Number <- c()
    Sample <- c()
    while (r > sumCNM(sample.number,q - 1) & r<sumCNM(sample.number,q) + 1){
      Number <- c(Number, as.integer(allIntersect[2*r-1][[1]]))
      Sample <- c(Sample, paste(sample.combination[[q]][,(r-sumCNM(sample.number,q - 1))],collapse = "∩"))
      r=r+1
    }
    sample.order.q[1:length(Sample),1] <- Sample
    sample.order.q[1:length(Number),2] <- Number
    sample.order.q <- sample.order.q[order(sample.order.q$Number,decreasing = TRUE),]
    sample.order.q <- sample.order.q[which(sample.order.q$Number!=0),]
    sample.order.q$Sample <- factor(sample.order.q$Sample,levels=sample.order.q$Sample)
    order_sample_w <- as.character(sample.order.q$Sample)
    final.frame1$Type <- factor(final.frame1$Type,levels = unique(final.frame1$Type))
    final.frame1$Sample <- factor(final.frame1$Sample,levels = order_sample_w)
    final.frame1$Number <- as.integer(final.frame1$Number)
    all.final.frame <- rbind(final.frame1,all.final.frame)
    q = q + 1
  }
  #process data which is going to be painted with points and lines
  point.line.frame <- data.frame(0)
  colnames(point.line.frame) <- c("combinations")
  point.line.frame$sample <- c(0)
  i=1
  rows.num = 1
  while (i > 0 & i < length(levels(all.final.frame$Sample)) + 1){
    sample.number.i <- length(strsplit(as.character(levels(all.final.frame$Sample)[i]),split = "∩")[[1]])
    point.line.frame[rows.num:(rows.num+sample.number.i - 1),1] <- as.character(levels(all.final.frame$Sample)[i])
    point.line.frame[rows.num:(rows.num+sample.number.i - 1),2] <- strsplit(as.character(levels(all.final.frame$Sample))[i],split = "∩")[[1]]
    rows.num <- rows.num+sample.number.i
    i=i+1
  }
  all.final.frame$Type  <-  factor(all.final.frame$Type, levels = unique(all.final.frame$Type))
  #draw stack plots
  key.point <- ggplot(all.final.frame)+
    aes(x = Sample,y = Number,fill = Type)+
    geom_bar(stat = "identity",position = "stack",width = 0.7)+
    labs(x = "",width = 1.0)+
    labs(y = "Mutation number")
  if (show.num == "TRUE") {
    key.point <- key.point +
      geom_text(aes(label = Number), size = 3, colour = 'black',
                hjust = .5, position = position_stack(vjust=0.5))
  }
  key.point <- key.point +
    theme(axis.ticks.x  = element_blank(),
          panel.border =element_blank(),
          axis.text.x =element_blank(),
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=14),
          plot.margin = unit(c(0.08,0.2,0,0.1),"inches"),
          legend.spacing  = unit(c(0.09,0,0,0),"inches"),
          legend.key.width  = unit(0.2,"inches"))+
    scale_fill_manual(values =c( "#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2",
                                 "#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2","#91D1C2B2",
                                 "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"))+
    scale_y_continuous(expand = c(0,0))
  # draw point-line plot
  point.line.frame$combinations <- factor(point.line.frame$combinations,levels = levels(all.final.frame$Sample))
  point.line.plot <- ggplot(point.line.frame)+
    aes(x=combinations,y=sample)
  if(length(levels(all.final.frame$Sample))<21){
    point.line.plot <- point.line.plot + geom_point(size=3.5)
  }
  else{
    point.line.plot <- point.line.plot + geom_point(size=2.5)
  }
  point.line.plot <- point.line.plot+
    geom_path(mapping = aes(group=combinations),inherit.aes = TRUE)+
    labs(x="",width=1.0)+
    labs(y="")+
    theme(panel.grid = element_blank(),
          panel.border =element_blank(),
          axis.text.x =element_blank(),
          axis.title.y = element_text(size=14),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = unit(c(0.15,0,0,0.76),"inches"))+
    scale_y_discrete(position = "right")
  # Keep the width of two plots the same
  point.line.plot  <-  ggplot_gtable(ggplot_build(point.line.plot))
  bar.plot  <-  ggplot_gtable(ggplot_build(key.point)) 
  point.line.plot$widths  <-  bar.plot$widths
  # put pictures together
  if(length(levels(all.final.frame$Sample)) < 21){
    gg <- ggdraw()+draw_plot(bar.plot,0,0.3,1,0.7)+draw_plot(point.line.plot,0,0,1,0.35)
    if(savePlot){
      ggsave(paste("mutSharedPrivate.",patientID,".pdf",sep = ""),width = 10,height = 10,plot = gg)
    }
  } 
  else {
    if (length(levels(all.final.frame$Sample)) > 60) {
      gg <- ggdraw()+draw_plot(bar.plot,0,0.3,1,0.7)+draw_plot(point.line.plot,0,0,1,0.35)
      if(savePlot){
        ggsave(paste("mutSharedPrivate.",patientID,".pdf",sep = ""),width = 13,height = 10,plot = gg)
      }
    } 
    else {
      gg <- ggdraw()+draw_plot(bar.plot,0,0.3,1,0.7)+draw_plot(point.line.plot,0,0,1,0.35)
      if(savePlot){
        ggsave(paste("mutSharedPrivate.",patientID,".pdf",sep = ""),width = 10,height = 10,plot = gg)
      }
    }
  }  
  return(gg)
}
        