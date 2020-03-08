#' mutSharedPrivate
#' 
#' @description  Use R code to find the intersect mutations and their types in several samples of one patient
#' @param maf Maf object return from read.Maf()
#' @param show.num a logic parameter to determine whether to show the number of each mutations in the stack plot
#'
#' @importFrom tidyr unite
#' @importFrom dplyr select
#' @examples
#' maf.File <- system.file("extdata/maf", "HCC6046.maf", package = "Meskit")
#' sampleInfo.File <- system.file("extdata", "HCC6046.sampleInfo.txt", package = "Meskit")
#' ccf.cluster.File <- system.file("extdata/ccf", "HCC6046.cluster.tsv", package = "Meskit")
#' ccf.loci.File <- system.file("extdata/ccf", "HCC6046.loci.tsv", package = "Meskit")
#' maf <- readMaf(mafFile = maf.File, sampleInfo = sampleInfo.File, 
#'                 ccfClusterTsvFile = ccf.cluster.File, 
#'                 ccfLociTsvFile = ccf.loci.File,
#'                 refBuild = "hg19")
#' mutSharedPrivate(maf)
#' @export mutSharedPrivate

mutSharedPrivate <- function(maf, show.num = FALSE){
    maf.dat <- maf@data
    patient.ids <- unique(maf.dat$Patient_ID)
    msp.list <- maf.dat %>% group_by(Patient_ID) %>% 
                group_map(~doMutSharedPrivate(.x,show.num = show.num), keep = TRUE) %>%
                rlang::set_names(patient.ids)
    return(msp.list)
}

doMutSharedPrivate <- function(df, show.num){
    df$Hugo_Symbol <- as.character(df$Hugo_Symbol)
    while(TRUE){
        numPos <-  unlist(lapply(df$Hugo_Symbol, function(x){
            len <- length(strsplit(x,split = ",")[[1]])
            return(len)
        }))
        dfNew <- df[numPos >= 2, ]
        dfNew$Hugo_Symbol <- as.character(unlist(
            lapply(dfNew$Hugo_Symbol, function(x){
                s <- strsplit(x, ",")[[1]]
                if(length(s) > 2){
                    return(paste(s[-1],collapse = ","))
                }
                else{
                    return(s[2])
                }
            })
        ))
        df$Hugo_Symbol <- as.character(unlist(
            lapply(df$Hugo_Symbol, function(x){
                return(strsplit(x, ",")[[1]][1])
            })
        ))
        df <- rbind(df,dfNew)
        if(length(which(numPos>=2)) == 0){
            break
        }
    }
    df <- dplyr::select(df,c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification, Reference_Allele,Tumor_Seq_Allele2,Tumor_Sample_Barcode))
    df <- na.omit(df)
    #organize dataframe
    df$Tumor_Sample_Barcode <- factor(df$Tumor_Sample_Barcode,levels = unique(df$Tumor_Sample_Barcode))
    df$Variant_Classification <- factor(df$Variant_Classification,levels = unique(df$Variant_Classification))
    df <- df[order(df$Tumor_Sample_Barcode), ]
    df$Reference_Allele <- as.character(df$Reference_Allele)
    df$Tumor_Seq_Allele2 <- as.character(df$Tumor_Seq_Allele2)
    sample.combination <- lapply(1:length(levels(df$Tumor_Sample_Barcode)),function(x){
        combn(levels(df$Tumor_Sample_Barcode),x)
    })
    neededData <- select(tidyr::unite(df, "neededData", Hugo_Symbol, Chromosome, Start_Position, End_Position ,Reference_Allele,Tumor_Seq_Allele2, sep = "_"), neededData)
    df <- cbind(df,neededData)
    finalFrame <- data.frame(Sample = "x",Type = "y",Number = 0)
    finalFrame$Sample <- as.vector(finalFrame$Sample)
    finalFrame$Type <- as.vector(finalFrame$Type)
    finalFrame$Number <- as.vector(finalFrame$Number)
    pointLineFrame <- data.frame(combinations = "x",sample = "y")
    pointLineFrame$combinations <- as.vector(pointLineFrame$combinations)
    pointLineFrame$sample <- as.vector(pointLineFrame$sample)
    i = 1
    while(i <= length(sample.combination)){
        isData <- apply(sample.combination[[i]],2,function(x){
            if(length(x) == 1){
                neededData <- df[which(df$Tumor_Sample_Barcode == x),]$neededData
                neededData <- unlist(lapply(neededData, function(j){
                    sampleGet <- unique(df[which(df$neededData == j),]$Tumor_Sample_Barcode)
                    if(length(sampleGet) > 1){
                        return(NULL)
                    }
                    else{
                        return(j)
                    }
                }))
                isNeededData <- neededData
            }
            else{
                neededData <- lapply(x,function(s){
                    return(df[which(df$Tumor_Sample_Barcode == s),]$neededData)
                })          
                isNeededData <- Reduce(intersect, neededData)
            }
            isDf <- df[which(df$Tumor_Sample_Barcode %in% x & df$neededData %in% isNeededData),]
            if(length(isNeededData) != 0){
                isType <- as.character(unique(isDf$Variant_Classification) ) 
                numList <-  unlist(lapply(isType, function(g){
                    num <- length(which(isDf$Variant_Classification == g)) 
                    return(num)
                }))
                sample <- paste(x,collapse = "∩")
                for(i in 1:length(isType)){
                    if(numList[i] == 0){
                        next
                    }
                    else{
                        type <- isType[i]
                        num <- numList[i]/length(x)
                        row <- c(sample, type, num, x)
                        finalFrame <<- rbind(finalFrame, row)
                    }
                }
                for(i in 1:length(x)){
                    row <- c(sample , x[i])
                    pointLineFrame <<- rbind(pointLineFrame, row) 
                }
            }
            
        })
        i <- i +1
    }
    finalFrame <- finalFrame[-1,]
    pointLineFrame <- pointLineFrame[-1,]
    finalFrame$Sample  <-  factor(finalFrame$Sample, levels = rev(unique(finalFrame$Sample)))
    finalFrame$Type  <-  factor(finalFrame$Type, levels = sort(unique(finalFrame$Type)))
    finalFrame$Number <- as.integer(finalFrame$Number)
    sample <- rev(unique(finalFrame$Sample))
    sampleNum <- unlist(lapply(sample,function(x){
        length(strsplit(as.character(x) ,"∩")[[1]])
    }))
    dataNum <- unlist(lapply(sample,function(x){
        sum(finalFrame[finalFrame$Sample == x,]$Number)
    }))
    sampleOrder <- c()
    pl <- 0
    for(i in unique(sampleNum)){
        posList <- order(dataNum[which(sampleNum == i)],decreasing = TRUE)
        posList <- posList + pl
        pl <- length(posList) + pl
        sampleOrder <- append(sampleOrder, posList) 
    }
    finalFrame$Sample  <-  factor(finalFrame$Sample, levels = levels(finalFrame$Sample)[sampleOrder])
    pointLineFrame$combinations <- factor(pointLineFrame$combinations, levels = levels(finalFrame$Sample))
    #draw stack plots
    keyPoint <- ggplot(finalFrame)+
        aes(x = Sample,y = Number,fill = Type)+
        geom_bar(stat = "identity",position = "stack",width = 0.7)+
        labs(x = "",width = 1.0)+
        labs(y = "Mutation number")
    if (show.num == "TRUE") {
        keyPoint <- keyPoint +
            geom_text(aes(label = Number), size = 3, colour = 'black',
                      hjust = .5, position = position_stack(vjust=0.5))
    }
    keyPoint <- keyPoint +
        theme(axis.ticks.x  = element_blank(),
              panel.border =element_blank(),
              axis.text.x =element_blank(),
              panel.grid.major =element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line.y  = element_line(colour = "black", size = 0.7),
              axis.line.x = element_blank(),
              axis.title.y = element_text(size=14),
              axis.text.y = element_text(size=12, colour = "black"),
              plot.margin = unit(c(0.08,0.2,0,0.1),"inches"),
              legend.spacing  = unit(c(0.09,0,0,0),"inches"),
              legend.key.width  = unit(0.2, "inches"),
              legend.text = element_text(size=12, colour = "black"))+
        scale_fill_manual(values =c( "#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2",
                                     "#8491B4B2","#91D1C2B2","#DC0000B2","#7E6148B2","#91D1C2B2",
                                     "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"))+
        scale_y_continuous(expand = c(0,0))
    # draw point-line plot
    pointLinePlot <- ggplot(pointLineFrame)+
        aes(x=combinations,y=sample)
    if(length(levels(finalFrame$Sample))<21){
        pointLinePlot <- pointLinePlot + geom_point(size=3.5)
    }else{
        pointLinePlot <- pointLinePlot + geom_point(size=2.5)
    }
    pointLinePlot <- pointLinePlot+
        geom_path(mapping = aes(group=combinations),inherit.aes = TRUE)+
        labs(x="",width=1.0)+
        labs(y="")+
        theme(panel.grid = element_blank(),
              panel.border =element_blank(),
              axis.text.x =element_blank(),
              axis.title.y = element_text(size=14),
              axis.text.y = element_text(size=12, colour = "black"),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.x = element_blank(),
              axis.line.y = element_blank(),
              plot.margin = unit(c(0.15,0,0,0.76),"inches"))+
        scale_y_discrete(position = "right")
    # Keep the width of two plots the same
    pointLinePlot  <-  ggplot_gtable(ggplot_build(pointLinePlot))
    barPlot  <-  ggplot_gtable(ggplot_build(keyPoint)) 
    pointLinePlot$widths  <-  barPlot$widths
    # put pictures together
    if(length(levels(finalFrame$Sample)) < 21){
        gg <- ggdraw()+draw_plot(barPlot,0,0.3,1,0.7)+draw_plot(pointLinePlot,0,0,1,0.35)
    }else{
        if (length(levels(finalFrame$Sample)) > 60) {
            gg <- ggdraw()+draw_plot(barPlot,0,0.3,1,0.7)+draw_plot(pointLinePlot,0,0,1,0.35)
        }else{
            gg <- ggdraw()+draw_plot(barPlot,0,0.3,1,0.7)+draw_plot(pointLinePlot,0,0,1,0.35)
        }
    }  
    return(gg)
}