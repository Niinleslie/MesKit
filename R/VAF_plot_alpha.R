#' A VAF plot painter
#' @description Read maf file as data.frame, cluster data with maftools and finally draw variant allele frequency(VAF)
#' frequency distribution curve with ggplot2 as well as ggridges. We could use different parameters to control output 
#' images from different samples or conclude all samples' VAF information in one image.
#' 
#' @import ggplot2 maftools ggridges ggsci dplyr
#' 
#' @param maf_file specify a maf document/directory as the input of the function
#' @param sample_option specify single/all sample names (Tumor_Sample_Barcodes, tsb). Default "OFA".
#' @param theme_option select a coloring scheme from ggsci. Default "aaas".
#' @param file_format choose an output file format accessable for ggsave. Default "png".
#' @return Images of selected samples' VAF
#'
#' @examples
#' \dontrun{
#' VAF_plot(maf_file, sample_option = "OFA", theme_option = "aaas") # draw a VAF image that contains all samples' VAF distribution curves with different themes.
#' VAF_plot(maf_file, sample_option = "All") # draw VAF images for every sample respectively
#' VAF_plot(maf_file, sample_option = "tsb1", file_format = "pdf") # draw a VAF image for sample tsb1 and save as a pdf file.
#' 
#'}


# import pkgs
library(ggplot2)
library(maftools)
library(ggridges)
library(ggsci)
library(dplyr)

# source MATH_score function
setwd("/home/ninomoriaty/R_Project/MesKit/R")
source("MATH_Score.R")

############ Major function ############
VAF_plot <-function(maf, sample_option = "OFA", theme_option = "aaas", file_format = "png", show.MATH = T)
{
  # read .maf file
  maf_input <- maf@data[,-c("patient", "lesion", "time")]
  laml <- read.maf(maf_input)
  # extract vaf info
  n <- length(maf_input$Hugo_Symbol)
  vaf_input_mt <- data.frame(maf_input$Hugo_Symbol, maf_input$VAF, maf_input$Tumor_Sample_Barcode)
  colnames(vaf_input_mt) <- c("Hugo_Symbol", "VAF", "Samples")
  cluster_all <- data.frame()
  # extract all tumor sample barcode
  tsb_ls <- data.frame(unique(vaf_input_mt$Samples))
  colnames(tsb_ls) <- c("samples")
  

  
  # sample options
  if (sample_option == "All")
  {
  # print all samples respectively
  for (counter_mt in 1:length(tsb_ls[,1]))
    {
    for (sample_name_mt in tsb_ls)
      {
      # calculate MATH_score
      if (show.MATH){
        MATH.score <- MATH_score(maf_input, c(as.character(sample_name_mt)[counter_mt]))
        MATH.score <- MATH.score[which(MATH.score$Tumor_Sample_Barcode == as.character(sample_name_mt)[counter_mt]), ]$MATH_score
      } else {
        MATH.score <- NA
      }
      sample_mt <- vaf_input_mt[which(vaf_input_mt$Samples %in% as.character(sample_name_mt)[counter_mt]),]
      cluster_mt = inferHeterogeneity(maf = laml, tsb = as.character(sample_mt[1,3]), vafCol = 'VAF', useSyn = TRUE)$"clusterData"
      colnames(cluster_mt)[6] = "VAF"
      # print VAF pictures for all samples
      pic <- VAF_draw(cluster_mt, theme_option, as.character(sample_name_mt)[counter_mt], MATH.score)
      ggsave(pic, filename =  paste(as.character(sample_name_mt)[counter_mt], "_VAF_Cluster", ".", file_format,sep=""), width = 12, height = 9)
      }
    }
  } else if (sample_option == "OFA")
  {
    # one pic for all sample
    patientID <- maf@patientID
    # collect all samples' cluster results
    for (counter_mt in 1:length(tsb_ls[,1]))
      {
      for (sample_name_mt in tsb_ls)
        {
        # calculate MATH_score
        if (show.MATH){
          MATH.score <- MATH_score(maf_input, c(sample_name_mt))
          MATH.score <- MATH.score[which(MATH.score$Tumor_Sample_Barcode == "ITH MATH score"), ]$MATH_score
        } else {
          MATH.score <- NA
        }
        sample_mt <- vaf_input_mt[which(vaf_input_mt$Samples %in% as.character(sample_name_mt)[counter_mt]),]
        cluster_mt_cha <- paste("cluster_mt_", counter_mt," <- inferHeterogeneity(maf = laml, tsb = as.character(sample_mt[1,3]), vafCol = \'VAF\', useSyn = TRUE)$\"clusterData\"", sep ="")
        eval(parse(text = cluster_mt_cha))
        cluster_mt_cha <- paste("colnames(cluster_mt_", counter_mt, ")[6] = \"VAF\"",sep ="")
        eval(parse(text = cluster_mt_cha))
        cluster_mt_cha <- paste("cluster_all <- rbind(cluster_all, cluster_mt_", counter_mt, ")",sep ="")
        eval(parse(text = cluster_mt_cha))
        }
      }
    colnames(cluster_all)[6] = "VAF"
    
    pic <- eval(parse(text = VAF_OFA(cluster_all, theme_option, tsb_ls, sample_option, MATH.score)))
    ggsave(pic, filename =  paste(patientID, "_VAF_Cluster",".", file_format, sep=""), width = 12, height = 9)
  } else 
  {
  # specific sample
    # calculate MATH_score
    if (show.MATH){
      MATH.score <- MATH_score(maf_input, c(sample_option))
      MATH.score <- MATH.score[which(MATH.score$Tumor_Sample_Barcode == sample_option), ]$MATH_score
    } else {
      MATH.score <- NA
    }
    sample_mt <- vaf_input_mt[which(vaf_input_mt$Samples %in% sample_option),]
    cluster_mt = inferHeterogeneity(maf = laml, tsb = as.character(sample_mt[1,3]), vafCol = 'VAF', useSyn = TRUE)$"clusterData"
    colnames(cluster_mt)[6] = "VAF"
    
    pic <- VAF_draw(cluster_mt, theme_option, sample_option, MATH.score)
    ggsave(pic, filename =  paste(sample_option,"_VAF_Cluster",".", file_format,sep=""), width = 12, height = 9)
  }
}

############ General Toolbox ############
# VAF main draw vlines for all sample_option
VAF_vline <- function(cluster_mt, pic, tsb_ls, sample_option, tsb, ingredients = NULL)
{
  # data prepare
  VAF_vline_cha <- ""
  cluster_ls <- unique(cluster_mt$cluster)
  # density information of the curve for a tsb
  density_info <- data.frame(layer_data(pic))
  
  # OFA specific: get scaling ratio
  if (typeof(ingredients) == "list")
  {
    iscale <- ingredients$iscale[1]
    scale <- ingredients$scale[1]
  }
  
  # Obtain vline Coordinate(x, xend, y, yend)
  for (cluster_name in cluster_ls)
  {
    x_end <- max(cluster_mt[which(cluster_mt$cluster == cluster_name)]$VAF)
    x_end_alter <- density_info$x[which.min(abs(outer(density_info$x,x_end,FUN="-")))]
    y_end <- density_info$y[which(density_info$x == x_end_alter)]
    if (sample_option == "OFA")
    {
      # Scale and draw lines
      density <- density_info$density[which(density_info$x == x_end_alter)]
      VAF_vline_cha <- paste(VAF_vline_cha, "geom_segment(data = cluster_mt_", which(tsb_ls == tsb)," ,aes(x = ", x_end_alter,", xend = ", x_end_alter, ", y = ", which(tsb_ls == tsb),", yend = ", which(tsb_ls == tsb) + density*iscale*scale,"), size = 0.3, colour=\"#00AAFF\", linetype=\"dashed\") + ",sep="")
    }else
    {
      VAF_vline_cha <- paste(VAF_vline_cha, "geom_segment(aes(x = ", x_end_alter,", xend = ", x_end_alter, ", y = 0, yend = ", y_end,"), size = 0.3, colour=\"#00AAFF\", linetype=\"dashed\") + ",sep="")
    }
  }
  VAF_vline_cha
}

############ Functions for specific sample_option: "All","tsb" ############
# VAF painter
VAF_draw <- function(cluster_mt, theme_option, sample_option, MATH.score){
  # A draft for density infomation(density_info) of ggplot
  picv <- ggplot(cluster_mt, aes(x = VAF)) + geom_line(size = 1, colour = "#00AAFF", stat = "density")
  if (is.na(MATH.score)){
    # generate character/string for ggplot and paint the picture
    VAF_draw_cha = paste("ggplot(cluster_mt, aes(x = VAF)) + 
                       theme_bw() + 
                       theme(title=element_text(size = 18), text = element_text(size = 18), panel.grid=element_blank(),panel.border=element_blank(), axis.line=element_line(size=0.25)) + 
                       geom_line(size = 1, colour = \"#00AAFF\", stat = \"density\") + geom_rug(aes(y = 0, colour = cluster), sides = \"b\") + ", 
                         VAF_vline(cluster_mt, picv, tsb_ls, sample_option),
                         "scale_color_", theme_option, "() + scale_fill_", theme_option, "()", sep="")
    eval(parse(text = VAF_draw_cha))
  } else {
    VAF_draw_cha = paste("ggplot(cluster_mt, aes(x = VAF)) + 
                       theme_bw() + 
                       theme(plot.title = element_text(size=18, hjust=1, vjust=0.5, face='bold'), title=element_text(size = 18), text = element_text(size = 18), panel.grid=element_blank(),panel.border=element_blank(), axis.line=element_line(size=0.25)) + 
                       ggtitle(\"MATH Score: ", as.character(MATH.score), "\") + 
                       geom_line(size = 1, colour = \"#00AAFF\", stat = \"density\") + geom_rug(aes(y = 0, colour = cluster), sides = \"b\") + ", 
                         VAF_vline(cluster_mt, picv, tsb_ls, sample_option),
                         "scale_color_", theme_option, "() + scale_fill_", theme_option, "()", sep="")
    eval(parse(text = VAF_draw_cha))
  }
  }

############ Functions for specific sample_option: "OFA" ############
# VAF drawer for OFA: generate character/string for follow-up painting with ggplot. 
VAF_OFA <- function(cluster_all, theme_option, tsb_ls, sample_option, MATH.score){
  if (is.na(MATH.score)){
    VAF_ofa_cha = paste("ggplot(cluster_all, aes(x=VAF, y=Tumor_Sample_Barcode)) +
                      theme_bw() + 
                      theme(title=element_text(size = 18), text = element_text(size = 18), panel.grid=element_blank(),panel.border=element_blank(), axis.line=element_line(size=0.25)) + 
                      geom_point(aes(x=VAF, y=Tumor_Sample_Barcode, color = cluster), alpha = 0.5) +
                      geom_density_ridges(color = \"#00AAFF\", fill = \"whitesmoke\", calc_ecdf = TRUE, alpha = 0.5) + ",
                        VAF_vline_ofa(cluster_all, tsb_ls, sample_option), 
                        "scale_color_", theme_option, "() + scale_fill_", theme_option, "()", sep="")
  } else {
    VAF_ofa_cha = paste("ggplot(cluster_all, aes(x=VAF, y=Tumor_Sample_Barcode)) +
                      theme_bw() + 
                      theme(plot.title = element_text(size=18, hjust=1, vjust=0.5, face='bold'), title=element_text(size = 18), text = element_text(size = 18), panel.grid=element_blank(),panel.border=element_blank(), axis.line=element_line(size=0.25)) + 
                      ggtitle(\"MATH Score: ", as.character(MATH.score), "\") + 
                      geom_point(aes(x=VAF, y=Tumor_Sample_Barcode, color = cluster), alpha = 0.5) +
                      geom_density_ridges(color = \"#00AAFF\", fill = \"whitesmoke\", calc_ecdf = TRUE, alpha = 0.5) + ",
                        VAF_vline_ofa(cluster_all, tsb_ls, sample_option), 
                        "scale_color_", theme_option, "() + scale_fill_", theme_option, "()", sep="")
  }
  VAF_ofa_cha
}

# VAF draw vlines for ofa
VAF_vline_ofa <- function(cluster_all, tsb_ls, sample_option)
{
  # data prepare
  VAF_vline_ofa <- ""
  # density information of the curve for all tsbs
  gr <- ggplot(cluster_all, aes(x = VAF, y = Tumor_Sample_Barcode)) + geom_density_ridges() 
  ingredients <- ggplot_build(gr) %>% purrr::pluck("data", 1)
  
  for (tsb in tsb_ls$samples)
  {
    # renew tsb's data
    VAF_vline_cha <- ""
    x_end_ls <- data.frame()
    y_end_ls <- data.frame()
    cluster_mt <- cluster_all[which(cluster_all$Tumor_Sample_Barcode == tsb)]
    colnames(cluster_mt)[6] = "VAF"
    cluster_ls <- unique(cluster_mt$cluster)
    # A draft for density infomation(density_info) of ggplot
    picv <- ggplot(cluster_mt, aes(x = VAF)) + geom_line(size = 1, colour = "#00AAFF", stat = "density")
    density_info <- data.frame(layer_data(picv))
    # collect vlines for a single tsb
    VAF_vline_cha <- VAF_vline(cluster_mt, picv, tsb_ls, sample_option, tsb, ingredients)
    # collect vlines for all tsbs
    VAF_vline_ofa <- paste(VAF_vline_ofa, VAF_vline_cha)
  }
  VAF_vline_ofa
}


########## Directory #######
# maf_file1 = "/home/ninomoriaty/R_Project/data/maf/311252.maf"
# maf_file2 = "/home/ninomoriaty/R_Project/data/maf/313544.maf"
# maf_file3 = "/home/ninomoriaty/R_Project/data/maf/313935.maf"
# maf_file4 = "/home/ninomoriaty/R_Project/data/maf/313953.maf"
# maf_file5 = "/home/ninomoriaty/R_Project/data/maf/314007.maf"
# maf_file6 = "/home/ninomoriaty/R_Project/data/maf/314069.maf"
# maf_file7 = "/home/ninomoriaty/R_Project/data/maf/314155.maf"
# maf_file_ls = c(maf_file1, maf_file2, maf_file3, maf_file4, maf_file5, maf_file6, maf_file7)
# for (counter in maf_file_ls){
#   VAF_plot(counter, "OFA")
#   VAF_plot(counter, "All")
# }
