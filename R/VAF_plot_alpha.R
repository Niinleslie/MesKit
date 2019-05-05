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

############ Major function ############
VAF_plot <-function(maf_file, sample_option = "OFA", theme_option = "aaas", file_format = "png")
{
  # read .maf file
  maf_input <- read.table(maf_file, header = TRUE, fill = TRUE, sep = '\t', quote = "")
  laml <- read.maf(maf=maf_file)
  # extract vaf info
  samples <- data.frame(maf_input[,ncol(maf_input)])
  cluster_all <- data.frame()
  n <- length(maf_input[,1])
  vaf_input_mt <- data.frame(maf_input[,1], maf_input[,ncol(maf_input)-3], samples)
  colnames(vaf_input_mt) <- c("Hugo_Symbol", "VAF", "Samples")
  # extract all tumor sample barcode
  tsb_ls <- as.data.frame(as.data.frame(table(samples))["samples"][which(as.data.frame(table(samples))["samples"]$samples != ""),])
  colnames(tsb_ls) <- c("samples")
  
  # sample options
  if (sample_option == "All")
  {
  # all samples and output respectively
  for (counter_mt in 1:length(tsb_ls[,1]))
    {
    for (sample_name_mt in tsb_ls)
      {
      sample_mt <- vaf_input_mt[which(vaf_input_mt$Samples %in% as.character(sample_name_mt)[counter_mt]),]
      cluster_mt = inferHeterogeneity(maf = laml, tsb = as.character(sample_mt[1,3]), vafCol = 'VAF', useSyn = TRUE)$"clusterData"
      colnames(cluster_mt)[6] = "VAF"
      pic <- VAF_draw(cluster_mt, theme_option, sample_option)
      ggsave(pic, filename =  paste(as.character(sample_name_mt)[counter_mt], "_VAF_Cluster", ".", file_format,sep=""), width = 12, height = 9)
      }
    }
  } else if (sample_option == "OFA")
  {
    # one pic for all sample
    patientID = strsplit(as.character(tsb_ls[1,1]), "-")[[1]][1]
    # collect all samples' cluster results
    for (counter_mt in 1:length(tsb_ls[,1]))
      {
      for (sample_name_mt in tsb_ls)
        {
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
    
    pic <- eval(parse(text = VAF_OFA(cluster_all, theme_option, tsb_ls, sample_option)))
    ggsave(pic, filename =  paste(patientID, "_VAF_Cluster",".", file_format, sep=""), width = 12, height = 9)
  } else 
  {
  # specific sample
    sample_mt <- vaf_input_mt[which(vaf_input_mt$Samples %in% sample_option),]
    cluster_mt = inferHeterogeneity(maf = laml, tsb = as.character(sample_mt[1,3]), vafCol = 'VAF', useSyn = TRUE)$"clusterData"
    colnames(cluster_mt)[6] = "VAF"
    
    pic <- VAF_draw(cluster_mt, theme_option, sample_option)
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
      VAF_vline_cha <- paste(VAF_vline_cha, "geom_segment(data = cluster_mt_", which(tsb_ls == tsb)," ,aes(x = ", x_end_alter,", xend = ", x_end_alter, ", y = ", which(tsb_ls == tsb),", yend = ", which(tsb_ls == tsb) + density*iscale*scale,"), size = 0.3, colour=\"cadetblue3\", linetype=\"dashed\") + ",sep="")
    }else
    {
      VAF_vline_cha <- paste(VAF_vline_cha, "geom_segment(aes(x = ", x_end_alter,", xend = ", x_end_alter, ", y = 0, yend = ", y_end,"), size = 0.3, colour=\"cadetblue3\", linetype=\"dashed\") + ",sep="")
    }
  }
  VAF_vline_cha
}

############ Functions for specific sample_option: "All","tsb" ############
# VAF painter
VAF_draw <- function(cluster_mt, theme_option, sample_option)
  {
  # A draft for density infomation(density_info) of ggplot
  picv <- ggplot(cluster_mt, aes(x = VAF)) + geom_line(size = 1, colour = "cadetblue3", stat = "density")
  # generate character/string for ggplot and paint the picture
  VAF_draw_cha = paste("ggplot(cluster_mt, aes(x = VAF)) + 
                       theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(), axis.line=element_line(size=0.25)) + 
                       geom_line(size = 1, colour = \"cadetblue3\", stat = \"density\") + geom_rug(aes(y = 0, colour = cluster), sides = \"b\") + ", 
                       VAF_vline(cluster_mt, picv, tsb_ls, sample_option),
                       "scale_color_", theme_option, "() + scale_fill_", theme_option, "()", sep="")
  eval(parse(text = VAF_draw_cha))
  }

############ Functions for specific sample_option: "OFA" ############
# VAF drawer for OFA: generate character/string for follow-up painting with ggplot. 
VAF_OFA <- function(cluster_all, theme_option, tsb_ls, sample_option)
{
  VAF_ofa_cha = paste("ggplot(cluster_all, aes(x=VAF, y=Tumor_Sample_Barcode)) +
                      theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(), axis.line=element_line(size=0.25)) + 
                      geom_point(aes(x=VAF, y=Tumor_Sample_Barcode, color = cluster), alpha = 0.5) +
                      geom_density_ridges(color = \"cadetblue3\", fill = \"whitesmoke\", calc_ecdf = TRUE, alpha = 0.5) + ",
                      VAF_vline_ofa(cluster_all, tsb_ls, sample_option), 
                      "scale_color_", theme_option, "() + scale_fill_", theme_option, "()", sep="")
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
    picv <- ggplot(cluster_mt, aes(x = VAF)) + geom_line(size = 1, colour = "cadetblue3", stat = "density")
    density_info <- data.frame(layer_data(picv))
    # collect vlines for a single tsb
    VAF_vline_cha <- VAF_vline(cluster_mt, picv, tsb_ls, sample_option, tsb, ingredients)
    # collect vlines for all tsbs
    VAF_vline_ofa <- paste(VAF_vline_ofa, VAF_vline_cha)
  }
  VAF_vline_ofa
}


########## Directory #######
# maf_dir = "/home/ninomoriaty/R_Project/patients_snv_indel.imputed.maf"
# maf_file = "/home/ninomoriaty/R_Project/311252_snv_indel.imputed.maf"
# sample_option = "311252-S"
# theme_option = "aaas"
# file_format = "png"

