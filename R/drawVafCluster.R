
drawVAFCombine <- function(subdata, xlab){
    
   id <- unique(subdata$ID)
   patient <- unique(subdata$Patient_ID)

   # build color vector for later use
   color_scale <- c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", 
                    "#008280FF", "#BB0021FF", "#5F559BFF", "#A20056FF", 
                    "#808180FF", "#1B1919FF")
   
   names(color_scale) <- c(seq_len(9), "outlier")
   
   ## initialize variable in ggplot for biocheck error
   VAF <- NULL
   cluster <- NULL
   ## generate  plot and specific titles for minifigures
   scaleFUN <- function(x) sprintf("%.1f", x)
   p <- ggplot(subdata, aes(x = V)) + 
         theme_bw() +
         theme(legend.position='right', 
               plot.title=element_text(size=13.5,hjust = 0,vjust = 0.5,face = "bold"),
               panel.grid=element_blank(), 
               panel.border=element_blank(),
               axis.line=element_line(size=0.7), 
               axis.title=element_text(size=13),
               axis.text=element_text(size=12, colour = "black"))+
          geom_line(size=1, colour="#00C0EB", stat="density") + 
          geom_point(aes(y=0, colour=cluster), alpha=0.5) + 
          # drawVAFCombineVline(subdata) + 
          # geom_rug(aes(y=0, colour=cluster), sides="b") + 
          scale_colour_manual(values=color_scale) + 
          ggtitle(paste0(patient,": ", id)) + 
          labs(y = "Density", colour = "Cluster",x = xlab) +
            scale_y_continuous(labels=scaleFUN)
   
   return(p)
}


## Draw vlines for all plotOption
drawVAFCombineVline <- function(subdata){

   ## initialize variable in ggplot for biocheck error
   VAF <- NULL
   cluster <- NULL
   x <- NULL
   xend <- NULL
   y <- NULL
   yend <- NULL

   ## A draft for density infomation(density_info) of ggplot
   picv <- ggplot(subdata, aes(x=V)) +
        geom_line(size=1, colour="#00C0EB", stat="density")
   ## density information of the curve for a tsb
   densityInfo <- data.frame(layer_data(picv))

   df_vline <- data.frame()
   cluster_list <- unique(subdata$cluster)
   ## Obtain vline Coordinate(x, xend, y, yend)
   for (cluster_name in cluster_list){
      x_end <- max(subdata[which(
         subdata$cluster == cluster_name), ]$VAF)
      x_end_alter <- densityInfo$x[which.min(
         abs(outer(densityInfo$x,x_end,FUN="-")))]
      y <- 0
      y_end <- densityInfo$y[which(
         densityInfo$x == x_end_alter)]
      sub <- data.frame(x = x_end_alter, xend = x_end_alter, y = y, yend = y_end)
      df_vline <- rbind(df_vline, sub)
   }
   vline <- geom_segment(aes(x= x, xend= xend, y= y, yend= yend),
                         data = df_vline,
                        size=0.5, colour= "grey", linetype= "dashed")
   return(vline)
}

# ## Functions for specific plotOption: "compare"
# ## VAF painter for OFA
# drawVAFCompare <- function(maf_data,withinTumor){
#     min.vaf <- min(maf_data$VAF)
#     max.vaf <- max(maf_data$VAF)
#     patient <- unique(maf_data$Patient_ID)
#     
#     # build color vector for later use
#     color_scale <- c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", 
#                      "#008280FF", "#BB0021FF", "#5F559BFF", "#A20056FF", 
#                      "#808180FF", "#1B1919FF")
#     names(color_scale) <- c(seq_len(9), "outlier")
#     pic <- ggplot(maf_data,aes(x=VAF, y=ID)) +
#            theme_bw() + 
#             theme(plot.title=element_text(size=16, hjust=0.5, vjust=0.5, face='bold'), 
#                   panel.grid=element_blank(), 
#                   panel.border=element_blank(),
#                   axis.title=element_text(size=13), 
#                   axis.text=element_text(size=12, colour = "black"),
#                   axis.line=element_line(size=0.7)) +  
#             ggtitle(paste0("VAF clustering of ", patient)) +  
#             ggridges::geom_density_ridges(fill="whitesmoke",calc_ecdf=TRUE, alpha=0.5) + 
#             geom_point(aes(x=VAF, y=ID, color=cluster),alpha=0.5, show.legend=TRUE) +  
#             ggridges::geom_density_ridges(color="#00C0EB",fill=NA, calc_ecdf=TRUE, alpha=0.5, size=1)+
#             drawVAFCompareVline(maf_data)+ 
#             scale_colour_manual(values=color_scale) +  
#             scale_x_continuous(limits = c(min.vaf,max.vaf))
#     
#     if(withinTumor){
#         pic <- pic + labs(y = "Tumor")
#     }else{
#         pic <- pic + labs(y = "Sample")
#     }
#     return(pic)
# }
# 
# 
# ## VAF draw vlines for ofa
# drawVAFCompareVline <- function(maf_data)
# {
#    ## density information of the curve 
#    gr <- ggplot(maf_data, aes(x=VAF, y=ID)) + 
#       ggridges::geom_density_ridges() 
#    ingredients <- ggplot_build(gr) %>% purrr::pluck("data", 1)
#    iscale <- ingredients$iscale[1]
#    scale <- ingredients$scale[1]
#       
#    df_vline <- data.frame()
#    for (id in unique(maf_data$ID)){
#       subdata <- maf_data[ID == id]
#       cluster_list <- unique(subdata$cluster)
#       ## A draft for density infomation(density_info) of ggplot
#       picv <- ggplot(subdata, aes(x=VAF)) + 
#          geom_line(size=1, colour="#00C0EB", stat="density")
#       densityInfo <- data.frame(layer_data(picv))
#       for (cluster_name in cluster_list){
#           x_end <- max(subdata[cluster == cluster_name]$VAF)
#           x_end_alter <- densityInfo$x[which.min(
#               abs(outer(densityInfo$x,x_end,FUN="-")))]
#           y <- which(id == unique(maf_data$ID))
#           density <- densityInfo$density[which(
#               densityInfo$x == x_end_alter)]
#           y_end <- y+ density*iscale*scale
#           sub <- data.frame(x = x_end_alter, xend = x_end_alter, y = y, yend = y_end)
#           df_vline <- rbind(df_vline, sub)
#       }
#    }
#    vline <- geom_segment(aes(x= x, xend= xend, y= y, yend= yend),
#                          data = df_vline,
#                          size=0.5, colour= "grey", linetype= "dashed")
#    return(vline)
# }