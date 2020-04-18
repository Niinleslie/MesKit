
testPowerLaw <- function(
   df, 
   min.vaf, 
   max.vaf, 
   min.depth,
   R2.threshold,
   min.mut.count,
   plot,
   withinType){
   
   subPowerLaw <- function(
      sample.df,
      min.vaf,
      max.vaf,
      min.depth,
      R2.threshold,
      min.mut.count,
      plot,
      withinType){
      
      if(withinType){
         type.id <- unique(sample.df$Tumor_Type)
         sample.df <- sample.df %>% 
            dplyr::mutate(VAF_adj = Type_Average_VAF_adj)
      }
      else{
         sample.id <- unique(sample.df$Tumor_Sample_Barcode)
      }
      sample.df <- sample.df %>%
         dplyr::filter(
            Clonal_Status == "Subclonal" & 
               Ref_allele_depth + Alt_allele_depth > min.depth &
               VAF_adj > min.vaf &
               VAF_adj < max.vaf &
               !is.na(VAF_adj)
         )  
      if(nrow(sample.df) < min.mut.count){
         R2.out = data.frame()
         vaf.plot  = NA
         if(withinType){
            warning(paste0("Sample ", type.id, ": There is no enough eligible mutations can be used."))
         }else{
            warning(paste0("Sample ", sample.id, ": There is no enough eligible mutations can be used."))
         }
      }else{
         vaf <- sample.df$VAF_adj
         # max.vaf <- max(sample.df$VAF_adj)
         # min.vaf <- min(sample.df$VAF_adj)
         # 
         ## get cumulative distribution 
         
         # breaks = ceiling((max.vaf - min.vaf)/0.005)
         # vafCount =  hist(1/sample.df$VAF_adj, breaks = breaks, plot = F)
         # vafCumsum <- data.frame(inv_f = vafCount$breaks, count = c(0,cumsum(vafCount$counts)))
         
         breaks <- seq(max.vaf, min.vaf, -0.005)
         mut.count <- sapply(breaks,function(x,vaf){sum(vaf > x)},vaf = vaf)
         vafCumsum <- data.frame(count = mut.count, f = breaks)
         vafCumsum$inv_f <- 1/vafCumsum$f - 1/max.vaf
         vafCumsum$n_count <- vafCumsum$count/max(vafCumsum)
         vafCumsum$t_count <- vafCumsum$inv_f/(1/min.vaf - 1/max.vaf)
         
         ## area of theoretical curve
         theoryA <- integrate(approxfun(vafCumsum$inv_f,vafCumsum$t_count),
                              min(vafCumsum$inv_f),
                              max(vafCumsum$inv_f),stop.on.error = F)$value
         # area of emprical curve
         dataA <- integrate(approxfun(vafCumsum$inv_f,vafCumsum$n_count),
                            min(vafCumsum$inv_f),
                            max(vafCumsum$inv_f),stop.on.error = F)$value
         # Take absolute difference between the two
         area <- abs(theoryA - dataA)
         # Normalize so that metric is invariant to chosen limits
         area<- area / (1 / min.vaf - 1 / max.vaf)
         
         
         ## calculate mean distance
         meandist <- mean(abs(vafCumsum$n_count - vafCumsum$t_count))
         
         ## calculate kolmogorovdist 
         n = length(vaf)
         cdfs <- 1 - ((1/sort(vaf) - 1/max.vaf) /(1/min.vaf - 1/max.vaf))
         dp <- max((1:n) / n - cdfs)
         dn <- - min((0:(n-1)) / n - cdfs)
         kolmogorovdist  <- max(c(dn, dp))
         
         ## R squared
         lmModel <- lm(vafCumsum$count ~ vafCumsum$inv_f + 0)
         lmLine = summary(lmModel)
         R2 = lmLine$adj.r.squared
         
         
         if(withinType){
            R2.out <- data.frame(
               Patient = as.character(unique(sample.df$Patient_ID)),
               Tumor_Type = type.id,
               Eligible_Mut_Count = nrow(sample.df ),
               Area = area,
               Kolmogorov_Distance = kolmogorovdist,
               Mean_Distance = meandist,
               R2 = R2, 
               Type = dplyr::if_else(
                  R2 >= R2.threshold,
                  "neutral",
                  "non-neutral") 
            )
         }else{
            R2.out <- data.frame(
               Patient = as.character(unique(sample.df$Patient_ID)),
               Tumor_Sample_Barcode = sample.id,
               Eligible_Mut_Count = nrow(sample.df ),
               Area = area,
               Kolmogorov_Distance = kolmogorovdist,
               Mean_Distance = meandist,
               R2 = R2, 
               Type = dplyr::if_else(
                  R2 >= R2.threshold,
                  "neutral",
                  "non-neutral") 
            ) 
         }
         
         vaf.plot <- NA 
         if(plot){
            
            Arealabel <- as.character(paste0("italic(Area) == ", round(area,4)))
            KDlabel <- as.character(paste0("italic(Kolmogorov_Distance) ==", round(kolmogorovdist,4) ))
            Mdlabel <- as.character(paste0("italic(Mean_Distance) == ", round(meandist,4) ))
            R2label <- as.character(paste0("italic(R)^2 == ", round(R2,4) ))
            
            x.min <- min(vafCumsum$f)
            x.max <- max(vafCumsum$f)
            x.breaks <- seq(x.min,x.max,(x.max-x.min)/2)
            x.breaks.pos <- 1/x.breaks - 1/max.vaf
            x.breaks.label <- paste("1/", round(x.breaks,2),sep="")
            y.min <- min(vafCumsum$count)
            y.max <- max(vafCumsum$count,
                         (max(vafCumsum$inv_f)*lmModel$coefficients[1]))
            vaf.plot <- ggplot(data = vafCumsum, mapping = aes(x = inv_f, y = count)) +
               geom_point()+
               geom_smooth(method=lm,formula = y ~ x + 0, color="red",se = FALSE)+
               theme(panel.grid =element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_blank(),
                     axis.ticks = element_line(size = 1),
                     axis.title = element_text(size = 13,face = "bold",colour = "black"),
                     axis.text = element_text(size = 10,face = "bold",colour = "black"),
                     axis.ticks.length = unit(.25, "cm"))+
               geom_segment(aes(x = min(inv_f),xend= max(inv_f), y=-Inf,yend=-Inf), size = 1)+
               geom_segment(aes(y = y.min ,yend = y.max,x=-Inf,xend=-Inf), size = 1.5)+
               scale_x_continuous(breaks = x.breaks.pos,
                                  labels = x.breaks.label)+
               scale_y_continuous(breaks = seq(y.min,y.max,(y.max-y.min)/4),
                                  labels = round(seq(y.min,y.max,(y.max-y.min)/4)))+
               # labels = c(round(y.min),
               #            round(y.min+(y.max-y.min)/4),
               #            round(y.min+(y.max-y.min)*2/4),
               #            round(y.max)))+
               labs(title= unique(sample.df$Tumor_Sample_Barcode),
                    x="Inverse allelic frequency 1/vaf",
                    y="Cumulative number of SSNVs")+
               annotate("text",
                        x = x.min,
                        y = y.max,
                        label = Arealabel,
                        size = 4,
                        fontface = "bold",
                        parse = TRUE,
                        hjust = 0)+
               annotate("text",
                        x = x.min,
                        y = y.max*0.95,
                        label = KDlabel,
                        size = 4,
                        fontface = "bold",
                        parse = TRUE,
                        hjust = 0)+
               annotate("text",
                        x = x.min,
                        y = y.max*0.9,
                        label = Mdlabel,
                        size = 4,
                        fontface = "bold",
                        parse = TRUE,
                        hjust = 0)+
               annotate("text",
                        x = x.min,
                        y = y.max*0.85,
                        label = R2label,
                        size = 4,
                        fontface = "bold",
                        parse = TRUE,
                        hjust = 0)
            
            if(withinType){
               vaf.plot <- vaf.plot +  
                  labs(title= unique(sample.df$Tumor_Type),
                       x="Inverse allelic frequency 1/vaf",
                       y="Cumulative number of SSNVs")
            }
            else{
               vaf.plot <- vaf.plot + labs(title= unique(sample.df$Tumor_Sample_Barcode),
                                           x="Inverse allelic frequency 1/vaf",
                                           y="Cumulative number of SSNVs")
            }
         }
         
      }
      
      return(list(model.fitting.out = R2.out, model.fitting.plot = vaf.plot))		
   }
   
   if(withinType){
      patient.R2 <- df %>%
         dplyr::group_by(Tumor_Type) %>%
         dplyr::group_map(~subPowerLaw(.,
                                       min.vaf,
                                       max.vaf,
                                       min.depth,
                                       R2.threshold,
                                       plot,
                                       min.mut.count = min.mut.count,
                                       withinType = withinType), 
                          keep = TRUE) %>%
         rlang::set_names(unique(df$Tumor_Type))
   }else{
      patient.R2 <- df %>%
         dplyr::group_by(Tumor_Sample_Barcode) %>%
         dplyr::group_map(~subPowerLaw(.,
                                       min.vaf,
                                       max.vaf,
                                       min.depth,
                                       R2.threshold,
                                       plot,
                                       min.mut.count = min.mut.count,
                                       withinType = withinType), 
                          keep = TRUE) %>%
         rlang::set_names(unique(df$Tumor_Sample_Barcode))
   }
   
   
   
   return(patient.R2)       
   #return(list(data.frame = neutralTest.out, plot.list = patient.plot))    	
}
