plotPowerLaw <- function(vafCumsum, test.df, id, max.vaf, lmModel, patient){
 area <- test.df$Area
 kolmogorovdist = test.df$Kolmogorov_Distance
 meandist = test.df$Mean_Distance 
 R2 = test.df$R2
 Arealabel <- as.character(paste0("italic(Area) == ", round(area,4)))
 KDlabel <- as.character(paste0("italic(Kolmogorov_Distance) ==", round(kolmogorovdist,4) ))
 Mdlabel <- as.character(paste0("italic(Mean_Distance) == ", round(meandist,4) ))
 R2label <- as.character(paste0("italic(R)^2 == ", round(R2,4) ))
 
 x.min <- min(vafCumsum$f)
 x.max <- max(vafCumsum$f)
 x.breaks <- seq(x.min, x.max, (x.max-x.min)/2)
 x.breaks.pos <- 1/x.breaks - 1/max.vaf
 x.breaks.label <- paste("1/", round(x.breaks,2),sep="")
 y.min <-  0
 y.max <- max(vafCumsum$count,
              (max(vafCumsum$inv_f)*lmModel$coefficients[1]))
 p <- NA
 p <- ggplot(data = vafCumsum, mapping = aes(x = inv_f, y = count)) +
  geom_point()+
  geom_smooth(method=lm,formula = y ~ x + 0, color="red",se = FALSE)+
  theme(panel.grid =element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 1),
        axis.title = element_text(size = 13,face = "bold",colour = "black"),
        axis.text = element_text(size = 10,face = "bold",colour = "black"),
        plot.title = element_text(size = 13.5,face = "bold"),
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
  labs(title= paste0(patient, ":", id) ,
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
 
 return(p)
}

