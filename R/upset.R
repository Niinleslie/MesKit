#load the packages
library("ggsci")
library("ggplot2")
library("cowplot")

# Clean the working environment
rm(list=ls())

#'Use R code to find the intersect mutations and their types in several samples of one patient
#'
#' @param dataframe a dataframe transferred from a maf file,which contains information going to be analyzed
#' @param patientID the ID of a patient whose cancer information is going to be analyzed
#' @param text_show a logic parameter to determine whether to show the number of each mutations in the stack plot
#'
#' @return a stack plot
#' @export
#' @examples
#' supset_analysis(maf,311525,TRUE)

maf<-read.table(Arg[1],sep = "\t",header = TRUE)

paramater2<-as.character(Arg[2])

for (i in 1:length(maf$Hugo_Symbol)){
  
  maf$Patient[i]<-strsplit(as.character(maf$Tumor_Sample_Barcode[i]),split = "-")[[1]][1]
  
}

maf$Patient<-as.factor(maf$Patient)


# get the information of the chosen patient

upset_analysis<-function(dataframe=NULL,patientID=NULL,text_show=FALSE){
  
  maf.used_patient_n<-dataframe[which (dataframe$Patient==as.character(patientID)),]
  
  primitive_length<-length(maf.used_patient_n$Hugo_Symbol)
    
  upsetanalysis<-function(df=NULL,patientID=NULL,text_show = FALSE){
    
    primitive_length<-length(df$Hugo_Symbol)
    
    
    #separate by quote
    
    df$Hugo_Symbol<-as.character(df$Hugo_Symbol)
    
    new_gene<-c()
    for (n in 1:length(df$Hugo_Symbol)){
    
      new_gene<-c(new_gene,na.omit(strsplit(df$Hugo_Symbol[n],split = ",")[[1]][2:100]))
    
    }
    
    n=1
    
    while (n>0 & n<primitive_length+1){
    
      gene_number_n<-as.numeric(length(strsplit(as.character(df$Hugo_Symbol[n]),split = ",")[[1]]))
    
      if (gene_number_n>1){
    
        df[n,1]<-strsplit(df$Hugo_Symbol[n],split = ",")[[1]][1]
    
        while (gene_number_n>1){
    
          l=length(df$Hugo_Symbol)
    
          df$Hugo_Symbol<-as.character(df$Hugo_Symbol)
    
          df[l+1,2:9]<-df[n,2:9]
    
          df[l+1,1]<-as.character(strsplit(df$Hugo_Symbol[n],split = ",")[[1]][gene_number_n])
    
          gene_number_n=gene_number_n-1
    
          }
    
        }
    
      n=n+1
    }
    
    new_gene_number<-length(df$Hugo_Symbol)-primitive_length
    
    for (n in 1:new_gene_number){
    
      df[primitive_length+n,1]<-new_gene[n]
    
    }
    
    
    #discrete the different data
    
    for (n in 1:length(df$Hugo_Symbol)){
    
      if (length(strsplit(as.character(df[n,5]),split = "\t")[[1]])>2){
    
        df[n,]<-NA
    
      }
      
    }
    
    df<-df[,c(1:5,7,9,15,16)]
    
    df<-na.omit(df)
    
    
    #organize dataframe
    
    df$Tumor_Sample_Barcode<-factor(df$Tumor_Sample_Barcode,levels = unique(df$Tumor_Sample_Barcode))
    
    df$Variant_Classification<-factor(df$Variant_Classification,levels = unique(df$Variant_Classification))
    
    df<-df[order(df$Tumor_Sample_Barcode),]
    
    
    #get the intersection genes
    
    combination<-function(vector){
    
      n<-length(vector)
   
      num<-0
    
      mycycle<-1
    
      for (i in 1:n){
    
        num=num+choose(n,i)
    
      }
    
      result<-list()
    
      for (j in 1:n){
    
        oneresult<-list(combn(vector,j))
    
        result=c(result,oneresult)
    
        mycycle=mycycle+1
    
        if (mycycle==num)
          break
    
        }
    
      return(result)
    }
    
    df$Reference_Allele<-as.character(df$Reference_Allele)
    
    df$Tumor_Seq_Allele2<-as.character(df$Tumor_Seq_Allele2)
    
    sample_combination<-combination(c(levels(df$Tumor_Sample_Barcode)))
    
    for (n in 1:length(df$Hugo_Symbol)){
    
      df[n,10]<-paste(c(as.character(df[n,][c(1:4,6,7)])),collapse = "_")
    
    }
    
    colnames(df)[10]<-c("needed_data")
    
    #get intersection
    
    all_intersect<-list()
    
    get_intersect<-function(sample,information){
      
      combination_number<-length(sample)/length(sample[,1])
      
      sample_number<-length(sample[,1])
      
      n<-1
      
      while (n>0 & n<combination_number+1){
      
        to_be_intersected_n<-list()
      
        a=1
      
        while (a<sample_number+1){
      
          to_be_intersected_n<-c(to_be_intersected_n,list(information$needed_data[which (information$Tumor_Sample_Barcode==sample[,n][a])]))
      
          a<-a+1
      
        }
      
        result_n<-Reduce(intersect,to_be_intersected_n)
      
        all_intersect<<-c(all_intersect,list(length(result_n)))
      
        all_intersect<<-c(all_intersect,list(result_n))
      
        n<-n+1
      
        }
    }
    
    process<-lapply(sample_combination,get_intersect,df)
    
    
    #get the dataframe only includes genes
    
    o<-1
    
    all_intersect_only<-list()
    
    while(o>0 & o<length(all_intersect)+1){
    
      if (o%%2==0){
    
        all_intersect_only<-c(all_intersect_only,list(all_intersect[o][[1]]))
    
      }
    
      o<-o+1
    }
    
    
    #create the final dataframe
    
    all_final_frame<-data.frame(0,0,0)
    
    colnames(all_final_frame)<-c("Sample","Type","Number")
    
    all_final_frame<-all_final_frame[-1,]
    
    
    #get the first dataframe with one sample
    
    sample_number<-length(levels(df$Tumor_Sample_Barcode))
    
    single_represent<-c(1:length(levels(df$Tumor_Sample_Barcode)))

    final_frame0<-data.frame(0,0,0)
    
    colnames(final_frame0)<-c("Sample","Type","Number")
    
    final_frame0<-final_frame0[-1,]
    
    n=1
    
    while(n>0 & n<sample_number+1){
    
      sub_n<-setdiff(single_represent,single_represent[n])

      seperate_n<-setdiff(unlist(all_intersect_only[n]),unlist(all_intersect_only[sub_n]))
    
      final_seperate_n<-c()
    
      a=1
    
      while (a>0 & a<length(seperate_n)+1){
    
        final_seperate_n<-c(final_seperate_n,as.character(unique(df$Variant_Classification[which (df$needed_data==seperate_n[a])])))
    
        a=a+1
    
      }
      
    
      try_seperate_n<-data.frame(0,0,0)
    
      colnames(try_seperate_n)<-c("Sample","Type","Number")
    
      try_seperate_n[1:length(names(table(final_seperate_n))),2]<-names(table(final_seperate_n))
    
      try_seperate_n[1:length(names(table(final_seperate_n))),3]<-as.character(table(final_seperate_n))
    
      try_seperate_n$Sample<-sample_combination[[1]][n]
    
      final_frame0<-rbind(final_frame0,try_seperate_n)
    
      n=n+1
    
    }
    
    final_frame0$Type<-factor(final_frame0$Type,levels = unique(final_frame0$Type))
    
    final_frame0$Sample<-factor(final_frame0$Sample,levels = unique(final_frame0$Sample))
    
    final_frame0$Number<-as.integer(final_frame0$Number)
    
    
    #create orders
    
    sample_order_0<-data.frame(0,0)
    
    colnames(sample_order_0)<-c("Sample","Number")
    
    sample_order_0<-sample_order_0[-1,]
    
    sample_order_0[1:sample_number,1]<-levels(df$Tumor_Sample_Barcode)[1:sample_number]
    
    for (n in 1:sample_number){
    
      sample_order_0[n,2]<-sum(as.integer(final_frame0$Number[which(final_frame0$Sample==levels(final_frame0$Sample)[n])]))
    
    }
    
    sample_order_0<-sample_order_0[order(sample_order_0$Number,decreasing =TRUE),]
    
    sample_order_0$Sample<-factor(sample_order_0$Sample,levels=sample_order_0$Sample)
    
    order_sample_0<-as.character(sample_order_0$Sample)
    
    
    #organize the final dataframe with one sample
    
    final_frame0$Type<-factor(final_frame0$Type,levels = unique(final_frame0$Type))
    
    final_frame0$Sample<-factor(final_frame0$Sample,levels = order_sample_0)
    
    final_frame0$Number<-as.integer(final_frame0$Number)
    
    all_final_frame<-rbind(all_final_frame,final_frame0)
    
    
    #get other dataframes
    
    CNM<-function(a,b){
    
      return(factorial(a)/(factorial(b)*factorial(a-b)))
    
    }
    
    
    sum_CNM<-function(a,b){
    
      sum_final=0
    
      while(b>0){
    
        sum_final<-sum_final+CNM(a,b)
    
        b=b-1
    
      }
    
      return(sum_final)
    
    }
    
    
    q<-2
    
    while(q>1 & q<sample_number+1){
    
      p<-sum_CNM(sample_number,(q-1))+1
    
      final_frame1<-data.frame(0,0,0)
    
      colnames(final_frame1)<-c("Sample","Type","Number")
    
      final_frame1<-final_frame1[-1,]
    
      while(p>sum_CNM(sample_number,q-1) & p<sum_CNM(sample_number,q)+1){
    
        n=1
    
        types_n<-c()
    
        while (n>0 & n<(length(all_intersect_only[p][[1]])+1)){
    
          types_n<-c(types_n,as.character(unique(df$Variant_Classification[which (df$needed_data==all_intersect_only[p][[1]][n])])))
    
          n=n+1
    
        }
    
        if(length(types_n)!=0){
    
          try_n<-data.frame(0,0,0)
          
          colnames(try_n)<-c("Sample","Type","Number")
          
          try_n[1:length(names(table(types_n))),2]<-names(table(types_n))
          
          try_n[1:length(names(table(types_n))),3]<-as.character(table(types_n))
          
          try_n$Sample[1:length(names(table(types_n)))]<-paste(sample_combination[[q]][,(p-sum_CNM(sample_number,q-1))],collapse = "∩")
          
          try_n$Type<-factor(try_n$Type,levels = levels(df$Variant_Classification))
          
          final_frame1<-rbind(final_frame1,try_n)
    
        }
          p=p+1
    
      }
      
    
      r<-sum_CNM(sample_number,q-1)+1
    
      sample_order_q<-data.frame(0,0)
    
      colnames(sample_order_q)<-c("Sample","Number")
    
      sample_order_q<-sample_order_q[-1,]
    
      Number<-c()
    
      Sample<-c()
    
      while (r>sum_CNM(sample_number,q-1) & r<sum_CNM(sample_number,q)+1){
    
        Number<-c(Number,as.integer(all_intersect[2*r-1][[1]]))
    
        Sample<-c(Sample,paste(sample_combination[[q]][,(r-sum_CNM(sample_number,q-1))],collapse = "∩"))
    
        r=r+1
    
      }
    
      sample_order_q[1:length(Sample),1]<-Sample
    
      sample_order_q[1:length(Number),2]<-Number
    
      sample_order_q<-sample_order_q[order(sample_order_q$Number,decreasing = TRUE),]
      
      sample_order_q<-sample_order_q[which(sample_order_q$Number!=0),]
    
      sample_order_q$Sample<-factor(sample_order_q$Sample,levels=sample_order_q$Sample)
    
      order_sample_w<-as.character(sample_order_q$Sample)
    
      final_frame1$Type<-factor(final_frame1$Type,levels = unique(final_frame1$Type))
    
      final_frame1$Sample<-factor(final_frame1$Sample,levels = order_sample_w)
    
      final_frame1$Number<-as.integer(final_frame1$Number)
    
      all_final_frame<-rbind(final_frame1,all_final_frame)
    
      q=q+1
    
    }
    
    
    #process data which is going to be painted with points and lines
    
    point_line_frame<-data.frame(0)
    
    colnames(point_line_frame)<-c("combinations")
    
    point_line_frame$sample<-c(0)
    
    i=1
    
    row_numbers=1
    
    while (i>0 & i<length(levels(all_final_frame$Sample))+1){
    
      sample_numbers_i<-length(strsplit(as.character(levels(all_final_frame$Sample)[i]),split = "∩")[[1]])
    
      point_line_frame[row_numbers:(row_numbers+sample_numbers_i-1),1]<-as.character(levels(all_final_frame$Sample)[i])
    
      point_line_frame[row_numbers:(row_numbers+sample_numbers_i-1),2]<-strsplit(as.character(levels(all_final_frame$Sample))[i],split = "∩")[[1]]
    
      row_numbers<-row_numbers+sample_numbers_i
    
      i=i+1
    
    }
    
    
    #draw stack plots
    
    key_point<-ggplot(all_final_frame)+
    
      aes(x=Sample,y=Number,fill=Type)+
    
      geom_bar(stat = "identity",position = "stack",width = 0.7)+
    
      labs(x="",width = 1.0)+
      
      labs(y="Mutation number")
    
    if (text_show == "TRUE") {
    
      key_point<-key_point+
   
      geom_text(aes(label = Number), size = 3, colour = 'black',
      
                hjust = .5, position = position_stack(vjust=0.5))
    
    }
    
    key_point <- key_point+
    
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
    
      scale_fill_npg()+
    
      scale_y_continuous(expand = c(0,0))
    
    
    #draw point-line plot
    
    point_line_frame$combinations<-factor(point_line_frame$combinations,levels = levels(all_final_frame$Sample))
    
    point_line_plot<-ggplot(point_line_frame)+
      
      aes(x=combinations,y=sample)
    
    if(length(levels(all_final_frame$Sample))<21){
      
      
      point_line_plot<-point_line_plot+geom_point(size=3.5)
      
    }else{
      
      point_line_plot<-point_line_plot+geom_point(size=2.5)
      
    }
    
    point_line_plot<-point_line_plot+
      
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
    
    
    
    #put pictures together
    
    if(length(levels(all_final_frame$Sample)) < 21){
      
      pdf(paste("upset_analysis.",patientID,".pdf",sep = ""),width = 10,height = 10)
      
      gg<-ggdraw()+draw_plot(key_point,0,0.3,1,0.7)+draw_plot(point_line_plot,0,0,0.9,0.35)
      
    } else {
      
      if (length(levels(all_final_frame$Sample)) > 60) {
      
        pdf(paste("upset_analysis.",patientID,".pdf",sep = ""),width = 13,height = 10)
      
        gg<-ggdraw()+draw_plot(key_point,0,0.3,1,0.7)+draw_plot(point_line_plot,0,0,0.93,0.35)
      
      } else {
        
          pdf(paste("upset_analysis.",patientID,".pdf",sep = ""),width = 10,height = 10)
        
          gg<-ggdraw()+draw_plot(key_point,0,0.3,1,0.7)+draw_plot(point_line_plot,0,0,0.9,0.35)
      
        }
    
    }  
    
    print(gg)
    
    dev.off()
    
  }
  
  upsetanalysis(maf.used_patient_n,patientID = patientID,text_show = paramater2)
  
}