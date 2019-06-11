#load packages
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(pathview)
library(topGO)

options(stringsAsFactors = FALSE)
options(bitmapType='cairo')
options(warn = -1)


# Clean the working environment
rm(list=ls())


#' Use R code to do the GO analysis.
#'
#' @param gene_data a dataframe transferred from a maf file,which contains information going to be analyzed
#' @param patientID the ID of a patient whose cancer information is going to be analyzed
#' @param type the analyze type of GO,one of "BP","CC" and "MF". Default type="ALL“
#' @param pval cutoff value of pvalue. Default pval=0.01
#' @param qval qvalue cutoff. Default qval=0.05
#'
#' @return
#' @export a .xls file contains and two pictures of the result of GO analysis, may also get a cnetplot.
#'
#' @examples 
#' GO_analysis(maf,311252,BP,0.01,0.05)
#' 

GO_analysis <- function(gene_data , patientID , type = "ALL", pval = 0.01, qval = 0.05){
  
  #name the patients whit their IDs
  
  for (i in 1:length(gene_data$Hugo_Symbol)){
    
    gene_data$Patient[i]<-strsplit(as.character(gene_data$Tumor_Sample_Barcode[i]),split = "-")[[1]][1]
    
  }
  
  gene_data$Patient<-as.factor(gene_data$Patient)
  
  
  #get the information of certain patient
  
  df<-gene_data[which(gene_data$Patient==as.character(patientID)),]
  
  primitive_length<-length(df$Hugo_Symbol)
  
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
  
  new_gene_number<-length(df$Hugo_Symbol)- primitive_length
  
  for (n in 1:new_gene_number){
    
    df[primitive_length+n,1]<-new_gene[n]
    
  }
  
  data_genes = unique(df[,1])
  
  trans = bitr(data_genes, fromType="SYMBOL", toType=c("ENTREZID", "GENENAME"), OrgDb="org.Hs.eg.db")
  
  
  # GO analysis
  
  GO.analysis = function(genes = data_genes, type = type, pval = pval, qval = qval){
    
    ego <- enrichGO(gene          = data_genes,
                    
                    OrgDb         = org.Hs.eg.db,
                    
                    keyType       = 'SYMBOL',
                    
                    ont           = type,
                    
                    pAdjustMethod = 'BH',
                    
                    pvalueCutoff  = pval,
                    
                    qvalueCutoff  = qval)
    
    write.table(ego@result, file = paste("GO_", patientID,"_",type ,"_enrich.xls" ,sep = ""),
                
                sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    
    if (!is.null(ego) && nrow(ego@result) > 0){
      
      str_length = max(na.omit(nchar(ego@result$Description[1:nrow(ego)])))
      
      str_height = nrow(ego)
      
      if (str_height >15){
        
        pdf(paste( patientID, type, sep = ""), width = 3+(str_length)/10, height = str_height/3)
        
      }else{
        
        pdf(paste( patientID, type, sep = ""), width = 3+(str_length)/10, height = 5)
        
      }
      
      print(barplot(ego, showCategory = nrow(ego)))
      
      print(dotplot(ego, showCategory = nrow(ego)))
      
      dev.off()
      
      if(nrow(ego)>0){
      
        pdf(paste(patientID, type, sep = ""), width = 15, height = 15)
        
        #enrichMap(ego)
        
        print(cnetplot(ego, categorySize="pvalue",showCategory = nrow(ego)))
    
        dev.off()
      }
      
    }
    
  }
  
  GO.analysis(genes = data_genes,type = type,pval = pval,qval = qval)
  
}
