#load packages
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(pathview)
library(topGO)

#options(stringsAsFactors = FALSE)
#options(bitmapType='cairo')
options(warn = -1)


# Clean the working environment
rm(list=ls())


#' USE R code to do the pathway analysis.
#'
#' @param gene_data a dataframe transferred from a maf file,which contains information going to be analyzed
#' @param patientID the ID of a patient whose cancer information is going to be processed
#' @param analysis_type the analysis type, one of "KEGG","Reactome", and "Both"
#'
#' @return
#' @export a .xls file contains all the information of the result and two pictures showing the cluster result of KEGG analysis; 
#'         a .xls file contains all the information of the reactome result and two pictures showing the cluster result of it;
#' @examples
#' Pathway_analysis(maf,311252,"Both")


#Pathway analysis

Pathway_analysis<-function( gene_data , patientID , analysis_type = "Both"){
  
  #name the patients whit their IDs
  
  for (i in 1:length(gene_data$Hugo_Symbol)){
    
    gene_data$Patient[i]<-strsplit(as.character(gene_data$Tumor_Sample_Barcode[i]),split = "-")[[1]][1]
    
  }
  
  gene_data$Patient<-as.factor(gene_data$Patient)
  
  
  
  #get the information of certain patient
  
  df<-gene_data[which(as.character(gene_data$Patient)==patientID),]
  
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
  
  
  
  # KEGG analysis
  
  KEGG_analysis<-function(){
    
    kegg <- enrichKEGG(gene       = trans$ENTREZID,
                       
                       organism     = 'hsa',
                       
                       pvalueCutoff = 0.05)
    
    
    if (nrow(kegg@result) > 0){
      
      str_length = max(nchar(kegg@result$Description))
      
      str_height = nrow(kegg)
      
      write.table(kegg@result, file =  paste("KEGG_enrich_",patientID,".xls",sep = ""),
                  
                  sep = "\t", quote = FALSE, row.names = T, col.names = NA)
      
      if (str_height > 15){
        
        pdf(paste("KEGG_",patientID,".pdf",sep = ""), width = 3+(str_length)/10, height = str_height/3)
        
      }else{
        
        pdf(paste("KEGG_",patientID,".pdf",sep = ""), width = 3+(str_length)/10, height = 5)
        
      }
      
      print(barplot(kegg, showCategory = nrow(kegg)))
      
      print(dotplot(kegg, showCategory = nrow(kegg)))
      
      dev.off()
      
      
      for (i in 1:nrow(kegg@result)){
        
        #browseKEGG(kegg, 'hsa04392')
        
        pathway.id = kegg@result$ID[i]
        
        pathview(gene.data  = data_genes,
                 
                 pathway.id = pathway.id,
                 
                 gene.idtype ="KEGG",
                 
                 out.suffix = "enrich",
                 
                 species    = "hsa",
                 
                 kegg.native = TRUE)
        
        
        
        pathview(gene.data   = data_genes,
                 
                 pathway.id  = pathway.id,
                 
                 gene.idtype = "KEGG",
                 
                 out.suffix  = "enrich",
                 
                 species     = "hsa",
                 
                 kegg.native = F,
                 
                 same.layer  = F)
        
        system(paste("mv ", pathway.id, "* ",  sep = ""))
        
      }
      
    }
    
  }
  
  
  # Reactome analysis
  
  Reactome_analysis<-function(){
  
    reactome_ana <- enrichPathway(gene         = trans$ENTREZID,
                                  
                                  organism     = 'human',
                                  
                                  pvalueCutoff = 0.05)
    
    
    write.table(reactome_ana@result, file = paste("Reactome_",patientID,"_enrich.xls"),
                
                sep = "\t", quote = FALSE, row.names = T, col.names = NA)
    
    if (nrow(reactome_ana@result) > 0){
      
      str_length = max(nchar(reactome_ana@result$Description))
      
      str_height = nrow(reactome_ana)
      
      if (str_height > 15){
        
        pdf(paste("Reactome_enrich_",patientID,".pdf"), width = 3+(str_length)/10, height = str_height/3)
        
      }else{
        
        pdf(paste("Reactome_enrich_",patientID,".pdf"), width = 3+(str_length)/10, height = 5)
        
      }
      
      print(barplot(reactome_ana, showCategory = nrow(reactome_ana)))
      
      print(dotplot(reactome_ana, showCategory = nrow(reactome_ana)))
      
      dev.off()
      
      #pdf("Reactome_enrich_plot2.pdf", width = 25, height = 25)
      
      #enrichMap(reactome_ana)
      
      #dev.off()
      
    }

  }
  
  if(as.character(analysis_type)=="Both"){
    
    KEGG_analysis()
    
    Reactome_analysis()
    
  }else{
    
    if(as.character(analysis_type)=="KEGG"){
      
      KEGG_analysis()
      
    }else{
      
      Reactome_analysis()
      
    }
    
  }
  
}
