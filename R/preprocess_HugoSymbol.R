preprocess_HugoSymbol <- function(mafData){
   
   remove.all_none.dat <- mafData %>% 
      dplyr::rowwise() %>% 
      ## remove NONE
      dplyr::filter(!all(strsplit(Hugo_Symbol,",|;")[[1]] == "NONE")) %>% 
      as.data.table()
      
   h <- lapply(remove.all_none.dat$Hugo_Symbol,function(x){
      s <- strsplit(x,",|;")[[1]]
      s <- s[s!= "NONE"]
      if(length(s) == 1){
         return(s)
      }else if(all(stringr::str_detect(s,"^LOC"))){
         s <- s[1]
      }else if(any(stringr::str_detect(s,"^LOC"))){
         s <- s[!stringr::str_detect(s,"^LOC")][1]
      }else{
         s <- s[1]
      }
      return(s)
   })
   
   remove.all_none.dat$Hugo_Symbol <- unlist(h)
   mafData <- as.data.frame(remove.all_none.dat) 
   
   return(mafData)
}

# mafData <- data.table::fread(
#    file = mafFile,
#    quote = "",
#    header = TRUE,
#    data.table = TRUE,
#    fill = TRUE,
#    sep = '\t',
#    skip = "Hugo_Symbol",
#    stringsAsFactors = FALSE
# )
# sort.maf <- sortHugo_Symbol(mafData)
# write.table(sort.maf, file = "HCC6046.maf", sep = "\t",quote = F,row.names = F)
