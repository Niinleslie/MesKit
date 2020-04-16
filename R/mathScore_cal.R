## MATH Caculation
calMATH <- function(VAF){
   VAF = VAF[!is.na(VAF)] 
   
   MAD <- 1.4826*median(abs(VAF - median(VAF)))
   MATH <- 100 * MAD / median(VAF)
   return(round(MATH, digits=3))
}
