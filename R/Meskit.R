#' Run the default MesKit app for analysis locally 
#' 
#' \code{runMesKit} run MesKit locally 
#' @author Mengni Liu
#' @export 
runMesKit<-function(){ 
  shiny::runApp(system.file("shiny", package = "MesKit"))
} 