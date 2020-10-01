#' Run the default MesKit app for analysis locally 
#' 
#' \code{runMesKit} run MesKit locally 
#' @author Mengni Liu
#' @examples
#' \donttest{
#' runMesKit()
#' }
#' @return a shiny app window
#' @export 
runMesKit<-function(){ 
  shiny::runApp(system.file("shiny", package = "MesKit"))
} 