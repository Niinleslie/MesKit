<<<<<<< HEAD
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
=======
#' Run the default MesKit app for analysis locally 
#' 
#' \code{runMesKit} run MesKit locally 
#' @author Mengni Liu
#' 
#' @return a shiny app window
#' @export 
runMesKit<-function(){ 
  shiny::runApp(system.file("shiny", package = "MesKit"))
>>>>>>> 32b9b431882b8ee2872c4b5b86c6bc237b4f3e12
} 