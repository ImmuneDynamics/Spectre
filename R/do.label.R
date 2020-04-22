#' do.label
#' 
#' @export

do.label <- function(dat,
                     label.name,
                     
                     marker.name,
                     marker.threshold,
                     marker.value){
  
  message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")
  
  if(marker.threshold == ">"){
    dat[[label.name]] <- ifelse(dat[[marker.name]] > marker.value, "1", "0")
  }
  
  if(marker.threshold == ">="){
    dat[[label.name]] <- ifelse(dat[[marker.name]] >= marker.value, "1", "0")
  }
  
  if(marker.threshold == "=="){
    dat[[label.name]] <- ifelse(dat[[marker.name]] == marker.value, "1", "0")
  }
  
  if(marker.threshold == "<="){
    dat[[label.name]] <- ifelse(dat[[marker.name]] <= marker.value, "1", "0")
  }
  
  if(marker.threshold == "<"){
    dat[[label.name]] <- ifelse(dat[[marker.name]] < marker.value, "1", "0")
  }
  
  return(dat)
}