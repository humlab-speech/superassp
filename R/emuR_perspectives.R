
enable_formantOverlay <- function(emuDBhandle,perspective){
  perspectiveNames <- list_perspectives(emuDBhandle)$name
  trackNames <- list_ssffTrackDefinitions(emuDBhandle)$name
  
  #Stop processing if the perspective is not defined in the database
  if(! perspective %in% perspectiveNames) {stop("The perspective  ",perspective," is not defined in the database ", emuDBhandle$dbName,"!")}
  
  #Stop processing if a track FORMANTS is not defined in the database
  if(! "FORMANTS" %in% trackNames) {stop("In order to enable formant overlays, a track named 'FORMANTS' must be defined in the database !")}
  
  which(grepl(perspective,perspectiveNames)) -> perspid
  dbConfig = emuR:::load_DBconfig(ae)
  
  dbConfig$EMUwebAppConfig$perspectives[[perspid]]$signalCanvases$assign[[1]] <- list("signalCanvasName"="SPEC","ssffTrackName"="FORMANTS")
  res <- emuR:::store_DBconfig(emuDBhandle,dbConfig = dbConfig)
  return(res)
}

set_specOverlay <- function(emuDBhandle,perspective,trackname){
  perspectiveNames <- list_perspectives(emuDBhandle)$name
  trackNames <- list_ssffTrackDefinitions(emuDBhandle)$name
  
  #Stop processing if the perspective is not defined in the database
  if(! perspective %in% perspectiveNames) {stop("The perspective  ",perspective," is not defined in the database ", emuDBhandle$dbName,"!")}
  
  #Stop processing if a track FORMANTS is not defined in the database
  if(! trackname %in% trackNames) {stop("In order to enable an overlay on the spectrogram, the track must already be defined in the database.")}
  
  which(grepl(perspective,perspectiveNames)) -> perspid
  dbConfig = emuR:::load_DBconfig(ae)
  
  dbConfig$EMUwebAppConfig$perspectives[[perspid]]$signalCanvases$assign[[1]] <- list("signalCanvasName"="SPEC","ssffTrackName"=trackname)
  res <- emuR:::store_DBconfig(emuDBhandle,dbConfig = dbConfig)
  return(res)
}


set_overlayTrack <- function(emuDBhandle,perspective,trackname, overlay.on="SPEC",overwrite=FALSE){
  
  #Ensure that the overlay place is valid
  overlay.on <- toupper(overlay.on)
  if(! overlay.on  %in% c("SPEC","OSCI")) {stop("You can only specify an overlay on the SPEC and OSCI displays!")}
  
  perspectiveNames <- list_perspectives(emuDBhandle)$name
  trackNames <- list_ssffTrackDefinitions(emuDBhandle)$name
  
  #Stop processing if the perspective is not defined in the database
  if(! perspective %in% perspectiveNames) {stop("The perspective  ",perspective," is not defined in the database ", emuDBhandle$dbName,"!")}
  
  #Stop processing if the track is not defined in the database
  if(! trackname %in% trackNames) {stop("The track  ",trackname," is not defined in the database ", emuDBhandle$dbName,"!")}
  
  which(grepl(perspective,perspectiveNames)) -> perspid
  dbConfig = emuR:::load_DBconfig(ae)
  
  #Check and stop processing if an overlay is alrady set 
  overlay <- dbConfig$EMUwebAppConfig$perspectives[[perspid]]$signalCanvases$assign
  
  if(length(overlay) > 0 ){
    for(ov in 1:length(overlay)){
      if(overlay[[ov]]$signalCanvasName == overlay.on ){
        
        if(! overwrite) {stop("Cannot set an overlay on ", overlay.on, " as one is already defined in the database.\nPlease set overwrite=TRUE if you wish to overwrite the previous setting.")}
        
        dbConfig$EMUwebAppConfig$perspectives[[perspid]]$signalCanvases$assign[[ov]] <- list("signalCanvasName"=overlay.on,"ssffTrackName"= trackname)
      }else{
        #In this case, existing overlay settings do not exist
        dbConfig$EMUwebAppConfig$perspectives[[perspid]]$signalCanvases$assign[[ov+1]] <- list("signalCanvasName"=overlay.on,"ssffTrackName"= trackname)
      }
    }
  }
  
  res <- emuR:::store_DBconfig(emuDBhandle,dbConfig = dbConfig)
  return(res)
}
