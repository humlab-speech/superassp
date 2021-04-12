
vat <- function(listOfFiles,beginTime=0,endTime=0){
  vatpath <- "/Users/frkkan96//Documents/MATLAB/VoiceAnalysisToolbox/"
  
  inFile <- tempfile(fileext = ".wav")
  outFile <- gsub("wav$","mat",inFile)
  file.copy(listOfFiles[1],inFile)
  
  scriptline1 <- paste0("cd ",vatpath,";soundfile = '",inFile,"'")
  scriptline2 <- paste0("[measures_vector, measures_names, f0] = local_voice_analysis(soundfile,",beginTime,",",endTime,")")
  outscript <- paste0("save('",outFile,"','measures_vector','measures_names','f0')")
  
  matlabr::run_matlab_code(c(scriptline1,scriptline2,outscript),endlines = TRUE,paths_to_add = vatpath,verbose = FALSE,figure_windows=FALSE)
  
  if(!file.exists(outFile)){
    stop("There was ")
  }
  matDat <- R.matlab::readMat(outFile)
  out <- data.frame(measures_names=as.vector(unlist(matDat$measures.names)),
                    values=t(matDat$measures.vector)) %>%
    tidyr::pivot_wider(names_from = "measures_names",values_from="values")
  
  return(out)
}