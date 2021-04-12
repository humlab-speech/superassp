

#' 
#' 
#' library(reticulate)
#' 
#' use_condaenv("opensauce_env",required = TRUE)
#' spw <- import("scipy.io.wavfile")
#' pr <- import("pyreaper")
#' np <- import("numpy", convert = TRUE)
#' 
#' w <- spw$read("/Users/frkkan96/Desktop/aaa.wav")
#' wd <- w[[2]]
#' wn <- r_to_py(wd, convert=TRUE)
#' out <- pr$reaper(wn,w[[1]])
#' 
#' wrassp::read.AsspDataObj("/Users/frkkan96/Desktop/aaa.wav") -> r
#' np_array(as.vector(unlist(r$audio)),dtype = "int16_t") -> ad
#' sr <- attr(r,"sampleRate")
#' out <- pr$reaper(r_to_py(ad),5000)