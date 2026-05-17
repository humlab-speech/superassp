# NOTE: Uses superassp:::get_extension(), :::get_definedtracks(), :::knownLossless() (internal API)
library(readr)
library(testthat)


# nonSSFFFunctions <-c("rfcana","afdiff","affilter")
# 
# 
# testFile <- file.path("..","signalfiles","msajc003.wav")
# 
# wrassp_funs <- setdiff(names(wrassp::wrasspOutputInfos),nonSSFFFunctions)
# 
# for(f in wrassp_funs){
#   test_that(paste("Confirm that",f,"can becalled to generate valid SSFF file in superassp"),{
#     ssff <- do.call(f,list(testFile,toFile=FALSE))
#     
#     ext <- superassp:::get_extension(f)
#     tracks <- superassp:::get_definedtracks(f)
#     
#     expect_true(base::setequal(names(ssff),tracks))
#     
#     tf <- tempfile(fileext = ext)
#     wrassp::write.AsspDataObj(ssff,file = tf)
#     tfRead <- wrassp::read.AsspDataObj(fname=tf)
#     
#     expect_true(base::setequal(names(tfRead),tracks))
#     
#   })
# }



wrassp_funs <- c("trk_acf","trk_rms","trk_pitch_ksv","trk_pitch_mhs","trk_formant_forest","trk_zcr","trk_lpc","trk_lar","trk_arf","trk_rfc","trk_cepstrum","trk_css_spectrum","trk_lps_spectrum")

#wrassp_funs <- c("trk_lps_spectrum")

knownLossless <- superassp:::knownLossless()

testFiles <- normalizePath(list.files(
  system.file("samples", "sustained", package = "superassp"),
  full.names = TRUE
))

for(f in wrassp_funs){
  for(testFile in testFiles){
    ext <- tools::file_ext(testFile)
    test_that(paste("Confirm that '",f,"' can handle files with extension '",ext,"'", sep=""),{
      if( ! ext %in% knownLossless ){
        expect_warning(ssff <- do.call(f,list(listOfFiles=testFile,toFile=FALSE)))
      }else{
        # C-level primary path handles native formats silently; av path handles
        # modern formats (flac, aiff, …) — no specific conversion message required.
        ssff <- do.call(f,list(testFile,toFile=FALSE))
      }

      
      ext <- superassp:::get_extension(f)
      tracks <- superassp:::get_definedtracks(f)
      
      expect_true(base::setequal(names(ssff),tracks))
      
      #This checks that the generated signal track is valid / can be read in
      tf <- tempfile(fileext = ext)
      wrassp::write.AsspDataObj(ssff,file = tf)
      tfRead <- wrassp::read.AsspDataObj(fname=tf)
      
      expect_true(base::setequal(names(tfRead),tracks))
      unlink(tf,force = FALSE)
      
    })
  }

}