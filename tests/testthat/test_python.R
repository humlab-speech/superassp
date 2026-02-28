# NOTE: Uses superassp:::get_extension(), :::get_definedtracks() (internal API)
library(readr)
library(testthat)
library(superassp)
library(reticulate)

testFile <- testthat::test_path("..", "signalfiles", "msajc003.wav")
# crepe not working


python_funs <- c("trk_pyin","trk_seenc","trk_yin","trk_swipe","trk_reaper","trk_rapt","reaper_pm","trk_harvest","trk_dio")
# 137:"trk_crepe","trk_kaldi_pitch",


for(f in python_funs){
  test_that(paste("Confirm that",f,"can generate valid SSFF files"),{
    ssff <- do.call(f,list(listOfFiles=testFile,toFile=FALSE))
    
    ext <- superassp:::get_extension(f)
    tracks <- superassp:::get_definedtracks(f)
    
    expect_true(base::setequal(names(ssff),tracks))
    
    tf <- tempfile(fileext = ext)
    wrassp::write.AsspDataObj(ssff,file = tf)
    tfRead <- wrassp::read.AsspDataObj(fname=tf)
    #cat(paste(names(tfRead),collapse = "-"))
    #cat(paste(tracks,collapse = "-"))
    expect_true(base::setequal(names(tfRead),tracks))
    Sys.sleep(0.5)
  })
}
