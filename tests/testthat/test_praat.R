# NOTE: Uses superassp:::get_extension(), :::get_definedtracks() (internal API)
library(readr)
library(testthat)
library(superassp)

testFile <- testthat::test_path("..", "signalfiles", "msajc003.wav")


praat_funs <- c("trk_formant", "trk_intensity","trk_praatsauce","trk_pitch_cc")



for(f in praat_funs){
  test_that(paste("Confirm that",f,"can generate valid SSFF files"),{
    ssff <- do.call(f,list(testFile,toFile=FALSE))
    
    ext <- superassp:::get_extension(f)
    tracks <- superassp:::get_definedtracks(f)
    
    expect_true(base::setequal(names(ssff),tracks))
    
    tf <- tempfile(fileext = ext)
    wrassp::write.AsspDataObj(ssff,file = tf)
    tfRead <- wrassp::read.AsspDataObj(fname=tf)
    
    expect_true(base::setequal(names(tfRead),tracks))
    
  })
}

