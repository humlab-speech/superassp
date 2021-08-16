context("Testing praat_formant_burg")
library(readr)
library(testthat)

testFile <- file.path("..","signalfiles","msajc003.wav")

praat_funs <- c("praat_formant_burg", "praat_intensity","praat_moments","praat_sauce")

for(f in praat_funs){
  test_that(paste("Confirm that",f,"can generate valid SSFF files"),{
    ssff <- do.call(f,list(testFile,toFile=FALSE))
    
    ext <- superassp::get_extension(f)
    tracks <- superassp::get_definedtracks(f)
    
    expect_true(base::setequal(names(ssff),tracks))
    
    tf <- tempfile(fileext = ext)
    wrassp::write.AsspDataObj(ssff,file = tf)
    tfRead <- wrassp::read.AsspDataObj(fname=tf)
    
    expect_true(base::setequal(names(tfRead),tracks))
    
  })
}

