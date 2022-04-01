context("Testing that wrassp signal processing functions work in superassp")
library(readr)
library(testthat)

nonSSFFFunctions <-c("rfcana","afdiff","affilter")


testFile <- file.path("..","signalfiles","msajc003.wav")

wrassp_funs <- setdiff(names(wrassp::wrasspOutputInfos),nonSSFFFunctions)

for(f in wrassp_funs){
  test_that(paste("Confirm that",f,"can becalled to generate valid SSFF file in superassp"),{
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

