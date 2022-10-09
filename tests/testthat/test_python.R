context("Testing Python signal processing functions")
library(readr)
library(testthat)
library(superassp)
library(reticulate)

testFile <- file.path("..","signalfiles","msajc003.wav")

python_funs <- c("swipe","reaper","rapt","reaper_pm","yin","pyin","crepe","harvest","dio","yaapt")
#python_funs <- "reaper_pm"

for(f in python_funs){
  test_that(paste("Confirm that",f,"can generate valid SSFF files"),{
    ssff <- do.call(f,list(listOfFiles=testFile,toFile=FALSE))
    
    ext <- superassp::get_extension(f)
    tracks <- superassp::get_definedtracks(f)
    
    expect_true(base::setequal(names(ssff),tracks))
    
    tf <- tempfile(fileext = ext)
    wrassp::write.AsspDataObj(ssff,file = tf)
    tfRead <- wrassp::read.AsspDataObj(fname=tf)
    #cat(paste(names(tfRead),collapse = "-"))
    #cat(paste(tracks,collapse = "-"))
    expect_true(base::setequal(names(tfRead),tracks))
    
  })
}
