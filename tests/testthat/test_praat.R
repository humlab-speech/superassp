context("Testing praat_formant_burg")
library(readr)
library(testthat)

testFile <- file.path("..","signalfiles","msajc003.wav")

test_that("Praat_formant_burg can generate signal files",{
  ssff <- praat_formant_burg(testFile,toFile=FALSE)

  expect_true(base::setequal(names(ssff),c("fm","bw")))
  
  
  
  
})
