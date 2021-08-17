
library(testthat)
library(superassp)

testFile <- file.path("..","signalfiles","msajc003.wav")



wrassp::forest(testFile,toFile=FALSE) -> inSSFF 

settings <- expand.grid(lag=1,order=1:4)

for(r in 1:nrow(settings)){
  #currLag <- settings[r,"lag"]
  currOrder<- settings[r,"order"]
  test_that(paste0("Test that differentiation order=",currOrder," works"),{
    
    difftrack(inSSFF,order=currOrder,toFile=FALSE) -> outSSFF
    
    expect_equal(length(outSSFF$bw),length(inSSFF$bw) )
    
    }
  )
}
