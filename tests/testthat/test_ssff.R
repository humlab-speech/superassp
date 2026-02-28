
library(testthat)
library(superassp)

testFile <- testthat::test_path("..", "signalfiles", "msajc003.wav")



wrassp::trk_forest(testFile,toFile=FALSE) -> inSSFF 

settings <- expand.grid(lag=1,order=1:4)

for(r in 1:nrow(settings)){
  #currLag <- settings[r,"lag"]
  currOrder<- settings[r,"order"]
  test_that(paste0("Test that differentiation order=",currOrder," works"),{
    
    difftrack(inSSFF,order=currOrder,toFile=FALSE) -> outSSFF
    dlab <- paste(rep("d",currOrder),collapse = "")
    
    for(tr in names(inSSFF)){
      otr <- paste0(dlab,tr)
      expect_equal(length(outSSFF[[otr]]),length(inSSFF[[tr]]) )
      
    }

    
    }
  )
}
