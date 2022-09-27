library(testthat)
library(superassp)

signalfiles <- c("../signalfiles/generated/vowel14s.wav","../signalfiles/generated/vowel14s_stereo.wav")

slicefunctions <- c("praat_voice_report","praat_voice_tremor")


for(fun in slicefunctions){
  
  test_that(paste0("Check that slice generating function",fun," can be applied to the a (portion of) a long mono file"),{
  
    sfile <-  signalfiles[1]
      
    out <- do.call(fun,list(listOfFiles=sfile))
    expect_type(out,"list")
    #just a portion
    out <- do.call(fun,list(listOfFiles=sfile,beginTime=1,endTime=2))
    expect_type(out,"list")
    
  })
  
  test_that(paste0("Check that slice generating function",fun," can be applied to a (portion of) a long stereo file"),{
    
    sfile <-  signalfiles[2]
    
    out <- do.call(fun,list(listOfFiles=sfile))
    expect_type(out,"list")
    #just a portion
    out <- do.call(fun,list(listOfFiles=sfile,beginTime=1.0,endTime=2.0))
    expect_type(out,"list")
  })
  
}