context("Testing praat_formant_burg")
library(readr)


path2demoData = file.path(tempdir(),"emuR_demoData")
unlink(path2demoData, recursive = TRUE)

emuR::create_emuRdemoData()

ae <- emuR::load_emuDB(file.path(path2demoData,"ae_emuDB"))


test_that("Formant tracks computed by praat can be loaded by emuR::serve()",{
  
  emuR::list_files(ae,"wav") -> wavs
  superassp:::praat_formant_burg(wavs$absolute_file_path,toFile = TRUE)
  
  
  list_files(ae,"fms") -> formants
  filesExistsBoolean <- file.exists(formants$absolute_file_path)
  ssffLst <- c()
  #Check that the generated files are still valid SSFF files
  for(ssff_file in formants$absolute_file_path){
    wrassp::read.AsspDataObj(ssff_file) -> currSSFF
    ssffLst <- c(ssffLst,is.AsspDataObj(currSSFF))
  }
  #Now, the state of existing and being a valid SSFF file should
  # be the same.
  cat(ssffLst)
  cat(filesExistsBoolean)
  expect_equal(ssffLst,filesExistsBoolean)
}
          )







#Cleanup
#DBI::dbDisconnect(ae$connection)
#rm(ae)
#unlink(path2demoData,recursive = TRUE)