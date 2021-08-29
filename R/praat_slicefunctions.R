#' Compute the Acoustic Voice Quality Index (AVQI)
#'
#' This function computes the Acoustic Voice Quality Index (AVQI) from a set of
#' continuous speech and sustained vowel samples. Praat is used to compute the
#' AVQI and the result is therefore identical to the code published in
#' \insertCite{Latoszek.2019.10.1002/lary.27350}{superassp}. The user may
#' provide multiple continuous speech and sustained vowel samples for the
#' speaker, and these will then be concatenates together before computing the
#' AVQI value for the speaker and in that recording session.
#'
#' If the user provides an `pdf.path`, a PDF of the AVQI analysis as implemented
#' in the Praat script of
#' \insertCite{Latoszek.2019.10.1002/lary.27350}{superassp} will be placed
#' there. The name of the PDF will be `speaker.ID`_`session.datetime`, or "NA"
#' if not provided. The default behavior is to not overwrite existing PDF files,
#' as this could result in loss of data. The user should make sure that an
#' appropriate `speaker.ID` and `session.datetime` values are provided to mark
#' the output appropriately.
#'
#' @param svDF A data.frame containing columns "absolute_file_path","start", and
#'   "end". Each row should contain the full path to a sound file, and "start"
#'   and "end" points of a sustained vowel in that sound file. The multiple
#'   sustained vowels will be concatenated together to compute the AVQI, and it
#'   is therefore important that the user that the start and end points are
#'   indeed inside of the sustained vowel as to not influence the results by
#'   introducing unvoiced frames.
#' @param csDF A data.frame containing columns "absolute_file_path","start", and
#'   "end". Each row should contain the full path to a sound file, and "start"
#'   and "end" points of a portion of continuous speech in that sound file.
#'   Multiple portions of speech will be concatenated together to compute the
#'   AVQI, and it is therefore important that the user that the start and end
#'   points are indeed inside of the portion of produced speech as to not
#'   influence the results by introducing unvoiced frames.
#' @param speaker.name The name of the speaker. Only used for in produced PDF
#'   output files.
#' @param speaker.ID This will be used to identify the output returned list and
#'   in PDF output, and could therefore be either the ID of a speaker or of a
#'   Speaker + session compilation.
#' @param speaker.dob The date of birth of the speaker. If provided, the PDF
#'   will be marked with this information.
#' @param session.datetime The date and time when the recording was performed
#'   (as a string). If provided, the PDF will be marked with this information.
#' @param pdf.path This is where PDF output files will be stored. If not
#'   provided, no PDF files will be produced.
#' @param simple.output The AVQI Praat function can produce either a full page
#'   report of the voice profile, or a much more condensed version filling just
#'   a portion of the page. If `simple.output=TRUE`, the simplified version will
#'   be produced.
#' @param overwrite.pdfs Should existing PDF files be overwritten in the PDF
#'   output directory? Defaults to a safe behavior where older PDFs are not
#'   overwritten.
#' @param praat_path An explicit path to the Praat binary.
#' 
#' @return A list with the following fields: 
#' \describe{
#' \item{ID}{The speaker / speaker + session identifier of the output}
#' \item{CPPS}{Smoothed Cepstral Peak Prominence value}
#' \item{HNR}{An Harmonic-to-noise estimate}
#' \item{Shim_local}{A (local) Shimmer measurement (in \%)}
#' \item{Shim_local_DB}{A (local) Shimmer measurement, in decibels}
#' \item{LTAS_Slope}{The slope of the Long Time Average Spectrum (in dB)}
#' \item{LTAS_Tilt}{The Long Time Average Spectrum tilt (in dB)}
#' \item{AVQI}{Acoustic Voice Quality Index summarizing the measures above}
#' }
#' 
#' @references
#'  \insertAllCited{}

praat_avqi <- function(svDF,
                       csDF,
                       speaker.name=NULL,
                       speaker.ID=NULL,
                       speaker.dob=NULL,
                       session.datetime=NULL,
                       pdf.path=NULL,
                       simple.output=FALSE,
                       overwrite.pdfs=FALSE,
                       praat_path=NULL){
  
  
  requiredDFColumns <- c("absolute_file_path","start","end")
  
  listOfFiles <- unique(
    c(svDF$absolute_file_path,
      csDF$absolute_file_path)
  )
  
  if(! have_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
  }
  
  if(! setequal(listOfFiles,unique(svDF$absolute_file_path)) | ! setequal(listOfFiles,unique(csDF$absolute_file_path))){
    stop("The 'svDF' and 'csDF' may contain only times for files in the 'listOfFiles', and similarly must _both_ include at least one row for each file name in 'listOfFiles'.")
  }
  # 
  if(! requiredDFColumns %in% names(svDF) || ! requiredDFColumns %in% names(csDF)  ){
    stop("The 'svDF' and 'csDF' structures must both contain columns named ",paste(requiredDFColumns,collapse=",",sep=""),".")
  }
  
  praat_script <- ifelse(dir.exists("inst"), ## This means that we are developing
                         file.path("inst","praat","AVQI301.praat"),
                         file.path(system.file(package = "superassp",mustWork = TRUE),"praat","AVQI301.praat"))
  
  #return(praat_script)
  #praat_script <- "/Users/frkkan96/Documents/src/superassp/inst/praat/AVQI203.praat"
  avqi <- tjm.praat::wrap_praat_script(praat_location = get_praat(),
                                       script_code_to_run = readLines(praat_script)
                                       ,return="last-argument")
  
  #Check that all files exists before we begin
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)){
    filedNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ",paste(filedNotExists, collapse = ", "))
  }
  
  # Set up a (CLEAN) directory for interchange
  avqiDir <- file.path(tempdir(check=TRUE),"avqtemp")
  unlink(avqiDir,recursive = TRUE,force=FALSE,expand=FALSE)
  dir.create(avqiDir)
  
  
  #The empty vector of file names that should be returned
  outListOfFiles <- c()
  
  
  #Copy Sustained Vowel portions from the file
  
  #Pre-generate names of output files
  svDF$OutFileName <- file.path(avqiDir,paste0("sv",1:nrow(svDF),".wav"))
  
  for(r in 1:nrow(svDF)){
    #Here we finally copy out all signal file content into separate files
    currSound <- wrassp::read.AsspDataObj(fname=svDF[r,"absolute_file_path"],begin=svDF[r,"start"],end=svDF[r,"end"])
    wrassp::write.AsspDataObj(currSound,file=svDF[r,"OutFileName"])
  }
  
  #Copy Continous Speech portions from the file
  
  #Pre-generate names of output files
  csDF$OutFileName <- file.path(avqiDir,paste0("cs",1:nrow(csDF),".wav"))
  
  for(r in 1:nrow(csDF)){
    #Here we finally copy out all signal file content into separate files
    currSound <- wrassp::read.AsspDataObj(fname=csDF[r,"absolute_file_path"],begin=csDF[r,"start"],end=csDF[r,"end"])
    wrassp::write.AsspDataObj(currSound,file=csDF[r,"OutFileName"])
  }
  
  # Now we are all set up to run the Praat script
  # boolean Illustrated_version 1
  # sentence name_patient Fredrik Karlsson
  # sentence Date_of_birth 1975-12-31
  # sentence Assessment_date 2021-12-31
  # sentence Input_directory /Users/frkkan96/Documents/src/superassp/tests/signalfiles/AVQI/input
  # boolean Generate_PDF_files 1
  # sentence Speaker_ID 1
  # sentence Output_directory /Users/frkkan96/Documents/src/superassp/tests/signalfiles/AVQI/output
  # sentence Output_file /Users/frkkan96/Documents/src/superassp/tests/signalfiles/AVQI/output/avqi.csv
  
  outAVQITabFile <- avqi(ifelse(simple.output,1,0),
                         ifelse(! is.null(speaker.name),speaker.name,"NA"),
                         ifelse(! is.null(speaker.dob),speaker.dob,"NA"),
                         ifelse(! is.null(session.datetime),session.datetime,"NA"),
                         avqiDir,
                         ifelse(! is.null(pdf.path),1,0),
                         ifelse(! is.null(speaker.ID),speaker.ID,"NA"),
                         avqiDir,
                         file.path(avqiDir,"avqi.csv")
  )
  
  inTable <- read.csv(file=outAVQITabFile
                      ,header=TRUE
                      ,na.strings =c("--undefined--","NA"),
                      sep = ",")
  assertthat::are_equal(nrow(inTable),1)
  
  # Now copy PDF files to the pdf.output path
  if(!is.null(pdf.path)){
    pdfFiles <- list.files(avqiDir,pattern=".*[.]pdf",full.names=TRUE)
    for(currFile in pdfFiles){
      file.copy(from=currFile,to=pdf.path,overwrite=overwrite.pdfs)
    }
  }
  return(as.list(inTable))  
}

# FOR INTERACTIVE TESTING
# df <- data.frame("absolute_file_path"="~/Desktop/kaa_yw_pb.wav","start"=0,"end"=3)
# df2 <- data.frame("absolute_file_path"=c("~/Desktop/kaa_yw_pb.wav","~/Desktop/kaa_yw_pb.wav"),"start"=c(0,0),"end"=c(3,3))
# 
# praat_avqi(df, df,pdf.path="~/Desktop/test/",simple.output = TRUE) -> outTab