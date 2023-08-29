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
#' @param min.sv The minimal continuous vowel duration required to make accurate measurements (in milliseconds). 
#' If the total duration of sustained vowels in the `svDF` tibble is smaller than this, the function will 
#'  produce an error and quit processing. Defaults to 1000 ms (1 second) and should rarely be shorter than that. 
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
#' \item{Shim_local}{A (local) Shimmer measurement (in %)}
#' \item{Shim_local_DB}{A (local) Shimmer measurement, in decibels}
#' \item{LTAS_Slope}{The slope of the Long Time Average Spectrum (in dB)}
#' \item{LTAS_Tilt}{The Long Time Average Spectrum tilt (in dB)}
#' \item{AVQI}{Acoustic Voice Quality Index summarizing the measures above}
#' }
#' 
#' @export
#' 
#' @references
#'  \insertAllCited{}

praat_avqi <- function(svDF,
                       csDF,
                       min.sv=1000,
                       speaker.name=NULL,
                       speaker.ID=speaker.name,
                       speaker.dob=NULL,
                       session.datetime=NULL,
                       pdf.path=NULL,
                       simple.output=FALSE,
                       overwrite.pdfs=FALSE,
                       praat_path=NULL){
  
  
  requiredDFColumns <- c("listOfFiles","start","end")
  
  listOfFiles <- unique(
    c(svDF$listOfFiles,
      csDF$listOfFiles)
  )
  
  if(! have_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
  }
  
  if(! all(requiredDFColumns %in% names(svDF))|| ! all(requiredDFColumns %in% names(csDF))  ){
    stop("The 'svDF' and 'csDF' structures must both contain columns named ",paste(requiredDFColumns,collapse=",",sep=""),".")
  }

  totalSVdur <- sum(svDF$end - svDF$start)

  if( totalSVdur < min.sv){
    stop("The total sustained vowel duration is less than the threshold length min.sv!")
  }
    
  praat_dsp_directory <- make_dsp_environment()
  
  praat_script <- ifelse(PRAAT_DEVEL== TRUE,
                         file.path("inst","praat","AVQI301.praat"),
                         file.path(system.file(package = "superassp",mustWork = TRUE),"praat","AVQI301.praat"))
  
  #return(praat_script)
  #praat_script <- "/Users/frkkan96/Documents/src/superassp/inst/praat/AVQI301.praat"
  avqi <- cs_wrap_praat_script(praat_location = get_praat(),
                                       script_code_to_run = readLines(praat_script),
                               directory=praat_dsp_directory
                                       ,return="last-argument")
  
  #Check that all files exists before we begin
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)){
    filedNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ",paste(filedNotExists, collapse = ", "))
  }
  
  ## Set up a (CLEAN) directory for interchange
  #avqiDir <- file.path(tempdir(check=TRUE),"avqtemp")
  #unlink(avqiDir,recursive = TRUE,force=FALSE,expand=FALSE)
  #dir.create(avqiDir)
  
  
  #The empty vector of file names that should be returned
  outListOfFiles <- c()
  
  
  #Copy Sustained Vowel portions from the file
  
  #Pre-generate names of output files
  svDF$OutFileName <- file.path(praat_dsp_directory,paste0("sv",1:nrow(svDF),".wav"))
  
  for(r in 1:nrow(svDF)){
    #Here we finally copy out all signal file content into separate files
    currSound <- read.AsspDataObj(fname=as.character(svDF[[r,"listOfFiles"]]),
                                          begin=as.numeric(svDF[[r,"start"]]) /1000, #Segment lists returned by emuR::query have times in milliseconds
                                          end=as.numeric(svDF[[r,"end"]]) / 1000, #Segment lists returned by emuR::query have times in milliseconds
                                          samples=FALSE)
    write.AsspDataObj(currSound,file=svDF[[r,"OutFileName"]])
  }
  
  #Copy Continous Speech portions from the file
  
  #Pre-generate names of output files
  csDF$OutFileName <- file.path(praat_dsp_directory,paste0("cs",1:nrow(csDF),".wav"))
  
  for(r in 1:nrow(csDF)){
    #Here we finally copy out all signal file content into separate files
    currSound <- read.AsspDataObj(fname=csDF[[r,"listOfFiles"]],
                                          begin=csDF[[r,"start"]]/1000, #Segment lists returned by emuR::query have times in milliseconds
                                          end=csDF[[r,"end"]]/1000, #Segment lists returned by emuR::query have times in milliseconds
                                          samples=FALSE)
    write.AsspDataObj(currSound,file=csDF[[r,"OutFileName"]])
  }
  
  #AVQI 3.01
  #boolean Simple_version 1 
  #sentence name_patient Fredrik Karlsson 
  #sentence Date_of_birth 1975-12-31 
  #sentence Assessment_date 2021-12-31 
  #sentence Input_directory ../../tests/signalfiles/AVQI/input 
  #boolean Generate_PDF_files 1 
  #sentence Speaker_ID "1" 
  #sentence Output_directory /Users/frkkan96/Documents/src/superassp/tests/signalfiles/AVQI/output 
  #sentence Output_file /Users/frkkan96/Documents/src/superassp/tests/signalfiles/AVQI/output/avqi.csv
  
  outAVQITabFile <- avqi(ifelse(simple.output,1,0),
                         ifelse(! is.null(speaker.name),speaker.name,"NA"),
                         ifelse(! is.null(speaker.dob),speaker.dob,"NA"),
                         ifelse(! is.null(session.datetime),session.datetime,"NA"),
                         praat_dsp_directory,
                         ifelse(! is.null(pdf.path),1,0),
                         as.character(speaker.ID),
                         praat_dsp_directory,
                         file.path(praat_dsp_directory,"avqi.csv")
  )

  
  inTable <- read.csv(file=outAVQITabFile
                      ,header=TRUE
                      ,na.strings =c("--undefined--","NA"),
                      sep = ",")
  
  assertthat::are_equal(nrow(inTable),1)
  
  # Now copy PDF files to the pdf.output path
  if(!is.null(pdf.path)){
    pdfFiles <- list.files(praat_dsp_directory,pattern=".*[.]pdf",full.names=TRUE)
    for(currFile in pdfFiles){
      file.copy(from=currFile,to=pdf.path,overwrite=overwrite.pdfs)
    }
  }
  clear_dsp_environment(praat_dsp_directory )
  
  logger::log_trace("Computed an AVQI value from ",nrow(svDF)," sustained vowels and ",nrow(csDF), "continuous speech utterances.")
  
  return(as.list(inTable))  
}
attr(praat_avqi,"outputType") <-  c("list")
attr(praat_avqi,"tracks") <-  c("AVQI_VERSION","Speaker","ID","CPPS", "HNR","Shim_local","Shim_local_DB","LTAS_Slope","LTAS_Tilt","AVQI")
attr(praat_avqi,"ext") <-  c("avqi") 

#Interactive testing
# 
# sv <- data.frame(listOfFiles=c(
#   "tests/signalfiles/AVQI/input/sv1.wav",
#   "tests/signalfiles/AVQI/input/sv2.wav",
#   "tests/signalfiles/AVQI/input/sv3.wav",
#   "tests/signalfiles/AVQI/input/sv4.wav"),
#   start=rep(0.0633328955584327 *1000,4),
#   end=rep(2.8305593763398864*1000,4)
#   )
# 
# cs <- data.frame(listOfFiles=c(
#   "tests/signalfiles/AVQI/input/cs1.wav",
#   "tests/signalfiles/AVQI/input/cs2.wav",
#   "tests/signalfiles/AVQI/input/cs3.wav",
#   "tests/signalfiles/AVQI/input/cs4.wav"),
#   start=rep(0.08250407973624065 *1000,4),
#   end=rep(3.738389187054969 *1000,4)
# )
# 
# praat_avqi(sv,cs,speaker.name="Fredrik Karlsson",speaker.ID=1,speaker.dob="1975-01-14",session.datetime = date(), pdf.path = "/Users/frkkan96/Desktop/",simple.output = TRUE) -> avqi_out



#' Compute the components of a Praat Voice report
#' 
#' The Praat program defines a voice report containing a range of fundamental properties of a voice sample. The most common application of the voice report is on a sustained vowel. This function computes the report from a sectio of a recording using Praat, and returns the voice measures as a list. The function also enable the user to mark just a part of the sustained vowel for analysis using an offset and a subsample length. In this scenario, the user specifies the start and end times (`beginTime` and `endTime`, respectively) of the sustained vowel. Then, the user specifies a `selectionOffset`, which is the number of seconds into the vowel where extraction for analysis will start. Finally, the user specifies a `selectionLength`, which is the (maximum) length of the extracted part. This means that if the user has a sustained vowel staring 1 second into the signal and extends for 2 seconds (very short), and the user asks for a 2 second extraction starting 0.5 s into the vowel, what will actually be analysed is a portion from 1.5s to 3s (a 1.5s signal, and not the 2s that the user asked for). This behavior is there so that the user is not inadvertently adding parts that are not part of a sustaind vowel production. The user may of course always choose to disregard measurements that were based on a too short sample.
#' 
#' 
#' which may be advantageous in cases where it may be suspected that 
#'
#' @param listOfFiles The full path of the sound file.
#' @param beginTime The time point (in s) in the sound file where the sustained vowel starts. If `NULL`, the start of the sound file will also be viewed as the start of the sustained vowel production.
#' @param endTime The time point (in s) in the sound file where the sustained vowel ends. If `NULL`, everyting up until the end of the sound file will be considered part of the sustained vowel.
#' @param selectionOffset An optional offset to be added to the time of the sustained vowel production when determining where the start of the extracted portion of the vowel.
#' @param selectionLength An optional (maximal) length of the selection. 
#' @param windowShape The window shape used for extracting the vowel. May be one of "rectangular", "triangular", "parabolic", "Hanning", "Hamming", "Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4", "Gaussian5", "Kaiser1", and "Kaiser2".
#' @param relativeWidth The relative width of the window used for extracting the vowel portion.
#' @param minF The minimum pitch (f~0~) to be considered.
#' @param maxF The maximum pitch (f~0~) to be considered.
#' @param max_period_factor The larges possible differences between consecutive intervals that will be used in computing jitter. Please consult the Praat manual for further information.
#' @param max_ampl_factor The larges possible differences between consecutive intervals that will be used in computing schimmer Please consult the Praat manual for further information.
#' @param silence_threshold The silence threshold. Please consult the Praat manual for further information.
#' @param voicing_threshold The voicing  threshold. Please consult the Praat manual for further information.
#' @param octave_cost The octave cost. Please consult the Praat manual for further information.
#' @param octave_jump_cost The octave jump cost. Please consult the Praat manual for further information.
#' @param voiced_unvoiced_cost  The cost for voiced to unvoiced change detection. Please consult the Praat manual for further information.
#' @param praat_path An optional explicit path to the Praat binary. Not usually required.
#'
#' @return A list of voice parameters:
#' \describe{
#' \item{Median pitch}{The median pitch (f~0~) of the sample (in Hz)}
#' \item{Mean pitch}{The mean pitch (f~0~) of the sample (in Hz)}
#' \item{Standard deviation}{The standard deviation of pitch (f~0~, in Hz) of the sample.}
#' \item{Minimum pitch}{The lowest pitch (f~0~) detected (in Hz)}
#' \item{Maximum pitch}{The highest pitch (f~0~) detected (in Hz)}
#' \item{Number of pulses}{The number of pulses detected}
#' \item{Number of periods}{The number of periods detected}
#' \item{Mean period}{The average period length}
#' \item{Standard deviation of period}{The standard deviation of period length}
#' \item{Fraction of locally unvoiced frames}{The fraction of frames detected as unvoiced in the sample.}
#' \item{Number of voice breaks}{Number of voice breaks}
#' \item{Degree of voice breaks}{The number of voice breaks in relation to the number of frames}
#' \item{Jitter (local)}{The average absolute difference between consequtive periods, divided by the average period (in %). See the Praat manual for more information.}
#' \item{Jitter (local, absolute)}{The average absolute difference between consequtive periods, in seconds. See the Praat manual for more information.}
#' \item{Jitter (rap)}{The three point Relative Average Pertubation: the average absolute difference between a period and the three point local average, divided by the average period (in %).}
#' \item{Jitter (ppq5)}{The five point Relative Average Pertubation: the average absolute difference between a period and the five point local average, divided by the average period (in %).}
#' \item{Jitter (ddp)}{The average absolute difference between consequtive differences between periods, divided by the average period (in %).}
#' \item{Shimmer (local)}{The average absolute difference between amplitudes of consequtive periods, divided by the average amplitude (in %). }
#' \item{Shimmer (local, dB)}{The average absolute difference between amplitudes of consequtive periods (in dB).}
#' \item{Shimmer (apq3)}{The three point Amplitude Pertubation Quotient: the average absolute difference between the amplitude of a period and the three point local average, divided by the average amplitude (in %).}
#' \item{Shimmer (apq5)}{The five point Amplitude Pertubation Quotient: the average absolute difference between the amplitude of a period and the five point local average, divided by the average amplitude (in %).}
#' \item{Shimmer (apq11)}{The 11 point Amplitude Pertubation Quotient: the average absolute difference between the amplitude of a period and the 11 point local average, divided by the average amplitude (in %).}
#' \item{Shimmer (dda)}{The average absolute difference between consequtive differences between amplitudes of consequtive periods, divided by the average period (in %).}
#' \item{Mean autocorrelation}{The average autocorrelation of the signal.}
#' \item{Mean noise-to-harmonics ratio}{The average NHR of the voice sample.}
#' \item{Mean harmonics-to-noise ratio}{The average HNR of the voice sample.}
#' }
#' @export
#'
#'

praat_voice_report <- function(listOfFiles,
                               beginTime=NULL,
                               endTime=NULL,
                               selectionOffset=NULL,
                               selectionLength=NULL,
                               windowShape="Gaussian1",
                               relativeWidth=1.0,
                               minF=75,
                               maxF=600,
                               max_period_factor=1.3,
                               max_ampl_factor=1.6,
                               silence_threshold=0.03,
                               voicing_threshold=0.45,
                               octave_cost=0.01,
                               octave_jump_cost=0.35,
                               voiced_unvoiced_cost=0.14,
                               praat_path=NULL){
  
  
  if(! have_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
  }
  
  #Make sure a valid windows hape is provided
  if(!windowShape %in% c("rectangular", "triangular", "parabolic", "Hanning", "Hamming", "Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4", "Gaussian5", "Kaiser1","Kaiser2")){
    stop("Invalid window shape. Permitted values are  \"rectangular\", \"triangular\", \"parabolic\", \"Hanning\", \"Hamming\", \"Gaussian1\", \"Gaussian2\", \"Gaussian3\", \"Gaussian4\", \"Gaussian5\", \"Kaiser1\", and \"Kaiser2\"")
  }
  
  salt <- paste(beginTime, endTime,collapse = "-")
  
  praat_dsp_directory <- make_dsp_environment()
  

  praat_script <- ifelse(PRAAT_DEVEL== TRUE,
                         file.path("inst","praat","praat_voice_report.praat"),
                         file.path(system.file(package = "superassp",mustWork = TRUE),"praat","praat_voice_report.praat"))
  
  voice_report <- cs_wrap_praat_script(praat_location = get_praat(),
                                       script_code_to_run = readLines(praat_script),
                                       directory=praat_dsp_directory
                                       ,return="last-argument")
  origSoundFile <- normalizePath(listOfFiles)
  
  soundFile <- tempfile(fileext = ".wav")
  R.utils::createLink(soundFile,origSoundFile)

  outputfile <- tempfile(fileext = ".csv")
  
  # Now we are all set up to run the Praat script
  # sentence SoundFile /Users/frkkan96/Desktop/aaa.wav
  # real StartTime 0.0
  # real EndTime 0.0
  # real SelectionOffset 0.0
  # real SelectionLength 2.0
  # real Minimum_f0 75.0
  # real Maximum_f0 600
  # real Maximum_period_factor 1.3
  # real Maximum_amplitude_factor 1.6
  # real Silence_threshold 0.03
  # real Voicing_threshold 0.45
  # real Octave_cost 0.01
  # real Octave_jump_cost 0.35
  # real Voiced/unvoiced_cost 0.14
  # word WindowType Gaussian1
  # real WindowWidth 1.0
  # sentence OutFile /Users/frkkan96/Desktop/aaa.csv

  outVRtabFile <- voice_report(soundFile,
                         ifelse(is.null(beginTime),0.0,beginTime),
                         ifelse(is.null(endTime),0.0,endTime), 
                         ifelse(is.null(selectionOffset),0.0,selectionOffset),
                         ifelse(is.null(selectionLength),0.0,selectionLength),
                         minF,
                         maxF,
                         max_period_factor,
                         max_ampl_factor,
                         silence_threshold,
                         voicing_threshold,
                         octave_cost,
                         octave_jump_cost,
                         voiced_unvoiced_cost,
                         windowShape,
                         relativeWidth,
                         outputfile)
  

  inTable <- read.csv(file=outVRtabFile
                      ,header=TRUE
                      ,na.strings =c("--undefined--","NA"),
                      sep = ";",
                      check.names = FALSE)
  
  assertthat::are_equal(nrow(inTable),1)
  
  
  clear_dsp_environment(praat_dsp_directory )
  
  return(as.list(inTable[,-(1:4)])) 
}

attr(praat_voice_report,"outputType") <-  c("list")
attr(praat_voice_report,"ext") <-  c("pvr") 
attr(praat_voice_report,"tracks") <- c("Start Time","End Time","Selection start","Selection end","Median pitch","Mean pitch","Standard deviation","Minimum pitch","Maximum pitch","Number of pulses","Number of periods","Mean period","Standard deviation of period","Fraction of locally unvoiced frames","Number of voice breaks","Degree of voice breaks","Jitter (local)","Jitter (local, absolute)","Jitter (rap)","Jitter (ppq5)","Jitter (ddp)","Shimmer (local)","Shimmer (local, dB)","Shimmer (apq3)","Shimmer (apq5)","Shimmer (apq11)","Shimmer (dda)","Mean autocorrelation","Mean noise-to-harmonics ratio","Mean harmonics-to-noise ratio")


#' Compute the Dysphonia Severity Index
#'
#' This function computes the Dysphonia Severity Index (DSI)
#' \insertCite{Wuyts:2000vb}{superassp} as implemented into Praat by
#' \insertCite{Maryn.2017.10.1016/j.jvoice.2017.01.002;textual}{superassp}. The
#' user is asked to supply data.frames indicating the "absolute_file_path" as
#' well as "start", and "end" times for parts of sound files that together form
#' the basis for the DSI computation for a speaker.
#'
#' The user has to indicate at least one sample in which the participant speaks
#' as softly as possible (`softDF`), at least one sample from which the maximum
#' f~0~ could be deduced (`highpitchDF`), at least one sample where a single
#' vowel is maximally prolonged (`maxprolongedDF`). The user can also provide a
#' sample of a maximally stable vowel (`stableDF`), but if such a sample is not
#' provided then the (`maxprolongedDF`) sample will be reused instead for the
#' computation of the DSI sub-component Jitter. The user may submit multiple
#' sound samples for all these sets of acoustic inputs, and the sounds files
#' will then be combined before the DSI sub components are computed. The largest
#' Maximum performance time will be used.
#'
#'
#' @param softDF A data.frame containing the file path and start and end times
#'   of samples indicating the softest voice production of the speaker.
#' @param highpitchDF A data.frame containing the file path and start and end
#'   times of samples from which the highest f0 of the speaker can be computed.
#' @param maxprolongedDF A data.frame containing the file path and start and end
#'   times of samples indicating maximally prolonged vowels of the speaker.
#' @param stableDF An optional data.frame containing the file path and start and
#'   end times of samples indicating a stable sustained vowel. If not provided,
#'   the `maxprolongedDF` will be used instead.
#' @param use.calibration Should a calibration factor be added to measured
#'   intensities before computing DSI?
#' @param db.calibration The number of dB:s to be added to measured intensity
#'   values before computing DSI.
#' @param speaker.name An optional name of the speaker. Only used when PDF
#'   output is produced.
#' @param speaker.ID An ID of the speaker. This will be included in the list
#'   output.
#' @param speaker.dob An optional date of birth of the speaker. Only used when
#'   PDF output is produced.
#' @param session.datetime A string indicating the session date (and time). Only
#'   used when PDF output is produced.
#' @param pdf.path A path where a PDF report file will be produced.
#' @param overwrite.pdfs Should existing PDF files be overwritten in the PDF
#'   output directory? Defaults to a safe behavior where older PDFs are not
#'   overwritten.
#' @param praat_path An explicit path to the Praat binary.
#'
#' @return A list with the following fields: 
#' \describe{ 
#' \item{ID}{The speaker speaker + session identifier of the output}
#'   \item{Maximum.phonation.time}{The speakers maximum phonation time.}
#'   \item{Softest.intensity.of.voiced.speech}{The intensity of the participants
#'   softest voice (in dB)} 
#'   \item{Maximum.fundamental.frequency}{The maximum f~0~
#'   frequency (in Hz)} 
#'   \item{Jitter.ppq5}{The five-point Period Pertubation
#'   Quotient measurement (in %)} 
#'   \item{Dysphonia.Severity.Index}{The resulting
#'   Dysphonia Severity Index } 
#'   }
#'
#'   A meta-analysis in which 1330 healthy participants (aged 17.3−94) were included indicated an average DSI of 3.05 in healthy speakers, with a conﬁdence level  between
#'   2.13 and 3.98 \insertCite{Sobol.2020.10.1016/j.jvoice.2020.04.010}{superassp}. 
#'
#'
#' @export
#'
#' @references 
#' \insertAllCited{}
#' 

praat_dsi <- function(softDF,
                      highpitchDF,
                      maxprolongedDF,
                      stableDF=NULL,
                      use.calibration=FALSE,
                      db.calibration=10,
                      speaker.name=NULL,
                      speaker.ID=NULL,
                      speaker.dob=NULL,
                      session.datetime=NULL,
                      pdf.path=NULL,
                      overwrite.pdfs=FALSE,
                      praat_path=NULL){
  
  
  requiredDFColumns <- c("absolute_file_path","start","end")
  # Reuse max prolonged vowels as stable productions.
  
  if(is.null(stableDF) ) {
    stableDF <- maxprolongedDF
  }
  if(nrow(softDF) < 1 || nrow(highpitchDF) < 1 || nrow(maxprolongedDF) < 1 || nrow(stableDF) < 1 ){
    stop("All of the supplied dataframes indicating samples to be analysed must be non-empty.")
  }
  
  listOfFiles <- unique(
    c(softDF$absolute_file_path,
      highpitchDF$absolute_file_path,
      maxprolongedDF$absolute_file_path,
      stableDF$absolute_file_path)
  )
  
  if(! have_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
  }
  

  if(! requiredDFColumns %in% names(softDF) || 
     ! requiredDFColumns %in% names(highpitchDF)  || 
     ! requiredDFColumns %in% names(maxprolongedDF) || 
     ! requiredDFColumns %in% names(stableDF)   
  ){
    stop("All dataframes must both contain columns named ",paste(requiredDFColumns,collapse=",",sep=""),".")
  }
  
  praat_dsp_directory <- make_dsp_environment()

  
  praat_script <- ifelse(PRAAT_DEVEL== TRUE,
                         file.path("inst","praat","DSI201.praat"),
                         file.path(system.file(package = "superassp",mustWork = TRUE),"praat","DSI201.praat"))
  

  dsi <- cs_wrap_praat_script(praat_location = get_praat(),
                            script_code_to_run = readLines(praat_script),
                            directory=praat_dsp_directory,      
                            return="last-argument")
  
  #Check that all files exists before we begin
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)){
    filedNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ",paste(filedNotExists, collapse = ", "))
  }
  
  # Set up a (CLEAN) directory for interchange
  dsiDir <- file.path(tempdir(check=TRUE),"dsitemp")
  unlink(dsiDir,recursive = TRUE,force=FALSE,expand=FALSE)
  dir.create(dsiDir)
  
  
  #The empty vector of file names that should be returned
  outListOfFiles <- c()

  #Copy soft voice samples
  
  #Pre-generate names of output files
  softDF$OutFileName <- file.path(dsiDir,paste0("im",1:nrow(softDF),".wav"))
  
  for(r in 1:nrow(softDF)){
    #Here we finally copy out all signal file content into separate files
    currSound <- read.AsspDataObj(fname=softDF[r,"absolute_file_path"],begin=softDF[r,"start"],end=softDF[r,"end"])
    write.AsspDataObj(currSound,file=softDF[r,"OutFileName"])
  }
  
  #Copy highest pitch samples
  
  #Pre-generate names of output files
  highpitchDF$OutFileName <- file.path(dsiDir,paste0("fh",1:nrow(highpitchDF),".wav"))
  
  for(r in 1:nrow(highpitchDF)){
    #Here we finally copy out all signal file content into separate files
    currSound <- read.AsspDataObj(fname=highpitchDF[r,"absolute_file_path"],begin=highpitchDF[r,"start"],end=highpitchDF[r,"end"])
    write.AsspDataObj(currSound,file=highpitchDF[r,"OutFileName"])
  }
  
  #Copy maximally prolonged vowel samples
  
  #Pre-generate names of output files
  maxprolongedDF$OutFileName <- file.path(dsiDir,paste0("mpt",1:nrow(maxprolongedDF),".wav"))
  
  for(r in 1:nrow(maxprolongedDF)){
    #Here we finally copy out all signal file content into separate files
    currSound <- read.AsspDataObj(fname=maxprolongedDF[r,"absolute_file_path"],begin=maxprolongedDF[r,"start"],end=maxprolongedDF[r,"end"])
    write.AsspDataObj(currSound,file=maxprolongedDF[r,"OutFileName"])
  }
  
  #Copy stable vowel samples
  
  #Pre-generate names of output files
  stableDF$OutFileName <- file.path(dsiDir,paste0("ppq",1:nrow(stableDF),".wav"))
  
  for(r in 1:nrow(stableDF)){
    #Here we finally copy out all signal file content into separate files
    currSound <- read.AsspDataObj(fname=stableDF[r,"absolute_file_path"],begin=stableDF[r,"start"],end=stableDF[r,"end"])
    write.AsspDataObj(currSound,file=stableDF[r,"OutFileName"])
  }
  
  
  # Now we are all set up to run the Praat script
  # boolean Apply_calibration 1
  # # sentence Calibration_factor 1*self+10
  # sentence Name_of_patient Fredrik Karlsson
  # sentence Date_of_birth 1975-12-31
  # sentence Assessment_date 2021-12-31
  # sentence Input_directory ../../tests/signalfiles/DSI/input
  # #/Users/frkkan96/Documents/src/superassp/tests/signalfiles/DSI/input
  # boolean Generate_PDF_files 1
  # sentence Speaker_ID 1
  # sentence Output_directory /Users/frkkan96/Documents/src/superassp/tests/signalfiles/DSI/output
  # sentence Output_file /Users/frkkan96/Documents/src/superassp/tests/signalfiles/DSI/output/dsi.csv
  
  outDSITabFile <- dsi(ifelse(use.calibration,1,0),
                       db.calibration,
                       ifelse(! is.null(speaker.name),shQuote(speaker.name),"NA"),
                       ifelse(! is.null(speaker.dob),shQuote(speaker.dob),"NA"),
                       ifelse(! is.null(session.datetime),shQuote(session.datetime),"NA"),
                       dsiDir,
                       ifelse(! is.null(pdf.path),1,0),
                       ifelse(! is.null(speaker.ID),speaker.ID,"NA"),
                       dsiDir,
                       file.path(dsiDir,"dsi.csv")
  )
  
  inTable <- read.csv(file=outDSITabFile
                      ,header=TRUE
                      ,na.strings =c("--undefined--","NA"),
                      sep = ",")
  
  assertthat::are_equal(nrow(inTable),1)
  
  # Now copy PDF files to the pdf.output path
  if(!is.null(pdf.path)){
    pdfFiles <- list.files(dsiDir,pattern=".*[.]pdf",full.names=TRUE)
    for(currFile in pdfFiles){
      file.copy(from=currFile,to=pdf.path,overwrite=overwrite.pdfs)
    }
  }
  
  clear_dsp_environment(  praat_dsp_directory )
  
  return(as.list(inTable))  
}
attr(praat_dsi,"outputType") <-  c("list")
attr(praat_dsi,"tracks") <-  c("ID","Maximum.phonation.time", "Softest.intensity.of.voiced.speech","Maximum.fundamental.frequency","Jitter.ppq5","Dysphonia.Severity.Index")
attr(praat_dsi,"ext") <-  c("dsi") 

#' Computes 18 measures of (vocal) tremor
#'
#' This function calls Praat and applies the "tremor 3.05" package
#' \insertCite{tremor305}{superassp} to compute measures of
#' voice tremor \insertCite{Bruckl2017}{superassp} on a
#' sustained vowel. The user may compute the measures for the entire
#' sample or identify a portion of the sound file where the sustained is
#' produced (`beginTime` and `endTime` arguments). The user may additionally
#' specify an offset (`selectionOffset`) from the `beginTime` where measurements
#' should start, as well as a (maximum) selection length. Using these arguments
#' the user may mark where all sustained vowels are produced in a database, and
#' apply this function across only a 2 second sample  from 1 second into the
#' vowel, and thus avoid vowel initial effects in phonation affecting measurements.
#'
#' @param listOfFiles The full path of the sound file.
#' @param beginTime The time point (in s) in the sound file where the sustained vowel starts. If `NULL`, the start of the sound file will also be viewed as the start of the sustained vowel production.
#' @param endTime The time point (in s) in the sound file where the sustained vowel ends. If `NULL`, everyting up until the end of the sound file will be considered part of the sustained vowel.
#' @param selectionOffset An optional offset to be added to the time of the sustained vowel production when determining where the start of the extracted portion of the vowel.
#' @param selectionLength An optional (maximal) length of the selection. 
#' @param windowShape The window shape used for extracting the vowel. May be one of "rectangular", "triangular", "parabolic", "Hanning", "Hamming", "Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4", "Gaussian5", "Kaiser1", and "Kaiser2".
#' @param relativeWidth The relative width of the window used for extracting the vowel portion.
#' @param minF The minimum pitch (f~0~) to be considered.
#' @param maxF The maximum pitch (f~0~) to be considered.
#' @param windowShift The number of milliseconds between analysis window centers.
#' @param max_period_factor The larges possible differences between consecutive intervals that will be used in computing jitter. Please consult the Praat manual for further information.
#' @param max_ampl_factor The larges possible differences between consecutive intervals that will be used in computing schimmer Please consult the Praat manual for further information.
#' @param silence_threshold The silence threshold. Please consult the Praat manual for further information.
#' @param voicing_threshold The voicing  threshold. Please consult the Praat manual for further information.
#' @param octave_cost The octave cost. Please consult the Praat manual for further information.
#' @param octave_jump_cost The octave jump cost. Please consult the Praat manual for further information.
#' @param voiced_unvoiced_cost  The cost for voiced to unvoiced change detection. Please consult the Praat manual for further information.
#' @param na.zero Should measurements of voice tremor that could not be determined be returned as a zero value, or as an `NA` (the default).
#' @param praat_path An optional explicit path to the Praat binary. Not usually required.
#'
#' @return A list of voice tremor measurements :
#' \describe{
#' \item{FCoM}{frequency contour magnitude}
#' \item{FTrC}{(maximum) frequency tremor cyclicality}
#' \item{FMon}{number of frequency modulations above thresholds}
#' \item{FTrF [Hz]}{(strongest) frequency tremor frequency}
#' \item{FTrI [%]}{frequency tremor intensity index}
#' \item{FTrP}{frequency tremor power index}
#' \item{FTrCIP}{frequency tremor cyclicality intensity product}
#' \item{FTrPS}{frequency tremor product sum}
#' \item{FCoHNR[dB]}{frequency contour harmonicity-to-noise ratio}
#' \item{ACoM}{amplitude contour magnitude}
#' \item{ATrC}{(maximum) amplitude tremor cyclicality}
#' \item{AMoN}{number of amplitude modulations above thresholds}
#' \item{ATrF [Hz]}{(strongest) amplitude tremor frequency}
#' \item{ATrI [%]}{amplitude tremor intensity index}
#' \item{ATrP}{amplitude tremor power index}
#' \item{ATrCIP}{amplitude tremor cyclicality intensity product }
#' \item{ATrPS}{amplitude tremor product sum}
#' \item{ACoHNR[dB]}{amplitude contour harmonicity-to-noise ratio}
#' }
#' @export
#'
#' @references 
#' \insertAllCited{}
#'

praat_voice_tremor <- function(listOfFiles,
                               beginTime=NULL,
                               endTime=NULL,
                               selectionOffset=NULL,
                               selectionLength=3.0,
                               windowShape="Gaussian1",
                               relativeWidth=1.0,
                               minF=75,
                               maxF=600,
                               windowShift=5,
                               max_period_factor=1.3,
                               max_ampl_factor=1.6,
                               silence_threshold=0.03,
                               voicing_threshold=0.45,
                               octave_cost=0.01,
                               octave_jump_cost=0.35,
                               voiced_unvoiced_cost=0.14,
                               envelope.amplitude=TRUE,
                               min.tremor.hz=1.5,
                               max.tremor.hz=15,
                               contour.magnitude.threshold=0.01,
                               tremor.cyclicality.threshold=0.15,
                               freq.tremor.octave.cost=0.01,
                               ampl.tremor.octave.cost=0.01,
                               na.zero=FALSE,
                               praat_path=NULL){


  if(! have_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
  }
  
  #Make sure a valid window shape is provided
  if(!windowShape %in% c("rectangular", "triangular", "parabolic", "Hanning", "Hamming", "Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4", "Gaussian5", "Kaiser1","Kaiser2")){
    stop("Invalid window shape. Permitted values are  \"rectangular\", \"triangular\", \"parabolic\", \"Hanning\", \"Hamming\", \"Gaussian1\", \"Gaussian2\", \"Gaussian3\", \"Gaussian4\", \"Gaussian5\", \"Kaiser1\", and \"Kaiser2\"")
  }
  
  
  praat_script <- ifelse(PRAAT_DEVEL== TRUE,
                         file.path("inst","praat","tremor3.05","console_tremor305.praat"),
                         file.path(system.file(package = "superassp",mustWork = TRUE),"praat","tremor3.05","console_tremor305.praat"))
  
  additional_script_names <- c("amptrem.praat","analysisinout.praat","freqtrem.praat","getCyclicality.praat","runinout.praat","singleruninout.praat","tremIntIndex.praat","tremProdSum.praat") 
  
  # Set up a (CLEAN) directory for additional scripts
  proceduresInDir <-ifelse(PRAAT_DEVEL== TRUE,
                         file.path("inst","praat","tremor3.05","procedures"),
                         file.path(system.file(package = "superassp",mustWork = TRUE),"praat","tremor3.05","procedures")) 
  
  praat_dsp_directory <- make_dsp_environment()

  

  #Copy additional files
  copied <- file.copy(paste0(proceduresInDir,.Platform$file.sep),praat_dsp_directory,overwrite = TRUE,recursive = TRUE)

  
  voice_tremor <- cs_wrap_praat_script(praat_location = get_praat(),
                                               script_code_to_run = readLines(praat_script),
                                               directory=  praat_dsp_directory,
                                               return="last-argument")
  origSoundFile <- normalizePath(listOfFiles)
  

  soundFile <- tempfile(fileext = ".wav")
  R.utils::createLink(soundFile,origSoundFile)
  
  outputfile <- tempfile(fileext = ".csv")
  
  # Now we are all set up to run the Praat script
  # real StartTime_(s) 0.0
  # real EndTime_(s) 0.0
  # real SelectionOffset 0.0
  # real SelectionLength 2.0
  # word WindowType Gaussian1
  # real WindowWidth 1.0
  # positive Analysis_time_step_(s) 0.015
  # comment Arguments for initial pitch extraction
  # positive Minimal_pitch_(Hz) 60
  # positive Maximal_pitch_(Hz) 350
  # positive Silence_threshold 0.03
  # positive Voicing_threshold 0.3
  # positive Octave_cost 0.01
  # positive Octave-jump_cost 0.35
  # positive Voiced_/_unvoiced_cost 0.14
  # comment Arguments for tremor extraction from contours
  # optionmenu Amplitude_extraction_method 2
  # option Integral [RMS per pitch period]
  # option Envelope [To AmplitudeTier (period)]
  # positive Minimal_tremor_frequency_(Hz) 1.5
  # positive Maximal_tremor_frequency_(Hz) 15
  # positive Contour_magnitude_threshold 0.01
  # positive Tremor_cyclicality_threshold 0.15
  # positive Frequency_tremor_octave_cost 0.01
  # positive Amplitude_tremor_octave_cost 0.01
  # optionmenu Output_of_indeterminate_values 2
  # option indeterminate values are replaced by zeros
  # option indeterminate values are -- undefined --
  #   sentence Path_of_sound_to_be_analyzed /Users/frkkan96/Desktop/aaa.wav
  # sentence Path_and_name_of_result_csv /Users/frkkan96/Desktop/aaa.csv
  
  outVTtabFile <- voice_tremor(ifelse(is.null(beginTime),0.0,beginTime),
                               ifelse(is.null(endTime),0.0,endTime), 
                               ifelse(is.null(selectionOffset),0.0,selectionOffset),
                               ifelse(is.null(selectionLength),0.0,selectionLength),
                               windowShape,
                               relativeWidth,
                               windowShift/1000, #Praat wants seconds
                               minF,
                               maxF,
                               silence_threshold,
                               voicing_threshold,
                               octave_cost,
                               octave_jump_cost,
                               voiced_unvoiced_cost,
                               ifelse(envelope.amplitude,"Envelope [To AmplitudeTier (period)]","Integral [RMS per pitch period]"), 
                               min.tremor.hz,
                               max.tremor.hz,
                               contour.magnitude.threshold,
                               tremor.cyclicality.threshold,
                               freq.tremor.octave.cost,
                               ampl.tremor.octave.cost,
                               ifelse(na.zero,"indeterminate values are replaced by zeros","indeterminate values are -- undefined --"),
                               soundFile,
                               outputfile)
  
  
  inTable <- read.csv(file=outVTtabFile
                      ,header=TRUE
                      ,na.strings =c("--undefined--","NA"),
                      sep = ",",
                      check.names = FALSE)
  

  assertthat::are_equal(nrow(inTable),1)
  
  clear_dsp_environment(  praat_dsp_directory )
  
  return(as.list(inTable[,-(1:4)]))
}
attr(praat_voice_tremor,"outputType") <-  c("list")
attr(praat_voice_tremor,"ext") <-  c("pvt") 
attr(praat_voice_tremor,"tracks") <- c("FCoM","FTrC","FMon","FTrF [Hz]","FTrI [%]","FTrP","FTrCIP","FTrPS","FCoHNR[dB]","ACoM","ATrC","AMoN","ATrF [Hz]","ATrI [%]","ATrP","ATrCIP","ATrPS","ACoHNR[dB]")



# FOR INTERACTIVE TESTING
#df <- data.frame("absolute_file_path"="~/Desktop/kaa_yw_pb.wav","start"=0,"end"=3)
#df2 <- data.frame("absolute_file_path"=c("~/Desktop/kaa_yw_pb.wav","~/Desktop/kaa_yw_pb.wav"),"start"=c(0,0),"end"=c(3,3))
#  
#praat_avqi(df, df,pdf.path="~/Desktop/test/",simple.output = TRUE) -> out
#praat_dsi(df, df,df2,pdf.path="~/Desktop/test/") -> out
# praat_dsi(, df,df2,pdf.path="~/Desktop/test/") -> out
#praat_voice_report("~/Desktop/aaa_sample.wav") -> out
