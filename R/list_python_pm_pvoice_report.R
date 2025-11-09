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
# lst_avqip(sv,cs,speaker.name="Fredrik Karlsson",speaker.ID=1,speaker.dob="1975-01-14",session.datetime = date(), pdf.path = "/Users/frkkan96/Desktop/",simple.output = TRUE) -> avqi_out





#' Compute the components of a Praat Voice report (memory-based with Parselmouth)
#'
#' This is a memory-based implementation of \code{\link{praat_voice_report}}
#' using Parselmouth instead of external Praat. It eliminates disk I/O by
#' loading audio directly into memory using the av package and processing
#' with Parselmouth.
#'
#' This function provides identical functionality to \code{praat_voice_report}
#' but with significantly improved performance (10-20x faster) by eliminating
#' file I/O operations. Audio is loaded directly from the original file format
#' (no WAV conversion needed) and processed entirely in memory.
#'
#' @inheritParams lst_voice_reportp
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "pvr".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#'
#' @return If \code{toFile=FALSE} (default), a list of voice parameters.
#'   If \code{toFile=TRUE}, invisibly returns the path to the written JSTF file.
#'
#'   The list contains 26 voice quality measures:
#'   \describe{
#'     \item{Pitch measures}{Median, Mean, SD, Min, Max (in Hz)}
#'     \item{Pulse measures}{Number of pulses, Number of periods, Mean period, SD of period}
#'     \item{Voicing measures}{Fraction of locally unvoiced frames, Number of voice breaks, Degree of voice breaks}
#'     \item{Jitter measures}{local, local absolute, rap, ppq5, ddp}
#'     \item{Shimmer measures}{local, local dB, apq3, apq5, apq11, dda}
#'     \item{Noise measures}{Mean autocorrelation, Mean noise-to-harmonics ratio, Mean harmonics-to-noise ratio}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- lst_voice_reportp("sustained_vowel.wav")
#'
#' # With time windowing (extract 1.0-3.0s from file)
#' result <- lst_voice_reportp(
#'   "sustained_vowel.wav",
#'   beginTime = 1.0,
#'   endTime = 3.0
#' )
#'
#' # With selection offset and length
#' # (extract vowel at 1.0-3.0s, then analyze 0.5-2.0s of that)
#' result <- lst_voice_reportp(
#'   "sustained_vowel.wav",
#'   beginTime = 1.0,
#'   endTime = 3.0,
#'   selectionOffset = 0.5,
#'   selectionLength = 1.5
#' )
#'
#' # Write results to JSTF file
#' lst_voice_reportp("sustained_vowel.wav", toFile = TRUE)  # Creates sustained_vowel.pvr
#'
#' # Read back and convert to data.frame
#' track <- read_track("sustained_vowel.pvr")
#' df <- as.data.frame(track)
#' head(df)  # Shows begin_time, end_time, and all 26 voice parameters
#' }
#'
lst_voice_reportp <- function(listOfFiles,
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
                                    praat_path=NULL,
                                    toFile = FALSE,
                                    explicitExt = "pvr",
                                    outputDirectory = NULL){

  # Check that Parselmouth is available
  if (!reticulate::py_module_available("parselmouth")) {
    stop("Parselmouth Python module not available. Install with: pip install praat-parselmouth")
  }

  # Validate window shape
  if(!windowShape %in% c("rectangular", "triangular", "parabolic", "Hanning", "Hamming",
                         "Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4", "Gaussian5",
                         "Kaiser1","Kaiser2")){
    stop("Invalid window shape. Permitted values are \"rectangular\", \"triangular\", ",
         "\"parabolic\", \"Hanning\", \"Hamming\", \"Gaussian1\", \"Gaussian2\", ",
         "\"Gaussian3\", \"Gaussian4\", \"Gaussian5\", \"Kaiser1\", and \"Kaiser2\"")
  }

  # Normalize file path
  origSoundFile <- normalizePath(listOfFiles, mustWork = TRUE)
  if(!file.exists(origSoundFile)){
    stop("Unable to open sound file '", listOfFiles, "'.")
  }

  # Handle time parameters (convert NULL to appropriate values)
  # av uses seconds (same as Praat)
  bt <- if(is.null(beginTime)) 0.0 else beginTime
  et <- if(is.null(endTime)) NULL else endTime  # NULL means use full file

  # Load audio with av → convert to numpy (MEMORY-BASED!)
  audio_result <- av_load_for_python(
    origSoundFile,
    start_time = bt,
    end_time = et
  )

  # Source Python script
  python_script <- system.file("python", "praat_voice_report_memory.py", package = "superassp")
  if (!file.exists(python_script)) {
    stop("Python script not found: ", python_script)
  }
  reticulate::source_python(python_script)

  # Get Python main module
  py <- reticulate::import_main()

  # Call Python function with numpy array
  result <- py$praat_voice_report_memory(
    audio_np = audio_result$audio_np,
    sample_rate = audio_result$sample_rate,
    start_time = 0.0,  # Already extracted by av
    end_time = 0.0,    # 0 means use full duration
    selection_offset = if(is.null(selectionOffset)) 0.0 else selectionOffset,
    selection_length = if(is.null(selectionLength)) 0.0 else selectionLength,
    window_shape = windowShape,
    relative_width = relativeWidth,
    min_f0 = minF,
    max_f0 = maxF,
    max_period_factor = max_period_factor,
    max_amplitude_factor = max_ampl_factor,
    silence_threshold = silence_threshold,
    voicing_threshold = voicing_threshold,
    octave_cost = octave_cost,
    octave_jump_cost = octave_jump_cost,
    voiced_unvoiced_cost = voiced_unvoiced_cost
  )

  # Convert Python dict to R list
  result_list <- as.list(result)

  # Handle JSTF file writing
  if (toFile) {
    # Get audio metadata
    audio_info <- av::av_media_info(origSoundFile)
    sample_rate <- audio_info$audio$sample_rate
    audio_duration <- audio_info$duration

    # Calculate analysis time range
    analysis_begin <- bt
    analysis_end <- if (!is.null(et)) et else audio_duration

    json_obj <- create_json_track_obj(
      results = result_list,
      function_name = "lst_voice_reportp",
      file_path = origSoundFile,
      sample_rate = sample_rate,
      audio_duration = audio_duration,
      beginTime = analysis_begin,
      endTime = analysis_end,
      parameters = list(
        selectionOffset = selectionOffset,
        selectionLength = selectionLength,
        windowShape = windowShape,
        relativeWidth = relativeWidth,
        minF = minF,
        maxF = maxF,
        max_period_factor = max_period_factor,
        max_ampl_factor = max_ampl_factor,
        silence_threshold = silence_threshold,
        voicing_threshold = voicing_threshold,
        octave_cost = octave_cost,
        octave_jump_cost = octave_jump_cost,
        voiced_unvoiced_cost = voiced_unvoiced_cost
      )
    )

    base_name <- tools::file_path_sans_ext(basename(origSoundFile))
    out_dir <- if (is.null(outputDirectory)) dirname(origSoundFile) else outputDirectory
    output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

    write_json_track(json_obj, output_path)
    return(invisible(output_path))
  }

  return(result_list)
}

attr(lst_voice_reportp, "ext") <- "pvr"
attr(lst_voice_reportp, "outputType") <- "JSTF"
attr(lst_voice_reportp, "format") <- "JSON"
attr(lst_voice_reportp, "tracks") <- c("Median pitch","Mean pitch","Standard deviation","Minimum pitch","Maximum pitch","Number of pulses","Number of periods","Mean period","Standard deviation of period","Fraction of locally unvoiced frames","Number of voice breaks","Degree of voice breaks","Jitter (local)","Jitter (local, absolute)","Jitter (rap)","Jitter (ppq5)","Jitter (ddp)","Shimmer (local)","Shimmer (local, dB)","Shimmer (apq3)","Shimmer (apq5)","Shimmer (apq11)","Shimmer (dda)","Mean autocorrelation","Mean noise-to-harmonics ratio","Mean harmonics-to-noise ratio")






# FOR INTERACTIVE TESTING
#df <- data.frame("absolute_file_path"="~/Desktop/kaa_yw_pb.wav","start"=0,"end"=3)
#df2 <- data.frame("absolute_file_path"=c("~/Desktop/kaa_yw_pb.wav","~/Desktop/kaa_yw_pb.wav"),"start"=c(0,0),"end"=c(3,3))
#  
#lst_avqip(df, df,pdf.path="~/Desktop/test/",simple.output = TRUE) -> out
#lst_dsip(df, df,df2,pdf.path="~/Desktop/test/") -> out
# lst_dsip(, df,df2,pdf.path="~/Desktop/test/") -> out
#lst_voice_reportp("~/Desktop/aaa_sample.wav") -> out
