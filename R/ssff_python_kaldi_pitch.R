
#' Estimate pitch using the Kaldi modifies version of RAPT
#'
#' @section Deprecation Warning:
#' **This function requires torchaudio < 2.9** as `torchaudio.functional.compute_kaldi_pitch`
#' has been removed in torchaudio 2.9+. Consider using alternative pitch trackers:
#' [crepe()], [pyin()], [swipe()], [rapt()], or [reaper()].
#'
#' For more information, see: https://github.com/pytorch/audio/issues/3902
#'
#' @description
#' The algorithm used is a version of the [RAPT][rapt] algorithm
#' that considers voicing also in voiceless frames and conputes a
#' Normalized Cross Correlation Function (NCCF) that can be used to
#' estimate the probability of voicing \insertCite{Ghahremani.2014.10.1109/icassp.2014.6854049}{superassp}.
#'
#' The function calls the [torchaudio](https://github.com/pytorch/audio) \insertCite{yang2021torchaudio}{superassp} library to do the pitch estimates and therefore
#' relies on it being present in a properly set up python environment to work. Please refer to the [torchaudio manual](https://pytorch.org/audio/main/generated/torchaudio.functional.compute_kaldi_pitch.html) for further information.
#'
#'
#' @inheritParams rapt
#' @param soft_min_f0 (float, optional) – Minimum f0, applied in soft way, must not exceed min-f0 (default: 10.0)
#' @param penalty_factor Cost factor for fO change. (default: 0.1)
#' @param lowpass_cutoff Cutoff frequency for LowPass filter (Hz) (default: 1000)
#' @param resample_frequency Frequency that we down-sample the signal to. Must be more than twice `lowpass_cutoff`. (default: 4000)
#' @param delta_pitch  Smallest relative change in pitch that our algorithm measures. (default: 0.005)
#' @param nccf_ballast Increasing this factor reduces NCCF for quiet frames (default: 7000)
#' @param lowpass_filter_width  Integer that determines filter width of lowpass filter, more gives sharper filter. (default: 1)
#' @param psample_filter_width  Integer that determines filter width when upsampling NCCF. (default: 5)
#' @param max_frames_latency  Maximum number of frames of latency that we allow pitch tracking to introduce into the feature processing (affects output only if `frames_per_chunk` > 0 and `simulate_first_pass_online`=`TRUE`) (default: 0)
#' @param frames_per_chunk The number of frames used for energy normalization. (default: 0)
#' @param simulate_first_pass_online  If true, the function will output features that correspond to what an online decoder would see in the first pass of decoding – not the final version of the features, which is the default. (default: `FALSE`) Relevant if `frames_per_chunk > 0`.
#' @param recompute_frame  Only relevant for compatibility with online pitch extraction. A non-critical parameter; the frame at which we recompute some of the forward pointers, after revising our estimate of the signal energy. Relevant if `frames_per_chunk > 0`. (default: 500)
#' @param snip_edges  If this is set to false, the incomplete frames near the ending edge won’t be snipped, so that the number of frames is the file size divided by the `windowShift`. This makes different types of features give the same number of frames. (default: True)
#'
#' @return An SSFF track object containing two tracks (f0 and nccf) that are
#'   either returned (toFile == FALSE) or stored on disk.
#'
#' @seealso rapt
#'
#' @references \insertAllCited{}
#' @export


kaldi_pitch <- function(listOfFiles,
                 beginTime=0,
                 endTime=0,
                 windowShift=5,
                 windowSize=25,
                 minF=70,
                 maxF=200,
                 softMinF0 = 10.0,
                 voiced_voiceless_cost = 0.10,
                 owpass_cutoff =1000.0,
                 resample_frequency = 4000.0,
                 deltaChange = 0.005,
                 nccfBallast = 7000,
                 lowpass_cutoff =1000,
                 lowpass_filter_width = 1,
                 upsample_filter_width = 5,
                 max_frames_latency = 0,
                 frames_per_chunk = 0,
                 simulate_first_pass_online = FALSE,
                 recompute_frame=500,
                 snip_edges = TRUE,
                 explicitExt="kap",
                 outputDirectory=NULL,
                 toFile=TRUE,
                 conda.env=NULL){

  if(!is.null(conda.env) && !conda.env %in%  reticulate::conda_list()$name){
    stop("The conda environment ",conda.env, " does not exist.\n Please ")
  }

  if(length(listOfFiles) > 1 & ! toFile){
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
  }

  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles,
      beginTime = beginTime,
      endTime=endTime
    )
  },error=function(e){stop("The beginTime and endTime must either be a single value or the same length as listOfFiles")})
  #Check that all files exists before we begin
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)){
    filedNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ",paste(filedNotExists, collapse = ", "))
  }
  #The empty vector of file names that should be returned
  outListOfFiles <- c()

  for(i in 1:nrow(fileBeginEnd)){
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"],mustWork = TRUE)

    beginTime <- fileBeginEnd[i, "beginTime"]
    endTime <- fileBeginEnd[i, "endTime"]

    # Initialize Python environment
    py <- reticulate::import_main()

    py$soundFile <- reticulate::r_to_py(origSoundFile)
    py$windowShift <- reticulate::r_to_py(as.double(windowShift))
    py$windowSize <- reticulate::r_to_py(as.double(windowSize))
    py$maxF <- reticulate::r_to_py(as.double(maxF))
    py$minF <- reticulate::r_to_py(as.double(minF))
    py$beginTime <- reticulate::r_to_py(as.double(beginTime))
    py$endTime <- reticulate::r_to_py(as.double(endTime))
    py$softMinF0 <- reticulate::r_to_py(as.double(softMinF0))
    py$voiced_voiceless_cost <- reticulate::r_to_py(as.double(voiced_voiceless_cost))
    py$lowpass_cutoff <- reticulate::r_to_py(as.double(lowpass_cutoff))
    py$resample_frequency <- reticulate::r_to_py(as.double(resample_frequency))
    py$deltaChange <- reticulate::r_to_py(as.double(deltaChange))
    py$nccfBallast <- reticulate::r_to_py(as.double(nccfBallast))
    py$lowpass_filter_width <- reticulate::r_to_py(as.integer(lowpass_filter_width))
    py$upsample_filter_width <- reticulate::r_to_py(as.integer(upsample_filter_width))
    py$max_frames_latency <- reticulate::r_to_py(as.integer(max_frames_latency))
    py$frames_per_chunk <- reticulate::r_to_py(as.integer(frames_per_chunk))
    py$simulate_first_pass_online <- reticulate::r_to_py(as.logical(simulate_first_pass_online))
    py$recompute_frame <- reticulate::r_to_py(as.integer(recompute_frame))
    py$snip_edges <- reticulate::r_to_py(as.logical(snip_edges))


    # Check if compute_kaldi_pitch is available (deprecated in torchaudio 2.9+)
    tryCatch({
      reticulate::py_run_string("import torch\
import gc\
import torchaudio\
import torchaudio.functional as F \
import torchaudio.transforms as T \
import math \
\
# Check if compute_kaldi_pitch is available\
if not hasattr(F, 'compute_kaldi_pitch'):\
    raise AttributeError('torchaudio.functional.compute_kaldi_pitch has been removed in torchaudio 2.9+. Please use an alternative pitch tracker (crepe, pyin, swipe, rapt, or reaper) or downgrade torchaudio: pip install torchaudio<2.9')\
\
metadata = torchaudio.info(soundFile) \
\
if  beginTime > 0 and endTime > beginTime : \
	startSample = math.floor(beginTime * metadata.sample_rate) \
else: \
	startSample = 0 \
\
if  endTime > 0 and endTime > beginTime :\
	endSample = math.ceil(endTime * metadata.sample_rate)\
	nSamples = endSample - startSample \
else: \
	nSamples = -1 \
\
SPEECH_WAVEFORM, SAMPLE_RATE = torchaudio.load(soundFile,frame_offset=startSample, num_frames= nSamples) \
\
pitch_feature = F.compute_kaldi_pitch(waveform=SPEECH_WAVEFORM,\
	sample_rate=SAMPLE_RATE, \
	frame_length= windowSize, \
	frame_shift= windowShift, \
	min_f0= minF, \
	max_f0= maxF, \
	soft_min_f0= softMinF0, \
	penalty_factor= voiced_voiceless_cost, \
	lowpass_cutoff= lowpass_cutoff, \
	resample_frequency= resample_frequency, \
	delta_pitch= deltaChange, \
	nccf_ballast= nccfBallast, \
	lowpass_filter_width= lowpass_filter_width, \
	upsample_filter_width= upsample_filter_width, \
	max_frames_latency= max_frames_latency, \
	frames_per_chunk= frames_per_chunk, \
	simulate_first_pass_online= simulate_first_pass_online, \
	recompute_frame= recompute_frame, \
	snip_edges=snip_edges)\
pitch, nccf = pitch_feature[..., 1], pitch_feature[..., 0] \
nppitch = pitch.numpy()\
npnccf = nccf.numpy()\
end_time = SPEECH_WAVEFORM.shape[1] / SAMPLE_RATE\
del pitch\
del nccf\
del pitch_feature\
del SPEECH_WAVEFORM\
gc.collect()")
    }, error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("compute_kaldi_pitch.*removed|no attribute.*compute_kaldi_pitch", msg, ignore.case = TRUE)) {
        stop(
          "kaldi_pitch() is unavailable: torchaudio.functional.compute_kaldi_pitch has been removed in torchaudio 2.9+.\n\n",
          "Options:\n",
          "  1. Use alternative pitch trackers: crepe(), pyin(), swipe(), rapt(), or reaper()\n",
          "  2. Downgrade torchaudio: pip install 'torchaudio<2.9'\n\n",
          "For more information, see: https://github.com/pytorch/audio/issues/3902",
          call. = FALSE
        )
      } else {
        stop(msg, call. = FALSE)
      }
    })


    inTable <- data.frame( "f0" = as.vector(py$nppitch),
                           "nccf"=as.vector(py$npnccf))

    startTime = windowShift

    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16", "REAL32")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.

    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate

    attr(outDataObj, "origFreq") <-  as.numeric(py$SAMPLE_RATE)
    startTime <- 1/sampleRate
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"

    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2) # == binary

    # Cross-correlation track
    f0Table <- inTable %>%
      dplyr::select(f0) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))


    nof0Values <- nrow(f0Table)
    names(f0Table) <- NULL
    outDataObj = addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")

    # Normalized Cross Correlation Function
    nccfTable <- inTable %>%
      dplyr::select(nccf) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))

    noNCCFValues <- nrow(nccfTable)
    names(nccfTable) <- NULL
    outDataObj = addTrack(outDataObj, "nccf", as.matrix(nccfTable[,1]), "REAL32")



    #return(outDataObj)
    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R

    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){

      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))

      missing_f0_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$f0))
      missing_nccf_vals = matrix(0,
                                  nrow = nr_of_missing_samples,
                                  ncol = ncol(outDataObj$nccf))

      # prepend values
      outDataObj$f0 = rbind(missing_f0_vals, outDataObj$f0)
      outDataObj$nccf = rbind(missing_nccf_vals, outDataObj$nccf)


      # fix start time
      attr(outDataObj, "startTime") = startTime - nr_of_missing_samples * (1/sampleRate)
    }

    assertthat::assert_that(is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by the kaldi_pitch function is invalid.")

    ssff_file <- sub("wav$",explicitExt,origSoundFile)
    if(!is.null(outputDirectory)){
      ssff_file <- file.path(outputDirectory,basename(ssff_file))
    }

    attr(outDataObj,"filePath") <- as.character(ssff_file)
    if(toFile){
      write.AsspDataObj(dobj=outDataObj,file=ssff_file)
      #Here we can be sure that the list is a valid SSFF object, so the
      # so we add TRUE to the out vector
      outListOfFiles <- c(listOfFiles,TRUE)
    }
  }
  if(toFile){
    return(length(outListOfFiles))
  }else{
    return(outDataObj)
  }

}



attr(kaldi_pitch,"ext") <-  c("kap")
attr(kaldi_pitch,"tracks") <-  c("f0","nccf")
attr(kaldi_pitch,"outputType") <-  c("SSFF")


