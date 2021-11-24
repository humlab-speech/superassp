#' Formant estimation (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::forest
#' @inherit wrassp::forest return references description details sections seealso
#' @seealso [wrassp::forest]
forest = wrassp::forest
attr(forest,"ext") <-  wrassp::wrasspOutputInfos[["forest"]]$ext 
attr(forest,"tracks") <-  wrassp::wrasspOutputInfos[["forest"]]$tracks
attr(forest,"outputType") <-  wrassp::wrasspOutputInfos[["forest"]]$outputType


#' Analysis of short-term autocorrelation function (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::acfana
#' @inherit wrassp::acfana return references description details sections seealso
#' @seealso [wrassp::acfana]
acfana = wrassp::acfana
attr(acfana,"ext") <-  wrassp::wrasspOutputInfos[["acfana"]]$ext 
attr(acfana,"tracks") <-  wrassp::wrasspOutputInfos[["acfana"]]$tracks
attr(acfana,"outputType") <-  wrassp::wrasspOutputInfos[["acfana"]]$outputType

#' Computes the first difference of the signal (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::afdiff
#' @inherit wrassp::afdiff return references description details sections seealso
#' @seealso [wrassp::afdiff]
afdiff = wrassp::afdiff
attr(afdiff,"ext") <-  wrassp::wrasspOutputInfos[["afdiff"]]$ext 
attr(afdiff,"tracks") <-  wrassp::wrasspOutputInfos[["afdiff"]]$tracks
attr(afdiff,"outputType") <-  wrassp::wrasspOutputInfos[["afdiff"]]$outputType


#' Filters the audio signal (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::affilter
#' @inherit wrassp::affilter return references description details sections seealso
#' @seealso [wrassp::affilter]
affilter = wrassp::affilter
attr(affilter,"ext") <-  wrassp::wrasspOutputInfos[["affilter"]]$ext 
attr(affilter,"tracks") <-  wrassp::wrasspOutputInfos[["affilter"]]$tracks
attr(affilter,"outputType") <-  wrassp::wrasspOutputInfos[["affilter"]]$outputType

#' Short-term cepstral analysis (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::cepstrum
#' @inherit wrassp::cepstrum return references description details sections seealso
#' @seealso [wrassp::cepstrum]
cepstrum = wrassp::cepstrum
attr(cepstrum,"ext") <-  wrassp::wrasspOutputInfos[["cepstrum"]]$ext 
attr(cepstrum,"tracks") <-  wrassp::wrasspOutputInfos[["cepstrum"]]$tracks
attr(cepstrum,"outputType") <-  wrassp::wrasspOutputInfos[["cepstrum"]]$outputType

#' Cepstral smoothed version of 'dftSpectrum' (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::cssSpectrum
#' @inherit wrassp::cssSpectrum return references description details sections seealso
#' @seealso [wrassp::cssSpectrum]
cssSpectrum = wrassp::cssSpectrum
attr(cssSpectrum,"ext") <-  wrassp::wrasspOutputInfos[["cssSpectrum"]]$ext 
attr(cssSpectrum,"tracks") <-  wrassp::wrasspOutputInfos[["cssSpectrum"]]$tracks
attr(cssSpectrum,"outputType") <-  wrassp::wrasspOutputInfos[["cssSpectrum"]]$outputType

#' Short-term DFT spectral analysis (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::dftSpectrum
#' @inherit wrassp::dftSpectrum return references description details sections seealso
#' @seealso [wrassp::dftSpectrum]
dftSpectrum = wrassp::dftSpectrum
attr(dftSpectrum,"ext") <-  wrassp::wrasspOutputInfos[["dftSpectrum"]]$ext 
attr(dftSpectrum,"tracks") <-  wrassp::wrasspOutputInfos[["dftSpectrum"]]$tracks
attr(dftSpectrum,"outputType") <-  wrassp::wrasspOutputInfos[["dftSpectrum"]]$outputType

#' F0 analysis of the signal (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::ksvF0
#' @inherit wrassp::ksvF0 return references description details sections seealso
#' @seealso [wrassp::ksvF0]
ksvF0 = wrassp::ksvF0
attr(ksvF0,"ext") <-  wrassp::wrasspOutputInfos[["ksvF0"]]$ext 
attr(ksvF0,"tracks") <-  wrassp::wrasspOutputInfos[["ksvF0"]]$tracks
attr(ksvF0,"outputType") <-  wrassp::wrasspOutputInfos[["ksvF0"]]$outputType

#' Linear Predictive smoothed version of 'dftSpectrum' (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::lpsSpectrum
#' @inherit wrassp::lpsSpectrum return references description details sections seealso
#' @seealso [wrassp::lpsSpectrum]
lpsSpectrum = wrassp::lpsSpectrum
attr(lpsSpectrum,"ext") <-  wrassp::wrasspOutputInfos[["lpsSpectrum"]]$ext 
attr(lpsSpectrum,"tracks") <-  wrassp::wrasspOutputInfos[["lpsSpectrum"]]$tracks
attr(lpsSpectrum,"outputType") <-  wrassp::wrasspOutputInfos[["lpsSpectrum"]]$outputType

#' Pitch analysis of the speech signal using Michel's (M)odified (H)armonic (S)ieve algorithm (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::mhsF0
#' @inherit wrassp::mhsF0 return references description details sections seealso
#' @seealso [wrassp::mhsF0]
mhsF0 = wrassp::mhsF0
attr(mhsF0,"ext") <-  wrassp::wrasspOutputInfos[["mhsF0"]]$ext 
attr(mhsF0,"tracks") <-  wrassp::wrasspOutputInfos[["mhsF0"]]$tracks
attr(mhsF0,"outputType") <-  wrassp::wrasspOutputInfos[["mhsF0"]]$outputType


#' Linear Prediction analysis (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::rfcana
#' @inherit wrassp::rfcana return references description details sections seealso
#' @seealso [wrassp::rfcana]
rfcana = wrassp::rfcana
attr(rfcana,"ext") <-  wrassp::wrasspOutputInfos[["rfcana"]]$ext 
attr(rfcana,"tracks") <-  wrassp::wrasspOutputInfos[["rfcana"]]$tracks
attr(rfcana,"outputType") <-  wrassp::wrasspOutputInfos[["rfcana"]]$outputType

#' Analysis of short-term Root Mean Square amplitude (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::rmsana
#' @inherit wrassp::rmsana return references description details sections seealso
#' @seealso [wrassp::rmsana]
rmsana = wrassp::rmsana
attr(rmsana,"ext") <-  wrassp::wrasspOutputInfos[["rmsana"]]$ext 
attr(rmsana,"tracks") <-  wrassp::wrasspOutputInfos[["rmsana"]]$tracks
attr(rmsana,"outputType") <-  wrassp::wrasspOutputInfos[["rmsana"]]$outputType

#' Analysis of the averages of the short-term positive and negative zero-crossing rates (From the 'wrassp' package)
#' @export
#' @inheritParams wrassp::zcrana
#' @inherit wrassp::zcrana return references description details sections seealso
#' @seealso [wrassp::zcrana]
zcrana = wrassp::zcrana
attr(zcrana,"ext") <-  wrassp::wrasspOutputInfos[["zcrana"]]$ext 
attr(zcrana,"tracks") <-  wrassp::wrasspOutputInfos[["zcrana"]]$tracks
attr(zcrana,"outputType") <-  wrassp::wrasspOutputInfos[["zcrana"]]$outputType
