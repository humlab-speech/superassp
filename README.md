
# Rationale

The `superassp` package provides access to an efficient, unified, and consistent collection of digital speech processing (DSP) routines of value to speech researchers.

Each function has either a 'trk_' or a 'lst_' prefix to the name, which indicates the output type as being either a track resulting from continous windowed processing of a signal, or as an R list of values. In addition, all DSP functions have attributes attached to them that also divulge the track or list field names the user can expect in the output. Each function also has an assocaiated suggested file extension that, if used consistently, ensure that applications of multiple DSP routines to the same speech recording does not risk overwriting each other. 

We aim to provide a consistent naming of formal arguments to functions so that for the argument determining for instance the time step between analysis intervals and analysis window size (both millisecond scale, always) are identically named regardless of the original implmentation's naming scheme. Further, we harmonize how portions of a signal are specified (begin and end times in seconds, always). All DSP functions return a in-memory representation by default, but when writing to disk is requested, periodically samples tracks are written in the Simple Speech Signal File Format (SSFF). Lists are stored in a JSON-based format.



This package can be seen as the succuessor of the "Advanced Speech Signal Processor" (libassp) (and `wrassp` packages) , and incorporates the libassp DSP functions as efficiently as possible. However, `superassp` also inporporates routines from several other code bases such as [Speech Signal Processing Toolkit](https://sp-tk.sourceforge.net) (SPTK), [Edinburgh Speech Tools Library](https://www.cstr.ed.ac.uk/projects/speech_tools/manual-1.2.0/) (ESTK), the [openSMILE](https://github.com/audeering/opensmile) audil feature extractor package, [The Snack Sound Toolkit](https://github.com/scottypitcher/tcl-snack), and smaller specialised libraries (sometimes through a python call). 
Routines that use the functionality of Praat to do the signal processing use the `pladdrr` R package to do so.

All wrapper functions support:

- Any media format via the [av](https://github.com/ropensci/av) package (WAV, MP3, MP4, MKV, AVI, etc.) if the file format is not natively supported by the routine.
- Output to SSFF files (`toFile = TRUE`) or in-memory `AsspDataObj` (`toFile = FALSE`)
- Automatic parallel processing for batch operations



## Installation

The package requires the Praat program to be installed in the user's PATH (or in '/Applications' on Mac OS).

Then simply install the package using
```r
install.packages("devtools") # If not installed already
devtools::install_github("humlab-speech/superassp",dependencies = "Imports")
```

## Quick Start: Pitch Tracking Examples

The SPTK C++ wrapper functions (`trk_rapt`, `trk_swipe`, `trk_reaper`, `trk_dio`, `trk_harvest`) provide the easiest way to extract F0 from any media file:

```r
library(superassp)
wfile <- system.file("samples","sustained","a1.wav",package="superassp")

# Extract F0 from a WAV file
f0_data <- trk_rapt(wfile, toFile = FALSE)

# Extract F0 from video (audio automatically extracted)
f0_data <- trk_rapt(system.file("samples","sustained","a1.mp4",package="superassp"), toFile = FALSE, minF = 75, maxF = 300)

f0_data <- trk_swipe(wfile, toFile = FALSE, minF = 75, maxF = 300)

f0_data <- trk_dio(wfile, toFile = FALSE)

f0_data <- trk_harvest(wfile, toFile = FALSE)

f0_data <- trk_crepe(wfile, toFile = FALSE)

f0_data <- trk_yin(wfile, toFile = FALSE)

f0_data <- trk_pyin(wfile, toFile = FALSE)

f0_data <- trk_pitch_cc(wfile, toFile = FALSE)



mb <- microbenchmark::microbenchmark(
  "SWIPE" =trk_swipe(wfile, toFile = FALSE, minF = 75, maxF = 300),
  "DIO" = trk_dio(wfile, toFile = FALSE),
  "Harvest" trk_harvest(wfile, toFile = FALSE),
"CREPE"<- trk_crepe(wfile, toFile = FALSE),
"YIN"= trk_yin(wfile, toFile = FALSE),
"pYIN"= trk_pyin(wfile, toFile = FALSE),
 "Praat Pitch"=trk_pitch_cc(wfile, toFile = FALSE)
,times=1L)





```

The f~o~ tracking functions all support:

- Time windowing with `beginTime` and `endTime`
- Custom f~0~ range with `minF` and `maxF`
- Frame shift control with `windowShift` (milliseconds)
- Voicing threshold adjustment with `voicing_threshold`

```r 

```

## Formant tracking

The `superassp` package also makes several methods for formant identification and quantification (in terms of frequencies and bandwidths) available.



```r
library(superassp)
wfile <- system.file("samples","sustained","a1.wav",package="superassp")

# Extract formant information from from a WAV file using the libassp `forest` function
fm_data <- trk_forest(system.file("samples","sustained","a1.wav",package="superassp"), toFile = FALSE)

# Extract formant information from video (audio automatically extracted)
fm_data <- trk_forest(system.file("samples","sustained","a1.mp4",package="superassp"), toFile = FALSE)

# Use Praat's Burg algorithm through the `pladdrr`package
fm_data <- trk_formant(system.file("samples","sustained","a1.wav",package="superassp"), toFile = FALSE)


fm_data <- trk_deepformant(system.file("samples","sustained","a1.wav",package="superassp"), toFile = FALSE)


```



# Other packages

This package was heavilly inspired by the [wrassp](https://github.com/IPS-LMU/wrassp) package, that the import of libassp C code and R code from that package is acknowledged.

Other code bases on wich this package was built:

 - [Speech Signal Processing Toolkit](https://sp-tk.sourceforge.net)
 - [Edinburgh Speech Tools Library](https://www.cstr.ed.ac.uk/projects/speech_tools/manual-1.2.0/)
 - [openSMILE](https://github.com/audeering/opensmile)
 - [The Snack Sound Toolkit](https://github.com/scottypitcher/tcl-snack)
 - [av](https://github.com/ropensci/av)