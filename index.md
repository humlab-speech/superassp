# Rationale

[![R-CMD-check](https://github.com/humlab-speech/superassp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/humlab-speech/superassp/actions/workflows/R-CMD-check.yaml)
[![lintr](https://github.com/humlab-speech/superassp/actions/workflows/lintr.yml/badge.svg)](https://github.com/humlab-speech/superassp/actions/workflows/lintr.yml)
[![lint-changed-files](https://github.com/humlab-speech/superassp/actions/workflows/lint-changed-files.yaml/badge.svg)](https://github.com/humlab-speech/superassp/actions/workflows/lint-changed-files.yaml)
[![test-coverage](https://github.com/humlab-speech/superassp/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/humlab-speech/superassp/actions/workflows/test-coverage.yaml)
[![pkgdown](https://github.com/humlab-speech/superassp/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/humlab-speech/superassp/actions/workflows/pkgdown.yaml)

The `superassp` package provides access to an efficient, unified, and
consistent collection of digital speech processing (DSP) routines of
value to speech researchers.

Each function has either a ‘trk\_’ or a ‘lst\_’ prefix to the name,
which indicates the output type as being either a track resulting from
continuous windowed processing of a signal, or as an R list of values.
In addition, all DSP functions have attributes attached to them that
also divulge the track or list field names the user can expect in the
output. Each function also has an associated suggested file extension
that, if used consistently, ensure that applications of multiple DSP
routines to the same speech recording does not risk overwriting each
other.

We aim to provide a consistent naming of formal arguments to functions
so that for the argument determining for instance the time step between
analysis intervals and analysis window size (both millisecond scale,
always) are identically named regardless of the original implmentation’s
naming scheme. Further, we harmonize how portions of a signal are
specified (begin and end times in seconds, always). All DSP functions
return a in-memory representation by default, but when writing to disk
is requested, periodically samples tracks are written in the Simple
Speech Signal File Format (SSFF). Lists are stored in a JSON-based
format.

This package can be seen as the successor of the “Advanced Speech Signal
Processor” (libassp) (and `wrassp` packages) , and incorporates the
libassp DSP functions as efficiently as possible. However, `superassp`
also incorporates routines from several other code bases such as [Speech
Signal Processing Toolkit](https://sp-tk.sourceforge.net) (SPTK),
[Edinburgh Speech Tools
Library](https://www.cstr.ed.ac.uk/projects/speech_tools/manual-1.2.0/)
(ESTK), the [openSMILE](https://github.com/audeering/opensmile) audio
feature extractor package, [The Snack Sound
Toolkit](https://github.com/scottypitcher/tcl-snack), and smaller
specialised libraries. Routines that use the functionality of Praat to
do the signal processing use the `pladdrr` R package to do so.

All wrapper functions support:

- Any media format via the [av](https://github.com/ropensci/av) package
  (WAV, MP3, MP4, MKV, AVI, etc.) if the file format is not natively
  supported by the routine.
- Output to SSFF files (`toFile = TRUE`) or in-memory `AsspDataObj`
  (`toFile = FALSE`)
- Automatic parallel processing for batch operations

## Continuous Integration and Validation

This repository uses a staged CI setup:

- **Fast PR feedback**: `R-CMD-check`, `lintr`, and `lint-changed-files`
- **Deeper main-branch validation**: `test-coverage`, `pkgcheck`, and
  `pkgdown`

Recommended protected-branch required checks:

- `R-CMD-check`
- `lintr`
- `lint-changed-files`

Coverage and pkgcheck are quality signals for ongoing package health and
can be kept advisory or made required as branch policy evolves.

## Installation

All functions are fully self-contained (bundled C/C++ libraries; no
external installs, including for Praat-backed functions such as
`trk_formant_burg`, `trk_praatsauce`, `lst_pharyngeal`,
`lst_voice_tremor`, `lst_voice_report`, and other `pladdrr`-based
wrappers). Praat’s own C++ source is vendored directly into the
`pladdrr` dependency and compiled into it — no separate Praat
installation is required.

Install the package using

``` r

install.packages("devtools") # If not installed already
devtools::install_github("humlab-speech/superassp",dependencies = "Imports")
```

## Quick Start: Pitch Tracking Examples

The SPTK C++ wrapper functions (`trk_pitch_rapt`, `trk_pitch_swipe`,
`trk_pitch_reaper`, `trk_pitch_dio`, `trk_pitch_harvest`) provide the
easiest way to extract F0 from any media file:

``` r
library(superassp)
wfile <- system.file("samples","sustained","a1.wav",package="superassp")

# Extract F0 from a WAV file
f0_data <- trk_pitch_rapt(wfile, toFile = FALSE)

# Extract F0 from video (audio automatically extracted)
f0_data <- trk_pitch_rapt(system.file("samples","sustained","a1.mp4",package="superassp"), toFile = FALSE, minF = 75, maxF = 300)

f0_data <- trk_pitch_swipe(wfile, toFile = FALSE, minF = 75, maxF = 300)

f0_data <- trk_pitch_dio(wfile, toFile = FALSE)

f0_data <- trk_pitch_harvest(wfile, toFile = FALSE)

f0_data <- trk_pitch_crepe(wfile, toFile = FALSE)

f0_data <- trk_pitch_yin(wfile, toFile = FALSE)

f0_data <- trk_pitch_pyin(wfile, toFile = FALSE)

f0_data <- trk_pitch_cc(wfile, toFile = FALSE)



mb <- microbenchmark::microbenchmark(
  "SWIPE" =trk_pitch_swipe(wfile, toFile = FALSE, minF = 75, maxF = 300),
  "DIO" = trk_pitch_dio(wfile, toFile = FALSE),
  "Harvest" trk_pitch_harvest(wfile, toFile = FALSE),
"CREPE"<- trk_pitch_crepe(wfile, toFile = FALSE),
"YIN"= trk_pitch_yin(wfile, toFile = FALSE),
"pYIN"= trk_pitch_pyin(wfile, toFile = FALSE),
 "Praat Pitch"=trk_pitch_cc(wfile, toFile = FALSE)
,times=1L)


```

The f_(o) tracking functions all support:

- Time windowing with `beginTime` and `endTime`
- Custom f₀ range with `minF` and `maxF`
- Frame shift control with `windowShift` (milliseconds)
- Voicing threshold adjustment with `voicing_threshold`

``` r
```

## Formant tracking

The `superassp` package also makes several methods for formant
identification and quantification (in terms of frequencies and
bandwidths) available.

``` r

library(superassp)
wfile <- system.file("samples","sustained","a1.wav",package="superassp")

# Extract formant information from from a WAV file using the libassp `forest` function
fm_data <- trk_formant_forest(system.file("samples","sustained","a1.wav",package="superassp"), toFile = FALSE)

# Extract formant information from video (audio automatically extracted)
fm_data <- trk_formant_forest(system.file("samples","sustained","a1.mp4",package="superassp"), toFile = FALSE)

# Use Praat's Burg algorithm through the `pladdrr`package
fm_data <- trk_formant_burg(system.file("samples","sustained","a1.wav",package="superassp"), toFile = FALSE)


fm_data <- trk_deepformant(system.file("samples","sustained","a1.wav",package="superassp"), toFile = FALSE)
```

## Choosing an algorithm

All pitch trackers emit a fundamental-frequency track (0 marks unvoiced
frames); they differ in method, speed, and robustness:

| Function | Method / origin | When to reach for it |
|----|----|----|
| `trk_pitch_rapt` | RAPT autocorrelation (SPTK) | Fast, robust general-purpose default |
| `trk_pitch_swipe` | SWIPE′ (SPTK) | Accurate, good on noisy speech |
| `trk_pitch_reaper` | REAPER (Google) | Also yields voicing / GCI information |
| `trk_pitch_yin` / `trk_pitch_pyin` | YIN / probabilistic YIN | Classic tracker; pYIN gives smoother voicing |
| `trk_pitch_crepe` | CREPE deep neural net (ONNX) | Most accurate; heavier (downloads a model) |
| `trk_pitch_srh` | Summation of Residual Harmonics | Strong under additive noise |
| `trk_pitch_dio` / `trk_pitch_harvest` | DIO / Harvest (WORLD) | Fast (DIO) vs. accurate (Harvest) |
| `trk_pitch_ac` / `trk_pitch_cc` / `trk_pitch_shs` / `trk_pitch_spinet` | Praat | Praat-faithful behaviour |
| `trk_pitch_ksv` / `trk_pitch_mhs` | ASSP (bundled C) | No external dependencies |

Formant trackers:

| Function | Method / origin | When to reach for it |
|----|----|----|
| `trk_formant_burg` | Praat Burg LPC (pladdrr) | Praat-faithful default |
| `trk_formant_forest` | ASSP forest tracker | Bundled C, no external deps |
| `trk_formant_snack` | Snack ESPS-style LPC | Familiar to Snack/WaveSurfer users |
| `trk_formant_cgdzp` | Complex-group-delay (COVAREP) | Robust bandwidth estimates |
| `trk_formant_tvwlp` | Time-varying weighted LP | Better on rapid transitions |
| `trk_formant_deepformants` / `trk_formant_formantnet` | Deep neural nets (ONNX) | Highest accuracy; heavier |

Start with `trk_pitch_rapt` and `trk_formant_burg`/`trk_formant_forest`;
move to the neural or SWIPE′/SRH options when accuracy or
noise-robustness demands it. See
[`vignette("getting_started", package = "superassp")`](https://humlab-speech.github.io/superassp/articles/getting_started.md)
for a fuller tour.

# Other packages

This package was heavilly inspired by the
[wrassp](https://github.com/IPS-LMU/wrassp) package, that the import of
libassp C code and R code from that package is acknowledged.

Other code bases on wich this package was built:

- [Speech Signal Processing Toolkit](https://sp-tk.sourceforge.net)
- [Edinburgh Speech Tools
  Library](https://www.cstr.ed.ac.uk/projects/speech_tools/manual-1.2.0/)
- [openSMILE](https://github.com/audeering/opensmile)
- [The Snack Sound Toolkit](https://github.com/scottypitcher/tcl-snack)
- [av](https://github.com/ropensci/av)
