# Track LP-smoothed spectrum

Computes a short-term spectral envelope smoothed by linear prediction
using the *libassp* C library (Scheffers 2012) . LP smoothing produces a
parametric model-based envelope; prefer `trk_lps_spectrum` over
`trk_css_spectrum` when the LP order is known and explicit control over
the number of resonances is desired.

## Usage

``` r
trk_lps_spectrum(
  listOfFiles,
  beginTime = 0,
  centerTime = FALSE,
  endTime = 0,
  resolution = 40,
  fftLength = 0,
  windowSize = 20,
  windowShift = 5,
  window = "BLACKMAN",
  order = 0,
  preemphasis = -0.95,
  deemphasize = TRUE,
  toFile = TRUE,
  explicitExt = "lps",
  outputDirectory = NULL,
  assertLossless = NULL,
  logToFile = FALSE,
  keepConverted = FALSE,
  convertOverwrites = FALSE,
  verbose = TRUE
)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- beginTime:

  Numeric. Start of analysis window in seconds. Default 0 (file start).

- centerTime:

  Numeric or logical. Single-frame analysis time point in seconds;
  overrides `beginTime`, `endTime`, and `windowShift`. Default `FALSE`.

- endTime:

  Numeric. End of analysis window in seconds. Default 0 (file end).

- resolution:

  Numeric. Target FFT frequency resolution in Hz. Default 40.0.

- fftLength:

  Integer. Explicit FFT length in points; overrides `resolution`.
  Default 0 (use `resolution`).

- windowSize:

  Numeric. Analysis window size in milliseconds. Default 20 ms.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate (1000 /
  windowShift Hz). Default 5 ms.

- window:

  Character. Analysis window function type. Default `"BLACKMAN"`. See
  [AsspWindowTypes](https://humlab-speech.github.io/superassp/reference/AsspWindowTypes.md).

- order:

  Integer. LP prediction order; 0 defaults to sample rate in kHz + 3.
  Default 0.

- preemphasis:

  Numeric. Pre-emphasis factor applied before LP analysis. Default
  -0.95.

- deemphasize:

  Logical. Undo the spectral tilt introduced by pre-emphasis when
  computing the output spectrum. Default `TRUE`.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written (invisibly). If `FALSE`, return an `AsspDataObj`. Default
  `TRUE`.

- explicitExt:

  Character. Output file extension. Default `"lps"`.

- outputDirectory:

  Character. Directory for output files. `NULL` (default) writes
  alongside the input file.

- assertLossless:

  Character vector of additional file extensions to treat as losslessly
  encoded.

- logToFile:

  Logical. Write processing log to a file in `outputDirectory` rather
  than the console. Default `FALSE`.

- keepConverted:

  Logical. Retain intermediate transcoded files. Default `FALSE`.

- convertOverwrites:

  Logical. Allow transcoding to overwrite existing files. Default
  `FALSE`.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `LPS[dB]`:

  REAL32, dB, n_frames x (FFT_length/2 + 1) columns. LP-smoothed
  spectral amplitude from 0 Hz to the Nyquist rate.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

`order = 0` sets the LP order to sample rate in kHz + 3. When
`deemphasize = TRUE` (default), the output spectrum is corrected for the
pre-emphasis spectral tilt applied before LP analysis.

## See also

[wrassp::lpsSpectrum](https://rdrr.io/pkg/wrassp/man/lpsSpectrum.html)

[AsspWindowTypes](https://humlab-speech.github.io/superassp/reference/AsspWindowTypes.md)

[av::av_audio_convert](https://docs.ropensci.org/av//reference/encoding.html)

## Author

Raphael Winkelmann

Lasse Bombien

Fredrik Nylén

## Examples

``` r
# get path to audio file
path2wav <- list.files(
   system.file("samples", "sustained", package = "superassp"),
   pattern = glob2rx("a1.wav"), full.names = TRUE)

# calculate linear prediction smoothed spectrum
res <- trk_lps_spectrum(path2wav, toFile=FALSE)
#> Applying `method(trk_lps_spectrum, class_character)()` to 1 recording
resolution <- attr(res,"origFreq") / ncol(res[[1]])

# plot spectral values at midpoint of signal
plot(y=res[["CSS[dB]"]][400,],
    x=seq(1,ncol(res[[1]]),1)* resolution,
    type='l',
    xlab='Frequency (Hz)',
    ylab='Amplitude (dB)')

```
