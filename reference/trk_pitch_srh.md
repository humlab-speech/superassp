# Track fundamental frequency using the Summation of Residual Harmonics (SRH)

Extracts F0 and a voiced/unvoiced decision using SRH (Drugman & Alwan
2011), a two-pass harmonic-summation pitch estimator operating on the
LPC residual. SRH is robust in noisy conditions and produces an
integrated VAD decision. Audio is resampled to 16 kHz internally. The
fixed 10 ms hop differs from RAPT/SWIPE, which honour the `windowShift`
parameter.

## Usage

``` r
trk_pitch_srh(listOfFiles, beginTime = 0, endTime = 0, minF = 50, maxF = 500, toFile = TRUE, explicitExt = "srh", outputDirectory = NULL, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- beginTime:

  Start time for the extracted portion in seconds. Default: NULL
  (beginning of signal). Note: uses `beginTime`/`endTime` (seconds)
  matching DSP function conventions, unlike
  [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  which uses `begin`/`end`.

- endTime:

  The end time of the section of the sound files that should be analysed
  (in seconds). Use 0 for end of file.

- minF:

  Numeric. Minimum F0 in Hz for the internal pitch estimator. Lower
  values allow lower-pitched voices but may increase false positives.
  Default 40.0 Hz.

- maxF:

  Numeric. Maximum F0 in Hz to treat as voiced. Default 400 Hz (speech).
  Must be \<= 2093.75 Hz (model maximum; C7). For music, use 2093.75.

- toFile:

  Logical. If `TRUE`, write SSFF output files and return the count
  written. If `FALSE`, return an `AsspDataObj` (single file only).
  Default `TRUE`.

- explicitExt:

  By default, a character "d" will be prepended to the file name suffix
  when writing the output to file. The user can also specify an explicit
  extension which will be used instead.

- outputDirectory:

  The directory where the slice file should be stored. If not defiled
  (NULL), the sparse slice file will placed in the same folder as the
  media file.

- verbose:

  Logical. Show a progress bar (sequential path) or a progress-aware
  parallel apply (`pbapply`/`pbmcapply`, if installed).

## Value

If `toFile = FALSE`: an `AsspDataObj` with tracks:

- `f0`:

  REAL32, fundamental frequency in Hz, n_frames × 1. Zero indicates
  unvoiced frames.

- `vad`:

  REAL32, voiced/unvoiced decision (0 = unvoiced, 1 = voiced), n_frames
  × 1.

Frame rate: 100 Hz (fixed 10 ms hop). If `toFile = TRUE`: integer count
of files written, returned invisibly.

## References

(Drugman and Alwan 2011)

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract F0 using SRH
trk_pitch_srh("recording.wav")

# Process with custom F0 range
trk_pitch_srh("speech.wav", minF = 80, maxF = 300)
} # }
```
