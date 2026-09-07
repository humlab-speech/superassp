# Track short-term RMS amplitude

Computes the short-term root-mean-square (RMS) amplitude of audio
signals using the *libassp* C library (Scheffers 2012) . By default,
output is expressed in dB (short-term power). A fast, low-memory energy
measure suitable as a voicing or loudness feature.

## Usage

``` r
trk_rms(listOfFiles, beginTime = 0, centerTime = FALSE, endTime = 0, windowShift = 5, windowSize = 20, effectiveLength = TRUE, linear = FALSE, window = "HAMMING", toFile = TRUE, explicitExt = "rms", outputDirectory = NULL, assertLossless = NULL, logToFile = FALSE, convertOverwrites = FALSE, keepConverted = FALSE, verbose = TRUE)
```

## Arguments

- listOfFiles:

  Character vector of audio file paths. Any format supported by av is
  accepted; non-native inputs are transcoded automatically.

- linear:

  Logical. If `TRUE`, return linear RMS amplitude instead of dB. Default
  `FALSE` (dB scale).

- beginTime:

  Start time for the extracted portion in seconds. Default: NULL
  (beginning of signal). Note: uses `beginTime`/`endTime` (seconds)
  matching DSP function conventions, unlike
  [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  which uses `begin`/`end`.

- centerTime:

  Numeric or logical. Single-frame analysis time point in seconds;
  overrides `beginTime`, `endTime`, and `windowShift`. Default `FALSE`.

- endTime:

  The end time of the section of the sound files that should be analysed
  (in seconds). Use 0 for end of file.

- windowShift:

  Numeric. Frame shift in milliseconds; sets output frame rate
  (`1000 / windowShift` Hz). Default 5.0 ms (200 Hz). Must be strictly
  less than 32 ms (the 512-sample analysis window at 16 kHz). Values
  other than the training default (5 ms) may slightly reduce accuracy.

- windowSize:

  Numeric. Smoothing filter window size in milliseconds, applied to both
  median (periodicity) and mean (F0) post-processing filters. Default 15
  ms.

- effectiveLength:

  Logical. Make window size effective rather than exact. Default
  `FALSE`.

- window:

  Character. Analysis window function type. Default `"BLACKMAN"`. See
  [AsspWindowTypes](https://humlab-speech.github.io/superassp/reference/AsspWindowTypes.md)
  for supported types.

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

- assertLossless:

  Character vector of additional file extensions to treat as losslessly
  encoded.

- logToFile:

  Logical. Write processing log to a file in `outputDirectory` rather
  than the console. Default `FALSE`.

- convertOverwrites:

  Logical. Allow transcoding to overwrite existing files. Default
  `FALSE`.

- keepConverted:

  Logical. Retain intermediate transcoded files. Default `FALSE`.

- verbose:

  Logical. Show a progress bar (sequential path) or a progress-aware
  parallel apply (`pbapply`/`pbmcapply`, if installed).

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `RMS[dB]`:

  REAL32, dB (or linear amplitude if `linear = TRUE`), n_frames x 1
  column. Short-term RMS amplitude per frame.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

Window function defaults to HAMMING (not BLACKMAN as in other libassp
functions). Set `linear = TRUE` for raw amplitude values without log
conversion.

## See also

[wrassp::rmsana](https://rdrr.io/pkg/wrassp/man/rmsana.html)
