# Track short-term cepstral coefficients

Computes the short-term real cepstrum of audio signals via FFT using the
*libassp* C library (Scheffers 2012; Oppenheim and Schafer 2004;
Childers et al. 1977) . Output coefficients index quefrency from 0 to
FFT_length/2 in steps of 0.5 ms (at the default 40 Hz resolution).
Useful for pitch-period detection and spectral tilt estimation.

## Usage

``` r
trk_cepstrum(
  listOfFiles,
  beginTime = 0,
  centerTime = FALSE,
  endTime = 0,
  resolution = 40,
  fftLength = 0,
  windowShift = 5,
  window = "BLACKMAN",
  toFile = TRUE,
  explicitExt = "cep",
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

- resolution:

  Numeric. Target FFT frequency resolution in Hz; the FFT length is set
  to the smallest power-of-2 meeting this target. Default 40.0.

- fftLength:

  Integer. Explicit FFT length in points; overrides `resolution`.
  Default 0 (use `resolution`).

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

- keepConverted:

  Logical. Retain intermediate transcoded files. Default `FALSE`.

- convertOverwrites:

  Logical. Allow transcoding to overwrite existing files. Default
  `FALSE`.

- verbose:

  Logical. Show a progress bar (sequential path) or a progress-aware
  parallel apply (`pbapply`/`pbmcapply`, if installed).

## Value

If `toFile = FALSE`: an `AsspDataObj` with track:

- `C[dB]`:

  REAL32, dB amplitude, n_frames x (FFT_length/2 + 1) columns. Each
  column corresponds to a quefrency of col_index × 0.5 ms.

Frame rate: `1000 / windowShift` Hz (default 200 Hz). If
`toFile = TRUE`: integer count of files written, returned invisibly.

## Details

The number of coefficients per frame equals FFT_length / 2 + 1
(one-sided, not mirrored). Use `fftLength` to fix the transform size;
otherwise `resolution` governs it. Pre-emphasis is not applied by this
function; apply `trk_afdiff` first if needed.

## References

Childers DG, Skinner DP, Kemerait RC (1977). “The cepstrum: A guide to
processing.” *Proceedings of the IEEE*, **65**(10), 1428–1443. ISSN
0018-9219.
[doi:10.1109/proc.1977.10747](https://doi.org/10.1109/proc.1977.10747)
.  
  
Oppenheim AV, Schafer RW (2004). “From frequency to quefrency: a history
of the cepstrum.” *IEEE Signal Processing Magazine*, **21**(5), 95–106.
ISSN 1053-5888.
[doi:10.1109/msp.2004.1328092](https://doi.org/10.1109/msp.2004.1328092)
.  
  
Scheffers M (2012). “Advanced Speech Signal Processor.”
<https://sourceforge.net/projects/libassp/files/libassp/>.

## See also

[wrassp::cepstrum](https://rdrr.io/pkg/wrassp/man/cepstrum.html)

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

# calulate cepstrum
res <- trk_cepstrum(path2wav, toFile=FALSE)
#> Applying `method(trk_cepstrum, class_character)()` to 1 recording

# plot cepstral values at midpoint of signal
plot(y=res[["C[dB]"]][400,],
    x=seq(1,ncol(res[["C[dB]"]])),
    type='l',
    xlab='Quefrency (ms)',
    ylab='Amplitude (dB)')

```
